#= ============ =#
#  ActorBlanket  #
#= ============ =#

import NNeutronics
using NLopt
using JuMP

mutable struct ActorBlanket <: ReactorAbstractActor
    dd::IMAS.dd
    par::ParametersActor
    blanket_multiplier::Real
    thermal_power_extraction_efficiency::Real
    objective::Symbol
end

function ParametersActor(::Type{Val{:ActorBlanket}})
    par = ParametersActor(nothing)
    par.blanket_multiplier = Entry(Real, "", "Neutron thermal power multiplier in blanket"; default = 1.2)
    par.thermal_power_extraction_efficiency = Entry(
        Real,
        "",
        "Fraction of thermal power that is carried out by the coolant at the blanket interface, rather than being lost in the surrounding strutures.";
        default=1.0
    )
    par.objective = Switch([:leakage, :TBR], "", "Blanket layers optimization objective"; default = :leakage)
    return par
end

"""
    ActorBlanket(dd::IMAS.dd, act::ParametersAllActors; kw...)

Blanket actor

!!! note 
    Stores data in `dd.blanket`
"""
function ActorBlanket(dd::IMAS.dd, act::ParametersAllActors; kw...)
    par = act.ActorBlanket(kw...)
    actor = ActorBlanket(dd, par)
    step(actor, act)
    finalize(actor)
    return actor
end

function ActorBlanket(dd::IMAS.dd, par::ParametersActor; kw...)
    logging_actor_init(ActorBlanket)
    par = par(kw...)
    return ActorBlanket(dd, par, par.blanket_multiplier, par.thermal_power_extraction_efficiency, par.objective)
end

function _step(actor::ActorBlanket, act)
    dd = actor.dd
    empty!(dd.blanket)
    blanket = IMAS.get_build(dd.build, type=_blanket_, fs=_hfs_, raise_error_on_missing=false)
    if blanket === missing
        @warn "No blanket present for ActorBlanket to do anything"
        return actor
    end

    dd = actor.dd
    eqt = dd.equilibrium.time_slice[]
    nnt = dd.neutronics.time_slice[]

    total_power_neutrons = IMAS.fusion_power(dd.core_profiles.profiles_1d[]) .* 4 / 5
    total_power_neutrons = sum(nnt.wall_loading.power)
    total_neutrons_per_second = total_power_neutrons / 2.259E-12  # 14.1 MeV = 2.259E-12 J
    total_power_radiated = 0.0 # IMAS.radiative_power(dd.core_profiles.profiles_1d[])

    fwr = dd.neutronics.first_wall.r
    fwz = dd.neutronics.first_wall.z

    # for all blanket modules
    blankets = [structure for structure in dd.build.structure if structure.type == Int(_blanket_)]
    resize!(dd.blanket.module, length(blankets))
    modules_effective_thickness = []
    modules_wall_loading_power = []
    Li6 = 0.0
    blanket_model_1d = NNeutronics.Blanket()
    for (istructure, structure) in enumerate(blankets)
        bm = dd.blanket.module[istructure]
        bm.name = structure.name

        # get the relevant blanket layers
        if sum(structure.outline.r) / length(structure.outline.r) < eqt.boundary.geometric_axis.r && length(blankets) > 1 # HFS
            fs = _hfs_
        else
            fs = _lfs_
        end
        d1 = IMAS.get_build(dd.build, type=_wall_, fs=fs, raise_error_on_missing=false)
        d2 = IMAS.get_build(dd.build, type=_blanket_, fs=fs)
        d3 = IMAS.get_build(dd.build, type=_shield_, fs=fs, return_only_one=false, raise_error_on_missing=false)
        if d3 !== missing
            # if there are multiple shields we choose the one closest to the plasma
            if fs == _hfs_
                d3 = d3[end]
            else
                d3 = d3[1]
            end
        end

        # assign blanket module layers (designed to handle missing wall and/or missing shield)
        resize!(bm.layer, 3)
        for (kl, dl) in enumerate(bm.layer)
            dl.name = "dummy layer $kl"
            dl.midplane_thickness = 0.0
            dl.material = "vacuum"
        end
        for (kl, dl) in enumerate([d1, d2, d3])
            if dl !== missing
                bm.layer[kl].name = dl.name
                bm.layer[kl].midplane_thickness = dl.thickness
                bm.layer[kl].material = dl.material
            end
        end

        # optimize layer thickensses
        total_leakage, d1_thickness, d2_thickness, d3_thickness, Li6, TBR = optimize_layers(blanket_model_1d, bm.layer[1].midplane_thickness, bm.layer[2].midplane_thickness, bm.layer[3].midplane_thickness, actor.objective)
        bm.layer[1].midplane_thickness, bm.layer[2].midplane_thickness, bm.layer[3].midplane_thickness = d1_thickness, d2_thickness, d3_thickness
        for (kl, dl) in enumerate([d1, d2, d3])
            if dl !== missing
                dl.name = bm.layer[kl].name
                dl.thickness = bm.layer[kl].midplane_thickness 
                dl.material = bm.layer[kl].material
            end
        end

        # identify first wall portion of the blanket module
        tmp = convex_hull(vcat(eqt.boundary.outline.r, structure.outline.r), vcat(eqt.boundary.outline.z, structure.outline.z); closed_polygon=true)
        index = findall(x -> x == 1, [IMAS.PolygonOps.inpolygon((r, z), tmp) for (r, z) in zip(fwr, fwz)])
        istart = argmin(abs.(fwz[index]))
        if fwz[index][istart+1] > fwz[index][istart]
            i = argmin(fwz[index])
            index = vcat(index[i:end], index[1:i-1])
        else
            i = argmin(fwz[index])
            index = vcat(index[i+1:end], index[1:i])
        end

        # calculate effective thickness of the layers based on the direction of the impinging neutron flux
        effective_thickness = zeros(length(index), 3)
        r_coords = zeros(4)
        z_coords = zeros(4)
        for (k, (r0, z0, fr, fz)) in enumerate(zip(fwr[index], fwz[index], nnt.wall_loading.flux_r[index], nnt.wall_loading.flux_z[index]))
            fn = norm(fr, fz)
            r_coords[1] = r0
            z_coords[1] = z0
            for (ilayer, layer) in enumerate([d1, d2, d3])
                if layer === missing
                    hit = []
                else
                    hit = IMAS.intersection([r0, fr / fn * 1000 + r0], [z0, fz / fn * 1000 + z0], layer.outline.r, layer.outline.z)
                end
                if isempty(hit)
                    r_coords[ilayer+1] = r_coords[ilayer]
                    z_coords[ilayer+1] = z_coords[ilayer]
                else
                    r_coords[ilayer+1] = hit[1][1]
                    z_coords[ilayer+1] = hit[1][2]
                end
            end
            effective_thickness[k, :] .= sqrt.(diff(r_coords) .^ 2.0 .+ diff(z_coords) .^ 2.0)
        end

        # assign totals in IDS
        total_neutron_capture_fraction = sum(nnt.wall_loading.power[index]) / total_power_neutrons
        # evaluate radiative_capture_fraction (could be done better, since ratiation source may not be coming from core)
        total_radiative_capture_fraction = total_neutron_capture_fraction
        resize!(bm.time_slice)
        bmt = bm.time_slice[]
        bmt.power_incident_neutrons = total_power_neutrons .* total_neutron_capture_fraction
        bmt.power_incident_radiated = total_power_radiated .* total_radiative_capture_fraction
        bmt.power_thermal_neutrons = bmt.power_incident_neutrons * actor.blanket_multiplier
        bmt.power_thermal_radiated = bmt.power_incident_radiated
        bmt.power_thermal_extracted = actor.thermal_power_extraction_efficiency * (bmt.power_thermal_neutrons + bmt.power_thermal_radiated)

        push!(modules_effective_thickness, effective_thickness)
        push!(modules_wall_loading_power, nnt.wall_loading.power[index])
    end

    # Optimize Li6/Li7 ratio to obtain target TBR
    function target_TBR(blanket_model, Li6, dd, modules_effective_thickness, modules_wall_loading_power, total_power_neutrons, target=nothing)
        total_tritium_breeding_ratio = 0.0
        for (ibm, bm) in enumerate(dd.blanket.module)
            bmt = bm.time_slice[]
            bm.layer[2].material = @sprintf("lithium-lead: Li6/7=%3.3f", Li6)
            bmt.tritium_breeding_ratio = 0.0
            module_wall_loading_power = sum(modules_wall_loading_power[ibm])
            for k in 1:length(modules_wall_loading_power[ibm])
                bmt.tritium_breeding_ratio += (NNeutronics.TBR(blanket_model, modules_effective_thickness[ibm][k, :]..., Li6) * modules_wall_loading_power[ibm][k] / module_wall_loading_power)
            end
            total_tritium_breeding_ratio += bmt.tritium_breeding_ratio * module_wall_loading_power / total_power_neutrons
        end
        if target === nothing
            return total_tritium_breeding_ratio
        else
            return (total_tritium_breeding_ratio - target)^2
        end
    end

    ActorCXbuild(dd, act)
    target_TBR(blanket_model_1d, Li6, dd, modules_effective_thickness, modules_wall_loading_power, total_power_neutrons)

    return actor
end

function optimize_layers(blanket_model, l1, l2, l3, objective = :leakage, target = 1.1)
    energy_grid = NNeutronics.energy_grid()
    modules_thickness = [l1, l2, l3]
    blanket_thickness = sum(modules_thickness)
    
    total_neutron_leakage(d1, d2, d3, Li6) = sum(NNeutronics.leakeage_energy(blanket_model, d1, d2, d3, Li6, energy_grid))
    total_TBR(d1, d2, d3, Li6) = NNeutronics.TBR(blanket_model, d1, d2, d3, Li6)
    function ∇f(g::AbstractVector{T}, x::T, y::T, z::T, a::T) where {T}
        g[1] = -exp(-x)
        g[2] = -exp(-y)
        g[3] = -exp(-z)
        g[4] = a
        return
    end

    function ∇g(g::AbstractVector{T}, x::T, y::T, z::T, a::T) where {T}
        g[1] = -exp(-x)
        g[2] = 1/y
        g[3] = 1/z
        g[4] = a
        return
    end

    model = Model(NLopt.Optimizer)
    empty!(model)
    set_optimizer_attribute(model, "algorithm", :LN_COBYLA)

    register(model, :total_neutron_leakage, 4, total_neutron_leakage, ∇f)
    register(model, :total_TBR, 4, total_TBR, ∇g)

    if l1 == 0.0
        @variable(model, d1 == 0.0)
    else
        @variable(model, d1 >= 0.0)
    end
    if l2 == 0.0
        @variable(model, d2 == 0.0)
    else
        @variable(model, d2 >= 0.0)
    end
    if l3 == 0.0
        @variable(model, d3 == 0.0)
    else
        @variable(model, d3 >= 0.0)
    end

    @variable(model, 0.0 <= Li6 <= 100.0)
    @NLconstraint(model, blanket_thickness == d1 + d2 + d3)

    if objective == :leakage
        @NLconstraint(model, target <= total_TBR(d1, d2, d3, Li6))
        @NLobjective(model, Min, total_neutron_leakage(d1, d2, d3, Li6))
        JuMP.optimize!(model)
        total_leakage = objective_value(model)
        l1, l2, l3, Li6 = (value(d1), value(d2), value(d3), value(Li6))
        TBR = total_TBR(l1, l2, l3, Li6)
    elseif objective == :TBR
        target = 0.25
        @NLconstraint(model, target >= total_neutron_leakage(d1, d2, d3, Li6))
        @NLobjective(model, Max, total_TBR(d1, d2, d3, Li6))
        JuMP.optimize!(model)
        TBR = objective_value(model)
        l1, l2, l3, Li6 = (value(d1), value(d2), value(d3), value(Li6))
        total_leakage = total_neutron_leakage(l1, l2, l3, Li6)
    end

    return total_leakage, l1, l2, l3, Li6, TBR
end