#= ============ =#
#  ActorBlanket  #
#= ============ =#

import NNeutronics
import NLopt
import JuMP

mutable struct ActorBlanket <: ReactorAbstractActor
    dd::IMAS.dd
    par::ParametersActor
    act::ParametersAllActors

    function ActorBlanket(dd::IMAS.dd, par::ParametersActor, act::ParametersAllActors; kw...)
        logging_actor_init(ActorBlanket)
        par = par(kw...)
        return new(dd, par, act)
    end
end

function ParametersActor(::Type{Val{:ActorBlanket}})
    par = ParametersActor(nothing)
    par.blanket_multiplier = Entry(Real, "", "Neutron thermal power multiplier in blanket"; default=1.2)
    par.thermal_power_extraction_efficiency = Entry(
        Real,
        "",
        "Fraction of thermal power that is carried out by the coolant at the blanket interface, rather than being lost in the surrounding strutures.";
        default=1.0
    )
    par.objective = Switch([:leakage, :TBR], "", "Blanket layers optimization objective"; default=:leakage)
    par.verbose = Entry(Bool, "", "verbose"; default=false)
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
    actor = ActorBlanket(dd, par, act)
    step(actor)
    finalize(actor)
    return actor
end

function _step(actor::ActorBlanket)
    dd = actor.dd
    empty!(dd.blanket)
    blanket = IMAS.get_build(dd.build, type=_blanket_, fs=_hfs_, raise_error_on_missing=false)
    if blanket === missing
        @warn "No blanket present for ActorBlanket to do anything"
        return actor
    end

    dd = actor.dd
    par = actor.par
    act = actor.act

    eqt = dd.equilibrium.time_slice[]
    nnt = dd.neutronics.time_slice[]

    total_power_neutrons = IMAS.fusion_power(dd.core_profiles.profiles_1d[]) .* 4 / 5
    total_power_neutrons = sum(nnt.wall_loading.power)
    total_neutrons_per_second = total_power_neutrons / 2.259E-12  # 14.1 MeV = 2.259E-12 J
    total_power_radiated = 0.0 # IMAS.radiative_power(dd.core_profiles.profiles_1d[])

    fwr = dd.neutronics.first_wall.r
    fwz = dd.neutronics.first_wall.z

    # all blanket modules
    blankets = [structure for structure in dd.build.structure if structure.type == Int(_blanket_)]
    resize!(dd.blanket.module, length(blankets))

    # Optimize midplane layers thicknesses that:
    # - achieve target TBR
    # - minimize leakage
    # - while fitting in current blanket thickness
    modules_Li6 = Float64[]
    blanket_model_1d = NNeutronics.Blanket()
    for (istructure, structure) in enumerate(blankets)
        bm = dd.blanket.module[istructure]
        bm.name = structure.name

        # get the relevant blanket layers, pre-optimization
        if sum(structure.outline.r) / length(structure.outline.r) < eqt.boundary.geometric_axis.r && length(blankets) > 1 # HFS
            fs = _hfs_
        else
            fs = _lfs_
        end
        d1 = IMAS.get_build(dd.build, type=_wall_, fs=fs, return_only_one=false, raise_error_on_missing=false)
        if d1 !== missing
            # if there are multiple walls we choose the one closest to the plasma
            if fs == _hfs_
                d1 = d1[end]
            else
                d1 = d1[1]
            end
        end
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

        # optimize layers thickeness
        res = optimize_layers(blanket_model_1d, bm.layer[1].midplane_thickness, bm.layer[2].midplane_thickness, bm.layer[3].midplane_thickness, dd.target.tritium_breeding_ratio, par.verbose)
        bm.layer[1].midplane_thickness, bm.layer[2].midplane_thickness, bm.layer[3].midplane_thickness, Li6, TBR, total_leakage = res
        push!(modules_Li6, Li6)
        for (kl, dl) in enumerate([d1, d2, d3])
            if dl !== missing
                dl.name = bm.layer[kl].name
                dl.thickness = bm.layer[kl].midplane_thickness
                dl.material = bm.layer[kl].material
            end
        end

    end

    # rebuild geometry
    ActorCXbuild(dd, act)

    # Evaluate TBR in 2D geometry
    modules_effective_thickness = []
    modules_wall_loading_power = []
    for (istructure, structure) in enumerate(blankets)
        bm = dd.blanket.module[istructure]

        # get the relevant optimized blanket layers
        if sum(structure.outline.r) / length(structure.outline.r) < eqt.boundary.geometric_axis.r && length(blankets) > 1 # HFS
            fs = _hfs_
        else
            fs = _lfs_
        end
        d1 = IMAS.get_build(dd.build, type=_wall_, fs=fs, return_only_one=false, raise_error_on_missing=false)
        if d1 !== missing
            # if there are multiple walls we choose the one closest to the plasma
            if fs == _hfs_
                d1 = d1[end]
            else
                d1 = d1[1]
            end
        end
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

        # calculate effective thickness of the layers based on the direction of the average impinging neutron flux
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
        bmt.power_thermal_neutrons = bmt.power_incident_neutrons * par.blanket_multiplier
        bmt.power_thermal_radiated = bmt.power_incident_radiated
        bmt.power_thermal_extracted = par.thermal_power_extraction_efficiency * (bmt.power_thermal_neutrons + bmt.power_thermal_radiated)

        push!(modules_effective_thickness, effective_thickness)
        push!(modules_wall_loading_power, nnt.wall_loading.power[index])
    end

    # Evaluate TBR in realistic geometry and optimize for single Li6 enrichment
    function target_TBR2D(blanket_model::NNeutronics.Blanket, Li6::Real, dd::IMAS.dd, modules_effective_thickness::Vector{<:Any}, modules_wall_loading_power::Vector{<:Any}, total_power_neutrons::Real, target::Float64=0.0)
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
        if target === 0.0
            return total_tritium_breeding_ratio
        else
            return (total_tritium_breeding_ratio - target)^2
        end
    end

    res = Optim.optimize(Li6 -> target_TBR2D(blanket_model_1d, Li6, dd, modules_effective_thickness, modules_wall_loading_power, total_power_neutrons, dd.target.tritium_breeding_ratio), 0.0, 100.0, Optim.GoldenSection(), rel_tol=1E-6)
    total_tritium_breeding_ratio = target_TBR2D(blanket_model_1d, res.minimizer, dd, modules_effective_thickness, modules_wall_loading_power, total_power_neutrons)

    @ddtime(dd.blanket.tritium_breeding_ratio = total_tritium_breeding_ratio)

    return actor
end

function optimize_layers(blanket_model::NNeutronics.Blanket, l1::Real, l2::Real, l3::Real, target_Li6::Real, verbose::Bool)
    energy_grid = NNeutronics.energy_grid()
    modules_thickness = [l1, l2, l3]
    blanket_thickness = sum(modules_thickness)

    model = JuMP.Model(NLopt.Optimizer)
    empty!(model)
    JuMP.set_optimizer_attribute(model, "algorithm", :LN_COBYLA)

    # leakage neutron model
    total_neutron_leakage(l1, l2, l3, Li6) = integrate(energy_grid, NNeutronics.leakeage_energy(blanket_model, l1, l2, l3, Li6, energy_grid))
    JuMP.register(model, :total_neutron_leakage, 4, total_neutron_leakage, autodiff=true)

    # TBR model
    JuMP_TBR(l1, l2, l3, Li6) = NNeutronics.TBR(blanket_model, l1, l2, l3, Li6)
    JuMP.register(model, :JuMP_TBR, 4, JuMP_TBR, autodiff=true)

    # variables and constraints
    if l1 == 0.0
        JuMP.@variable(model, d1 == 0.0)
    else
        JuMP.@variable(model, d1 >= 0.0)
    end
    if l2 == 0.0
        JuMP.@variable(model, d2 == 0.0)
    else
        JuMP.@variable(model, d2 >= 0.0)
    end
    if l3 == 0.0
        JuMP.@variable(model, d3 == 0.0)
    else
        JuMP.@variable(model, d3 >= 0.0)
    end
    JuMP.@variable(model, 0.0 <= Li6 <= 100.0)

    JuMP.@NLconstraint(model, blanket_thickness == d1 + d2 + d3)
    JuMP.@NLconstraint(model, JuMP_TBR(d1, d2, d3, Li6) >= target_Li6)
    JuMP.@NLobjective(model, Min, total_neutron_leakage(d1, d2, d3, Li6))
    JuMP.optimize!(model)

    l1, l2, l3, Li6 = (JuMP.value(d1), JuMP.value(d2), JuMP.value(d3), JuMP.value(Li6))
    total_leakage = total_neutron_leakage(l1, l2, l3, Li6)
    TBR = JuMP_TBR(l1, l2, l3, Li6)

    if verbose
        println(model)
        println(JuMP.solution_summary(model))
        println("l1=$l1  l2=$l2  l3=$l3  Li6=$Li6")
        println("TBR=$TBR")
        println("total_leakage=$total_leakage")
    end

    return l1, l2, l3, Li6, TBR, total_leakage
end