#= ============ =#
#  ActorBlanket  #
#= ============ =#

import NNeutronics

mutable struct ActorBlanket <: ReactorAbstractActor
    dd::IMAS.dd
    par::ParametersActor
    blanket_multiplier::Real
    thermal_power_extraction_efficiency::Real
end

function ParametersActor(::Type{Val{:ActorBlanket}})
    par = ParametersActor(nothing)
    par.blanket_multiplier = Entry(Real, "", "Neutron thermal power multiplier in blanket"; default = 1.2)
    par.optimize_parameter = Switch([:TBR, :leakage], "", "Blanket parameter to be optimized"; default = :TBR)
    par.optimize_variable = Switch([:d1, :d2, :d3, :Li6], "", "Blanket variable to be adjusted"; default = :Li6)
    par.thermal_power_extraction_efficiency = Entry(
        Real,
        "",
        "Fraction of thermal power that is carried out by the coolant at the blanket interface, rather than being lost in the surrounding strutures.";
        default=1.0
    )
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
    step(actor)
    finalize(actor)
    return actor
end

function ActorBlanket(dd::IMAS.dd, par::ParametersActor; kw...)
    logging_actor_init(ActorBlanket)
    par = par(kw...)
    return ActorBlanket(dd, par, par.blanket_multiplier, par.thermal_power_extraction_efficiency)
end

function _step(actor::ActorBlanket)
    @warn "currently only using GAMBL materials in blanket and first wall."
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

    blankets = [structure for structure in dd.build.structure if structure.type == Int(_blanket_)]
    resize!(dd.blanket.module, length(blankets))

    fwr = dd.neutronics.first_wall.r
    fwz = dd.neutronics.first_wall.z

    modules_effective_thickness = []
    modules_wall_loading_power = []

    for (istructure, structure) in enumerate(blankets)
        bm = dd.blanket.module[istructure]
        bm.name = structure.name

        # get the relevant blanket layers
        if sum(structure.outline.r) / length(structure.outline.r) < eqt.boundary.geometric_axis.r # HFS
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

    function target_leakage(blanket_model, Li6, dd, modules_effective_thickness, modules_wall_loading_power, total_power_neutrons, target=nothing)
        total_leakage = 0.0
        for (ibm, bm) in enumerate(dd.blanket.module)
            # bmt = bm.time_slice[]
            blank_mod_time_leakage = 0.0 # bmt.neutron_leakage = 0.0, need to add to data structure
            module_wall_loading_power = sum(modules_wall_loading_power[ibm])
            for k in 1:length(modules_wall_loading_power[ibm])
                blank_mod_time_leakage += (sum(NNeutronics.leakeage_energy(blanket_model, modules_effective_thickness[ibm][k, :]..., Li6)) * modules_wall_loading_power[ibm][k] / module_wall_loading_power)
            end
            total_leakage += blank_mod_time_leakage * module_wall_loading_power / total_power_neutrons
        end
        if target === nothing
            return total_leakage
        else
            return (total_leakage - target)^2
        end
    end

    blanket_model_1d = NNeutronics.Blanket()
    res = Optim.optimize(Li6 -> target_TBR(blanket_model_1d, Li6, dd, modules_effective_thickness, modules_wall_loading_power, total_power_neutrons, dd.target.tritium_breeding_ratio), 0.0, 100.0, Optim.GoldenSection(), rel_tol=1E-6)
    total_tritium_breeding_ratio = target_TBR(blanket_model_1d, res.minimizer, dd, modules_effective_thickness, modules_wall_loading_power, total_power_neutrons)
    Li6 = res.minimizer
    @ddtime(dd.blanket.tritium_breeding_ratio = total_tritium_breeding_ratio)
    total_leakage = target_leakage(blanket_model_1d, Li6, dd, modules_effective_thickness, modules_wall_loading_power, total_power_neutrons)
    println(total_leakage)
    return actor
end
