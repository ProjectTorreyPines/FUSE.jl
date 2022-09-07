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
    par.blanket_multiplier = Entry(Real, "", "Neutron thermal power multiplier in blanket"; default=1.2)
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
    dd = actor.dd

    eqt = dd.equilibrium.time_slice[]
    nnt = dd.neutronics.time_slice[]

    total_power_neutrons = IMAS.fusion_power(dd.core_profiles.profiles_1d[]) .* 4 / 5
    total_power_neutrons = sum(nnt.wall_loading.power)

    total_power_radiated = 0.0 # IMAS.radiative_power(dd.core_profiles.profiles_1d[])

    blankets = [structure for structure in dd.build.structure if structure.type == Int(_blanket_)]
    resize!(dd.blanket.module, length(blankets))

    for (k, structure) in enumerate(blankets)
        bm = dd.blanket.module[k]
        bm.name = structure.name

        # evaluate neutron_capture_fraction
        tmp = convex_hull(vcat(eqt.boundary.outline.r, structure.outline.r), vcat(eqt.boundary.outline.z, structure.outline.z); closed_polygon=true)
        index = findall(x -> x == 1, [IMAS.PolygonOps.inpolygon((r, z), tmp) for (r, z) in zip(dd.neutronics.first_wall.r, dd.neutronics.first_wall.z)])
        neutron_capture_fraction = sum(nnt.wall_loading.power[index]) / total_power_neutrons

        # evaluate radiative_capture_fraction (to be fixed)
        radiative_capture_fraction = 1.0 / length(dd.blanket.module)

        resize!(bm.time_slice)
        bmt = bm.time_slice[]

        bmt.power_incident_neutrons = total_power_neutrons .* neutron_capture_fraction
        bmt.power_incident_radiated = total_power_radiated .* radiative_capture_fraction

        bmt.power_thermal_neutrons = bmt.power_incident_neutrons * actor.blanket_multiplier
        bmt.power_thermal_radiated = bmt.power_incident_radiated

        bmt.power_thermal_extracted = actor.thermal_power_extraction_efficiency * (bmt.power_thermal_neutrons + bmt.power_thermal_radiated)

        # blanket layer structure (designed to handle missing wall and/or missing shield)
        resize!(bm.layer, 3)
        if sum(structure.outline.r) / length(structure.outline.r) < eqt.boundary.geometric_axis.r # HFS
            d1 = IMAS.get_build(dd.build, type=_wall_, fs=_hfs_, raise_error_on_missing=false)
            d2 = IMAS.get_build(dd.build, type=_blanket_, fs=_hfs_)
            d3 = IMAS.get_build(dd.build, type=_shield_, fs=_hfs_, return_only_one=false, raise_error_on_missing=false)
            if d3 !== missing
                d3 = d3[end]
            end
        else  # LFS
            d1 = IMAS.get_build(dd.build, type=_wall_, fs=_lfs_)
            d2 = IMAS.get_build(dd.build, type=_blanket_, fs=_lfs_)
            d3 = IMAS.get_build(dd.build, type=_shield_, fs=_lfs_, return_only_one=false, raise_error_on_missing=false)
            if d3 !== missing
                d3 = d3[1]
            end
        end
        for (kl, dl) in enumerate(bm.layer)
            dl.name = "dummy layer $kl"
            dl.thickness = 0.0
            dl.material = "vacuum"
        end
        for (kl, dl) in enumerate([d1, d2, d3])
            if dl !== missing
                for field in [:name, :thickness, :material]
                    setproperty!(bm.layer[kl], field, getproperty(dl, field))
                end
            end
        end
    end

    # Optimize Li6/Li7 ratio to obtain target TBR
    function target_TBR(blanket_model, Li6, dd, target=nothing)
        total_tritium_breeding_ratio = 0.0
        for bm in dd.blanket.module
            bmt = bm.time_slice[]
            bm.layer[2].material = @sprintf("lithium-lead: Li6/7=%3.3f", Li6)
            bmt.tritium_breeding_ratio = NNeutronics.TBR(blanket_model, [dl.thickness for dl in bm.layer]..., Li6)
            total_tritium_breeding_ratio += bmt.tritium_breeding_ratio * bmt.power_incident_neutrons / total_power_neutrons
        end
        if target === nothing
            return total_tritium_breeding_ratio
        else
            return (total_tritium_breeding_ratio - target)^2
        end
    end

    blanket_model_1d = NNeutronics.Blanket()
    res = Optim.optimize(Li6 -> target_TBR(blanket_model_1d, Li6, dd, dd.target.tritium_breeding_ratio), 0.0, 100.0, Optim.GoldenSection(), rel_tol=1E-6)
    total_tritium_breeding_ratio = target_TBR(blanket_model_1d, res.minimizer, dd)

    @ddtime(dd.blanket.tritium_breeding_ratio = total_tritium_breeding_ratio)

    return actor
end
