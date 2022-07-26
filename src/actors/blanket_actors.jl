#= ============ =#
#  ActorBlanket  #
#= ============ =#

Base.@kwdef mutable struct ActorBlanket <: ReactorAbstractActor
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
        default=1.0,
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
    par = par(kw...)
    return ActorBlanket(dd, par, par.blanket_multiplier, par.thermal_power_extraction_efficiency)
end

function step(actor::ActorBlanket)
    dd = actor.dd

    empty!(dd.blanket)
    eqt = dd.equilibrium.time_slice[]
    nnt = dd.neutronics.time_slice[]

    total_power_neutrons = IMAS.fusion_power(dd.core_profiles.profiles_1d[]) .* 4 / 5
    total_power_neutrons = sum(nnt.wall_loading.power)

    total_power_radiated = 0.0 # IMAS.radiative_power(dd.core_profiles.profiles_1d[])

    blankets = [structure for structure in dd.build.structure if structure.type == Int(_blanket_)]

    tritium_breeding_ratio = 0.0
    resize!(dd.blanket.module, length(blankets))
    for (k, structure) in enumerate(blankets)
        bm = dd.blanket.module[k]
        bm.name = structure.name

        # evaluate neutron_capture_fraction
        tmp = convex_hull(collect(zip(vcat(eqt.boundary.outline.r, structure.outline.r), vcat(eqt.boundary.outline.z, structure.outline.z))); closed_polygon=true)
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

        bmt.tritium_breeding_ratio = 1.0 # some function
        tritium_breeding_ratio += bmt.tritium_breeding_ratio * bmt.power_incident_neutrons
    end

    @ddtime(dd.blanket.tritium_breeding_ratio = tritium_breeding_ratio / total_power_neutrons)
end
