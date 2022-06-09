#= ============ =#
#  ActorBlanket  #
#= ============ =#

Base.@kwdef mutable struct ActorBlanket <: ReactorAbstractActor
    dd::IMAS.dd
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
    ActorBlanket(dd::IMAS.dd, act::ParametersActor; kw...)

Blanket actor

!!! note 
    Stores data in `dd.blanket`
"""
function ActorBlanket(dd::IMAS.dd, act::ParametersActor; kw...)
    par = act.ActorBlanket(kw...)
    actor = ActorBlanket(dd, par.blanket_multiplier, par.thermal_power_extraction_efficiency)
    step(actor)
    finalize(actor)
    return actor
end

function step(actor::ActorBlanket)
    dd = actor.dd

    empty!(dd.blanket)

    total_power_neutrons = IMAS.fusion_power(dd.core_profiles.profiles_1d[]) .* 4 / 5
    total_power_radiated = 0.0 # IMAS.radiative_power(dd.core_profiles.profiles_1d[])
    tritium_breeding_ratio = 0.0

    blankets = blanket_regions!(dd.build, dd.equilibrium.time_slice[])

    resize!(dd.blanket.module, length(blankets))
    for (k, structure) in enumerate(blankets)
        bm = dd.blanket.module[k]
        bm.name = structure.name

        neutron_capture_fraction = 1.0 / length(dd.blanket.module)
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
