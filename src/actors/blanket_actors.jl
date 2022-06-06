#= ============ =#
#  ActorBlanket  #
#= ============ =#

mutable struct ActorBlanket <: AbstractActor
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
    actor = ActorBlanket(dd; par.blanket_multiplier, par.thermal_power_extraction_efficiency)
    step(actor)
    finalize(actor)
    return actor
end

function step(actor::ActorBlanket)
    dd = actor.dd

    total_power_neutron = IMAS.fusion_power(dd.core_profiles.profiles_1d[]) .* 4.0
    total_power_radiative = 0.0 # IMAS.radiative_power(dd.core_profiles.profiles_1d[])

    tritium_breeding_ratio = 0.0

    for bm in dd.blanket.module
        neutron_capture_fraction = 1.0 / lenght(dd.blanket.module)
        radiative_capture_fraction = 1.0 / lenght(dd.blanket.module)

        bm.time_slice[].power_incident_neutron = total_power_neutron .* neutron_capture_fraction
        bm.time_slice[].power_incident_radiative = total_power_radiative .* radiative_capture_fraction

        bm.time_slice[].power_thermal_neutron = blanket_multiplier .* bm.time_slice[].power_incident_neutron
        bm.time_slice[].power_thermal_radiative = power_incident_radiative

        bm.time_slice[].power_thermal_extracted = thermal_power_extraction_efficiency * (power_thermal_neutron + power_thermal_radiative)

        bm.time_slice[].tritium_breeding_ratio = 1.0 # some function
        tritium_breeding_ratio += bm.time_slice[].tritium_breeding_ratio * bm.time_slice[].power_incident_neutron
    end

    @ddtime(bm.tritium_breeding_ratio = tritium_breeding_ratio / total_power_neutron)

end

function blanket_geometry(dd)
    resize(dd.blanket.module, 1)
    dd.blanket.module[1].name = "all"
end