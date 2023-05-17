#= ============== =#
#  ActorDivertors  #
#= ============== =#
Base.@kwdef mutable struct FUSEparameters__ActorDivertors{T} <: ParametersActor where {T<:Real}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    thermal_power_extraction_efficiency::Entry{T} = Entry(T, "-",
        "Fraction of thermal power that is carried out by the coolant at the divertor interface, rather than being lost in the surrounding strutures.";
        default=1.0)
end

mutable struct ActorDivertors <: ReactorAbstractActor
    dd::IMAS.dd
    par::FUSEparameters__ActorDivertors
    thermal_power_extraction_efficiency::Real
end

"""
    ActorSimpleDivertors(dd::IMAS.dd, act::ParametersAllActors; kw...)

Evaluates divertor loading and deposited power

!!! note 
    Stores data in `dd.divertors`
"""
function ActorDivertors(dd::IMAS.dd, act::ParametersAllActors; kw...)
    actor = ActorDivertors(dd, act.ActorDivertors; kw...)
    step(actor)
    finalize(actor)
    return actor
end

function ActorDivertors(dd::IMAS.dd, par::FUSEparameters__ActorDivertors; kw...)
    logging_actor_init(ActorDivertors)
    par = par(kw...)
    return ActorDivertors(dd, par, par.thermal_power_extraction_efficiency)
end

function _step(actor::ActorDivertors)
    dd = actor.dd

    empty!(dd.divertors)

    divertors = [structure for structure in dd.build.structure if structure.type == Int(_divertor_)]
    resize!(dd.divertors.divertor, length(divertors))

    total_power_incident = IMAS.total_power_source(IMAS.total_sources(dd.core_sources, dd.core_profiles.profiles_1d[]))
    dd.divertors.time = dd.core_sources.time
    for (k, structure) in enumerate(divertors)

        div = dd.divertors.divertor[k]
        div.name = structure.name

        power_incident = @ddtime(div.power_incident.data = total_power_incident ./ length(dd.divertors.divertor))

        power_neutrals = @ddtime(div.power_neutrals.data = 0.0)

        power_radiated = @ddtime(div.power_radiated.data = 0.0)

        power_recombination_neutrals = @ddtime(div.power_recombination_neutrals.data = 0.0)

        @ddtime(div.power_thermal_extracted.data = actor.thermal_power_extraction_efficiency * (
            power_incident +
            power_radiated +
            power_neutrals +
            power_recombination_neutrals))
    end

    return actor
end
