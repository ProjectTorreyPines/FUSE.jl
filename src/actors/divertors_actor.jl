#= ============== =#
#  ActorDivertors  #
#= ============== =#
Base.@kwdef mutable struct FUSEparameters__DivertorParameters <: ParametersActorSet 
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    inner = Entry(FUSEparameters__DivertorLegParameters, "-","inner divertor leg parameters")
    outer = Entry(FUSEparameters__DivertorLegParameters, "-","outer  divertor leg parameters")
    dome = Entry(FUSEparameters__DivertorDomeParameters, "-","divertor dome parameters")
end

Base.@kwdef mutable struct FUSEparameters__DivertorDomeParameters <: ParametersActorSet 
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
   
end

Base.@kwdef mutable struct FUSEparameters__DivertorLegParameters <: ParametersActorSet 
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_sett
    target = Entry(FUSEparameters__DivertorTargetParameters, "-","targe divertor parameters")
    baffle_cfr = Entry(FUSEparameters__DivertorBaffleParameters, "-","cfr divertor baffle parameters")
    baffle_pfr = Entry(FUSEparameters__DivertorBaffleParameters, "-","pfr divertor baffle parameters")
end

Base.@kwdef mutable struct FUSEparameters__DivertorBaffleParameters <: ParametersActorSet 
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_sett
end

Base.@kwdef mutable struct FUSEparameters__DivertorTargetParameters <: ParametersActorSet 
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_sett
end
#alternative dict based approach
DicDivertorParameters = Dict(:lower => Dict(:a=>1),:upper => Dict(:b=>2))

Base.@kwdef mutable struct FUSEparameters__ActorDivertors{T} <: ParametersActor 
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    lower = Entry(FUSEparameters__DivertorParameters, "-","lower divertor parameters")
    upper = Entry(FUSEparameters__DivertorParameters, "-","upper divertor parameters")
    thermal_power_extraction_efficiency = Entry(Real, "-",
        "Fraction of thermal power that is carried out by the coolant at the divertor interface, rather than being lost in the surrounding strutures.";
        default=1.0)
end

mutable struct ActorDivertors <: ReactorAbstractActor
    dd::IMAS.dd
    par::FUSEparameters__ActorDivertors
end

"""
    ActorSimpleDivertors(dd::IMAS.dd, act::ParametersAllActors; kw...)

Evaluates divertor loading and deposited power

!!! note 
    Stores data in `dd.divertors`
"""
function ActorDivertors(dd::IMAS.dd, act::ParametersAllActors; kw...)
    par = act.ActorDivertors(kw...)
    actor = ActorDivertors(dd, par)
    step(actor)
    finalize(actor)
    return actor
end

function ActorDivertors(dd::IMAS.dd, par::FUSEparameters__ActorDivertors; kw...)
    logging_actor_init(ActorDivertors)
    par = par(kw...)
    return ActorDivertors(dd, par)
end

function _step(actor::ActorDivertors)
    dd = actor.dd

    empty!(dd.divertors)

    divertors = [structure for structure in dd.build.structure if structure.type == Int(_divertor_)]
    resize!(dd.divertors.divertor, length(divertors))

    total_power_incident = IMAS.total_power_source(IMAS.total_sources(dd.core_sources, dd.core_profiles.profiles_1d[]))
    dd.divertors.time = div_time = dd.core_sources.time
    for (k, structure) in enumerate(divertors)

        div = dd.divertors.divertor[k]
        div.name = structure.name

        div.power_incident.time = dd.divertors.time
        div.power_incident.data = ones(length(div_time)) .* total_power_incident ./ length(dd.divertors.divertor)

        div.power_neutrals.time = div_time
        div.power_neutrals.data = zeros(length(div_time))

        div.power_radiated.time = div_time
        div.power_radiated.data = divertor_radiated(dd)

        div.power_recombination_neutrals.time = div_time
        div.power_recombination_neutrals.data = divertor_recombination_neutrals(dd)

        div.power_thermal_extracted.time = div_time
        div.power_thermal_extracted.data = actor.thermal_power_extraction_efficiency .* (
            div.power_incident.data .+
            div.power_radiated.data .+
            div.power_recombination_neutrals.data)
    end
end

"""
    divertor_radiated(dd::IMAS.dd)

dummy divertor_radiated now assumed at zero"""
function divertor_radiated(dd)
    return zeros(length(dd.divertors.time))
end

"""
    divertor_reflection(dd::IMAS.dd)

dummy divertor_reflection now assumed at zero"""
function divertor_reflection(dd)
    return zeros(length(dd.divertors.time))
end

"""
    divertor_recombination_neutrals(dd::IMAS.dd)

dummy divertor_recombination_neutrals power calculation now assumed at zero"""
function divertor_recombination_neutrals(dd)
    return zeros(length(dd.divertors.time))
end
