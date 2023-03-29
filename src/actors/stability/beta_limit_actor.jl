#= ============== =#
#  ActorBetaLimit  #
#= ============== =#

Base.@kwdef mutable struct FUSEparameters__ActorBetaLimit{T} <: ParametersActor where {T<:Real}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    model::Switch{Symbol} = Switch(Symbol, [:Li, :Basic, :None], "-", "Model for the limit calculation"; default=:None)
end

mutable struct ActorBetaLimit <: PlasmaAbstractActor
    dd::IMAS.dd
    par::FUSEparameters__ActorBetaLimit
    passed::Union{Nothing, Bool}
    value::Union{Nothing, Real}
end

"""
    ActorBetaLimit(dd::IMAS.dd, act::ParametersAllActors; kw...)


"""
function ActorBetaLimit(dd::IMAS.dd, act::ParametersAllActors; kw...)
    par = act.ActorBetaLimit(kw...)
    actor = ActorBetaLimit(dd, par)
    step(actor)
    finalize(actor)
    return actor
end

function ActorBetaLimit(dd::IMAS.dd, par::FUSEparameters__ActorBetaLimit; kw...)
    logging_actor_init(ActorBetaLimit)
    par = par(kw...)
    ActorBetaLimit(dd, par, nothing, nothing)
end

"""
    step(actor::ActorBetaLimit)

Runs ActorBetaLimit to evaluate the beta limit for the given equilibrium
"""
function _step(actor::ActorBetaLimit)
    dd = actor.dd
    eqt = dd.equilibrium.time_slice[]
    par = actor.par

    if par.model == :None 
        logging(Logging.Error, :actors, "ActorBetaLimit: limit check disabled")
    elseif par.model == :Basic
        #actor.val = model1!eqt)
        val = model1(eqt)
    elseif par.model == :Li
        #actor.val = model2(eqt)
    else
        error("ActorBetaLimit: model = $(par.model) is unknown")
    end

    actor.passed = check_pass(val)
    return actor
end

function check_pass(value)
    if value > 0
        return true
    else
        return false
    end
end

function model1(eqt)
    beta_normal = eqt.global_quantities.beta_normal
    target_value = 3.5
    actual_value = beta_normal
    value = actual_value - target_value
    return value
end #use par.value, or return value?

function model2(eqt)
    beta_normal = eqt.global_quantities.beta_normal
    plasma_inductance =  eqt.global_quantities.li_3
    target_value = 4.0
    actual_value = beta_normal / plasma_inductance
    value = actual_value - target_value
    return value
end


# function model_list()
#     print("here")
#     model_list = Vector([
#         :None,
#         :Basic,
#         :Li
#     ])
#     return model_list
# end
