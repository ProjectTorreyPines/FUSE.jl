#= ================= =#
#  ActorPlasmaLimits  #
#= ================= =#
Base.@kwdef mutable struct FUSEparameters__ActorPlasmaLimits{T<:Real} <: ParametersActor{T}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    _time::Float64 = NaN
    models::Entry{Vector{Symbol}} =
        Entry{Vector{Symbol}}("-", "Models used for checking plasma operational limits: $(supported_limit_models)"; default=deepcopy(default_limit_models))
    raise_on_breach::Entry{Bool} = Entry{Bool}("-", "Raise an error when one or more operational limits are breached"; default=true)
    #== display and debugging parameters ==#
    verbose::Entry{Bool} = act_common_parameters(; verbose=false)
end

mutable struct ActorPlasmaLimits{D,P} <: CompoundAbstractActor{D,P}
    dd::IMAS.dd{D}
    par::FUSEparameters__ActorPlasmaLimits{P}
    act::ParametersAllActors{P}
    function ActorPlasmaLimits(dd::IMAS.dd{D}, par::FUSEparameters__ActorPlasmaLimits{P}, act::ParametersAllActors{P}; kw...) where {D<:Real,P<:Real}
        logging_actor_init(ActorPlasmaLimits)
        par = par(kw...)
        return new{D,P}(dd, par, act)
    end
end

"""
    ActorPlasmaLimits(dd::IMAS.dd, act::ParametersAllActors; kw...)

Runs all the limit actors
"""
function ActorPlasmaLimits(dd::IMAS.dd, act::ParametersAllActors; kw...)
    actor = ActorPlasmaLimits(dd, act.ActorPlasmaLimits, act; kw...)
    step(actor)
    finalize(actor)
    return actor
end

"""
    _step(actor::ActorPlasmaLimits)

Runs through the selected plasma limits
"""
function _step(actor::ActorPlasmaLimits)
    dd = actor.dd
    par = actor.par
    act = actor.act

    # run all limit models
    for model_name in par.models
        func = try
            eval(model_name)
        catch e
            if typeof(e) <: UndefVarError
                error("Unknown model `$(repr(model_name))`. Supported models are $(supported_limit_models)")
            else
                rethrow(e)
            end
        end
        func(dd, act)
    end

    # identify which limits were breached
    failed_list = String[]
    passed_list = String[]
    for model in dd.limits.model
        txt = "$(Int(ceil(@ddtime(model.fraction)*100)))% of: $(model.identifier.description)"
        if Bool(@ddtime model.cleared)
            push!(passed_list, txt)
        else
            push!(failed_list, txt)
        end
    end

    # generate report text
    if !isempty(failed_list)
        failed = "\n\nSome stability rules exceed their limit threshold:\n $(join(failed_list, "\n "))"
    else
        failed = ""
    end
    if !isempty(passed_list)
        passed = "\n\nSome stability rules satisfy their limit threshold:\n $(join(passed_list, "\n "))"
    else
        passed = ""
    end
    txt = "act.ActorPlasmaLimits.models = $(par.models)$failed$passed"

    # error, warn, or print
    if !isempty(failed)
        if par.raise_on_breach
            error(txt)
        else
            @warn(txt)
        end
    elseif par.verbose
        println(txt)
    end

    return actor
end
