#= ================= =#
#  ActorPlasmaLimits  #
#= ================= =#
Base.@kwdef mutable struct FUSEparameters__ActorPlasmaLimits{T} <: ParametersActor where {T<:Real}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    error_on_exceeded_limits::Entry{Bool} = Entry(Bool, "-", "Error if limits are exceeded"; default=true)
    min_q95::Entry{T} = Entry(T, "-", "Minimum q95 (0.0 disables it)"; default=2.0)
    vertical_stability::Entry{Bool} = Entry(Bool, "-", "Vertical stability"; default=true)
    greenwald_fraction::Entry{T} = Entry(T, "-", "Greenwald fraction (0.0 disables it)"; default=0.0)
end

mutable struct ActorPlasmaLimits <: PlasmaAbstractActor
    dd::IMAS.dd
    par::FUSEparameters__ActorPlasmaLimits
    function ActorPlasmaLimits(dd::IMAS.dd, par::FUSEparameters__ActorPlasmaLimits; kw...)
        logging_actor_init(ActorPlasmaLimits)
        par = par(kw...)
        new(dd, par)
    end
end

"""
    ActorPlasmaLimits(dd::IMAS.dd, act::ParametersAllActors; kw...)

Tests if plasma operational limits are exceeded
"""
function ActorPlasmaLimits(dd::IMAS.dd, act::ParametersAllActors; kw...)
    par = act.ActorPlasmaLimits
    actor = ActorPlasmaLimits(dd, par; kw...)
    step(actor)
    finalize(actor)
    return actor
end

function add_limit_check!(exceeded_limits, what, value, f, limit)
    if f(value, limit)
        push!(exceeded_limits, "$what of $value $f limit of $limit")
    end
end

"""
    step(actor::ActorPlasmaLimits)

Tests if limits are exceeded
"""
function _step(actor::ActorPlasmaLimits)
    dd = actor.dd
    par = actor.par
    eqt = dd.equilibrium.time_slice[]
    cp1d = dd.core_profiles.profiles_1d[]

    exceeded_limits = String[]

    if par.min_q95 > 0.0
        value = abs(dd.equilibrium.time_slice[].global_quantities.q_95)
        limit = par.min_q95
        add_limit_check!(exceeded_limits, "q95", value, <, limit)
    end

    if par.vertical_stability
        value = dd.equilibrium.time_slice[].boundary.elongation
        limit = IMAS.elongation_limit(eqt)
        add_limit_check!(exceeded_limits, "Elongation", value, >, limit)
    end

    if par.greenwald_fraction != 0.0
        value = IMAS.greenwald_fraction(eqt, cp1d)
        limit = par.greenwald_fraction
        add_limit_check!(exceeded_limits, "Greenwald fraction", value, >, limit)
    end

    if !isempty(exceeded_limits)
        message = "Exceeded limits:"
        for txt in exceeded_limits
            message *= "\n* " * txt
        end

        if par.error_on_exceeded_limits
            error(message)
        else
            logging(Logging.Error, :actors, message)
        end
    end

    return actor
end
