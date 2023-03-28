#= ================= =#
#  ActorPlasmaLimits  #
#= ================= =#
Base.@kwdef mutable struct FUSEparameters__ActorPlasmaLimits{T} <: ParametersActor where {T<:Real}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    error_on_exceeded_limits::Entry{Bool} = Entry(Bool, "-", "Error if limits are exceeded"; default=true)
    limit_beta_model::Switch{Symbol} = Switch(Symbol, [:None, :Basic, :Li], "-", "Beta limit model to run"; default=:None)
    limit_beta_value::Entry{T} = Entry(T, "-", "Value for the beta limit"; default=0.0)
    limit_q95_model::Switch{Symbol} = Switch(Symbol, [:None, :Basic], "-", "Current limit model to run"; default=:None)
    limit_q95_value::Entry{T} = Entry(T, "-", "Value for the q95/current limit"; default=0.0)
    limit_density_model::Switch{Symbol} = Switch(Symbol, [:None, :Basic], "-", "Current limit model to run"; default=:None)
    limit_density_value::Entry{T} = Entry(T, "-", "Value for the qedge/current limit"; default=0.0)
    min_q95::Entry{T} = Entry(T, "-", "Minimum q95 (0.0 disables it)"; default=2.0)
    vertical_stability::Entry{Bool} = Entry(Bool, "-", "Vertical stability"; default=true)
    greenwald_fraction::Entry{T} = Entry(T, "-", "Greenwald fraction (0.0 disables it)"; default=0.0)
    beta_critical::Entry{T} = Entry(T, "-", "Critical normalized beta (0.0 disables it)"; default=0.0)
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
    par = act.ActorPlasmaLimits(kw...)
    actor = ActorPlasmaLimits(dd, par)
    step(actor)
    finalize(actor)
    return actor
end

function add_limit_check!(exceeded_limits, what, value, f, limit)
    if f(value, limit)
        push!(exceeded_limits, "$what of $value $f limit of $limit")
    end
end

function check_limit_beta!(exceeded_limits, dd, par)
    limit_name = "beta limit"

    if par.limit_beta_model == :None
        # Don't check limit

    elseif par.limit_beta_model == :Basic
        value = dd.equilibrium.time_slice[].global_quantities.beta_normal
        limit = par.limit_beta_value
        add_limit_check!(exceeded_limits, limit_name, value, >, limit)

    elseif par.limit_beta_model == :Li
        beta_normal = dd.equilibrium.time_slice[].global_quantities.beta_normal
        plasma_inductance =  dd.equilibrium.time_slice[].global_quantities.li_3
        value = beta_normal / plasma_inductance
        limit = par.limit_beta_value
        add_limit_check!(exceeded_limits, limit_name, value, >, limit)

    else
        error("$(par.limit_betan_model) is not implemented")

    end
end

function check_limit_q95!(exceeded_limits, dd, par)
    limit_name = "current limit"

    if par.limit_q95_model == :None
        # Don't check limit

    elseif par.limit_q95_model == :Basic
        value = abs(dd.equilibrium.time_slice[].global_quantities.q_95)
        limit = par.limit_q95_value
        add_limit_check!(exceeded_limits, "q95", value, <, limit)

    else
        error("$(par.limit_q95_model) is not implemented")

    end
end

function check_limit_density!(exceeded_limits, dd, par)
    limit_name = "density limit"

    eqt = dd.equilibrium.time_slice[]
    cp1d = dd.core_profiles.profiles_1d[]

    if par.limit_density_model == :None
        # Don't check limit

    elseif par.limit_density_model == :Basic
        value = IMAS.greenwald_fraction(eqt, cp1d)
        limit = par.limit_density_value
        add_limit_check!(exceeded_limits, limit_name, value, >, limit)

    else
        error("$(par.limit_density_model) is not implemented")

    end
end

"""
    step(actor::ActorPlasmaLimits)

Tests if limits are exceeded
"""
function _step(actor::ActorPlasmaLimits)
    dd = actor.dd
    par = actor.par

    exceeded_limits = String[]

    check_limit_beta!(exceeded_limits, dd, par)
    check_limit_q95!(exceeded_limits, dd, par)
    check_limit_density!(exceeded_limits, dd, par)

    if par.vertical_stability
        eqt = dd.equilibrium.time_slice[]
        value = dd.equilibrium.time_slice[].boundary.elongation
        limit = IMAS.elongation_limit(eqt)
        add_limit_check!(exceeded_limits, "Elongation", value, >, limit)
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
