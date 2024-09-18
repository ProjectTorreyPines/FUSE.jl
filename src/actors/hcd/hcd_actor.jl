#= ======== =#
#  ActorHCD  #
#= ======== =#
Base.@kwdef mutable struct FUSEparameters__ActorHCD{T<:Real} <: ParametersActor{T}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    _time::Float64 = NaN
    ec_model::Switch{Symbol} = Switch{Symbol}([:ECsimple, :none], "-", "EC source actor to run"; default=:ECsimple)
    ic_model::Switch{Symbol} = Switch{Symbol}([:ICsimple, :none], "-", "IC source actor to run"; default=:ICsimple)
    lh_model::Switch{Symbol} = Switch{Symbol}([:LHsimple, :none], "-", "LH source actor to run"; default=:LHsimple)
    nb_model::Switch{Symbol} = Switch{Symbol}([:NBsimple, :RABBIT, :none], "-", "NB source actor to run"; default=:NBsimple)
    pellet_model::Switch{Symbol} = Switch{Symbol}([:Pelletsimple, :none], "-", "Pellet source actor to run"; default=:Pelletsimple)
    power_scaling_cost_function::Entry{Function} = Entry{Function}("-", "EC, IC, LH, NB power optimization cost function, takes dd as input. Eg. dd -> (1.0 - IMAS.tau_e_thermal(dd) / IMAS.tau_e_h98(dd))")
end

mutable struct ActorHCD{D,P} <: CompoundAbstractActor{D,P}
    dd::IMAS.dd{D}
    par::FUSEparameters__ActorHCD{P}
    ec_actor::Union{Missing,ActorSimpleEC{D,P}}
    ic_actor::Union{Missing,ActorSimpleIC{D,P}}
    lh_actor::Union{Missing,ActorSimpleLH{D,P}}
    nb_actor::Union{Missing,ActorSimpleNB{D,P},ActorRABBIT{D,P}}
    pellet_actor::Union{Missing,ActorSimplePellet{D,P}}
end

"""
    ActorHCD(dd::IMAS.dd, act::ParametersAllActors; kw...)

Provides a common interface to run HCD actors
"""
function ActorHCD(dd::IMAS.dd, act::ParametersAllActors; kw...)
    actor = ActorHCD(dd, act.ActorHCD, act; kw...)
    step(actor)
    finalize(actor)
    return actor
end

function ActorHCD(dd::IMAS.dd, par::FUSEparameters__ActorHCD, act::ParametersAllActors; kw...)
    logging_actor_init(ActorHCD)
    par = par(kw...)
    if par.ec_model == :ECsimple
        ec_actor = ActorSimpleEC(dd, act.ActorSimpleEC)
    else
        ec_actor = missing
    end
    if par.ic_model == :ICsimple
        ic_actor = ActorSimpleIC(dd, act.ActorSimpleIC)
    else
        ic_actor = missing
    end
    if par.lh_model == :LHsimple
        lh_actor = ActorSimpleLH(dd, act.ActorSimpleLH)
    else
        lh_actor = missing
    end
    if par.nb_model == :NBsimple
        nb_actor = ActorSimpleNB(dd, act.ActorSimpleNB)
    elseif par.nb_model == :RABBIT
        nb_actor = ActorRABBIT(dd, act.ActorRABBIT)
    else
        nb_actor = missing
    end
    if par.pellet_model == :Pelletsimple
        pellet_actor = ActorSimplePellet(dd, act.ActorSimplePellet)
    else
        pellet_actor = missing
    end
    return ActorHCD(dd, par, ec_actor, ic_actor, lh_actor, nb_actor, pellet_actor)
end

"""
    _step(actor::ActorHCD)

Runs through the selected HCD actor's step
"""
function _step(actor::ActorHCD)
    dd = actor.dd
    par = actor.par

    # Call IMAS.sources!(dd) since the most would expect sources to be consistent when coming out of this actor
    IMAS.sources!(dd)

    if ismissing(par, :power_scaling_cost_function)
        _step_all_hcd(actor)

    else
        ps0 = deepcopy(dd.pulse_schedule)
        function scale_power_tau_cost(scale; dd, ps0, power_scaling_cost_function)
            scale_powers(dd.pulse_schedule, ps0, scale)
            _finalize(_step_all_hcd(actor))
            return power_scaling_cost_function(dd)^2
        end
        old_logging = actor_logging(dd, false)
        try
            res = Optim.optimize(scale -> scale_power_tau_cost(scale; dd, ps0, par.power_scaling_cost_function), 0.01, 100, Optim.GoldenSection(); rel_tol=1E-3)
            actor_logging(dd, old_logging)
            scale_power_tau_cost(res.minimizer; dd, ps0, par.power_scaling_cost_function)
        catch e
            actor_logging(dd, old_logging)
            rethrow(e)
        end
    end

    return actor
end

function scale_powers(ps, ps0, scale)
    for (beam, beam0) in zip(ps.ec.beam, ps0.ec.beam)
        beam.power_launched.reference .= beam0.power_launched.reference .* scale
    end
    for (antenna, antenna0) in zip(ps.ic.antenna, ps0.ic.antenna)
        antenna.power.reference .= antenna0.power.reference .* scale
    end
    for (antenna, antenna0) in zip(ps.lh.antenna, ps0.lh.antenna)
        antenna.power.reference .= antenna0.power.reference .* scale
    end
    for (unit, unit0) in zip(ps.nbi.unit, ps0.nbi.unit)
        unit.power.reference .= unit0.power.reference .* scale
    end
end

function _step_all_hcd(actor::ActorHCD)
    if actor.ec_actor !== missing
        step(actor.ec_actor)
    end
    if actor.ic_actor !== missing
        step(actor.ic_actor)
    end
    if actor.lh_actor !== missing
        step(actor.lh_actor)
    end
    if actor.nb_actor !== missing
        step(actor.nb_actor)
    end
    if actor.pellet_actor !== missing
        step(actor.pellet_actor)
    end
    return actor
end

"""
    _finalize(actor::ActorHCD)

Finalizes the selected CHD actor's finalize
"""
function _finalize(actor::ActorHCD)
    if actor.ec_actor !== missing
        finalize(actor.ec_actor)
    end
    if actor.ic_actor !== missing
        finalize(actor.ic_actor)
    end
    if actor.lh_actor !== missing
        finalize(actor.lh_actor)
    end
    if actor.nb_actor !== missing
        finalize(actor.nb_actor)
    end
    if actor.pellet_actor !== missing
        finalize(actor.pellet_actor)
    end
    return actor
end