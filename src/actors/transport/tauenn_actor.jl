import TAUENN

#= ================ =#
#     TAUENN actor   #
#= ================ =#
Base.@kwdef mutable struct FUSEparameters__ActorTauenn{T} <: ParametersActor where {T<:Real}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    error::Entry{T} = Entry{T}("-", "Target convergence error"; default=1E-2)
    eped_factor::Entry{T} = Entry{T}("-", "Scaling parameter for EPED-NN prediction"; default=1.0)
    rho_fluxmatch::Entry{T} = Entry{T}("-", "Radial location where flux-macthing is done"; default=0.6)
    T_shaping::Entry{T} = Entry{T}("-", "Shaping coefficient for the temperature profile"; default=1.8)
    T_ratio_pedestal::Entry{T} = Entry{T}("-", "Ion to electron temperature ratio in the pedestal"; default=1.0)
    transport_model::Switch{Symbol} = Switch{Symbol}([:tglfnn, :tglf, :h98y2, :ds03], "-", "Transport model"; default=:tglfnn)
    warn_nn_train_bounds::Entry{Bool} = Entry{Bool}("-", "Warn if EPED-NN / TGLF-NN training bounds are exceeded"; default=false)
    eped_only_powerlaw::Entry{Bool} = Entry{Bool}("-", "EPED-NN uses power-law pedestal fit (without NN correction)"; default=false)
    update_pedestal::Entry{Bool} = Entry{Bool}("-","update pedestal with eped_nn inside TAUENN" ;default=true)
    confinement_factor::Entry{T} = Entry{T}("-", "Confinement multiplier"; default=1.0)
    do_plot::Entry{Bool} = Entry{Bool}("-", "Plot"; default=false)
    verbose::Entry{Bool} = Entry{Bool}("-", "Verbose"; default=false)
end

mutable struct ActorTauenn{D,P} <: PlasmaAbstractActor
    dd::IMAS.dd{D}
    par::FUSEparameters__ActorTauenn{P}
    tauenn_parameters::TAUENN.TauennParameters
    tauenn_outputs::TAUENN.TauennOutputs
end

"""
    ActorTauenn(dd::IMAS.dd, act::ParametersAllActors; kw...)

Estimates the core-transport using TAUENN, which evolves the kinetic profiles according to heat and particle flux matching.

The pedestal is evolved using the EPED-NN model

!!! note 
    Stores data in `dd.core_profiles`
"""
function ActorTauenn(dd::IMAS.dd, act::ParametersAllActors; kw...)
    actor = ActorTauenn(dd, act.ActorTauenn; kw...)
    step(actor)
    finalize(actor)
    return actor
end

function ActorTauenn(dd::IMAS.dd, par::FUSEparameters__ActorTauenn; kw...)
    logging_actor_init(ActorTauenn)
    par = par(kw...)

    tauenn_parameters = TAUENN.TauennParameters(;
        par.error,
        par.eped_factor,
        par.rho_fluxmatch,
        par.T_shaping,
        par.T_ratio_pedestal,
        par.transport_model,
        par.confinement_factor,
        par.warn_nn_train_bounds,
        par.update_pedestal,
        par.eped_only_powerlaw)

    return ActorTauenn(dd, par, tauenn_parameters, TAUENN.TauennOutputs())
end

function _step(actor::ActorTauenn)
    dd = actor.dd
    par = actor.par
    if par.do_plot
        ps = plot(dd.core_sources; color=:gray)
        pp = plot(dd.core_profiles; color=:gray, label="")
    end
    actor.tauenn_outputs = TAUENN.tau_enn(dd, actor.tauenn_parameters; par.verbose)
    if par.do_plot
        display(plot!(ps, dd.core_sources))
        display(plot!(pp, dd.core_profiles))
    end
    if par.verbose
        display(actor.tauenn_parameters)
    end
    return actor
end