import TAUENN

#= ================ =#
#     TAUENN actor   #
#= ================ =#

mutable struct ActorTauenn <: PlasmaAbstractActor
    dd::IMAS.dd
    par::ParametersActor
    tauenn_parameters::TAUENN.TauennParameters
    tauenn_outputs::TAUENN.TauennOutputs
end

Base.@kwdef mutable struct FUSEparameters__ActorTauenn{T} <: ParametersActor where {T<:Real}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    error = Entry(Real, "", "Target convergence error"; default=1E-2)
    eped_factor = Entry(Real, "", "Scaling parameter for EPED-NN prediction"; default=1.0)
    rho_fluxmatch = Entry(Real, "", "Radial location where flux-macthing is done"; default=0.6)
    T_shaping = Entry(Real, "", "Shaping coefficient for the temperature profile"; default=1.8)
    temp_pedestal_ratio = Entry(Real, "", "Ion to electron temperature ratio in the pedestal"; default=1.0)
    transport_model = Switch(Symbol, [:tglfnn, :tglf, :h98y2, :ds03], "", "Transport model"; default=:tglfnn)
    warn_nn_train_bounds = Entry(Bool, "", "Warn if EPED-NN / TGLF-NN training bounds are exceeded"; default=false)
    update_pedestal = Entry(Bool, "","update pedestal with eped_nn inside TAUENN" ;default=true)
    confinement_factor = Entry(Real, "", "Confinement multiplier"; default=1.0)
    do_plot = Entry(Bool, "", "plot"; default=false)
    verbose = Entry(Bool, "", "verbose"; default=false)
end

"""
    ActorTauenn(dd::IMAS.dd, act::ParametersAllActors; kw...)

Estimates the core-transport using TAUENN, which evolves the kinetic profiles according to heat and particle flux matching.

The pedestal is evolved using the EPED-NN model

!!! note 
    Stores data in `dd.core_profiles`
"""
function ActorTauenn(dd::IMAS.dd, act::ParametersAllActors; kw...)
    par = act.ActorTauenn(kw...)
    actor = ActorTauenn(dd, par)
    step(actor)
    finalize(actor)
    return actor
end

function ActorTauenn(dd::IMAS.dd, par::ParametersActor; kw...)
    logging_actor_init(ActorTauenn)
    par = par(kw...)

    tauenn_parameters = TAUENN.TauennParameters(;
        par.error,
        par.eped_factor,
        par.rho_fluxmatch,
        par.T_shaping,
        par.temp_pedestal_ratio,
        par.transport_model,
        par.confinement_factor,
        par.warn_nn_train_bounds,
        par.update_pedestal)

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