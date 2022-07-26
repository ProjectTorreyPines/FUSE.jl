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

function ParametersActor(::Type{Val{:ActorTauenn}})
    par = ParametersActor(nothing)
    par.error = Entry(Real, "", "Target convergence error"; default=1E-2)
    par.eped_factor = Entry(Real, "", "Scaling parameter for EPED-NN prediction"; default=1.0)
    par.rho_fluxmatch = Entry(Real, "", "Radial location where flux-macthing is done"; default=0.6)
    par.T_shaping = Entry(Real, "", "Shaping coefficient for the temperature profile"; default=1.8)
    par.temp_pedestal_ratio = Entry(Real, "", "Ion to electron temperature ratio in the pedestal"; default=1.0)
    par.transport_model = Switch([:tglfnn, :tglf, :h98y2, :ds03], "", "Transport model"; default=:tglfnn)
    par.warn_nn_train_bounds = Entry(Bool, "", "Warn if EPED-NN / TGLF-NN training bounds are exceeded"; default=false)
    par.confinement_factor = Entry(Real, "", "Confinement multiplier"; default=1.0)
    par.do_plot = Entry(Bool, "", "plot"; default=false)
    par.verbose = Entry(Bool, "", "verbose"; default=false)
    return par
end

"""
    ActorTauenn(dd::IMAS.dd, act::ParametersAllActors; kw...)

This actor estimates the core-transport using Tauenn, which evolves the kinetic profiles according to heat and particle flux matching.

The pedestal in this actor is evolved using EPED-NN.

!!! note 
    Stores data in `dd.core_profiles`
"""
function ActorTauenn(dd::IMAS.dd, act::ParametersAllActors; kw...)
    par = act.ActorTauenn(kw...)
    actor = ActorTauenn(dd, par)
    if par.do_plot
        ps = plot(dd.core_sources; color=:gray)
        pp = plot(dd.core_profiles; color=:gray, label="")
    end
    step(actor)
    finalize(actor)
    if par.do_plot
        display(plot!(ps, dd.core_sources))
        display(plot!(pp, dd.core_profiles))
    end
    if par.verbose
        display(actor.tauenn_parameters)
    end
    return actor
end

function ActorTauenn(dd::IMAS.dd, par::ParametersActor; kw...)
    par = par(kw...)
    tauenn_parameters = TAUENN.TauennParameters(;
        par.error,
        par.eped_factor,
        par.rho_fluxmatch,
        par.T_shaping,
        par.temp_pedestal_ratio,
        par.transport_model,
        par.confinement_factor,
        par.warn_nn_train_bounds)
    return ActorTauenn(dd, par, tauenn_parameters, TAUENN.TauennOutputs())
end

function step(actor::ActorTauenn)
    actor.tauenn_outputs = TAUENN.tau_enn(actor.dd, actor.tauenn_parameters; actor.par.verbose)
    return actor
end