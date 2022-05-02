import TAUENN

#= ================ =#
#     TAUENN actor   #
#= ================ =#

mutable struct ActorTauenn <: ActorAbstract
    dd::IMAS.dd
    tauenn_parameters::TAUENN.TauennParameters
    tauenn_outputs::TAUENN.TauennOutputs
end

function ParametersActor(::Type{Val{:ActorTauenn}})
    par = ParametersActor(nothing)
    par.error = Entry(Real, "", "Error level"; default=1E-2)
    par.eped_factor = Entry(Real, "", "Scaling parameter for EPED-NN prediction"; default=1.0)
    par.rho_fluxmatch = Entry(Real, "", "radial location where flux-macthing is done"; default=0.6)
    par.T_shaping = Entry(Real, "", "shaping coefficient for the temperature profile"; default=1.8)
    par.temp_pedestal_ratio = Entry(Real, "", "ion to electron temperature ratio in the pedestal"; default=1.0)
    par.transport_model = Switch([:tglfnn, :h98y2, :tglf], "", "transport model"; default=:tglfnn)
    par.warn_nn_train_bounds = Entry(Bool, "", "warn if EPED-NN / TGLF-NN training bounds are exceeded"; default=false)
    par.do_plot = Entry(Bool, "", "plot"; default=false)
    par.verbose = Entry(Bool, "", "verbose"; default=false)
    return par
end

function ActorTauenn(dd::IMAS.dd, act::ParametersActor; kw...)
    par = act.ActorTauenn(kw...)
    if par.do_plot
        plot(dd.core_profiles; color=:gray, label="")
    end
    actor = ActorTauenn(dd;
        error=par.error,
        eped_factor=par.eped_factor,
        rho_fluxmatch=par.rho_fluxmatch,
        T_shaping=par.T_shaping,
        temp_pedestal_ratio=par.temp_pedestal_ratio,
        transport_model=par.transport_model,
        warn_nn_train_bounds=par.warn_nn_train_bounds)
    step(actor; verbose=par.verbose)
    finalize(actor)
    if par.do_plot
        display(plot!(dd.core_profiles))
    end
    if par.verbose
        display(actor.tauenn_parameters)
    end
    return dd
end

function ActorTauenn(dd::IMAS.dd; kw...)
    tauenn_parameters = TAUENN.TauennParameters()
    for key in keys(kw)
        setfield!(tauenn_parameters, key, kw[key])
    end
    return ActorTauenn(dd, tauenn_parameters, TAUENN.TauennOutputs())
end

function step(actor::ActorTauenn; verbose=false)
    actor.tauenn_outputs = TAUENN.tau_enn(actor.dd, actor.tauenn_parameters; verbose)
    return actor
end