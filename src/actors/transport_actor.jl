import TAUENN

#= ================ =#
#     TAUENN actor   #
#= ================ =#

mutable struct TauennActor <: AbstractActor
    dd::IMAS.dd
    tauenn_parameters::TAUENN.TauennParameters
    tauenn_outputs::TAUENN.TauennOutputs
end

function TauennActor(dd::IMAS.dd, par::Parameters; do_plot = false, verbose = false, kw...)
    if do_plot
        plot(dd.core_profiles; color = :gray, label = "")
    end
    actor = TauennActor(dd; kw...)
    step(actor; verbose = verbose)
    finalize(actor)
    if do_plot
        display(plot!(dd.core_profiles))
    end
    return dd
end

function TauennActor(dd::IMAS.dd; rho_fluxmatch = 0.6, eped_factor = 1.0, T_shaping = 1.8, temp_pedestal_ratio = 1.0, error = 1E-2, transport_model = :tglfnn, kw...)
    tauenn_parameters = TAUENN.TauennParameters()
    tauenn_parameters.eped_factor = eped_factor
    tauenn_parameters.rho_fluxmatch = rho_fluxmatch
    tauenn_parameters.T_shaping = T_shaping
    tauenn_parameters.temp_pedestal_ratio = temp_pedestal_ratio
    tauenn_parameters.transport_model = transport_model
    tauenn_parameters.error = error
    for key in keys(kw)
        setfield!(tauenn_parameters, key, kw[key])
    end
    return TauennActor(dd, tauenn_parameters, TAUENN.TauennOutputs())
end

function step(actor::TauennActor; verbose = false)
    actor.tauenn_outputs = TAUENN.tau_enn(actor.dd, actor.tauenn_parameters; verbose)
    return actor
end