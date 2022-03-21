import TAUENN

#= ================ =#
#     TAUENN actor   #
#= ================ =#

mutable struct TaueNNactor <: AbstractActor
    dd::IMAS.dd
    parameters::TAUENN.TauennParameters
    outputs::TAUENN.TauennOutputs
end

function TaueNNactor(dd::IMAS.dd; rho_fluxmatch = 0.6, eped_factor = 1.0, T_shaping = 1.8, temp_pedestal_ratio = 1.0, error = 1E-2, transport_model = :tglfnn, kw...)
    parameters = TAUENN.TauennParameters()
    parameters.eped_factor = eped_factor
    parameters.rho_fluxmatch = rho_fluxmatch
    parameters.T_shaping = T_shaping
    parameters.temp_pedestal_ratio = temp_pedestal_ratio
    parameters.transport_model = transport_model
    parameters.error = error
    for param in keys(kw)
        setfield!(parameters, param, kw[param])
    end
    return TaueNNactor(dd, parameters, TAUENN.TauennOutputs())
end

function step(actor::TaueNNactor; verbose = false)
    actor.outputs = TAUENN.tau_enn(actor.dd, actor.parameters; verbose)
    return actor
end