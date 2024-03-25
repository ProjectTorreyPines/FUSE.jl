import RABBIT
using Plots

#= =========== =#
#  ActorRABBIT  #
#= =========== =#
Base.@kwdef mutable struct FUSEparameters__ActorRABBIT{T<:Real} <: ParametersActorPlasma{T}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    _time::Float64 = NaN
end

mutable struct ActorRABBIT{D,P} <: SingleAbstractActor{D,P}
    dd::IMAS.dd{D}
    par::FUSEparameters__ActorRABBIT{P}
    outputs::Union{RABBIT.RABBIToutputs,Vector{<:RABBIT.RABBIToutputs}}
end

function ActorRABBIT(dd::IMAS.dd, par::FUSEparameters__ActorRABBIT; kw...)
    par = par(kw...)
    return ActorRABBIT(dd, par, RABBIT.RABBIToutputs[])
end

"""
    ActorRABBIT(dd::IMAS.dd, act::ParametersAllActors; kw...)

"""
function ActorRABBIT(dd::IMAS.dd, act::ParametersAllActors; kw...)
    actor = ActorRABBIT(dd, act.ActorRABBIT; kw...)
    step(actor)
    finalize(actor)
    return actor
end

function _step(actor::ActorRABBIT)
    dd = actor.dd
    
    all_inputs = RABBIT.FUSEtoRABBITinput(dd)

    powe_data, powi_data, rho_data, time_data = RABBIT.run_RABBIT(all_inputs; remove_inputs=true)
    output = RABBIT.RABBIToutputs()

    output.powe_data = powe_data
    output.powi_data = powi_data
    output.rho_data = rho_data 
    output.time_data = time_data 

    actor.outputs = output

    return actor

end

function _finalize(actor::ActorRABBIT)

    powe = reshape(actor.outputs.powe_data, length(actor.outputs.rho_data), length(actor.outputs.time_data))

    p = plot(actor.outputs.rho_data, powe[:,1])

    for i in 2:length(actor.outputs.time_data)
        plot!(p, actor.outputs.rho_data, powe[:,i])
    end
    display(p)

    xlabel!("rho")
    ylims!(0,1e5)
    display(ylabel!("Power density profile to electrons - W/m^3 "))
        
    return actor

end