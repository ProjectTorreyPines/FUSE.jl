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
    outputs::Union{RABBIT.RABBIToutput,Vector{<:RABBIT.RABBIToutput}}
end

function ActorRABBIT(dd::IMAS.dd, par::FUSEparameters__ActorRABBIT; kw...)
    par = par(kw...)
    return ActorRABBIT(dd, par, RABBIT.RABBIToutput[])
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

    powe_data, powi_data, jnbcd_data, bdep_data, torqdepo_data, rho_data, time_data = RABBIT.run_RABBIT(all_inputs; remove_inputs=true)
    output = RABBIT.RABBIToutput()

    output.powe_data = powe_data
    output.powi_data = powi_data
    output.jnbcd_data = jnbcd_data
    output.bdep_data = bdep_data
    output.torque_data = torqdepo_data
    output.rho_data = rho_data
    output.time_data = time_data

    actor.outputs = output

    return actor

end

function _finalize(actor::ActorRABBIT)
    dd = actor.dd
    cs = dd.core_sources
    output = actor.outputs

    num_t = length(output.time_data)
    source = resize!(cs.source, :nbi, "identifier.name" => "nbi"; wipe=true)

    if num_t > 2
        prof_1d = resize!(source.profiles_1d, num_t)
    else
        prof_1d = resize!(source.profiles_1d, 1)
    end

    source.profiles_1d[1].grid.rho_tor_norm = output.rho_data
    ion = resize!(source.profiles_1d[1].ion, 1)[1]

    source.profiles_1d[1].total_ion_energy = output.powi_data[:, 1]
    source.profiles_1d[1].electrons.energy = output.powe_data[:, 1]
    source.profiles_1d[1].electrons.particles = vec(sum(output.bdep_data[:, 1, :]; dims=2))
    source.profiles_1d[1].j_parallel = output.jnbcd_data[:, 1]
    source.profiles_1d[1].momentum_tor = vec(sum(output.torque_data[:, 1, :]; dims=2))

    if num_t > 2
        for i in 2:num_t
            source.profiles_1d[i].grid.rho_tor_norm = output.rho_data
            source.profiles_1d[i].total_ion_energy = output.powi_data[:, i]
            source.profiles_1d[i].electrons.energy = output.powe_data[:, i]
            source.profiles_1d[i].electrons.particles = vec(sum(output.bdep_data[:, i, :]; dims=2))
            source.profiles_1d[i].momentum_tor = vec(sum(output.torque_data[:, i, :]; dims=2))
        end
    end

    for nbu in dd.nbi.unit
        IMAS.ion_element!(ion, "H$(Int(floor(nbu.species.a)))"; fast=true)
    end

    return actor

end