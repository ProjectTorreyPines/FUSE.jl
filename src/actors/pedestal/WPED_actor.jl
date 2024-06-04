#= ============= =#
#  ActorWPED  #
#= ============= =#
Base.@kwdef mutable struct FUSEparameters__ActorWPED{T<:Real} <: ParametersActorPlasma{T}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    _time::Float64 = NaN
    #== actor parameters ==#
    ped_to_core_fraction::Entry{T} = Entry{T}("-", "ratio of pedestal stored energy to core stored energy"; default=0.3)
end

mutable struct ActorWPED{D,P} <: SingleAbstractActor{D,P}
    dd::IMAS.dd{D}
    par::FUSEparameters__ActorWPED{P}
    result::Nothing
end

"""
    ActorWPED(dd::IMAS.dd, act::ParametersAllActors; kw...)

Finds the temperature profile at the edge to match the ped_to_core_fraction of stored energy set in par.ped_to_core_fraction
"""
function ActorWPED(dd::IMAS.dd, act::ParametersAllActors; kw...)
    actor = ActorWPED(dd, act.ActorWPED; kw...)
    step(actor)
    finalize(actor)
    return actor
end


function ActorWPED(dd::IMAS.dd, par::FUSEparameters__ActorWPED; kw...)
    logging_actor_init(ActorWPED)
    par = par(kw...)
    return ActorWPED(dd, par, nothing)
end

"""
    _step(actor::ActorWPED)

Finds the temperature profile at the edge to match the ped_to_core_fraction of stored energy set in par.ped_to_core_fraction
"""
function _step(actor::ActorWPED{D,P}) where {D<:Real,P<:Real}
    dd = actor.dd
    par = actor.par

    cp1d = dd.core_profiles.profiles_1d[]

    rho_idx09 = argmin(abs.(cp1d.grid.rho_tor_norm .- 0.9))
    # Keep the Ti/Te ratio constant
    Ti_over_Te = cp1d.ion[1].temperature[rho_idx09] / cp1d.electrons.temperature[rho_idx09]
    res = Optim.optimize(x -> cost_t08(x, cp1d, par.ped_to_core_fraction, Ti_over_Te), 50, 5e3, Optim.GoldenSection(); rel_tol=1E-3)
    cost_t08(res.minimizer, cp1d, par.ped_to_core_fraction, Ti_over_Te)

    return actor
end

"""
    _finalize(actor::ActorWPED)

Writes results to dd.summary.local.pedestal and possibly updates core_profiles
"""
function _finalize(actor::ActorWPED)
    dd = actor.dd
    par = actor.par

    return actor
end


function cost_t08(T08, cp1d, ratio_wanted, Ti_over_Te)

    rho_idx09 = argmin(abs.(cp1d.grid.rho_tor_norm .- 0.9))

    Te = cp1d.electrons.temperature
    Tlin = LinRange(T08, Te[end], length(Te) - rho_idx09 + 1)
    cp1d.electrons.temperature[rho_idx09+1:end] .= Tlin[2:end]
    cp1d.electrons.temperature[1:rho_idx09] = cp1d.electrons.temperature[1:rho_idx09] .- cp1d.electrons.temperature[rho_idx09] .+ Tlin[1]

    Tlin = LinRange(T08, Te[end], length(Te) - rho_idx09 + 1) .* Ti_over_Te
    for ion in cp1d.ion
        ion.temperature[rho_idx09+1:end] .= Tlin[2:end]
        ion.temperature[1:rho_idx09] = ion.temperature[1:rho_idx09] .- ion.temperature[rho_idx09] .+ Tlin[1]
    end

    core, edge = core_and_edge_energy(cp1d)
    ratio = edge / core
    return ((ratio .- ratio_wanted) / ratio_wanted)^2
end


function core_and_edge_energy(cp1d)
    p = cp1d.pressure_thermal
    rho_idx09 = argmin(abs.(cp1d.grid.rho_tor_norm .- 0.9))
    p_edge = deepcopy(p)
    p_edge[1:rho_idx09] .= p[rho_idx09]
    p_core = p .- p_edge
    core_value = 3 / 2 .* IMAS.cumul_integrate(cp1d.grid.volume, p_core)
    edge_value = 3 / 2 .* IMAS.cumul_integrate(cp1d.grid.volume, p_edge)

    return core_value[end], edge_value[end]
end
