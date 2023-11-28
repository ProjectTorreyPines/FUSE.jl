using LinearAlgebra

#= ================== =#
#  ActorFluxCapacitor  #
#= ================== =#
Base.@kwdef mutable struct FUSEparameters__ActorFluxCapacitor{T} <: ParametersActor where {T<:Real}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
end

mutable struct ActorFluxCapacitor{D,P} <: PlasmaAbstractActor{D,P}
    dd::IMAS.dd{D}
    par::FUSEparameters__ActorFluxCapacitor{P}
    actor_fm::ActorFluxMatcher{D,P}
    q_history::Vector{Vector{D}}
    cp1d_history::Vector{IMAS.core_profiles__profiles_1d{D}}
end

"""
    ActorFluxCapacitor(dd::IMAS.dd, act::ParametersAllActors; kw...)

Optimizes q-profile to maximize fusion power, while keeping |q|>1
"""
function ActorFluxCapacitor(dd::IMAS.dd, act::ParametersAllActors; kw...)
    actor = ActorFluxCapacitor(dd, act.ActorFluxCapacitor, act; kw...)
    step(actor)
    finalize(actor)
    return actor
end

function ActorFluxCapacitor(dd::IMAS.dd{D}, par::FUSEparameters__ActorFluxCapacitor, act::ParametersAllActors; kw...) where {D<:Real}
    logging_actor_init(ActorFluxCapacitor)
    par = par(kw...)
    actor_fm = ActorFluxMatcher(dd, act.ActorFluxMatcher, act)
    cp1d_history = Vector{IMAS.core_profiles__profiles_1d{D}}()
    return ActorFluxCapacitor(dd, par, actor_fm, Vector{Vector{D}}(), cp1d_history)
end

"""
    _step(actor::ActorFluxCapacitor)

ActorFluxCapacitor step
"""
function _step(actor::ActorFluxCapacitor)
    dd = actor.dd
    par = actor.par

    eqt1d = dd.equilibrium.time_slice[].profiles_1d
    cp1d = dd.core_profiles.profiles_1d[]

    function cost(
        actor::ActorFluxCapacitor,
        eqt1d::IMAS.equilibrium__time_slice___profiles_1d,
        cp1d::IMAS.core_profiles__profiles_1d,
        zq::Vector{<:Real})

        # reset q profile
        eqt1d.q = deepcopy(actor.q_history[1])

        # reset cp1d
        # NOTE: Technically we should reset cp1d but we don't because so long as the FluxMatcher converges,
        # then the starting profiles going into FluxMatcher do not affect the inputs-outputs relationship that the optimizer sees.
        # Using the previous cp1d (ie. without reset) makes the FluxMatcher start from a closer solution, which improves convergence.
        # fill!(cp1d, deepcopy(actor.cp1d_history[1]))

        # update the q profile
        eq_gridpoints = [argmin(abs.(rho_x .- eqt1d.rho_tor_norm)) for rho_x in rho_transport]
        eq_rho_transport = eqt1d.rho_tor_norm[eq_gridpoints]
        eqt1d.q = IMAS.profile_from_z_transport(eqt1d.q, eqt1d.rho_tor_norm, eq_rho_transport, zq, 0.9)

        # run flux matcher
        finalize(step(actor.actor_fm))

        # save history
        push!(actor.q_history, deepcopy(eqt1d.q))
        push!(actor.cp1d_history, cp1d)

        # maximize fusion power
        Pfusion = IMAS.fusion_power(cp1d)
        error_fusion = actor.dd.requirements.power_electric_net / Pfusion

        # avoid dropping q below 1
        # cost function designed to be 1.0 at |q|=1 and 0.0 at x0
        x0 = 1.1
        α = 1.0
        qmin = minimum(abs.(eqt1d.q))
        cost_q_lt1 = (((exp.(x0 ./ qmin) ./ exp(1)) .- 1.0) ./ ((exp.(x0 ./ 1.0) ./ exp(1)) .- 1.0)) .^ α
        cost_q_lt1 *= (qmin < x0)

        total_cost = norm((error_fusion, cost_q_lt1))

        @show total_cost
        @show (error_fusion, Pfusion)
        @show (cost_q_lt1, qmin)

        return total_cost
    end

    # optimize q
    push!(actor.q_history, deepcopy(eqt1d.q))
    push!(actor.cp1d_history, deepcopy(cp1d))
    rho_transport = actor.actor_fm.par.rho_transport
    initial_zq = pack_q_profiles(eqt1d, rho_transport)
    res = Optim.optimize(x -> cost(actor, eqt1d, cp1d, x), initial_zq, Optim.NelderMead(), Optim.Options(; g_tol=1E-3))
    cost(actor, eqt1d, cp1d, res.minimizer)

    return actor
end

function pack_q_profiles(eqt1d::IMAS.equilibrium__time_slice___profiles_1d, rho_transport::AbstractVector{<:Real})
    eq_gridpoints = [argmin(abs.(rho_x .- eqt1d.rho_tor_norm)) for rho_x in rho_transport]
    return IMAS.calc_z(eqt1d.rho_tor_norm, eqt1d.q)[eq_gridpoints]
end