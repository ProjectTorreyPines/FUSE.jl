using LinearAlgebra

#= ================== =#
#  ActorFluxCapacitor  #
#= ================== =#
Base.@kwdef mutable struct FUSEparameters__ActorFluxCapacitor{T<:Real} <: ParametersActor{T}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    _time::Float64 = NaN
    rho_q_optimization::Entry{AbstractVector{T}} = Entry{AbstractVector{T}}("-", "ρ grid for q profile optimization (taken by default from FluxMatcher)")
end

mutable struct ActorFluxCapacitor{D,P} <: CompoundAbstractActor{D,P}
    dd::IMAS.dd{D}
    par::FUSEparameters__ActorFluxCapacitor{P}
    fluxmatcher::ActorFluxMatcher{D,P}
    qed::ActorQED{D,P}
    q_history::Vector{Vector{D}}
    cp1d_history::Vector{IMAS.core_profiles__profiles_1d{D}}
end

"""
    ActorFluxCapacitor(dd::IMAS.dd, act::ParametersAllActors; kw...)

Optimizes q-profile to maximize fusion power, while keeping |q|>1 and Jtor physical
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
    fluxmatcher = ActorFluxMatcher(dd, act.ActorFluxMatcher, act)
    qed = ActorQED(dd, act.ActorQED; Nt=0)
    finalize(step(qed))
    cp1d_history = Vector{IMAS.core_profiles__profiles_1d{D}}()
    return ActorFluxCapacitor(dd, par, fluxmatcher, qed, Vector{Vector{D}}(), cp1d_history)
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

    rho_top = try
        rho_ped = @ddtime(dd.summary.local.pedestal.position.rho_tor_norm)
        rho_top = 1.0 - (1.0 - rho_ped) * 2.0
    catch
        0.9
    end

    function cost(
        actor::ActorFluxCapacitor,
        eqt1d::IMAS.equilibrium__time_slice___profiles_1d,
        cp1d::IMAS.core_profiles__profiles_1d,
        rho_q_optimization::Vector{<:Real},
        zq::Vector{<:Real})

        # reset q profile
        eqt1d.q = deepcopy(actor.q_history[1])

        # reset cp1d
        # NOTE: Technically we should reset cp1d but we don't because so long as the FluxMatcher converges,
        # then the starting profiles going into FluxMatcher do not affect the inputs-outputs relationship that the optimizer sees.
        # Using the previous cp1d (ie. without reset) makes the FluxMatcher start from a closer solution, which improves convergence.
        # merge!(cp1d, actor.cp1d_history[1])

        # update the q profile
        eq_gridpoints = [argmin(abs.(rho_x .- eqt1d.rho_tor_norm)) for rho_x in rho_q_optimization]
        eqt1d.q = IMAS.profile_from_z_transport(eqt1d.q, eqt1d.rho_tor_norm, eqt1d.rho_tor_norm[eq_gridpoints], zq, rho_top)

        # run flux matcher
        finalize(step(actor.fluxmatcher))

        # evaluate corresponding toroidal current
        ι = QED.FE(eqt1d.rho_tor_norm, 1.0 ./ eqt1d.q)
        jt_R_new = QED.Jt_R(actor.qed.QO; ι, ρ=eqt1d.rho_tor_norm)

        # save history
        push!(actor.q_history, deepcopy(eqt1d.q))
        push!(actor.cp1d_history, cp1d)

        # maximize fusion power
        Pfusion = IMAS.fusion_power(cp1d)
        cost_low_fusion = actor.dd.requirements.power_electric_net / Pfusion

        # avoid dropping q below 1
        # cost function designed to be 1.0 at |q|=1 and 0.0 at x0
        x0 = 1.1
        α = 1.0
        qmin = minimum(abs.(eqt1d.q))
        cost_q_lt1 = (((exp.(x0 ./ qmin) ./ exp(1)) .- 1.0) ./ ((exp.(x0 ./ 1.0) ./ exp(1)) .- 1.0)) .^ α
        cost_q_lt1 *= (qmin < x0)

        # avoid current density below value at rho_top
        index = argmin(abs.(max(rho_top, rho_q_optimization[end]) .- eqt1d.rho_tor_norm))
        jt_R_new .*= sign(jt_R_new[index])
        cost_j_low = sum(abs.(jt_R_new[1:index-1][jt_R_new[1:index-1].<jt_R_new[index]])) / jt_R_new[index] * 100.0
        min_jt = minimum(jt_R_new[1:index-1])

        total_cost = norm((cost_low_fusion, cost_q_lt1, cost_j_low))

        @show total_cost
        @show (cost_j_low, min_jt)
        @show (cost_low_fusion, Pfusion)
        @show (cost_q_lt1, qmin)

        return total_cost
    end

    if ismissing(par, :rho_q_optimization)
        rho_q_optimization = collect(actor.fluxmatcher.par.rho_transport)
    else
        rho_q_optimization = collect(par.rho_q_optimization)
    end

    # optimize q
    push!(actor.q_history, deepcopy(eqt1d.q))
    push!(actor.cp1d_history, deepcopy(cp1d))
    initial_zq = pack_q_profiles(eqt1d, rho_q_optimization)
    res = Optim.optimize(x -> cost(actor, eqt1d, cp1d, rho_q_optimization, x), initial_zq, Optim.NelderMead(), Optim.Options(; g_tol=1E-3))
    cost(actor, eqt1d, cp1d, rho_q_optimization, res.minimizer)

    return actor
end

function pack_q_profiles(eqt1d::IMAS.equilibrium__time_slice___profiles_1d, rho_q_optimization::AbstractVector{<:Real})
    eq_gridpoints = [argmin(abs.(rho_x .- eqt1d.rho_tor_norm)) for rho_x in rho_q_optimization]
    return IMAS.calc_z(eqt1d.rho_tor_norm, eqt1d.q)[eq_gridpoints]
end