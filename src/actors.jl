#= =================== =#
#  init IMAS structures #
#= =================== =#

function init(ids::IMAS.IDS, time::Real)
    error("Function init() not defined for ids of type $(typeof(ids))")
end

function init(equilibrium::IMAS.equilibrium, time::Real=0.0;
        B0, R0, ϵ, δ, κ, beta_t, qstar,
        diverted::Bool=false,
        xpoint::Union{NTuple{2},Nothing}=(diverted ? (R0 * (1 - 1.1 * δ * ϵ), -R0 * 1.1 * κ * ϵ) : nothing),
        symmetric::Bool=(xpoint === nothing))
    time_index = get_time_index(equilibrium.time_slice, time)
    eqt = equilibrium.time_slice[time_index]
    eqt.profiles_1d.psi = [1.0]
    eqt.profiles_1d.r_outboard = [R0 + ϵ * R0]
    eqt.profiles_1d.r_inboard = [R0 - ϵ * R0]
    eqt.profiles_1d.elongation = [κ]
    eqt.profiles_1d.triangularity_lower = [δ]
    eqt.profiles_1d.triangularity_upper = [δ]
    set_field_time_array(equilibrium.vacuum_toroidal_field, :b0, time_index, B0)
    equilibrium.vacuum_toroidal_field.r0 = R0
    eqt.profiles_1d.q = [qstar]
    eqt.global_quantities.beta_tor = beta_t
    return equilibrium
end


#= ============ =#
#  AbstractActor #
#= ============ =#
abstract type AbstractActor end

"""
    Take a step with a given actor
"""
function Base.step(actor::AbstractActor)
    error("Function step() not defined for actor of type $(typeof(actor))")
end

"""
    store output data in IDS
"""
function finalize(actor::AbstractActor)
    error("Function finalize() not defined for actor of type $(typeof(actor))")
end

#= =========== =#
#  Equilibrium  #
#= =========== =#
using Equilibrium

abstract type EquilibriumActor <: AbstractActor end

mutable struct SolovevEquilibriumActor <: EquilibriumActor
    eq_in::IMAS.equilibrium
    time::Real
    S::SolovevEquilibrium
    eq_out::IMAS.equilibrium
end

function SolovevEquilibriumActor(equilibrium::IMAS.equilibrium, time::Real)
    time_index = get_time_index(equilibrium.time_slice, time)
    eqt = equilibrium.time_slice[time_index]
    a = eqt.profiles_1d.r_outboard[end] - eqt.profiles_1d.r_inboard[end]
    R0 = eqt.profiles_1d.r_outboard[end] + eqt.profiles_1d.r_inboard[end]
    κ = eqt.profiles_1d.elongation[end]
    δ = (eqt.profiles_1d.triangularity_lower[end] + eqt.profiles_1d.triangularity_upper[end]) / 2.0
    ϵ = a / R0
    B0 = abs(equilibrium.vacuum_toroidal_field.b0[time_index])
    B0_dir = Int(sign(B0))
    R0 = equilibrium.vacuum_toroidal_field.r0
    qstar = eqt.profiles_1d.q[end]
    Ip_dir = Int(sign(qstar) * B0_dir)
    alpha = 0.0
    S0 = solovev(B0, R0, ϵ, δ, κ, alpha, qstar, B0_dir=B0_dir, Ip_dir=Ip_dir)
    SolovevEquilibriumActor(equilibrium, time, S0, IMAS.equilibrium())
end

function Base.step(actor::SolovevEquilibriumActor; abs_error=1E-3, verbose=false)
    # non-linear optimization translating `alpha` into `beta_t`
    S0 = actor.S
    time_index = get_time_index(actor.eq_in.time_slice, actor.time)
    target_beta = actor.eq_in.time_slice[time_index].global_quantities.beta_tor

    function opti(x)
        S = solovev(S0.B0, S0.R0, S0.epsilon, S0.delta, S0.kappa, x[1], S0.qstar, B0_dir=S0.sigma_B0, Ip_dir=S0.sigma_Ip)
        v[1] = S
        precision = abs(S.beta_t - target_beta)
        cost = precision.^2
        return cost
    end
    
    x = [S0.alpha]
    v = Any[S0]
    for k in 1:1000
        g = ForwardDiff.gradient(opti, x, )
        x[1] -= (g[1] / abs(S0.beta_t))
        precision = abs(v[1].beta_t.value - target_beta)
        if verbose
            @printf("α=%3.3f β_t=%3.3e precision=%3.3e\n", x[1], v[1].beta_t.value, precision)
        end
        if precision < abs_error
            break
        end
    end
    
    actor.S = solovev(S0.B0, S0.R0, S0.epsilon, S0.delta, S0.kappa, x[1], S0.qstar, B0_dir=S0.sigma_B0, Ip_dir=S0.sigma_Ip)
end

function finalize(actor::SolovevEquilibriumActor)
    equilibrium = actor.eq_out
    time_index = get_time_index(equilibrium.time_slice, actor.time)
    eqt = equilibrium.time_slice[time_index]

    # eqt.profiles_1d.r_outboard =
    # eqt.profiles_1d.r_inboard = 
    # eqt.profiles_1d.elongation = 
    # eqt.profiles_1d.triangularity_upper = 
    # eqt.profiles_1d.triangularity_lower = 
    # equilibrium.vacuum_toroidal_field.b0 = 
    # equilibrium.vacuum_toroidal_field.r0 = 
    # ...
    return actor.eq_out
end