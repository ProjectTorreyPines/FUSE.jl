#= ============ =#
#  AbstractActor #
#= ============ =#

abstract type AbstractActor end

function step(actor::AbstractActor, ids::IDS, time_index::Integer)
    error("Function step() not defined for actor of type $(typeof(actor))")
end

function finalize(actor::AbstractActor, ids::IDS, time_index::Integer)
    error("Function finalize() not defined for actor of type $(typeof(actor))")
end

#= =========== =#
#  Equilibrium  #
#= =========== =#

using Equilibrium

abstract type EquilibriumActor <: AbstractActor end

struct SolovevEquilibriumActor <: EquilibriumActor
    S::SolovevEquilibrium
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

function SolovevEquilibriumActor(equilibrium::IMAS.equilibrium, time::Real; verbose=false)
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

    # non-linear optimization translating `alpha` into `beta_t`
    target_beta = equilibrium.time_slice[time_index].global_quantities.beta_tor
    alpha = -0.155 # initial guess
    S0 = solovev(B0, R0, ϵ, δ, κ, alpha, qstar, B0_dir=B0_dir, Ip_dir=Ip_dir)
    
    function opti(x)
        S = solovev(B0, R0, ϵ, δ, κ, x[1], qstar, B0_dir=B0_dir, Ip_dir=Ip_dir)
        cost = abs(S.beta_t - target_beta).^2
        if verbose
            @printf("α=%3.3f β_t=%3.3e cost=%3.3e\n",x[1].value,S.beta_t,cost.value)
        end
        if cost < (1E-3)^2
            throw(StopIteration())
        end
        return cost
    end
    
    x = [alpha]
    try
        for k in 1:1000
            g = ForwardDiff.gradient(opti, x)
            x[1] -= (g[1] / abs(S0.beta_t))
        end
    catch e
        if e isa StopIteration
            # ignore
        else
            rethrow()
        end
    end
    
    S1 = solovev(B0, R0, ϵ, δ, κ, x[1], qstar, B0_dir=B0_dir, Ip_dir=Ip_dir)
    
    # here sill need to do the work to return the data in IMAS format
    return SolovevEquilibriumActor(S1)
end