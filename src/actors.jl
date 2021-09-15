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

# initialize SolovevEquilibriumActor from an equilibrium IDS
function SolovevEquilibriumActor(equilibrium::IMAS.equilibrium, time_index::Integer)
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
    alpha = -0.155
    S = solovev(B0, R0, ϵ, δ, κ, alpha, qstar, B0_dir=B0_dir, Ip_dir=Ip_dir)
    return SolovevEquilibriumActor(S)
end

function init(eq::IMAS.equilibrium, time=0.0;
                B0, R0, ϵ, δ, κ, qstar,
                diverted::Bool=false,
                xpoint::Union{NTuple{2},Nothing}=(diverted ? (R0 * (1 - 1.1 * δ * ϵ), -R0 * 1.1 * κ * ϵ) : nothing),
                symmetric::Bool=(xpoint === nothing))
    time_index = get_time_index(eq.time_slice, time)
    eqt = eq.time_slice[time_index]
    eqt.profiles_1d.psi = [1.0]
    eqt.profiles_1d.r_outboard = [R0 + ϵ * R0]
    eqt.profiles_1d.r_inboard = [R0 - ϵ * R0]
    eqt.profiles_1d.elongation = [κ]
    eqt.profiles_1d.triangularity_lower = [δ]
    eqt.profiles_1d.triangularity_upper = [δ]
    set_field_time_array(eq.vacuum_toroidal_field, :b0, time_index, B0)
    eq.vacuum_toroidal_field.r0 = R0
    eqt.profiles_1d.q = [qstar]
    return eq
end

