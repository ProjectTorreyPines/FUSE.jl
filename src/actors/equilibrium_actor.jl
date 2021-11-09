#= ==== =#
#  init  #
#= ==== =#

"""
    init(eqt::IMAS.equilibrium__time_slice; B0::Real, R0::Real, ϵ::Real, δ::Real, κ::Real, beta_n::Real, ip::Real, x_point::Union{Vector, NTuple{2}, Bool}=false)

Initialize equilibrium IDS based on some basic Miller geometry parameters
"""
function init(eqt::IMAS.equilibrium__time_slice;
              B0::Real, R0::Real, ϵ::Real, δ::Real, κ::Real, beta_n::Real, ip::Real,
              x_point::Union{Vector,NTuple{2},Bool}=false)
    eqt.boundary.minor_radius = ϵ * R0
    eqt.boundary.geometric_axis.r = R0
    eqt.boundary.elongation = κ
    eqt.boundary.triangularity = δ
    eqt.boundary.triangularity = δ
    eqt.profiles_1d.psi = [1.0]
    eqt.profiles_1d.f = [B0 * R0]
    eqt.global_quantities.ip = ip
    eqt.global_quantities.beta_normal = beta_n
    if x_point === true
        x_point = (R0 * (1 - 1.1 * δ * ϵ), -R0 * 1.1 * κ * ϵ)
    end
    if isa(x_point, Union{Vector,Tuple})
        resize!(eqt.boundary.x_point, 1)
        eqt.boundary.x_point[1].r = x_point[1]
        eqt.boundary.x_point[1].z = x_point[2]
    end
    return eqt
end

#= ======= =#
#  Solovev  #
#= ======= =#

using Equilibrium
using Printf
import ForwardDiff
import Optim

mutable struct SolovevEquilibriumActor <: EquilibriumActor
    eq_in::IMAS.equilibrium__time_slice
    S::SolovevEquilibrium
    eq_out::IMAS.equilibrium__time_slice
end

"""
    function SolovevEquilibriumActor(eq_in::IMAS.equilibrium__time_slice, qstar=1.5, alpha=0.0, symmetric=true)

Constructor for the SolovevEquilibriumActor structure
“One size fits all” analytic solutions to the Grad–Shafranov equation
Phys. Plasmas 17, 032502 (2010); https://doi.org/10.1063/1.3328818

- qstar: Kink safety factor

- alpha: Constant affecting the pressure
"""
function SolovevEquilibriumActor(eq_in::IMAS.equilibrium__time_slice; qstar=1.5, alpha=0.0,
                                 symmetric=true) # symmetric should really be passed/detected through IMAS

    a = eq_in.boundary.minor_radius
    R0 = eq_in.boundary.geometric_axis.r
    κ = eq_in.boundary.elongation
    δ = eq_in.boundary.triangularity
    ϵ = a / R0
    B0 = eq_in.profiles_1d.f[end] / R0
    B0_dir = Int(sign(B0))
    Ip_dir = Int(sign(qstar) * B0_dir)

    if length(eq_in.boundary.x_point) > 0
        xpoint = (eq_in.boundary.x_point[1].r, eq_in.boundary.x_point[1].z)
    else
        xpoint = nothing
    end

    S0 = solovev(B0, R0, ϵ, δ, κ, alpha, qstar, B0_dir=B0_dir, Ip_dir=Ip_dir, symmetric=symmetric, xpoint=xpoint)

    eq_out = IMAS.equilibrium__time_slice()
    if ! is_missing(eq_in, :time)
        eq_out.time = eq_in.time
    end

    SolovevEquilibriumActor(eq_in, S0, eq_out)
end

"""
    IMAS2Equilibrium(eqt::IMAS.equilibrium__time_slice)

Convert IMAS.equilibrium__time_slice to Equilibrium.jl EFIT structure
"""
function IMAS2Equilibrium(eqt::IMAS.equilibrium__time_slice)
    dim1 = range(eqt.profiles_2d[1].grid.dim1[1], eqt.profiles_2d[1].grid.dim1[end], length=length(eqt.profiles_2d[1].grid.dim1))
    @assert collect(dim1) ≈ eqt.profiles_2d[1].grid.dim1
    dim2 = range(eqt.profiles_2d[1].grid.dim2[1], eqt.profiles_2d[1].grid.dim2[end], length=length(eqt.profiles_2d[1].grid.dim2))
    @assert collect(dim2) ≈ eqt.profiles_2d[1].grid.dim2
    psi = range(eqt.profiles_1d.psi[1], eqt.profiles_1d.psi[end], length=length(eqt.profiles_1d.psi))
    @assert collect(psi) ≈ eqt.profiles_1d.psi

    Equilibrium.efit(   Equilibrium.cocos(11), # COCOS
                        dim1, # Radius/R range
                        dim2, # Elevation/Z range
                        psi, # Polodial Flux range (polodial flux from magnetic axis)
                        eqt.profiles_2d[1].psi, # Polodial Flux on RZ grid (polodial flux from magnetic axis)
                        eqt.profiles_1d.f, # Polodial Current
                        eqt.profiles_1d.pressure, # Plasma pressure
                        eqt.profiles_1d.q, # Q profile
                        eqt.profiles_1d.psi .* 0, # Electric Potential
                        (eqt.global_quantities.magnetic_axis.r, eqt.global_quantities.magnetic_axis.z), # Magnetic Axis (raxis,zaxis)
                        Int(sign(eqt.profiles_1d.f[end]) * sign(eqt.global_quantities.ip)) # sign(dot(J,B))
                    )
end

#= == =#
# STEP #
#= == =#
"""
    Base.step(actor::SolovevEquilibriumActor; verbose=false)

Non-linear optimization to obtain a target `ip` and `beta_normal`
"""
function Base.step(actor::SolovevEquilibriumActor; verbose=false)
    S0 = actor.S

    target_ip = actor.eq_in.global_quantities.ip
    target_beta = actor.eq_in.global_quantities.beta_normal

    B0, R0, epsilon, delta, kappa, alpha, qstar, target_ip, target_beta = promote(S0.B0, S0.R0, S0.epsilon, S0.delta, S0.kappa, S0.alpha, S0.qstar, target_ip, target_beta)

    function opti(x)
        S = solovev(B0, R0, epsilon, delta, kappa, x[1], x[2], B0_dir=S0.sigma_B0, Ip_dir=S0.sigma_Ip, symmetric=true, xpoint=nothing)
        beta_cost = abs((Equilibrium.beta_n(S) - target_beta)^2 / target_beta)
        ip_cost = abs((Equilibrium.plasma_current(S) - target_ip)^2 / target_ip)
        return (beta_cost + ip_cost)^2
    end

    # having issues taking derivatives of SolovevEquilibriumActor and using Newton optimizer when S is diverted 
    if isa(B0, ForwardDiff.Dual)# || S0.diverted
        optim_method = Optim.NelderMead()
    else
        optim_method = Optim.Newton()
    end
    res = Optim.optimize(opti, [alpha, qstar], optim_method, Optim.Options(g_tol=1E-1); autodiff=:forward)
    
    if verbose
        println(res)
    end

    actor.S = solovev(B0, R0, epsilon, delta, kappa, res.minimizer[1], res.minimizer[2], B0_dir=S0.sigma_B0, Ip_dir=S0.sigma_Ip, symmetric=S0.symmetric, xpoint=S0.xpoint)
    return res
end

#= ====== =#
# FINALIZE #
#= ====== =#
"""
    finalize(actor::SolovevEquilibriumActor, n::Integer=129)::IMAS.equilibrium__time_slice

Store SolovevEquilibriumActor data in IMAS.equilibrium
"""
function finalize(actor::SolovevEquilibriumActor,
                  resolution::Integer=129,
                  rlims::NTuple{2,<:Real}=Equilibrium.limits(actor.S)[1],
                  zlims::NTuple{2,<:Real}=Equilibrium.limits(actor.S)[2])::IMAS.equilibrium__time_slice

    eqt = actor.eq_out
    eqt.profiles_1d.psi = collect(range(Equilibrium.psi_limits(actor.S)..., length=resolution))

    eqt.profiles_1d.pressure = Equilibrium.pressure(actor.S, eqt.profiles_1d.psi)
    eqt.profiles_1d.dpressure_dpsi = Equilibrium.pressure_gradient(actor.S, eqt.profiles_1d.psi)

    eqt.profiles_1d.f = Equilibrium.poloidal_current(actor.S, eqt.profiles_1d.psi)
    eqt.profiles_1d.f_df_dpsi = eqt.profiles_1d.f .* Equilibrium.poloidal_current_gradient(actor.S, eqt.profiles_1d.psi)

    eqt.global_quantities.magnetic_axis.r, eqt.global_quantities.magnetic_axis.z = Equilibrium.magnetic_axis(actor.S)

    resize!(eqt.profiles_2d, 1)
    eqt.profiles_2d[1].grid_type.index = 1
    eqt.profiles_2d[1].grid.dim1 = range(rlims[1], rlims[2], length=resolution)
    eqt.profiles_2d[1].grid.dim2 = range(zlims[1], zlims[2], length=Int(ceil(resolution * actor.S.kappa)))

    eqt.profiles_2d[1].psi = [actor.S(rr, zz) for rr in eqt.profiles_2d[1].grid.dim1, zz in eqt.profiles_2d[1].grid.dim2]

    eqt.profiles_2d[1].b_field_r = zeros(size(eqt.profiles_2d[1].psi)...)
    eqt.profiles_2d[1].b_field_tor = zeros(size(eqt.profiles_2d[1].psi)...)
    eqt.profiles_2d[1].b_field_z = zeros(size(eqt.profiles_2d[1].psi)...)
    for (kr, rr) in enumerate(eqt.profiles_2d[1].grid.dim1), (kz, zz) in enumerate(eqt.profiles_2d[1].grid.dim2)
        (eqt.profiles_2d[1].b_field_r[kr,kz], eqt.profiles_2d[1].b_field_tor[kr,kz], eqt.profiles_2d[1].b_field_z[kr,kz]) = Bfield(actor.S, rr, zz)
    end

    IMAS.flux_surfaces(eqt, actor.S.B0, actor.S.R0)

    return actor.eq_out
end