using Equilibrium
import ForwardDiff
import Optim

#= ============ =#
#  SolovevActor  #
#= ============ =#
mutable struct SolovevActor <: AbstractActor
    eq::IMAS.equilibrium
    S::SolovevEquilibrium
end

function ActorParameters(::Type{Val{:SolovevActor}})
    par = ActorParameters(nothing)
    par.ngrid = Entry(Integer, "", "ngrid"; default=129)
    par.volume = Entry(Real, "m³", "Scalar volume to match (optional)"; default=missing)
    par.area = Entry(Real, "m²", "Scalar area to match (optional)"; default=missing)
    par.verbose = Entry(Bool, "", "verbose"; default=false)
    return par
end

function SolovevActor(dd::IMAS.dd, act::ActorParameters; kw...)
    par = act.SolovevActor(kw...)
    actor = SolovevActor(dd.equilibrium)
    step(actor; par.verbose)
    finalize(actor, ngrid=par.ngrid, volume=evalmissing(par, :volume), area=evalmissing(par, :area))
end

"""
    function SolovevActor(eq::IMAS.equilibrium; qstar = 1.5, alpha = 0.0)

Constructor for the SolovevActor structure
“One size fits all” analytic solutions to the Grad–Shafranov equation
Phys. Plasmas 17, 032502 (2010); https://doi.org/10.1063/1.3328818

- qstar: Kink safety factor

- alpha: Constant affecting the pressure
"""
function SolovevActor(eq::IMAS.equilibrium; qstar=1.5, alpha=0.0)
    eqt = eq.time_slice[]
    a = eqt.boundary.minor_radius
    R0 = eqt.boundary.geometric_axis.r
    κ = eqt.boundary.elongation
    δ = eqt.boundary.triangularity
    ϵ = a / R0
    B0 = @ddtime eq.vacuum_toroidal_field.b0

    # check number of x_points
    if mod(length(eqt.boundary.x_point), 2) == 0
        symmetric = true
    else
        symmetric = false
    end

    if length(eqt.boundary.x_point) > 0
        x_point = (eqt.boundary.x_point[1].r, -abs(eqt.boundary.x_point[1].z))
    else
        x_point = nothing
    end
    S0 = Equilibrium.solovev(abs(B0), R0, ϵ, δ, κ, alpha, qstar, B0_dir=Int64(sign(B0)), Ip_dir=1, x_point=x_point, symmetric=symmetric)

    SolovevActor(eq, S0)
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

    Equilibrium.efit(Equilibrium.cocos(11), # COCOS
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

"""
    step(actor::SolovevActor; verbose=false)

Non-linear optimization to obtain a target `ip` and `beta_normal`
"""
function step(actor::SolovevActor; verbose=false)
    S0 = actor.S

    eqt = actor.eq.time_slice[]
    target_ip = abs(eqt.global_quantities.ip)
    target_beta = eqt.global_quantities.beta_normal

    B0, R0, epsilon, delta, kappa, alpha, qstar, target_ip, target_beta = promote(S0.B0, S0.R0, S0.epsilon, S0.delta, S0.kappa, S0.alpha, S0.qstar, target_ip, target_beta)

    function cost(x)
        # NOTE: Ip/Beta calculation is very much off in Equilibrium.jl for diverted plasmas because boundary calculation is wrong
        S = solovev(abs(B0), R0, epsilon, delta, kappa, x[1], x[2], B0_dir=sign(B0), Ip_dir=1, symmetric=true, x_point=nothing)
        beta_cost = (Equilibrium.beta_n(S) - target_beta) / target_beta
        ip_cost = (Equilibrium.plasma_current(S) - target_ip) / target_ip
        c = sqrt(beta_cost^2 + ip_cost^2)
        return c
    end

    res = Optim.optimize(cost, [alpha, qstar], Optim.NelderMead(), Optim.Options(g_tol=1E-3))

    if verbose
        println(res)
    end

    actor.S = solovev(abs(B0), R0, epsilon, delta, kappa, res.minimizer[1], res.minimizer[2], B0_dir=sign(B0), Ip_dir=1, symmetric=S0.symmetric, x_point=S0.x_point)

    # @show Equilibrium.beta_t(actor.S)
    # @show Equilibrium.beta_p(actor.S)
    # @show Equilibrium.beta_n(actor.S)
    # @show Equilibrium.plasma_current(actor.S)

    return res
end

"""
    function finalize(
        actor::SolovevActor;
        ngrid::Int = 129,
        rlims::NTuple{2,<:Real} = (maximum([actor.S.R0 * (1 - actor.S.epsilon * 2), 0.0]), actor.S.R0 * (1 + actor.S.epsilon * 2)),
        zlims::NTuple{2,<:Real} = (-actor.S.R0 * actor.S.epsilon * actor.S.kappa * 2, actor.S.R0 * actor.S.epsilon * actor.S.kappa * 2)
    )::IMAS.equilibrium__time_slice

Store SolovevActor data in IMAS.equilibrium format
"""
function finalize(
    actor::SolovevActor;
    ngrid::Int=129,
    volume::Union{Missing,Real}=missing,
    area::Union{Missing,Real}=missing,
    rlims::NTuple{2,<:Real}=(maximum([actor.S.R0 * (1 - actor.S.epsilon * 2), 0.0]), actor.S.R0 * (1 + actor.S.epsilon * 2)),
    zlims::NTuple{2,<:Real}=(-actor.S.R0 * actor.S.epsilon * actor.S.kappa * 1.7, actor.S.R0 * actor.S.epsilon * actor.S.kappa * 1.7)
)::IMAS.equilibrium__time_slice

    tc = transform_cocos(3, 11)

    eq = actor.eq
    eqt = eq.time_slice[]
    ip = eqt.global_quantities.ip
    sign_Ip = sign(ip)
    sign_Bt = sign(eqt.profiles_1d.f[end])

    Z0 = eqt.boundary.geometric_axis.z
    flip_z = 1.0
    if mod(length(eqt.boundary.x_point), 2) == 1 && eqt.boundary.x_point[1].z > Z0
        flip_z = -1.0
    end

    eq.vacuum_toroidal_field.r0 = actor.S.R0
    @ddtime eq.vacuum_toroidal_field.b0 = actor.S.B0 * sign_Bt

    empty!(eqt)

    eqt.global_quantities.ip = ip
    eqt.boundary.geometric_axis.r = actor.S.R0
    eqt.boundary.geometric_axis.z = Z0
    orig_psi = collect(range(Equilibrium.psi_limits(actor.S)..., length=ngrid))
    eqt.profiles_1d.psi = orig_psi * (tc["PSI"] * sign_Ip)

    eqt.profiles_1d.pressure = Equilibrium.pressure(actor.S, orig_psi)
    eqt.profiles_1d.dpressure_dpsi = Equilibrium.pressure_gradient(actor.S, orig_psi) / (tc["PSI"] * sign_Ip)

    eqt.profiles_1d.f = Equilibrium.poloidal_current(actor.S, orig_psi) * (tc["F"] * sign_Bt)
    eqt.profiles_1d.f_df_dpsi = Equilibrium.poloidal_current(actor.S, orig_psi) .* Equilibrium.poloidal_current_gradient(actor.S, orig_psi) * (tc["F_FPRIME"] * sign_Bt * sign_Ip)

    resize!(eqt.profiles_2d, 1)
    eqt.profiles_2d[1].grid_type.index = 1
    eqt.profiles_2d[1].grid.dim1 = range(rlims[1], rlims[2], length=ngrid)
    eqt.profiles_2d[1].grid.dim2 = range(zlims[1] + Z0, zlims[2] + Z0, length=Int(ceil(ngrid * actor.S.kappa)))

    eqt.profiles_2d[1].psi = [actor.S(rr, flip_z * (zz - Z0)) * (tc["PSI"] * sign_Ip) for rr in eqt.profiles_2d[1].grid.dim1, zz in eqt.profiles_2d[1].grid.dim2]
    IMAS.flux_surfaces(eqt)

    # correct equilibrium volume and area
    if !ismissing(volume)
        eqt.profiles_1d.volume .*= volume / eqt.profiles_1d.volume[end]
    end
    if !ismissing(area)
        eqt.profiles_1d.area .*= area / eqt.profiles_1d.area[end]
    end

    return eqt
end