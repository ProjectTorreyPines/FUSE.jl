import Equilibrium
import EFIT
import CHEASE
import ForwardDiff
import Optim

#= ============ =#
#  ActorSolovev  #
#= ============ =#
Base.@kwdef mutable struct ActorSolovev <: PlasmaAbstractActor
    eq::IMAS.equilibrium
    S::Equilibrium.SolovevEquilibrium
end

function ParametersActor(::Type{Val{:ActorSolovev}})
    par = ParametersActor(nothing)
    par.ngrid = Entry(Integer, "", "Grid size (for R, Z follows proportionally to plasma elongation)"; default=129)
    par.qstar = Entry(Real, "", "Initial guess of kink safety factor"; default=1.5)
    par.alpha = Entry(Real, "", "Initial guess of constant relating to beta regime"; default=0.0)
    par.volume = Entry(Real, "m³", "Scalar volume to match (optional)"; default=missing)
    par.area = Entry(Real, "m²", "Scalar area to match (optional)"; default=missing)
    par.verbose = Entry(Bool, "", "verbose"; default=false)
    return par
end

"""
    ActorSolovev(dd::IMAS.dd, act::ParametersActor; kw...)

Solovev equilibrium actor, based on:
“One size fits all” analytic solutions to the Grad–Shafranov equation
Phys. Plasmas 17, 032502 (2010); https://doi.org/10.1063/1.3328818
"""
function ActorSolovev(dd::IMAS.dd, act::ParametersActor; kw...)
    par = act.ActorSolovev(kw...)
    actor = ActorSolovev(dd.equilibrium)
    step(actor; par.verbose)
    finalize(actor; par.ngrid, volume=getproperty(par, :volume, missing), area=getproperty(par, :area, missing))
    # record optimized values of qstar and alpha in `act` for subsequent ActorSolovev calls
    par.qstar = actor.S.qstar
    par.alpha = actor.S.alpha
    return actor
end

"""
    function ActorSolovev(eq::IMAS.equilibrium; qstar = 1.5, alpha = 0.0)

Constructor for the ActorSolovev structure
“One size fits all” analytic solutions to the Grad–Shafranov equation
Phys. Plasmas 17, 032502 (2010); https://doi.org/10.1063/1.3328818

- qstar: Kink safety factor

- alpha: Constant affecting the pressure
"""
function ActorSolovev(eq::IMAS.equilibrium; qstar=1.5, alpha=0.0)
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

    ActorSolovev(eq, S0)
end

"""
    step(actor::ActorSolovev; verbose=false)

Non-linear optimization to obtain a target `ip` and `beta_normal`
"""
function step(actor::ActorSolovev; verbose=false)
    S0 = actor.S

    eqt = actor.eq.time_slice[]
    target_ip = abs(eqt.global_quantities.ip)
    target_beta = eqt.global_quantities.beta_normal

    B0, R0, epsilon, delta, kappa, alpha, qstar, target_ip, target_beta = promote(S0.B0, S0.R0, S0.epsilon, S0.delta, S0.kappa, S0.alpha, S0.qstar, target_ip, target_beta)

    function cost(x)
        # NOTE: Ip/Beta calculation is very much off in Equilibrium.jl for diverted plasmas because boundary calculation is wrong
        S = Equilibrium.solovev(abs(B0), R0, epsilon, delta, kappa, x[1], x[2], B0_dir=sign(B0), Ip_dir=1, symmetric=true, x_point=nothing)
        beta_cost = (Equilibrium.beta_n(S) - target_beta) / target_beta
        ip_cost = (Equilibrium.plasma_current(S) - target_ip) / target_ip
        c = beta_cost^2 + ip_cost^2
        return c
    end

    res = Optim.optimize(cost, [alpha, qstar], Optim.NelderMead(), Optim.Options(g_tol=1E-3))

    if verbose
        println(res)
    end

    actor.S = Equilibrium.solovev(abs(B0), R0, epsilon, delta, kappa, res.minimizer[1], res.minimizer[2], B0_dir=sign(B0), Ip_dir=1, symmetric=S0.symmetric, x_point=S0.x_point)

    # @show Equilibrium.beta_t(actor.S)
    # @show Equilibrium.beta_p(actor.S)
    # @show Equilibrium.beta_n(actor.S)
    # @show Equilibrium.plasma_current(actor.S)

    return res
end

"""
    finalize(
        actor::ActorSolovev;
        ngrid::Int=129,
        volume::Union{Missing,Real}=missing,
        area::Union{Missing,Real}=missing,
        rlims::NTuple{2,<:Real}=(maximum([actor.S.R0 * (1 - actor.S.epsilon * 2), 0.0]), actor.S.R0 * (1 + actor.S.epsilon * 2)),
        zlims::NTuple{2,<:Real}=(-actor.S.R0 * actor.S.epsilon * actor.S.kappa * 1.7, actor.S.R0 * actor.S.epsilon * actor.S.kappa * 1.7)
    )::IMAS.equilibrium__time_slice

Store ActorSolovev data in IMAS.equilibrium format
"""
function finalize(
    actor::ActorSolovev;
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
    gEQDSK2IMAS(GEQDSKFile::GEQDSKFile,eq::IMAS.equilibrium)

Convert IMAS.equilibrium__time_slice to Equilibrium.jl EFIT structure
"""
function gEQDSK2IMAS(g::EFIT.GEQDSKFile, eq::IMAS.equilibrium)

    tc = transform_cocos(1, 11) # chease output is cocos 1 , dd is cocos 11

    eqt = eq.time_slice[]
    eq1d = eqt.profiles_1d
    resize!(eqt.profiles_2d, 1)
    eq2d = eqt.profiles_2d[1]

    @ddtime(eq.vacuum_toroidal_field.b0 = g.bcentr)
    eq.vacuum_toroidal_field.r0 = g.rcentr

    eqt.global_quantities.magnetic_axis.r = g.rmaxis
    eqt.global_quantities.magnetic_axis.z = g.zmaxis
    eqt.global_quantities.ip = g.current

    eq1d.psi = g.psi .* tc["PSI"]
    eq1d.q = g.qpsi
    eq1d.pressure = g.pres
    eq1d.dpressure_dpsi = g.pprime .* tc["PPRIME"]
    eq1d.f = g.fpol .* tc["F"]
    eq1d.f_df_dpsi = g.ffprim .* tc["F_FPRIME"]

    eq2d.grid_type.index = 1
    eq2d.grid.dim1 = g.r
    eq2d.grid.dim2 = g.z
    eq2d.psi = g.psirz .* tc["PSI"]

    IMAS.flux_surfaces(eqt)
end

#= =========== =#
#  ActorCHEASE  #
#= =========== =#

# Defintion of the actor structure
Base.@kwdef mutable struct ActorCHEASE <: PlasmaAbstractActor
    dd::IMAS.dd
    j_tor_from::Symbol
    pressure_from::Symbol
    free_boundary::Bool
    chease::Union{Nothing,CHEASE.Chease}
end

# Definition of the `act` parameters relevant to the actor
function ParametersActor(::Type{Val{:ActorCHEASE}})
    par = ParametersActor(nothing)
    par.j_tor_from = Switch([:core_profiles, :equilibrium], "", "get j_tor from core_profiles or equilibrium"; default=:equilibrium)
    par.pressure_from = Switch([:core_profiles, :equilibrium], "", "get pressure from from core_profiles or equilibrium"; default=:equilibrium)
    par.free_boundary = Entry(Bool, "", "Convert fixed boundary equilibrium to free boundary one"; default=true)
    return par
end

# run actor with `dd` and `act` as arguments
"""
    ActorCHEASE(dd::IMAS.dd, act::ParametersActor; kw...)

This actor runs the Fixed boundary equilibrium solver CHEASE"""
function ActorCHEASE(dd::IMAS.dd, act::ParametersActor; kw...)
    par = act.ActorCHEASE(kw...)
    actor = ActorCHEASE(dd, par.j_tor_from, par.pressure_from, par.free_boundary, nothing)
    step(actor)
    finalize(actor)
    return actor
end

# define `step` function for this actor
function step(actor::ActorCHEASE)
    dd = actor.dd
    eqt = dd.equilibrium.time_slice[]
    eq1d = eqt.profiles_1d

    if actor.j_tor_from == :equilibrium && actor.pressure_from == :equilibrium
        j_tor = eq1d.j_tor
        pressure = eq1d.pressure
        psin = IMAS.norm01(eq1d.psi)
    elseif actor.j_tor_from == :equilibrium && actor.pressure_from == :core_profiles
        psin = IMAS.norm01(eq1d.psi)
        cp1d = dd.core_profiles.profiles_1d[]
        j_tor = eq1d.j_tor
        pressure = IMAS.interp1d(IMAS.norm01(cp1d.grid.psi), cp1d.pressure_thermal).(psin)
    elseif actor.j_tor_from == :core_profiles && actor.pressure_from == :equilibrium
        psin = IMAS.norm01(eq1d.psi)
        cp1d = dd.core_profiles.profiles_1d[]
        j_tor = IMAS.interp1d(IMAS.norm01(cp1d.grid.psi), cp1d.j_tor).(psin)
        pressure = eq1d.pressure
    else
        cp1d = dd.core_profiles.profiles_1d[]
        j_tor = cp1d.j_tor
        pressure = cp1d.pressure_thermal
        psin = IMAS.norm01(cp1d.grid.psi)
    end

    r_bound = eqt.boundary.outline.r
    z_bound = eqt.boundary.outline.z
    index = (z_bound .> minimum(z_bound) * 0.98) .& (z_bound .< maximum(z_bound) * 0.98)
    r_bound = r_bound[index]
    z_bound = z_bound[index]

    Bt_center = @ddtime(dd.equilibrium.vacuum_toroidal_field.b0)
    r_center = dd.equilibrium.vacuum_toroidal_field.r0

    r_geo = eqt.boundary.geometric_axis.r
    z_geo = eqt.boundary.geometric_axis.z
    Bt_geo = Bt_center * r_center / r_geo
    pressure_sep = pressure[end]

    ϵ = eqt.boundary.minor_radius / r_geo
    Ip = eqt.global_quantities.ip

    # Signs aren't conveyed properly 
    actor.chease = CHEASE.run_chease(ϵ, z_geo, pressure_sep, abs(Bt_geo), r_geo, abs(Ip), r_bound, z_bound, 82, psin, pressure, abs.(j_tor), clear_workdir=true)

    if actor.free_boundary
        # convert from fixed to free boundary equilibrium
        EQ = Equilibrium.efit(actor.chease.gfile, 1)
        psi_free_rz = VacuumFields.fixed2free(EQ, actor.chease.gfile.nbbbs)
        actor.chease.gfile.psirz = psi_free_rz
        EQ = Equilibrium.efit(actor.chease.gfile, 1)
        psi_b = Equilibrium.psi_boundary(EQ; r=EQ.r, z=EQ.z)
        psi_a = EQ.psi_rz(EQ.axis...)
        actor.chease.gfile.psirz = (psi_free_rz .- psi_a) * ((actor.chease.gfile.psi[end] - actor.chease.gfile.psi[1]) / (psi_b - psi_a)) .+ actor.chease.gfile.psi[1]
    end
end

# define `finalize` function for this actor
function finalize(actor::ActorCHEASE)
    gEQDSK2IMAS(actor.chease.gfile, actor.dd.equilibrium)
end