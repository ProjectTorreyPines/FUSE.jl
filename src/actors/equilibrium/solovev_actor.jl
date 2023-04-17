import MXHEquilibrium
import EFIT
import Optim

#= ============ =#
#  ActorSolovev  #
#= ============ =#
Base.@kwdef mutable struct FUSEparameters__ActorSolovev{T} <: ParametersActor where {T<:Real}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    ngrid::Entry{Int} = Entry(Int, "-", "Grid size (for R, Z follows proportionally to plasma elongation)"; default=129)
    qstar::Entry{T} = Entry(T, "-", "Initial guess of kink safety factor"; default=1.5)
    alpha::Entry{T} = Entry(T, "-", "Initial guess of constant relating to pressure"; default=0.0)
    volume::Entry{T} = Entry(T, "m³", "Scalar volume to match (optional)"; default=missing)
    area::Entry{T} = Entry(T, "m²", "Scalar area to match (optional)"; default=missing)
    verbose::Entry{Bool} = Entry(Bool, "-", "verbose"; default=false)
end

mutable struct ActorSolovev <: PlasmaAbstractActor
    dd::IMAS.dd
    par::FUSEparameters__ActorSolovev
    S::Union{Nothing,MXHEquilibrium.SolovevEquilibrium}
end

"""
    ActorSolovev(dd::IMAS.dd, act::ParametersAllActors; kw...)

Solovev equilibrium actor, based on:
“One size fits all” analytic solutions to the Grad–Shafranov equation
Phys. Plasmas 17, 032502 (2010); https://doi.org/10.1063/1.3328818
"""
function ActorSolovev(dd::IMAS.dd, act::ParametersAllActors; kw...)
    par = act.ActorSolovev
    actor = ActorSolovev(dd, par; kw...)
    step(actor)
    finalize(actor)
    # record optimized values of qstar and alpha in `act` for subsequent ActorSolovev calls
    par.qstar = actor.S.qstar
    par.alpha = actor.S.alpha
    return actor
end

function ActorSolovev(dd::IMAS.dd, par::FUSEparameters__ActorSolovev; kw...)
    logging_actor_init(ActorSolovev)
    par = par(kw...)
    return ActorSolovev(dd, par, nothing)
end

"""
    _step(actor::ActorSolovev)

Non-linear optimization to obtain a target `ip` and `pressure_core`
"""
function _step(actor::ActorSolovev)
    dd = actor.dd
    par = actor.par

    # initialize eqt from pulse_schedule and core_profiles
    prepare_eq(dd)
    eq = dd.equilibrium
    eqt = eq.time_slice[]

    # magnetic field
    B0 = @ddtime eq.vacuum_toroidal_field.b0
    target_ip = abs(eqt.global_quantities.ip)
    target_pressure_core = eqt.profiles_1d.pressure[1]

    # plasma shape as MXH
    pr, pz = eqt.boundary.outline.r, eqt.boundary.outline.z
    pr, pz = IMAS.resample_2d_path(pr, pz)
    pr, pz = IMAS.reorder_flux_surface!(pr, pz)
    mxh = IMAS.MXH(pr, pz, 4)
    Z0off = mxh.Z0 # Solovev has a bug for Z!=0.0
    mxh.Z0 -= Z0off
    plasma_shape = MXHEquilibrium.MillerExtendedHarmonicShape(mxh.R0, mxh.Z0, mxh.ϵ, mxh.κ, mxh.c0, mxh.c, mxh.s)
    # check number of x_points to infer symmetry
    if mod(length(eqt.boundary.x_point), 2) == 0
        symmetric = true
    else
        symmetric = false
    end

    # add x_point info
    if length(eqt.boundary.x_point) > 0
        x_point = (eqt.boundary.x_point[1].r, -abs(eqt.boundary.x_point[1].z) - Z0off)
    else
        x_point = nothing
    end

    # first run of Solovev
    S0 = MXHEquilibrium.solovev(abs(B0), plasma_shape, par.alpha, par.qstar; B0_dir=Int64(sign(B0)), Ip_dir=1, x_point, symmetric)

    # optimize `alpha` and `qstar` to get target_ip and target_pressure_core
    function cost(x)
        S = MXHEquilibrium.solovev(abs(B0), plasma_shape, x[1], x[2]; B0_dir=Int64(sign(B0)), Ip_dir=1, symmetric, x_point)
        psimag, psibry = MXHEquilibrium.psi_limits(S)
        pressure_cost = (MXHEquilibrium.pressure(S, psimag) - target_pressure_core) / target_pressure_core
        ip_cost = (MXHEquilibrium.plasma_current(S) - target_ip) / target_ip
        return sqrt(pressure_cost^2 + 100.0 * ip_cost^2)
    end
    res = Optim.optimize(cost, [S0.alpha, S0.qstar], Optim.NelderMead(), Optim.Options(g_tol=1E-3))
    if par.verbose
        println(res)
    end

    # final Solovev solution
    actor.S = MXHEquilibrium.solovev(abs(B0), plasma_shape, res.minimizer[1], res.minimizer[2]; B0_dir=Int64(sign(B0)), Ip_dir=1, symmetric, x_point)

    return actor
end

"""
    _finalize(actor::ActorSolovev)

Store ActorSolovev data in IMAS.equilibrium format
"""
function _finalize(actor::ActorSolovev)
    dd = actor.dd
    par = actor.par

    rlims = MXHEquilibrium.limits(actor.S.S; pad=0.3)[1]
    zlims = MXHEquilibrium.limits(actor.S.S; pad=0.3)[2]

    ngrid = par.ngrid
    tc = MXHEquilibrium.transform_cocos(3, 11)

    eq = dd.equilibrium
    eqt = eq.time_slice[]
    target_ip = eqt.global_quantities.ip
    sign_Ip = sign(target_ip)
    sign_Bt = sign(eqt.profiles_1d.f[end])
    target_psi_norm = getproperty(eqt.profiles_1d, :psi_norm, missing)
    target_pressure = getproperty(eqt.profiles_1d, :pressure, missing)
    target_j_tor = getproperty(eqt.profiles_1d, :j_tor, missing)

    Z0 = eqt.boundary.geometric_axis.z
    flip_z = 1.0
    if mod(length(eqt.boundary.x_point), 2) == 1 && eqt.boundary.x_point[1].z > Z0
        flip_z = -1.0
    end

    eq.vacuum_toroidal_field.r0 = actor.S.S.R0
    @ddtime eq.vacuum_toroidal_field.b0 = actor.S.B0 * sign_Bt

    empty!(eqt)

    eqt.global_quantities.ip = target_ip
    eqt.boundary.geometric_axis.r = actor.S.S.R0
    eqt.boundary.geometric_axis.z = Z0
    orig_psi = collect(range(MXHEquilibrium.psi_limits(actor.S)..., length=ngrid))
    eqt.profiles_1d.psi = orig_psi * (tc["PSI"] * sign_Ip)

    eqt.profiles_1d.pressure = MXHEquilibrium.pressure.(actor.S, orig_psi)
    eqt.profiles_1d.dpressure_dpsi = MXHEquilibrium.pressure_gradient.(actor.S, orig_psi) ./ (tc["PSI"] * sign_Ip)

    eqt.profiles_1d.f = MXHEquilibrium.poloidal_current.(actor.S, orig_psi) .* (tc["F"] * sign_Bt)
    eqt.profiles_1d.f_df_dpsi = MXHEquilibrium.poloidal_current.(actor.S, orig_psi) .* MXHEquilibrium.poloidal_current_gradient.(actor.S, orig_psi) .* (tc["F_FPRIME"] * sign_Bt * sign_Ip)

    resize!(eqt.profiles_2d, 1)
    eqt.profiles_2d[1].grid_type.index = 1
    eqt.profiles_2d[1].grid.dim1 = range(rlims[1], rlims[2], length=ngrid)
    eqt.profiles_2d[1].grid.dim2 = range(zlims[1] + Z0, zlims[2] + Z0, length=Int(ceil(ngrid * actor.S.S.κ)))
    eqt.profiles_2d[1].psi = [actor.S(rr, flip_z * (zz - Z0)) * (tc["PSI"] * sign_Ip) for rr in eqt.profiles_2d[1].grid.dim1, zz in eqt.profiles_2d[1].grid.dim2]

    IMAS.flux_surfaces(eqt)

    if true
        # force total plasma current to target_ip to avoid drifting after multiple calls of SolovevActor
        eqt.profiles_2d[1].psi = (eqt.profiles_2d[1].psi .- eqt.profiles_1d.psi[end]) .* (target_ip / eqt.global_quantities.ip) .+ eqt.profiles_1d.psi[end]
        eqt.profiles_1d.psi = (eqt.profiles_1d.psi .- eqt.profiles_1d.psi[end]) .* (target_ip / eqt.global_quantities.ip) .+ eqt.profiles_1d.psi[end]
        # match entry target_pressure and target_j_tor as if Solovev could do this
        if !ismissing(target_pressure)
            eqt.profiles_1d.pressure = IMAS.interp1d(target_psi_norm, target_pressure).(eqt.profiles_1d.psi_norm)
        end
        if !ismissing(target_j_tor)
            eqt.profiles_1d.j_tor = IMAS.interp1d(target_psi_norm, target_j_tor, :cubic).(eqt.profiles_1d.psi_norm)
        end
        IMAS.p_jtor_2_pprime_ffprim_f!(eqt.profiles_1d, actor.S.S.R0, actor.S.B0)
        IMAS.flux_surfaces(eqt)
    end

    # this is to fix the boundary back to its original input
    # since Solovev will not satisfy the original boundary
    # this is not needed since we the boundary info is now stored in pulse_schedule
    if false
        eqt.boundary.outline.r, eqt.boundary.outline.z = actor.mxh()
        eqt.boundary.elongation = actor.mxh.κ
        eqt.boundary.triangularity = sin(actor.mxh.s[1])
        eqt.boundary.squareness = -actor.mxh.s[2]
        eqt.boundary.minor_radius = actor.mxh.ϵ * actor.mxh.R0
        eqt.boundary.geometric_axis.r = actor.mxh.R0
        eqt.boundary.geometric_axis.z = actor.mxh.Z0
    end

    # correct equilibrium volume and area
    if !ismissing(par, :volume)
        eqt.profiles_1d.volume .*= par.volume / eqt.profiles_1d.volume[end]
    end
    if !ismissing(par, :area)
        eqt.profiles_1d.area .*= par.area / eqt.profiles_1d.area[end]
    end

    return actor
end
