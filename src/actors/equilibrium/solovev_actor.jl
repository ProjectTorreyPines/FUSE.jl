import MXHEquilibrium
import Optim

#= ============ =#
#  ActorSolovev  #
#= ============ =#
Base.@kwdef mutable struct FUSEparameters__ActorSolovev{T<:Real} <: ParametersActor{T}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    _time::Float64 = NaN
    #== actor parameters ==#
    ngrid::Entry{Int} = Entry{Int}("-", "Grid size (for R, Z follows proportionally to plasma elongation)"; default=129)
    qstar::Entry{T} = Entry{T}("-", "Initial guess of kink safety factor"; default=1.5)
    alpha::Entry{T} = Entry{T}("-", "Initial guess of constant relating to pressure"; default=0.0)
    #== data flow parameters ==#
    ip_from::Switch{Symbol} = switch_get_from(:ip)
    #== display and debugging parameters ==#
    verbose::Entry{Bool} = act_common_parameters(; verbose=false)
end

mutable struct ActorSolovev{D,P} <: SingleAbstractActor{D,P}
    dd::IMAS.dd{D}
    par::FUSEparameters__ActorSolovev{P}
    S::Union{Nothing,MXHEquilibrium.SolovevEquilibrium}
end

"""
    ActorSolovev(dd::IMAS.dd, act::ParametersAllActors; kw...)

Solovev equilibrium actor, based on:
“One size fits all” analytic solutions to the Grad–Shafranov equation
Phys. Plasmas 17, 032502 (2010); https://doi.org/10.1063/1.3328818
"""
function ActorSolovev(dd::IMAS.dd, act::ParametersAllActors; kw...)
    actor = ActorSolovev(dd, act.ActorSolovev; kw...)
    step(actor)
    finalize(actor)
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
    eq = dd.equilibrium
    eqt = eq.time_slice[]

    # magnetic field
    B0 = @ddtime eq.vacuum_toroidal_field.b0
    target_ip = eqt.global_quantities.ip
    target_pressure_core = eqt.profiles_1d.pressure[1]

    # plasma shape as MXH
    mxh = IMAS.MXH(eqt.boundary.outline.r, eqt.boundary.outline.z, 4)
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
    # display(plot(mxh))
    # @show abs(B0)
    # @show plasma_shape
    # @show par.alpha
    # @show par.qstar
    # @show x_point
    # @show symmetric
    S0 = MXHEquilibrium.solovev(abs(B0), plasma_shape, par.alpha, par.qstar; B0_dir=Int64(sign(B0)), Ip_dir=1, x_point, symmetric)

    # optimize `alpha` and `qstar` to get target_ip and target_pressure_core
    function cost(x)
        S = MXHEquilibrium.solovev(abs(B0), plasma_shape, x[1], x[2]; B0_dir=Int64(sign(B0)), Ip_dir=1, symmetric, x_point)
        psimag, psibry = MXHEquilibrium.psi_limits(S)
        pressure_cost = (MXHEquilibrium.pressure(S, psimag) - target_pressure_core) / target_pressure_core
        ip_cost = (MXHEquilibrium.plasma_current(S) - target_ip) / target_ip
        return sqrt(pressure_cost^2 + 100.0 * ip_cost^2)
    end
    res = Optim.optimize(cost, [S0.alpha, S0.qstar], Optim.NelderMead(), Optim.Options(; g_tol=1E-3))
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
    mxh_eq = actor.S

    eqt = dd.equilibrium.time_slice[]
    eqt1d = eqt.profiles_1d

    target_ip = eqt.global_quantities.ip
    target_psi_norm = getproperty(eqt1d, :psi_norm, missing)
    target_pressure = getproperty(eqt1d, :pressure, missing)
    target_j_tor = getproperty(eqt1d, :j_tor, missing)

    MXHEquilibrium_to_dd!(dd.equilibrium, mxh_eq, par.ngrid; cocos_in=3)

    # force total plasma current to target_ip to avoid drifting after multiple calls of SolovevActor
    eqt2d = findfirst(:rectangular, eqt.profiles_2d)
    eqt2d.psi = (eqt2d.psi .- eqt1d.psi[end]) .* (target_ip / eqt.global_quantities.ip) .+ eqt1d.psi[end]
    eqt1d.psi = (eqt1d.psi .- eqt1d.psi[end]) .* (target_ip / eqt.global_quantities.ip) .+ eqt1d.psi[end]
    # match entry target_pressure and target_j_tor as if Solovev could do this
    if !ismissing(target_pressure)
        eqt1d.pressure = IMAS.interp1d(target_psi_norm, target_pressure).(eqt1d.psi_norm)
    end
    if !ismissing(target_j_tor)
        eqt1d.j_tor = IMAS.interp1d(target_psi_norm, target_j_tor, :cubic).(eqt1d.psi_norm)
    end

    # pprime, ffprim, f from  p and j_tor
    dpressure_dpsi, f_df_dpsi, f = IMAS.p_jtor_2_pprime_ffprim_f(eqt1d, mxh_eq.S.R0, mxh_eq.B0)
    eqt1.dpressure_dpsi = dpressure_dpsi
    eqt1.f_df_dpsi = f_df_dpsi
    eqt1.f = f

    # record optimized values of qstar and alpha in `act` for subsequent calls to the same actor
    actor.par.qstar = actor.S.qstar
    actor.par.alpha = actor.S.alpha

    return actor
end

function MXHEquilibrium_to_dd!(eq::IMAS.equilibrium, mxh_eq::MXHEquilibrium.AbstractEquilibrium, ngrid::Int; cocos_in::Int=11)
    tc = MXHEquilibrium.transform_cocos(cocos_in, 11)

    rlims = MXHEquilibrium.limits(mxh_eq.S; pad=0.3)[1]
    zlims = MXHEquilibrium.limits(mxh_eq.S; pad=0.3)[2]

    eqt = eq.time_slice[]
    Ip = eqt.global_quantities.ip
    sign_Ip = sign(Ip)
    sign_Bt = sign(eqt.profiles_1d.f[end])

    Z0 = eqt.boundary.geometric_axis.z
    flip_z = 1.0
    if mod(length(eqt.boundary.x_point), 2) == 1 && eqt.boundary.x_point[1].z > Z0
        flip_z = -1.0
    end

    eq.vacuum_toroidal_field.r0 = mxh_eq.S.R0
    @ddtime eq.vacuum_toroidal_field.b0 = mxh_eq.B0 * sign_Bt

    t0 = eqt.time
    empty!(eqt)
    eqt.time = t0

    eqt.global_quantities.ip = Ip
    eqt.boundary.geometric_axis.r = mxh_eq.S.R0
    eqt.boundary.geometric_axis.z = Z0
    orig_psi = collect(range(MXHEquilibrium.psi_limits(mxh_eq)..., ngrid))
    eqt.profiles_1d.psi = orig_psi * (tc["PSI"] * sign_Ip)

    eqt.profiles_1d.pressure = MXHEquilibrium.pressure.(mxh_eq, orig_psi)
    eqt.profiles_1d.dpressure_dpsi = MXHEquilibrium.pressure_gradient.(mxh_eq, orig_psi) ./ (tc["PSI"] * sign_Ip)

    eqt.profiles_1d.f = MXHEquilibrium.poloidal_current.(mxh_eq, orig_psi) .* (tc["F"] * sign_Bt)
    eqt.profiles_1d.f_df_dpsi =
        MXHEquilibrium.poloidal_current.(mxh_eq, orig_psi) .* MXHEquilibrium.poloidal_current_gradient.(mxh_eq, orig_psi) .* (tc["F_FPRIME"] * sign_Bt * sign_Ip)

    eqt2d = resize!(eqt.profiles_2d, 1)[1]
    eqt2d.grid_type.index = 1
    eqt2d.grid.dim1 = range(rlims[1], rlims[2], ngrid)
    eqt2d.grid.dim2 = range(zlims[1] + Z0, zlims[2] + Z0, Int(ceil(ngrid * mxh_eq.S.κ)))
    eqt2d.psi = [mxh_eq(rr, flip_z * (zz - Z0)) * (tc["PSI"] * sign_Ip) for rr in eqt2d.grid.dim1, zz in eqt2d.grid.dim2]

    fw = IMAS.first_wall(IMAS.top_dd(eqt).wall)
    IMAS.flux_surfaces(eqt, fw.r, fw.z)

    return
end