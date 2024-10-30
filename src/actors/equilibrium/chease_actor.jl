import CHEASE

#= =========== =#
#  ActorCHEASE  #
#= =========== =#
Base.@kwdef mutable struct FUSEparameters__ActorCHEASE{T<:Real} <: ParametersActor{T}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    _time::Float64 = NaN
    #== actor parameters ==#
    free_boundary::Entry{Bool} = Entry{Bool}("-", "Convert fixed boundary equilibrium to free boundary one"; default=true)
    clear_workdir::Entry{Bool} = Entry{Bool}("-", "Clean the temporary workdir for CHEASE"; default=true)
    rescale_eq_to_ip::Entry{Bool} = Entry{Bool}("-", "Scale equilibrium to match Ip"; default=true)
end

mutable struct ActorCHEASE{D,P} <: SingleAbstractActor{D,P}
    dd::IMAS.dd{D}
    par::FUSEparameters__ActorCHEASE{P}
    chease::Union{Nothing,CHEASE.Chease}
end

"""
    ActorCHEASE(dd::IMAS.dd, act::ParametersAllActors; kw...)

Runs the Fixed boundary equilibrium solver CHEASE
"""
function ActorCHEASE(dd::IMAS.dd, act::ParametersAllActors; kw...)
    actor = ActorCHEASE(dd, act.ActorCHEASE; kw...)
    step(actor)
    finalize(actor)
    return actor
end

function ActorCHEASE(dd::IMAS.dd, par::FUSEparameters__ActorCHEASE; kw...)
    logging_actor_init(ActorCHEASE)
    par = par(kw...)
    return ActorCHEASE(dd, par, nothing)
end

"""
    _step(actor::ActorCHEASE)

Runs CHEASE on the r_z boundary, equilibrium pressure and equilibrium j_tor
"""
function _step(actor::ActorCHEASE)
    dd = actor.dd
    par = actor.par

    # initialize eqt from pulse_schedule and core_profiles
    eqt = dd.equilibrium.time_slice[]
    eqt1d = eqt.profiles_1d

    # boundary
    pr = eqt.boundary.outline.r
    pz = eqt.boundary.outline.z
    ab = sqrt((maximum(pr) - minimum(pr))^2 + (maximum(pz) - minimum(pz))^2) / 2.0
    pr, pz = limit_curvature(pr, pz, ab / 20.0)
    pr, pz = IMAS.resample_2d_path(pr, pz; n_points=2 * length(pr), method=:linear)

    # scalars
    Ip = eqt.global_quantities.ip
    Bt_center = eqt.global_quantities.vacuum_toroidal_field.b0
    r_center = eqt.global_quantities.vacuum_toroidal_field.r0
    r_geo = eqt.boundary.geometric_axis.r
    z_geo = eqt.boundary.geometric_axis.z
    Bt_geo = Bt_center * r_center / r_geo
    ϵ = eqt.boundary.minor_radius / r_geo

    # pressure and j_tor
    psin = eqt1d.psi_norm
    j_tor = [sign(j) == sign(Ip) ? j : 0.0 for j in eqt1d.j_tor]
    pressure = eqt1d.pressure
    rho_pol = sqrt.(psin)
    pressure_sep = pressure[end]

    # run and handle errors
    try
        actor.chease = CHEASE.run_chease(
            ϵ, z_geo, pressure_sep, Bt_geo,
            r_geo, Ip, pr, pz, 82,
            rho_pol, pressure, j_tor;
            par.rescale_eq_to_ip,
            par.clear_workdir)
    catch e
        display(plot(pr, pz; marker=:circle, aspect_ratio=:equal))
        display(plot(rho_pol, pressure; marker=:circle, xlabel="sqrt(ψ)", title="Pressure [Pa]"))
        display(plot(rho_pol, j_tor; marker=:circle, xlabel="sqrt(ψ)", title="Jtor [A]"))
        rethrow(e)
    end

    return actor
end

function _finalize(actor::ActorCHEASE)
    dd = actor.dd
    par = actor.par

    # convert from fixed to free boundary equilibrium
    if par.free_boundary
        eqt = dd.equilibrium.time_slice[]

        RA = actor.chease.gfile.rmaxis
        ZA = actor.chease.gfile.zmaxis

        EQ = MXHEquilibrium.efit(actor.chease.gfile, 1)
        ψbound = 0.0

        # Boundary control points
        iso_cps = VacuumFields.boundary_iso_control_points(EQ, 0.999)

        # Flux control points
        imid = argmax(eqt.boundary.outline.r)
        flux_cps = [VacuumFields.FluxControlPoint(r, z, ψbound, 1.0) for (r, z) in zip(eqt.boundary.outline.r[imid:imid], eqt.boundary.outline.z[imid:imid])]
        strike_weight = 1.0
        strike_cps = [VacuumFields.FluxControlPoint(strike_point.r, strike_point.z, ψbound, strike_weight) for strike_point in eqt.boundary.strike_point]
        append!(flux_cps, strike_cps)

        # Saddle control points
        saddle_weight = 1.0
        saddle_cps = [VacuumFields.SaddleControlPoint(x_point.r, x_point.z, saddle_weight) for x_point in eqt.boundary.x_point]

        # Coils locations
        if isempty(dd.pf_active.coil)
            coils = encircling_coils(eqt.boundary.outline.r, eqt.boundary.outline.z, RA, ZA, 8)
        else
            coils = VacuumFields.IMAS_pf_active__coils(dd; green_model=:quad, zero_currents=true)
        end

        # from fixed boundary to free boundary via VacuumFields
        psi_free_rz = VacuumFields.fixed2free(EQ, coils, EQ.r, EQ.z; iso_cps, flux_cps, saddle_cps, ψbound, λ_regularize=-1.0)
        actor.chease.gfile.psirz .= psi_free_rz'
    end

    # Convert gEQDSK data to IMAS
    try
        gEQDSK2IMAS(actor.chease.gfile, dd.equilibrium)
    catch e
        EQ = MXHEquilibrium.efit(actor.chease.gfile, 1)
        psi_b = MXHEquilibrium.psi_boundary(EQ; r=EQ.r, z=EQ.z)
        psi_a = EQ.psi_rz(EQ.axis...)
        delta = (psi_b - psi_a) / 10.0
        levels = range(psi_b - delta, psi_b + delta, 11)
        display(contour(EQ.r, EQ.z, transpose(actor.chease.gfile.psirz); levels, aspect_ratio=:equal, clim=(levels[1], levels[end])))
        rethrow(e)
    end

    if par.free_boundary
        pf_current_limits(dd.pf_active, dd.build)
    end

    return actor
end

"""
    gEQDSK2IMAS(GEQDSKFile::GEQDSKFile,eq::IMAS.equilibrium)

Convert IMAS.equilibrium__time_slice to MXHEquilibrium.jl EFIT structure
"""
function gEQDSK2IMAS(g::CHEASE.EFIT.GEQDSKFile, eq::IMAS.equilibrium)
    tc = MXHEquilibrium.transform_cocos(1, 11) # chease output is cocos 1 , dd is cocos 11

    eqt = eq.time_slice[]
    eqt1d = eqt.profiles_1d
    eqt2d = resize!(eqt.profiles_2d, 1)[1]

    @ddtime(eq.vacuum_toroidal_field.b0 = g.bcentr)
    eq.vacuum_toroidal_field.r0 = g.rcentr

    eqt.global_quantities.magnetic_axis.r = g.rmaxis
    eqt.boundary.geometric_axis.r = g.rcentr
    eqt.boundary.geometric_axis.z = g.zmid
    eqt.global_quantities.magnetic_axis.z = g.zmaxis
    eqt.global_quantities.ip = g.current

    eqt1d.psi = g.psi .* tc["PSI"]
    eqt1d.q = g.qpsi
    eqt1d.pressure = g.pres
    eqt1d.dpressure_dpsi = g.pprime .* tc["PPRIME"]
    eqt1d.f = g.fpol .* tc["F"]
    eqt1d.f_df_dpsi = g.ffprim .* tc["F_FPRIME"]

    eqt2d.grid_type.index = 1
    eqt2d.grid.dim1 = g.r
    eqt2d.grid.dim2 = g.z
    eqt2d.psi = g.psirz .* tc["PSI"]

    return nothing
end
