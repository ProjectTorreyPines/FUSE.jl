import CHEASE

#= =========== =#
#  ActorCHEASE  #
#= =========== =#
Base.@kwdef mutable struct FUSEparameters__ActorCHEASE{T} <: ParametersActor where {T<:Real}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    #== actor parameters ==#
    free_boundary::Entry{Bool} = Entry{Bool}("-", "Convert fixed boundary equilibrium to free boundary one"; default=true)
    clear_workdir::Entry{Bool} = Entry{Bool}("-", "Clean the temporary workdir for CHEASE"; default=true)
    rescale_eq_to_ip::Entry{Bool} = Entry{Bool}("-", "Scale equilibrium to match Ip"; default=true)
    #== data flow parameters ==#
    ip_from::Switch{Symbol} = switch_get_from(:ip)
end

mutable struct ActorCHEASE{D,P} <: PlasmaAbstractActor{D,P}
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
    eq1d = eqt.profiles_1d

    # boundary
    pr = eqt.boundary.outline.r
    pz = eqt.boundary.outline.z
    pr, pz = limit_curvature(pr, pz, (maximum(pr) - minimum(pr)) / 20.0)

    # scalars
    Ip = eqt.global_quantities.ip
    Bt_center = @ddtime(dd.equilibrium.vacuum_toroidal_field.b0)
    r_center = dd.equilibrium.vacuum_toroidal_field.r0
    r_geo = eqt.boundary.geometric_axis.r
    z_geo = eqt.boundary.geometric_axis.z
    Bt_geo = Bt_center * r_center / r_geo
    ϵ = eqt.boundary.minor_radius / r_geo

    # pressure and j_tor
    psin = eq1d.psi_norm
    j_tor = [sign(j) == sign(Ip) ? j : 0.0 for j in eq1d.j_tor]
    pressure = eq1d.pressure
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
        display(plot(pr, pz; marker=:dot, aspect_ratio=:equal))
        display(plot(rho_pol, pressure; marker=:dot, xlabel="sqrt(ψ)", title="Pressure [Pa]"))
        display(plot(rho_pol, j_tor; marker=:dot, xlabel="sqrt(ψ)", title="Jtor [A]"))
        rethrow(e)
    end

    return actor
end

function _finalize(actor::ActorCHEASE)
    dd = actor.dd
    par = actor.par

    # convert from fixed to free boundary equilibrium
    if par.free_boundary

        RA = actor.chease.gfile.rmaxis
        ZA = actor.chease.gfile.zmaxis

        EQ = MXHEquilibrium.efit(actor.chease.gfile, 1)
        psib = 0.0# MXHEquilibrium.psi_boundary(EQ; r=EQ.r, z=EQ.z)
        ψbound = 0.0

        eqt = dd.equilibrium.time_slice[]
        z_geo = eqt.boundary.geometric_axis.z

        # constraints for the private flux region
        z_geo = eqt.boundary.geometric_axis.z
        Rb, Zb = eqt.boundary.outline.r, eqt.boundary.outline.z
        upper_x_point = any(x_point.z > z_geo for x_point in eqt.boundary.x_point)
        lower_x_point = any(x_point.z < z_geo for x_point in eqt.boundary.x_point)
        fraction = 0.05
        Rx, Zx = free_boundary_private_flux_constraint(Rb, Zb; upper_x_point, lower_x_point, fraction, n_points=2)

        # Flux Control Points
        flux_cps = VacuumFields.FluxControlPoints(Rx, Zx, psib)
        append!(flux_cps, VacuumFields.boundary_control_points(EQ, 0.999, psib))
        append!(flux_cps, [VacuumFields.FluxControlPoint(eqt.boundary.x_point[1].r, eqt.boundary.x_point[1].z, psib)])

        # Saddle Control Points
        saddle_weight = length(flux_cps)
        saddle_cps = [VacuumFields.SaddleControlPoint(x_point.r, x_point.z, saddle_weight) for x_point in eqt.boundary.x_point]

        # Coils locations
        if isempty(dd.pf_active.coil)
            coils = encircling_coils(Rb, Zb, RA, ZA, 8)
        else
            coils = IMAS_pf_active__coils(dd; green_model=:simple)
        end

        #display(contour(EQ.r, EQ.z,actor.chease.gfile.psirz';aspect_ratio=:equal,levels=range(-7,15,100)))

        # from fixed boundary to free boundary via VacuumFields
        psi_free_rz = VacuumFields.fixed2free(EQ, coils, EQ.r, EQ.z; flux_cps, saddle_cps, ψbound=psib, λ_regularize=-1.0)
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
        IMAS.tweak_psi_to_match_psilcfs!(dd.equilibrium.time_slice[]; ψbound)
        pf_current_limits(dd.pf_active, dd.build)
    end

    return actor
end

"""
    gEQDSK2IMAS(GEQDSKFile::GEQDSKFile,eq::IMAS.equilibrium)

Convert IMAS.equilibrium__time_slice to MXHEquilibrium.jl EFIT structure
"""
function gEQDSK2IMAS(g::EFIT.GEQDSKFile, eq::IMAS.equilibrium)
    tc = MXHEquilibrium.transform_cocos(1, 11) # chease output is cocos 1 , dd is cocos 11

    eqt = eq.time_slice[]
    eq1d = eqt.profiles_1d
    eq2d = resize!(eqt.profiles_2d, 1)[1]

    @ddtime(eq.vacuum_toroidal_field.b0 = g.bcentr)
    eq.vacuum_toroidal_field.r0 = g.rcentr

    eqt.global_quantities.magnetic_axis.r = g.rmaxis
    eqt.boundary.geometric_axis.r = g.rcentr
    eqt.boundary.geometric_axis.z = g.zmid
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

    return nothing
end
