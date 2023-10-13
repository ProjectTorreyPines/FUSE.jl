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

mutable struct ActorCHEASE{D,P} <: PlasmaAbstractActor
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

    # uniformely distribute points on the boundary
    r_bound = eqt.boundary.outline.r
    z_bound = eqt.boundary.outline.z
    r_bound, z_bound = IMAS.resample_2d_path(r_bound, z_bound; n_points=201)

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
            r_geo, Ip, r_bound, z_bound, 82,
            rho_pol, pressure, j_tor;
            par.rescale_eq_to_ip,
            par.clear_workdir)
    catch e
        display(plot(r_bound, z_bound; marker=:dot, aspect_ratio=:equal))
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
        n_point_shot_boundary = 500
        eqt = dd.equilibrium.time_slice[]
        z_geo = eqt.boundary.geometric_axis.z
        r_bound = eqt.boundary.outline.r
        z_bound = eqt.boundary.outline.z
        # constraints for the private flux region
        upper_x_point = any(x_point.z > z_geo for x_point in eqt.boundary.x_point)
        lower_x_point = any(x_point.z < z_geo for x_point in eqt.boundary.x_point)
        fraction = 0.5
        Rx, Zx = free_boundary_private_flux_constraint(r_bound, z_bound; upper_x_point, lower_x_point, fraction, n_points=Int(ceil(fraction * n_point_shot_boundary)))

        # convert from fixed to free boundary equilibrium
        EQ = MXHEquilibrium.efit(actor.chease.gfile, 1)

        if true || isempty(dd.pf_active.coil)
            n_coils = 100
            psi_free_rz = VacuumFields.encircling_fixed2free(EQ, n_coils; Rx, Zx)
        else
            coils = IMAS_pf_active__coils(dd; green_model=:simple)
            psi_free_rz = VacuumFields.encircling_fixed2free(EQ, coils; Rx, Zx)
        end

        actor.chease.gfile.psirz = psi_free_rz
        # retrace the last closed flux surface (now with x-point) and scale psirz so to match original psi bounds
        EQ = MXHEquilibrium.efit(actor.chease.gfile, 1)
        psi_b = MXHEquilibrium.psi_boundary(EQ; r=EQ.r, z=EQ.z)
        psi_a = EQ.psi_rz(EQ.axis...)
        actor.chease.gfile.psirz = (psi_free_rz .- psi_a) * ((actor.chease.gfile.psi[end] - actor.chease.gfile.psi[1]) / (psi_b - psi_a)) .+ actor.chease.gfile.psi[1]
    end

    # Convert gEQDSK data to IMAS
    try
        gEQDSK2IMAS(actor.chease.gfile, dd.equilibrium)
    catch e
        EQ = MXHEquilibrium.efit(actor.chease.gfile, 1)
        psi_b = MXHEquilibrium.psi_boundary(EQ; r=EQ.r, z=EQ.z)
        psi_a = EQ.psi_rz(EQ.axis...)
        delta = (psi_b - psi_a) / 10.0
        levels = LinRange(psi_b - delta, psi_b + delta, 11)
        display(contour(EQ.r, EQ.z, transpose(actor.chease.gfile.psirz); levels, aspect_ratio=:equal, clim=(levels[1], levels[end])))
        rethrow(e)
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
    resize!(eqt.profiles_2d, 1)
    eq2d = eqt.profiles_2d[1]

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
    return eq2d.psi = g.psirz .* tc["PSI"]
end
