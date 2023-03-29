import CHEASE

#= =========== =#
#  ActorCHEASE  #
#= =========== =#
Base.@kwdef mutable struct FUSEparameters__ActorCHEASE{T} <: ParametersActor where {T<:Real}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    free_boundary::Entry{Bool} = Entry(Bool, "-", "Convert fixed boundary equilibrium to free boundary one"; default=true)
    clear_workdir::Entry{Bool} = Entry(Bool, "-", "Clean the temporary workdir for CHEASE"; default=true)
    rescale_eq_to_ip::Entry{Bool} = Entry(Bool, "-", "Scale equilibrium to match Ip"; default=true)
end

mutable struct ActorCHEASE <: PlasmaAbstractActor
    dd::IMAS.dd
    par::FUSEparameters__ActorCHEASE
    chease::Union{Nothing,CHEASE.Chease}
end

"""
    ActorCHEASE(dd::IMAS.dd, act::ParametersAllActors; kw...)

Runs the Fixed boundary equilibrium solver CHEASE
"""
function ActorCHEASE(dd::IMAS.dd, act::ParametersAllActors; kw...)
    par = act.ActorCHEASE(kw...)
    actor = ActorCHEASE(dd, par)
    step(actor)
    finalize(actor)
    return actor
end

function ActorCHEASE(dd::IMAS.dd, par::FUSEparameters__ActorCHEASE; kw...)
    logging_actor_init(ActorCHEASE)
    par = par(kw...)
    ActorCHEASE(dd, par, nothing)
end

"""
    prepare(dd::IMAS.dd, :ActorCHEASE, act::ParametersAllActors; kw...)

Prepare dd to run ActorCHEASE
* Copy pressure from core_profiles to equilibrium
* Copy j_parallel from core_profiles to equilibrium
"""
function prepare(dd::IMAS.dd, ::Type{Val{:ActorCHEASE}}, act::ParametersAllActors; kw...)
    eq1d = dd.equilibrium.time_slice[].profiles_1d
    cp1d = dd.core_profiles.profiles_1d[]
    eq1d.j_tor = IMAS.interp1d(cp1d.grid.psi_norm, cp1d.j_tor).(eq1d.psi_norm)
    eq1d.pressure = IMAS.interp1d(cp1d.grid.psi_norm, cp1d.pressure).(eq1d.psi_norm)
end

"""
    step(actor::ActorCHEASE)

Runs CHEASE on the r_z boundary, equilibrium pressure and equilibrium j_tor
"""
function _step(actor::ActorCHEASE)
    dd = actor.dd
    eqt = dd.equilibrium.time_slice[]
    eq1d = eqt.profiles_1d
    par = actor.par

    # remove points at high curvature points (ie. X-points)
    r_bound = eqt.boundary.outline.r
    z_bound = eqt.boundary.outline.z
    r_bound, z_bound = IMAS.resample_2d_path(r_bound, z_bound; n_points=201)
    index = abs.(IMAS.curvature(r_bound, z_bound)) .< 0.9
    r_bound = r_bound[index]
    z_bound = z_bound[index]

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
    j_tor = eq1d.j_tor
    pressure = eq1d.pressure
    rho_pol = sqrt.(psin)
    pressure_sep = pressure[end]

    # run and handle errors
    try
        actor.chease = CHEASE.run_chease(
            ϵ, z_geo, pressure_sep, Bt_geo,
            r_geo, Ip, r_bound, z_bound, 82,
            rho_pol, pressure, j_tor,
            rescale_eq_to_ip=par.rescale_eq_to_ip,
            clear_workdir=par.clear_workdir)
    catch
        display(plot(r_bound, z_bound; marker=:dot, aspect_ratio=:equal))
        display(plot(psin, pressure))
        display(plot(psin, abs.(j_tor)))
        rethrow()
    end

    # convert from fixed to free boundary equilibrium
    if par.free_boundary
        EQ = MXHEquilibrium.efit(actor.chease.gfile, 1)
        psi_free_rz = Float64.(VacuumFields.fixed2free(EQ, Int(ceil(length(r_bound)/2))))
        actor.chease.gfile.psirz = psi_free_rz
        # retrace the last closed flux surface (now with x-point) and scale psirz so to match original psi bounds
        EQ = MXHEquilibrium.efit(actor.chease.gfile, 1)
        psi_b = MXHEquilibrium.psi_boundary(EQ; r=EQ.r, z=EQ.z)
        psi_a = EQ.psi_rz(EQ.axis...)
        actor.chease.gfile.psirz = (psi_free_rz .- psi_a) * ((actor.chease.gfile.psi[end] - actor.chease.gfile.psi[1]) / (psi_b - psi_a)) .+ actor.chease.gfile.psi[1]
    end

    return actor
end

# finalize by converting gEQDSK data to IMAS
function finalize(actor::ActorCHEASE)
    try
        gEQDSK2IMAS(actor.chease.gfile, actor.dd.equilibrium)
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
    eq2d.psi = g.psirz .* tc["PSI"]

    IMAS.flux_surfaces(eqt)
end
