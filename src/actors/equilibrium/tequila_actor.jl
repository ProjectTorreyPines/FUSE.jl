import TEQUILA

#= =========== =#
#  ActorTEQUILA  #
#= =========== =#
Base.@kwdef mutable struct FUSEparameters__ActorTEQUILA{T} <: ParametersActor where {T<:Real}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    free_boundary::Entry{Bool} = Entry(Bool, "-", "Convert fixed boundary equilibrium to free boundary one"; default=true)
    number_of_radial_grid_points::Entry{Int} = Entry(Int, "-", "Number of TEQUILA radial grid points"; default=20)
    number_of_fourier_modes::Entry{Int} = Entry(Int, "-", "Number of modes for Fourier decomposition"; default=20)
    number_of_iterations::Entry{Int} = Entry(Int, "-", "Number of TEQUILA iterations"; default=5)
    do_plot::Entry{Bool} = Entry(Bool, "-", "Plot before and after actor"; default=false)
end

mutable struct ActorTEQUILA <: PlasmaAbstractActor
    dd::IMAS.dd
    par::FUSEparameters__ActorTEQUILA
    shot::Union{Nothing,TEQUILA.Shot}
end

"""
    ActorTEQUILA(dd::IMAS.dd, act::ParametersAllActors; kw...)

Runs the Fixed boundary equilibrium solver TEQUILA
"""
function ActorTEQUILA(dd::IMAS.dd, act::ParametersAllActors; kw...)
    par = act.ActorTEQUILA(kw...)
    actor = ActorTEQUILA(dd, par)
    step(actor)
    finalize(actor)
    return actor
end

function ActorTEQUILA(dd::IMAS.dd, par::FUSEparameters__ActorTEQUILA; kw...)
    logging_actor_init(ActorTEQUILA)
    par = par(kw...)
    ActorTEQUILA(dd, par, nothing)
end

"""
    prepare(dd::IMAS.dd, :ActorTEQUILA, act::ParametersAllActors; kw...)

Prepare dd to run ActorTEQUILA
* Copy pressure from core_profiles to equilibrium
* Copy j_parallel from core_profiles to equilibrium
"""
function prepare(dd::IMAS.dd, ::Type{Val{:ActorTEQUILA}}, act::ParametersAllActors; kw...)
    eq1d = dd.equilibrium.time_slice[].profiles_1d
    cp1d = dd.core_profiles.profiles_1d[]
    eq1d.j_tor = IMAS.interp1d(cp1d.grid.psi_norm, cp1d.j_tor).(eq1d.psi_norm)
    eq1d.pressure = IMAS.interp1d(cp1d.grid.psi_norm, cp1d.pressure).(eq1d.psi_norm)
end

"""
    step(actor::ActorTEQUILA)

Runs TEQUILA on the r_z boundary, equilibrium pressure and equilibrium j_tor
"""
function _step(actor::ActorTEQUILA)
    dd = actor.dd
    eqt = dd.equilibrium.time_slice[]
    eq1d = eqt.profiles_1d
    eq2d = eqt.profiles_2d
    par = actor.par

    # remove points at high curvature points (ie. X-points)
    r_bound = eqt.boundary.outline.r
    z_bound = eqt.boundary.outline.z
    r_bound, z_bound = IMAS.resample_2d_path(r_bound, z_bound; n_points=201)
    index = abs.(IMAS.curvature(r_bound, z_bound)) .< 0.9
    r_bound = r_bound[index]
    z_bound = z_bound[index]

    # To be used in next update:
    psin = eq1d.psi_norm
    j_tor = eq1d.j_tor
    pressure = eq1d.pressure
    rho_pol = sqrt.(psin)
    pressure_sep = pressure[end]

    mxh = IMAS.MXH(r_bound, z_bound, 5;optimize_fit=false)

    # TEQUILA shot
    shot = TEQUILA.Shot(par.number_of_radial_grid_points, par.number_of_fourier_modes, mxh)


    # will have to be done correctly (probalby by just passing j_tor and pressure array on rho_pol grid)
    dp_dψ = eq1d.dpressure_dpsi # should be changed to pressure and j_tor passing 
    f_df_dψ = eq1d.f_df_dpsi
    dp(x) = eq1d.dpressure_dpsi[20] #eq1d.dpressure_dpsi 
    fdf(x) = eq1d.f_df_dpsi[20]

    actor.shot = TEQUILA.solve(shot, dp, fdf, par.number_of_iterations)

    # convert from fixed to free boundary equilibrium (this is from chease but similar for teq)
    if par.free_boundary
        """
        EQ = MXHEquilibrium.efit(actor.chease.gfile, 1)
        psi_free_rz = VacuumFields.fixed2free(EQ, Int(ceil(length(r_bound)/2)))
        actor.chease.gfile.psirz = psi_free_rz
        # retrace the last closed flux surface (now with x-point) and scale psirz so to match original psi bounds
        EQ = MXHEquilibrium.efit(actor.chease.gfile, 1)
        psi_b = MXHEquilibrium.psi_boundary(EQ; r=EQ.r, z=EQ.z)
        psi_a = EQ.psi_rz(EQ.axis...)
        actor.chease.gfile.psirz = (psi_free_rz .- psi_a) * ((actor.chease.gfile.psi[end] - actor.chease.gfile.psi[1]) / (psi_b - psi_a)) .+ actor.chease.gfile.psi[1]
        """
    end

    if par.do_plot
        p=plot(dd.equilibrium;cx=true)
        for s in eachcol(actor.shot.surfaces[:,2:end])
            plot!(p,IMAS.MXH(s), color=:red)
        end
        display(p)
    end

    return actor
end

# finalize by converting TEQUILA shot to dd.equilibrium
function _finalize(actor::ActorTEQUILA)
    try
        tequila2imas(actor.shot, dd.equilibrium)
    catch e
        display(TEQUILA.plot_shot(actor.shot))
        rethrow(e)
    end
    return actor
end

function tequila2imas(shot::TEQUILA.Shot, eq::IMAS.equilibrium)
    # New rho-psi map
    # Flux surface trace and put important quantities back to imas

    # IMAS.flux_surfaces(eq.time_slice[])
end



# SAMPLE FUNCTION AS FOR CHEASE FOR INSPRIATON
"""
    gEQDSK2IMAS(GEQDSKFile::GEQDSKFile,eq::IMAS.equilibrium)

Convert IMAS.equilibrium__time_slice to MXHEquilibrium.jl EFIT structure
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
"""
