import TEQUILA

#= =========== =#
#  ActorTEQUILA  #
#= =========== =#
Base.@kwdef mutable struct FUSEparameters__ActorTEQUILA{T} <: ParametersActor where {T<:Real}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set

    #== actor parameters ==#
    free_boundary::Entry{Bool} = Entry{Bool}("-", "Convert fixed boundary equilibrium to free boundary one"; default=true)
    number_of_radial_grid_points::Entry{Int} = Entry{Int}("-", "Number of TEQUILA radial grid points"; default=21)
    number_of_fourier_modes::Entry{Int} = Entry{Int}("-", "Number of modes for Fourier decomposition"; default=20)
    number_of_MXH_harmonics::Entry{Int} = Entry{Int}("-", "Number of Fourier harmonics in MXH representation of flux surfaces"; default=10)
    number_of_iterations::Entry{Int} = Entry{Int}("-", "Number of TEQUILA iterations"; default=20)
    relax::Entry{Float64} = Entry{Float64}("-", "Relaxation on the Picard iterations"; default=1.0)
    tolerance::Entry{Float64} = Entry{Float64}("-", "Tolerance for terminating iterations"; default=1e-6)
    psi_norm_boundary_cutoff::Entry{Float64} = Entry{Float64}("-", "Cutoff psi_norm for determining boundary"; default=0.999)
    #== data flow parameters ==#
    ip_from::Switch{Union{Symbol,Missing}} = set_ip_from()
    #== display and debugging parameters ==#
    do_plot::Entry{Bool} = Entry{Bool}("-", "Plot before and after actor"; default=false)
    debug::Entry{Bool} = Entry{Bool}("-", "Print debug information withing TEQUILA solve"; default=false)
end

mutable struct ActorTEQUILA{D,P} <: PlasmaAbstractActor
    dd::IMAS.dd{D}
    par::FUSEparameters__ActorTEQUILA{P}
    shot::Union{Nothing,TEQUILA.Shot}
    ψbound::D
end

"""
    ActorTEQUILA(dd::IMAS.dd, act::ParametersAllActors; kw...)

Runs the Fixed boundary equilibrium solver TEQUILA
"""
function ActorTEQUILA(dd::IMAS.dd, act::ParametersAllActors; ip_from=:core_profiles, kw...)
    par = act.ActorTEQUILA(ip_from, kw...)
    actor = ActorTEQUILA(dd, par)
    step(actor)
    finalize(actor)
    return actor
end

function ActorTEQUILA(dd::IMAS.dd, par::FUSEparameters__ActorTEQUILA; kw...)
    logging_actor_init(ActorTEQUILA)
    par = par(kw...)
    return ActorTEQUILA(dd, par, nothing, 0.0)

end

"""
    _step(actor::ActorTEQUILA)

Runs TEQUILA on the r_z boundary, equilibrium pressure and equilibrium j_tor
"""
function _step(actor::ActorTEQUILA)
    dd = actor.dd
    par = actor.par

    par.do_plot && (p = plot(dd.equilibrium; cx=true, label="before"))

    eqt = dd.equilibrium.time_slice[]
    eq1d = eqt.profiles_1d

    # BCL 5/30/23: ψbound should be set time dependently, related to the flux swing of the OH coils
    #              For now setting to zero as initial eq1d.psi profile from prepare() can be nonsense
    actor.ψbound = 0.0 #IMAS.interp1d(eq1d.psi_norm, eq1d.psi)(par.psi_norm_boundary_cutoff)

    r_bound = eqt.boundary.outline.r
    z_bound = eqt.boundary.outline.z

    mxh = IMAS.MXH(r_bound, z_bound, par.number_of_MXH_harmonics; optimize_fit=true)

    rho_pol = sqrt.(eq1d.psi_norm)
    rho_pol[1] = 0.0
    P = TEQUILA.FE(rho_pol, eq1d.pressure)
    Jt = TEQUILA.FE(rho_pol, eq1d.j_tor)
    Pbnd = eq1d.pressure[end]
    Fbnd = eq1d.f[end]
    Ip = eqt.global_quantities.ip

    # TEQUILA shot
    shot = TEQUILA.Shot(par.number_of_radial_grid_points, par.number_of_fourier_modes, mxh; P, Jt, Pbnd, Fbnd, Ip_target=Ip)
    actor.shot = TEQUILA.solve(shot, par.number_of_iterations; tol=par.tolerance, par.debug, par.relax)

    return actor
end

# finalize by converting TEQUILA shot to dd.equilibrium
function _finalize(actor::ActorTEQUILA)
    try
        tequila2imas(actor.shot, actor.dd.equilibrium; actor.ψbound, actor.par.free_boundary)
    catch e
        display(plot(actor.shot))
        rethrow(e)
    end
    return actor
end

function tequila2imas(shot::TEQUILA.Shot, eq::IMAS.equilibrium; ψbound::Real=0.0, free_boundary::Bool=false)
    eqt = eq.time_slice[]
    eq1d = eqt.profiles_1d

    R0 = shot.R0fe(1.0)
    Z0 = shot.Z0fe(1.0)
    eq.vacuum_toroidal_field.r0 = R0
    @ddtime(eq.vacuum_toroidal_field.b0 = shot.Fbnd / R0)
    eqt.boundary.geometric_axis.r = R0
    eqt.boundary.geometric_axis.z = Z0

    Rax = shot.R0fe(0.0)
    Zax = shot.Z0fe(0.0)
    eqt.global_quantities.magnetic_axis.r = Rax
    eqt.global_quantities.magnetic_axis.z = Zax

    n_grid = 10 * length(shot.ρ)

    psit = shot.C[2:2:end, 1]
    psia = psit[1]
    psib = psit[end]
    eq1d.psi = range(psia, psib, n_grid)
    rhoi = TEQUILA.ρ.(Ref(shot), eq1d.psi)
    eq1d.pressure = MXHEquilibrium.pressure.(Ref(shot), eq1d.psi)
    eq1d.dpressure_dpsi = MXHEquilibrium.pressure_gradient.(Ref(shot), eq1d.psi)
    eq1d.f = TEQUILA.Fpol.(Ref(shot), rhoi)
    eq1d.f_df_dpsi = TEQUILA.Fpol_dFpol_dψ.(Ref(shot), rhoi)

    resize!(eqt.profiles_2d, 1)
    eq2d = eqt.profiles_2d[1]

    R0 = shot.surfaces[1, end]
    Z0 = shot.surfaces[2, end]
    ϵ = shot.surfaces[3, end]
    κ = shot.surfaces[4, end]
    Rdim = min(1.5 * R0 * ϵ, R0) # 50% bigger than the plasma, but a no bigger than R0
    Zdim = κ * Rdim
    Rgrid = range(R0 - Rdim, R0 + Rdim, n_grid)
    Zgrid = range(Z0 - Zdim, Z0 + Zdim, n_grid)

    eq2d.grid.dim1 = Rgrid
    eq2d.grid.dim2 = Zgrid
    eq2d.grid_type.index = 1
    eq2d.psi = zeros(n_grid, n_grid)

    if free_boundary
        # constraints for the private flux region
        Rb, Zb = TEQUILA.MXH(shot.surfaces[:, end])()
        upper_x_point = any(x_point.z > Z0 for x_point in eqt.boundary.x_point)
        lower_x_point = any(x_point.z < Z0 for x_point in eqt.boundary.x_point)
        Rx, Zx = free_boundary_private_flux_constraint(Rb, Zb; upper_x_point, lower_x_point, fraction=0.25, n_points=10)
        # convert from fixed to free boundary equilibrium
        eq2d.psi .= VacuumFields.fixed2free(shot, Int(ceil(length(Rb) / 2)), Rgrid, Zgrid; Rx, Zx, ψbound=psib)
        IMAS.tweak_psi_to_match_psilcfs!(eqt; ψbound)
    else
        Threads.@threads for (i, r) in enumerate(Rgrid)
            for (j, z) in enumerate(Zgrid)
                eq2d.psi[i, j] = shot(r, z) + ψbound
            end
        end
        eq1d.psi .+= ψbound
    end

end
