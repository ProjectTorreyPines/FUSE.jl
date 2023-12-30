import TEQUILA

#= =========== =#
#  ActorTEQUILA  #
#= =========== =#
Base.@kwdef mutable struct FUSEparameters__ActorTEQUILA{T} <: ParametersActor where {T<:Real}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    #== actor parameters ==#
    free_boundary::Entry{Bool} = Entry{Bool}("-", "Convert fixed boundary equilibrium to free boundary one"; default=true)
    number_of_radial_grid_points::Entry{Int} = Entry{Int}("-", "Number of TEQUILA radial grid points"; default=11)
    number_of_fourier_modes::Entry{Int} = Entry{Int}("-", "Number of modes for Fourier decomposition"; default=10)
    number_of_MXH_harmonics::Entry{Int} = Entry{Int}("-", "Number of Fourier harmonics in MXH representation of flux surfaces"; default=4)
    number_of_iterations::Entry{Int} = Entry{Int}("-", "Number of TEQUILA iterations"; default=20)
    relax::Entry{Float64} = Entry{Float64}("-", "Relaxation on the Picard iterations"; default=0.5)
    tolerance::Entry{Float64} = Entry{Float64}("-", "Tolerance for terminating iterations"; default=1e-3)
    #== data flow parameters ==#
    ip_from::Switch{Symbol} = switch_get_from(:ip)
    #== display and debugging parameters ==#
    do_plot::Entry{Bool} = Entry{Bool}("-", "Plot before and after actor"; default=false)
    debug::Entry{Bool} = Entry{Bool}("-", "Print debug information withing TEQUILA solve"; default=false)
end

mutable struct ActorTEQUILA{D,P} <: PlasmaAbstractActor{D,P}
    dd::IMAS.dd{D}
    par::FUSEparameters__ActorTEQUILA{P}
    shot::Union{Nothing,TEQUILA.Shot}
    ψbound::D
    old_boundary_outline_r::Vector{D}
    old_boundary_outline_z::Vector{D}
end

"""
    ActorTEQUILA(dd::IMAS.dd, act::ParametersAllActors; kw...)

Runs the Fixed boundary equilibrium solver TEQUILA
"""
function ActorTEQUILA(dd::IMAS.dd, act::ParametersAllActors; kw...)
    actor = ActorTEQUILA(dd, act.ActorTEQUILA; kw...)
    step(actor)
    finalize(actor)
    return actor
end

function ActorTEQUILA(dd::IMAS.dd{D}, par::FUSEparameters__ActorTEQUILA{P}; kw...) where {D<:Real,P<:Real}
    logging_actor_init(ActorTEQUILA)
    par = par(kw...)
    return ActorTEQUILA(dd, par, nothing, 0.0, D[], D[])
end

"""
    _step(actor::ActorTEQUILA)

Runs TEQUILA on the r_z boundary, equilibrium pressure and equilibrium j_tor
"""
function _step(actor::ActorTEQUILA)
    dd = actor.dd
    par = actor.par
    eqt = dd.equilibrium.time_slice[]
    eq1d = eqt.profiles_1d

    # BCL 5/30/23: ψbound should be set time dependently, related to the flux swing of the OH coils
    #              For now setting to zero as initial eq1d.psi profile from prepare() can be nonsense
    actor.ψbound = 0.0

    Ip_target = eqt.global_quantities.ip
    rho_pol = sqrt.(eq1d.psi_norm)
    rho_pol[1] = 0.0
    P = TEQUILA.FE(rho_pol, eq1d.pressure)
    Jt = TEQUILA.FE(rho_pol, [sign(j) == sign(Ip_target) ? j : 0.0 for j in eq1d.j_tor]) # don't allow for negative current
    Pbnd = eq1d.pressure[end]
    Fbnd = eq1d.f[end]

    # TEQUILA shot
    if actor.shot === nothing || actor.old_boundary_outline_r != eqt.boundary.outline.r || actor.old_boundary_outline_z != eqt.boundary.outline.z
        pr = eqt.boundary.outline.r
        pz = eqt.boundary.outline.z
        pr, pz = limit_curvature(pr, pz, (maximum(pr) - minimum(pr)) / 20.0)
        # kmin = argmin(pz)
        # @views pr[kmin] = sum(pr[kmin-1:kmin+1]) / 3.0
        # @views pz[kmin] = sum(pz[kmin-1:kmin+1]) / 3.0
        mxh = IMAS.MXH(pr, pz, par.number_of_MXH_harmonics; spline=true)
        # p=scatter(pr, pz)
        # plot!(p, mxh, lw=2)
        # display(p)
        actor.shot = TEQUILA.Shot(par.number_of_radial_grid_points, par.number_of_fourier_modes, mxh; P, Jt, Pbnd, Fbnd, Ip_target)
        solve_function = TEQUILA.solve
        concentric_first = true
    else
        # reuse flux surface information if boundary has not changed
        actor.shot = TEQUILA.Shot(actor.shot; P, Jt, Pbnd, Fbnd, Ip_target)
        solve_function = TEQUILA.solve!
        concentric_first = false
    end

    # solve
    try
        actor.shot = solve_function(actor.shot, par.number_of_iterations; tol=par.tolerance, par.debug, par.relax, concentric_first)
    catch e
        display(plot(eqt.boundary.outline.r, eqt.boundary.outline.z; marker=:dot, aspect_ratio=:equal))
        display(plot(rho_pol, eq1d.pressure; marker=:dot, xlabel="sqrt(ψ)", title="Pressure [Pa]"))
        display(plot(rho_pol, eq1d.j_tor; marker=:dot, xlabel="sqrt(ψ)", title="Jtor [A]"))
        rethrow(e)
    end

    actor.old_boundary_outline_r = eqt.boundary.outline.r
    actor.old_boundary_outline_z = eqt.boundary.outline.z

    return actor
end

# finalize by converting TEQUILA shot to dd.equilibrium
function _finalize(actor::ActorTEQUILA)
    try
        tequila2imas(actor.shot, actor.dd; actor.ψbound, actor.par.free_boundary)
    catch e
        display(plot(actor.shot))
        display(contour())
        rethrow(e)
    end
    return actor
end

function tequila2imas(shot::TEQUILA.Shot, dd::IMAS.dd; ψbound::Real=0.0, free_boundary::Bool=false)
    eq = dd.equilibrium
    eqt = eq.time_slice[]
    eq1d = eqt.profiles_1d

    R0 = shot.R0fe(1.0)
    Z0 = shot.Z0fe(1.0)
    eq.vacuum_toroidal_field.r0 = R0
    @ddtime(eq.vacuum_toroidal_field.b0 = shot.Fbnd / R0)
    eqt.boundary.geometric_axis.r = R0
    eqt.boundary.geometric_axis.z = Z0

    RA = shot.R0fe(0.0)
    ZA = shot.Z0fe(0.0)
    eqt.global_quantities.magnetic_axis.r = RA
    eqt.global_quantities.magnetic_axis.z = ZA

    R0 = shot.surfaces[1, end]
    Z0 = shot.surfaces[2, end]
    ϵ = shot.surfaces[3, end]
    κ = shot.surfaces[4, end]
    a = R0 * ϵ
    Rdim = min(1.5 * a, R0) # 50% bigger than the plasma, but a no bigger than R0
    Zdim = κ * 1.6 * a

    nψ_grid = 129
    nz_grid = 129
    nr_grid = Int(ceil(nz_grid * Rdim / Zdim))

    psit = shot.C[2:2:end, 1]
    psia = psit[1]
    psib = psit[end]
    eq1d.psi = range(psia, psib, nψ_grid)
    rhoi = TEQUILA.ρ.(Ref(shot), eq1d.psi)
    eq1d.pressure = MXHEquilibrium.pressure.(Ref(shot), eq1d.psi)
    eq1d.dpressure_dpsi = MXHEquilibrium.pressure_gradient.(Ref(shot), eq1d.psi)
    eq1d.f = TEQUILA.Fpol.(Ref(shot), rhoi)
    eq1d.f_df_dpsi = TEQUILA.Fpol_dFpol_dψ.(Ref(shot), rhoi)

    resize!(eqt.profiles_2d, 2)

    # MXH flux surface parametrization
    eq2d = eqt.profiles_2d[1]
    eq2d.grid.dim1 = shot.ρ
    eq2d.grid.dim2 = vcat([0.0, 0.0, 0.0, 0.0, 0.0], fill(1.0, (shot.M,)), fill(-1.0, (shot.M,)))
    eq2d.grid_type.index = 57 # inverse_rhotor_mxh : Flux surface type with radial label sqrt[Phi/pi/B0] (dim1), Phi being toroidal flux, and MXH coefficients R0, Z0, ϵ, κ, c0, c[...], s[...] (dim2)
    eq2d.psi = collect(shot.surfaces')

    # RZ
    Rgrid = range(R0 - Rdim, R0 + Rdim, nr_grid)
    Zgrid = range(Z0 - Zdim, Z0 + Zdim, nz_grid)
    eq2d = eqt.profiles_2d[2]
    eq2d.grid.dim1 = Rgrid
    eq2d.grid.dim2 = Zgrid
    eq2d.grid_type.index = 1
    eq2d.psi = fill(Inf, (length(eq2d.grid.dim1), length(eq2d.grid.dim2)))
    if free_boundary
        pr = eqt.boundary.outline.r
        pz = eqt.boundary.outline.z
        pr, pz = limit_curvature(pr, pz, (maximum(pr) - minimum(pr)) / 20.0)

        # Flux Control Points
        flux_cps = VacuumFields.boundary_control_points(shot, 0.999, psib)

        # Saddle Control Points
        saddle_weight = length(flux_cps) / length(eqt.boundary.x_point)
        saddle_cps = [VacuumFields.SaddleControlPoint(x_point.r, x_point.z, saddle_weight) for x_point in eqt.boundary.x_point]

        if isempty(dd.pf_active.coil)
            coils = encircling_coils(pr, pz, RA, ZA, 8)
        else
            coils = IMAS_pf_active__coils(dd; green_model=:simple)
        end

        psi_free_rz = VacuumFields.fixed2free(shot, coils, Rgrid, Zgrid; flux_cps, saddle_cps, ψbound=psib, λ_regularize=-1.0)
        eq2d.psi .= psi_free_rz'

        IMAS.tweak_psi_to_match_psilcfs!(eqt; ψbound)
        pf_current_limits(dd.pf_active, dd.build)
    else
        # to work with a closed boundary equilibrium for now we need
        # ψ outside of the CLFS to grow out until it touches the computation domain
        for (i, r) in collect(enumerate(Rgrid))
            for (j, z) in enumerate(Zgrid)
                for (r0, z0) in zip(eqt.boundary.outline.r, eqt.boundary.outline.z)
                    eq2d.psi[i, j] = min(eq2d.psi[i, j], sqrt((r - r0) .^ 2 + (z - z0) .^ 2) * (psib - psia) / 4)
                end
            end
        end
        Threads.@threads for (i, r) in collect(enumerate(Rgrid))
            for (j, z) in enumerate(Zgrid)
                value = shot(r, z)
                if value != 0.0
                    eq2d.psi[i, j] = value + ψbound
                end
            end
        end
        eq1d.psi .+= ψbound
    end

end
