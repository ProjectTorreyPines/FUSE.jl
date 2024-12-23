import TEQUILA

#= =========== =#
#  ActorTEQUILA  #
#= =========== =#
Base.@kwdef mutable struct FUSEparameters__ActorTEQUILA{T<:Real} <: ParametersActor{T}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    _time::Float64 = NaN
    #== actor parameters ==#
    free_boundary::Entry{Bool} = Entry{Bool}("-", "Convert fixed boundary equilibrium to free boundary one"; default=true)
    number_of_radial_grid_points::Entry{Int} = Entry{Int}("-", "Number of TEQUILA radial grid points"; default=31)
    number_of_fourier_modes::Entry{Int} = Entry{Int}("-", "Number of modes for Fourier decomposition"; default=8)
    number_of_MXH_harmonics::Entry{Int} = Entry{Int}("-", "Number of Fourier harmonics in MXH representation of flux surfaces"; default=4)
    number_of_iterations::Entry{Int} = Entry{Int}("-", "Number of TEQUILA iterations"; default=1000)
    relax::Entry{Float64} = Entry{Float64}("-", "Relaxation on the Picard iterations"; default=0.25, check=x -> @assert 0.0 <= x <= 1.0 "must be: 0.0 <= relax <= 1.0")
    tolerance::Entry{Float64} = Entry{Float64}("-", "Tolerance for terminating iterations"; default=1e-4)
    #== data flow parameters ==#
    fixed_grid::Switch{Symbol} = Switch{Symbol}([:poloidal, :toroidal], "-", "Fix P and Jt on this rho grid"; default=:toroidal)
    #== display and debugging parameters ==#
    do_plot::Entry{Bool} = act_common_parameters(; do_plot=false)
    debug::Entry{Bool} = Entry{Bool}("-", "Print debug information withing TEQUILA solve"; default=false)
    #== IMAS psi grid settings ==#
    R::Entry{Vector{Float64}} = Entry{Vector{Float64}}("m", "Psi R axis")
    Z::Entry{Vector{Float64}} = Entry{Vector{Float64}}("m", "Psi Z axis")
end

mutable struct ActorTEQUILA{D,P} <: CompoundAbstractActor{D,P}
    dd::IMAS.dd{D}
    par::FUSEparameters__ActorTEQUILA{P}
    act::ParametersAllActors{P}
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
    actor = ActorTEQUILA(dd, act.ActorTEQUILA, act; kw...)
    step(actor)
    finalize(actor)
    return actor
end

function ActorTEQUILA(dd::IMAS.dd{D}, par::FUSEparameters__ActorTEQUILA{P}, act::ParametersAllActors{P}; kw...) where {D<:Real,P<:Real}
    logging_actor_init(ActorTEQUILA)
    par = par(kw...)
    return ActorTEQUILA(dd, par, act, nothing, D(0.0), D[], D[])
end

"""
    _step(actor::ActorTEQUILA)

Runs TEQUILA on the r_z boundary, equilibrium pressure and equilibrium j_tor
"""
function _step(actor::ActorTEQUILA)
    dd = actor.dd
    D = eltype(dd)
    par = actor.par
    eqt = dd.equilibrium.time_slice[]
    eqt1d = eqt.profiles_1d

    # BCL 5/30/23: ψbound should be set time dependently, related to the flux swing of the OH coils
    #              For now setting to zero as initial eqt1d.psi profile from prepare() can be nonsense
    actor.ψbound = D(0.0)

    Ip_target = eqt.global_quantities.ip
    if par.fixed_grid === :poloidal
        rho = sqrt.(eqt1d.psi_norm)
        rho[1] = D(0.0)
        P = (TEQUILA.FE(rho, eqt1d.pressure), :poloidal)
        # don't allow current to change sign
        Jt = (TEQUILA.FE(rho, D[sign(j) == sign(Ip_target) ? j : D(0.0) for j in eqt1d.j_tor]), :poloidal)
        Pbnd = eqt1d.pressure[end]
    elseif par.fixed_grid === :toroidal
        rho = eqt1d.rho_tor_norm
        P = (TEQUILA.FE(rho, eqt1d.pressure), :toroidal)
        # don't allow current to change sign
        Jt = (TEQUILA.FE(rho, D[sign(j) == sign(Ip_target) ? j : D(0.0) for j in eqt1d.j_tor]), :toroidal)
        Pbnd = eqt1d.pressure[end]
    end

    Fbnd = eqt.global_quantities.vacuum_toroidal_field.b0 * eqt.global_quantities.vacuum_toroidal_field.r0

    # see if boundary has not changed
    same_boundary = false
    if actor.shot !== nothing
        if length(actor.old_boundary_outline_r) == length(eqt.boundary.outline.r)
            if (sum(abs.(actor.old_boundary_outline_r .- eqt.boundary.outline.r)) + sum(abs.(actor.old_boundary_outline_z .- eqt.boundary.outline.z))) /
               length(eqt.boundary.outline.r) < 1E-3
                same_boundary = true
            end
        end
    end

    # TEQUILA shot
    if actor.shot === nothing || !same_boundary
        pr = eqt.boundary.outline.r
        pz = eqt.boundary.outline.z
        ab = sqrt((maximum(pr) - minimum(pr))^2 + (maximum(pz) - minimum(pz))^2) / 2.0
        pr, pz = limit_curvature(pr, pz, ab / 20.0)
        pr, pz = IMAS.resample_2d_path(pr, pz; n_points=2 * length(pr), method=:linear)
        mxh = IMAS.MXH(pr, pz, par.number_of_MXH_harmonics; spline=true)
        actor.shot = TEQUILA.Shot(par.number_of_radial_grid_points, par.number_of_fourier_modes, mxh; P, Jt, Pbnd, Fbnd, Ip_target)
        solve_function = TEQUILA.solve
        concentric_first = true
        actor.old_boundary_outline_r = eqt.boundary.outline.r
        actor.old_boundary_outline_z = eqt.boundary.outline.z
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
        plot(eqt.boundary.outline.r, eqt.boundary.outline.z; marker=:dot, aspect_ratio=:equal)
        display(plot!(IMAS.MXH(actor.shot.surfaces[:, end])))
        display(plot(rho, eqt1d.pressure; marker=:dot, xlabel="ρ", title="Pressure [Pa]"))
        display(plot(rho, eqt1d.j_tor; marker=:dot, xlabel="ρ", title="Jtor [A]"))
        rethrow(e)
    end

    return actor
end

# finalize by converting TEQUILA shot to dd.equilibrium
function _finalize(actor::ActorTEQUILA)
    try
        tequila2imas(actor.shot, actor.dd, actor.par, actor.act; actor.ψbound)
    catch e
        display(plot(actor.shot))
        rethrow(e)
    end
    return actor
end

function tequila2imas(shot::TEQUILA.Shot, dd::IMAS.dd{D}, par::FUSEparameters__ActorTEQUILA, act::ParametersAllActors; ψbound::D) where {D<:Real}
    free_boundary = par.free_boundary
    eq = dd.equilibrium
    eqt = eq.time_slice[]
    eqt1d = eqt.profiles_1d

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

    nψ_grid = 129

    psit = shot.C[2:2:end, 1]
    psia = psit[1]
    psib = psit[end]
    eqt1d.psi = range(psia, psib, nψ_grid)
    rhoi = TEQUILA.ρ.(Ref(shot), eqt1d.psi)
    eqt1d.pressure = MXHEquilibrium.pressure.(Ref(shot), eqt1d.psi)
    eqt1d.dpressure_dpsi = MXHEquilibrium.pressure_gradient.(Ref(shot), eqt1d.psi)
    eqt1d.f = shot.F.(rhoi)
    eqt1d.f_df_dpsi = TEQUILA.Fpol_dFpol_dψ.(Ref(shot), rhoi; shot.invR, shot.invR2)

    resize!(eqt.profiles_2d, 2)

    # MXH flux surface parametrization
    eq2d = eqt.profiles_2d[1]
    eq2d.grid.dim1 = shot.ρ
    MXH_modes = (size(shot.surfaces, 1) - 5) ÷ 2
    eq2d.grid.dim2 = vcat(@SVector[0.0, 0.0, 0.0, 0.0, 0.0], range(1, MXH_modes), .-range(1, MXH_modes))
    eq2d.grid_type.index = 57 # inverse_rhotor_mxh : Flux surface type with radial label sqrt[Phi/pi/B0] (dim1), Phi being toroidal flux, and MXH coefficients R0, Z0, ϵ, κ, c0, c[...], s[...] (dim2)
    eq2d.psi = collect(shot.surfaces')

    # maximum expected size of divertor regions based on ActorCXbuild settings
    divertor_size = max(act.ActorCXbuild.divertor_hfs_size_fraction, act.ActorCXbuild.divertor_lfs_size_fraction)

    # RZ
    if ismissing(par, :Z)
        Zdim = κ * (1.1 + divertor_size) * a
        nz_grid = nψ_grid
        Zgrid = range(Z0 - Zdim, Z0 + Zdim, nz_grid)
    else
        Zgrid = par.Z
        Zdim = abs(-(extrema(Zgrid)...)) / 2
        nz_grid = length(Zgrid)
    end

    if ismissing(par, :R)
        Rdim = (1.1 + divertor_size) * a # divertor_size% bigger than the plasma, but a no bigger than R0
        nr_grid = Int(ceil(nz_grid * Rdim / Zdim))
        Rgrid = range(R0 - min(Rdim, R0), R0 + Rdim, nr_grid)
    else
        Rgrid = par.R
        Rdim = abs(-(extrema(Rgrid)...)) / 2
        nr_grid = length(Rgrid)
    end

    eq2d = eqt.profiles_2d[2]
    eq2d.grid.dim1 = Rgrid
    eq2d.grid.dim2 = Zgrid
    eq2d.grid_type.index = 1
    eq2d.psi = fill(Inf, (length(eq2d.grid.dim1), length(eq2d.grid.dim2)))

    if free_boundary
        # Boundary control points
        iso_cps = VacuumFields.boundary_iso_control_points(shot, 0.999)

        # Flux control points
        mag = VacuumFields.FluxControlPoint{D}(eqt.global_quantities.magnetic_axis.r, eqt.global_quantities.magnetic_axis.z, psia,  iso_cps[1].weight)
        flux_cps = VacuumFields.FluxControlPoint[mag]
        strike_weight = act.ActorPFactive.strike_points_weight / length(eqt.boundary.strike_point)
        strike_cps = [VacuumFields.FluxControlPoint{D}(strike_point.r, strike_point.z, ψbound, strike_weight) for strike_point in eqt.boundary.strike_point]
        append!(flux_cps, strike_cps)

        # Saddle control points
        saddle_weight = act.ActorPFactive.x_points_weight / length(eqt.boundary.x_point)
        saddle_cps = [VacuumFields.SaddleControlPoint{D}(x_point.r, x_point.z, saddle_weight) for x_point in eqt.boundary.x_point]
        push!(saddle_cps, VacuumFields.SaddleControlPoint{D}(eqt.global_quantities.magnetic_axis.r, eqt.global_quantities.magnetic_axis.z, iso_cps[1].weight))

        # Coils locations
        coils = VacuumFields.IMAS_pf_active__coils(dd; act.ActorPFactive.green_model, zero_currents=true)

        # from fixed boundary to free boundary via VacuumFields
        psi_free_rz = VacuumFields.fixed2free(shot, coils, Rgrid, Zgrid; iso_cps, flux_cps, saddle_cps, ψbound, λ_regularize=-1.0)
        eq2d.psi .= psi_free_rz'

        pf_current_limits(dd.pf_active, dd.build)

    else
        # to work with a closed boundary equilibrium for now we need
        # ψ outside of the CLFS to grow out until it touches the computation domain
        Threads.@threads for (i, r) in collect(enumerate(Rgrid))
            for (j, z) in enumerate(Zgrid)
                eq2d.psi[i, j] = shot(r, z; extrapolate=true) + ψbound
            end
        end
        eqt1d.psi .+= ψbound
    end

end
