import FRESCO

#= =========== =#
#  ActorFRESCO  #
#= =========== =#
@actor_parameters_struct ActorFRESCO{T} begin
    #== actor parameters ==#
    control::Switch{Symbol} = Switch{Symbol}([:vertical, :shape, :magnetics], "-", "Vertical control algorithm to be used"; default=:shape)
    number_of_iterations::Entry{Tuple{Int,Int}} = Entry{Tuple{Int,Int}}("-", "Number of outer and inner iterations"; default=(100, 3))
    relax::Entry{Float64} = Entry{Float64}("-", "Relaxation on the Picard iterations"; default=0.5)
    tolerance::Entry{Float64} = Entry{Float64}("-", "Tolerance for terminating iterations"; default=1e-4)
    fixed_grid::Switch{Symbol} = Switch{Symbol}([:psi_norm, :rho_tor_norm], "-", "Fix P and Jt on this grid"; default=:rho_tor_norm)
    nR::Entry{Int} = Entry{Int}("-", "Grid resolution along R"; default=129)
    nZ::Entry{Int} = Entry{Int}("-", "Grid resolution along Z"; default=129)
    vacuum_only::Entry{Bool} = Entry{Bool}("-", "Only compute the vacuum ψ of the dd.pf_active currents and append it as a new equilibrium time slice at dd.global_time (no Grad-Shafranov solve, pf_active currents untouched)"; default=false)
    vacuum_r0_b0_from::Switch{Symbol} = switch_get_from(:vacuum_r0_b0; default=:pulse_schedule)
    #== display and debugging parameters ==#
    do_plot::Entry{Bool} = act_common_parameters(; do_plot=false)
    debug::Entry{Bool} = Entry{Bool}("-", "Print debug information withing FRESCO solve"; default=false)
end

mutable struct ActorFRESCO{D,P} <: CompoundAbstractActor{D,P}
    dd::IMAS.DD{D}
    par::OverrideParameters{P,FUSEparameters__ActorFRESCO{P}}
    act::ParametersAllActors{P}
    canvas::Union{Nothing,FRESCO.Canvas}
    profile::Union{Nothing,FRESCO.PressureJt}
    converged::Bool
end

"""
    ActorFRESCO(dd::IMAS.DD, act::ParametersAllActors; kw...)

Solves the tokamak MHD equilibrium using the FRESCO free-boundary equilibrium solver.

FRESCO (Free-boundary Reactor Equilibrium Solver for Comprehensive Optimization) solves the
Grad-Shafranov equation for toroidal plasma equilibria with specified pressure and current profiles.

This actor performs:

  - Free-boundary equilibrium reconstruction using pressure and current profiles from core_profiles
  - Iterative solution of the Grad-Shafranov equation on a rectangular R-Z grid
  - Vertical control options (vertical position, shape control, or magnetics-based)
  - Updates to poloidal field coil currents for force balance
  - Conversion of equilibrium solution back to IMAS equilibrium format

The solver supports:

  - Configurable grid resolution (nR × nZ)
  - Multiple control algorithms for vertical stability
  - Active X-point specification
  - Green's function table reuse for computational efficiency
  - Relaxation and tolerance control for convergence
  - Vacuum-only mode (`vacuum_only=true`): appends the coils vacuum ψ at `dd.global_time` as a new
    equilibrium time slice, e.g. for breakdown studies. `ActorEquilibrium` uses it below `vacuum_ip_threshold`.

!!! note

    Reads data from `dd.equilibrium`, `dd.core_profiles`, `dd.pf_active`, and `dd.wall`,
    and updates `dd.equilibrium` and `dd.pf_active` with the solved equilibrium
"""
function ActorFRESCO(dd::IMAS.DD, act::ParametersAllActors; kw...)
    actor = ActorFRESCO(dd, act.ActorFRESCO, act; kw...)
    step(actor)
    finalize(actor)
    return actor
end

function ActorFRESCO(dd::IMAS.DD{D}, par::FUSEparameters__ActorFRESCO{P}, act::ParametersAllActors{P}; kw...) where {D<:Real,P<:Real}
    logging_actor_init(ActorFRESCO)
    par = OverrideParameters(par; kw...)
    return ActorFRESCO(dd, par, act, nothing, nothing, true)
end

"""
    _step(actor::ActorFRESCO)

Runs FRESCO on the r_z boundary, equilibrium pressure and equilibrium j_tor
"""
function _step(actor::ActorFRESCO{D,P}) where {D<:Real,P<:Real}
    dd = actor.dd
    par = actor.par
    act = actor.act

    # FRESCO requires wall information
    if IMAS.hasdata(dd.wall)
        fw_r, fw_z = IMAS.first_wall(dd.wall)
    elseif IMAS.hasdata(dd.pf_active) && !isempty(dd.equilibrium.time_slice) && findfirst(:rectangular, dd.equilibrium.time_slice[].profiles_2d) !== nothing
        fw_r, fw_z = IMAS.first_wall(dd.equilibrium.time_slice[], dd.pf_active)
    elseif IMAS.hasdata(dd.pf_active)
        fw_r, fw_z = IMAS.first_wall(dd.pf_active)
    else
        error("FRESCO needs at least pf_active.coils to be filled")
    end
    ΔR = maximum(fw_r) - minimum(fw_r)
    ΔZ = maximum(fw_z) - minimum(fw_z)
    Rs = range(max(0.01, minimum(fw_r) - ΔR / 20), maximum(fw_r) + ΔR / 20, par.nR)
    Zs = range(minimum(fw_z) - ΔZ / 20, maximum(fw_z) + ΔZ / 20, par.nZ)

    if par.vacuum_only
        return _step_vacuum(actor, Rs, Zs, fw_r, fw_z)
    end

    eqt = dd.equilibrium.time_slice[]

    # reuse green table if possible
    gt_kw = isnothing(actor.canvas) ? (;) : (; Green_table=actor.canvas.Green_table)

    # control points (NOTE: we purposely do not pass dd.pulse_schedule.position_control because we want the eqt_control_points to be as ActorEquilibrium set it up
    eqt_control_points = VacuumFields.equilibrium_control_points(eqt, IMAS.pulse_schedule__position_control{D}();
        act.ActorPFactive.boundary_weight,
        act.ActorPFactive.x_points_weight,
        act.ActorPFactive.strike_points_weight)
    mag_control_points = VacuumFields.magnetic_control_points(dd.magnetics;
        act.ActorPFactive.flux_loop_weight,
        act.ActorPFactive.magnetic_probe_weight)

    # Umbrella: a sufficiently degenerate input map (e.g. GSLite's committed
    # psi during/after a disruption) can make FRESCO throw during Canvas
    # construction or the solve itself (axis-guess BoundsError on a fresh
    # canvas with no previous axis, contour failures, ...). A throw here would
    # bypass the converged-guard in _finalize, so catch everything and treat
    # it exactly like a non-converged solve: keep the previous equilibrium
    # for this substep (restored by ActorEquilibrium) and move on.
    try
        actor.canvas = FRESCO.Canvas(
            dd,
            Rs,
            Zs;
            wall_r=fw_r,
            wall_z=fw_z,
            load_pf_passive=false,
            iso_cps=[eqt_control_points.iso_cps; mag_control_points.iso_cps],
            saddle_cps=eqt_control_points.saddle_cps,
            flux_cps=[eqt_control_points.flux_cps; mag_control_points.flux_cps],
            field_cps=mag_control_points.field_cps,
            gt_kw...
        )

        actor.profile = FRESCO.PressureJt(dd; grid=par.fixed_grid)

        status = FRESCO.solve!(actor.canvas, actor.profile, par.number_of_iterations...; par.relax, par.debug, par.control, par.tolerance)
        actor.converged = (status == 0)
    catch e
        isa(e, InterruptException) && rethrow(e)
        @warn "ActorFRESCO: solve threw — treating as non-converged, keeping previous equilibrium for this step" exception = e maxlog = 20
        actor.converged = false
    end

    return actor
end

"""
    _step_vacuum(actor::ActorFRESCO, Rs, Zs, fw_r, fw_z)

Vacuum-only mode (e.g. breakdown / field-null studies): evaluates the coils flux with FRESCO's
own `set_Ψvac!` using the `dd.pf_active` currents at `dd.global_time`.
The Canvas and its Green table are built once and reused as long as the grid is unchanged.
"""
function _step_vacuum(actor::ActorFRESCO{D,P}, Rs, Zs, fw_r, fw_z) where {D<:Real,P<:Real}
    dd = actor.dd
    canvas = actor.canvas
    if canvas === nothing || canvas.Rs != Rs || canvas.Zs != Zs
        coils = VacuumFields.MultiCoils(dd; load_pf_active=true, load_pf_passive=false)
        actor.canvas = canvas = FRESCO.Canvas(; Rs, Zs, Ip=zero(D), Fbnd=zero(D), coils, Rw=fw_r, Zw=fw_z, Rb_target=D[], Zb_target=D[])
    else
        for (icoil, mcoil) in zip(dd.pf_active.coil, canvas.coils)
            VacuumFields.set_current_per_turn!(mcoil, VacuumFields.current_per_turn(icoil))
        end
    end
    FRESCO.set_Ψvac!(canvas)
    actor.converged = true
    return actor
end

# vacuum-only mode: append the vacuum ψ as a new time slice at dd.global_time
function _finalize_vacuum(actor::ActorFRESCO)
    dd = actor.dd
    eq = dd.equilibrium
    canvas = actor.canvas

    # get r0 and b0 before adding the time slice in case vacuum_r0_b0_from=:equilibrium
    r0, b0 = IMAS.get_from(dd, Val(:vacuum_r0_b0), actor.par.vacuum_r0_b0_from)

    eqt = resize!(eq.time_slice)
    eqt.global_quantities.ip = 0.0
    eqt.global_quantities.free_boundary = 1
    eqt.global_quantities.vacuum_toroidal_field.r0 = r0
    eqt.global_quantities.vacuum_toroidal_field.b0 = b0
    @ddtime(eq.vacuum_toroidal_field.b0 = b0)
    eq.vacuum_toroidal_field.r0 = r0
    eq2d = resize!(eqt.profiles_2d, 1)[1]
    eq2d.grid_type.index = 1
    eq2d.grid.dim1 = canvas.Rs
    eq2d.grid.dim2 = canvas.Zs
    eq2d.psi = copy(canvas._Ψvac)
    return actor
end

# finalize by converting FRESCO canvas to dd.equilibrium and updating currents in dd.pf_active and dd.pf_passive
function _finalize(actor::ActorFRESCO)
    if actor.par.vacuum_only
        return _finalize_vacuum(actor)
    end
    if !actor.converged
        # A non-converged canvas carries an unreliable current distribution
        # (observed: recomputed ip at ~50% of target, inflated psi span whose
        # edge-anchored pressure reconstruction then ratchets betap/li).
        # Keep dd's existing equilibrium (the committed/previous slice) for
        # this substep instead of finalizing garbage.
        @warn "ActorFRESCO: solve did not converge — keeping previous equilibrium for this step" maxlog = 20
        return actor
    end
    canvas = actor.canvas
    profile = actor.profile
    dd = actor.dd
    eq = dd.equilibrium
    eqt = eq.time_slice[]
    eqt1d = eqt.profiles_1d
    eq2d = resize!(eqt.profiles_2d, 1)[1]

    eqt.global_quantities.magnetic_axis.r = canvas.Raxis
    eqt.global_quantities.magnetic_axis.z = canvas.Zaxis
    eqt.global_quantities.psi_boundary = canvas.Ψbnd
    eqt.global_quantities.psi_axis = canvas.Ψaxis

    Npsi = length(eqt1d.psi)
    eqt1d.psi = range(canvas.Ψaxis, canvas.Ψbnd, Npsi)
    eqt1d.dpressure_dpsi = FRESCO.Pprime.(Ref(canvas), Ref(profile), eqt1d.psi_norm)
    eqt1d.f_df_dpsi = FRESCO.FFprime.(Ref(canvas), Ref(profile), eqt1d.psi_norm)

    @ddtime(eq.vacuum_toroidal_field.b0 = eqt.global_quantities.vacuum_toroidal_field.b0)
    eq.vacuum_toroidal_field.r0 = eqt.global_quantities.vacuum_toroidal_field.r0

    fend = eqt.global_quantities.vacuum_toroidal_field.b0 * eqt.global_quantities.vacuum_toroidal_field.r0
    f2 = 2 * IMAS.cumtrapz(eqt1d.psi, eqt1d.f_df_dpsi)
    f2 .= f2 .- f2[end] .+ fend^2
    eqt1d.f = sign(fend) .* sqrt.(f2)

    pend = eqt1d.pressure[end]
    eqt1d.pressure = IMAS.cumtrapz(eqt1d.psi, eqt1d.dpressure_dpsi)
    eqt1d.pressure .+= pend .- eqt1d.pressure[end]

    eq2d.grid_type.index = 1
    eq2d.grid.dim1 = canvas.Rs
    eq2d.grid.dim2 = canvas.Zs
    FRESCO.update_interpolation!(canvas)
    eq2d.psi = [canvas._Ψitp(r, z) for r in eq2d.grid.dim1, z in eq2d.grid.dim2]

    # Set the currents in the pf_active and pf_passive
    for (icoil, mcoil) in zip(dd.pf_active.coil, actor.canvas.coils)
        VacuumFields.set_current_per_turn!(icoil, VacuumFields.current_per_turn(mcoil))
    end

    eqt.global_quantities.free_boundary = 1 # free boundary

    return actor
end
