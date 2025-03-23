#= ================ =#
#  ActorEquilibrium  #
#= ================ =#
Base.@kwdef mutable struct FUSEparameters__ActorEquilibrium{T<:Real} <: ParametersActor{T}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    _time::Float64 = NaN
    #== actor parameters ==#
    model::Switch{Symbol} = Switch{Symbol}([:TEQUILA, :FRESCO, :EGGO, :CHEASE, :replay, :none], "-", "Equilibrium actor to run"; default=:TEQUILA)
    symmetrize::Entry{Bool} = Entry{Bool}("-", "Force equilibrium up-down symmetry with respect to magnetic axis"; default=false)
    #== data flow parameters ==#
    j_p_from::Switch{Symbol} = Switch{Symbol}([:equilibrium, :core_profiles], "-", "Take j_tor and pressure profiles from this IDS"; default=:core_profiles)
    ip_from::Switch{Symbol} = switch_get_from(:ip)
    vacuum_r0_b0_from::Switch{Symbol} = switch_get_from(:vacuum_r0_b0; default=:pulse_schedule)
    #== display and debugging parameters ==#
    do_plot::Entry{Bool} = act_common_parameters(; do_plot=false)
end

mutable struct ActorEquilibrium{D,P} <: CompoundAbstractActor{D,P}
    dd::IMAS.dd{D}
    par::OverrideParameters{P,FUSEparameters__ActorEquilibrium{P}}
    act::ParametersAllActors{P}
    eq_actor::Union{Nothing,ActorTEQUILA{D,P},ActorFRESCO{D,P},ActorEGGO{D,P},ActorCHEASE{D,P},ActorReplay{D,P},ActorNoOperation{D,P}}
end

"""
    ActorEquilibrium(dd::IMAS.dd, act::ParametersAllActors; kw...)

Provides a common interface to run different equilibrium actors
"""
function ActorEquilibrium(dd::IMAS.dd, act::ParametersAllActors; kw...)
    actor = ActorEquilibrium(dd, act.ActorEquilibrium, act; kw...)
    step(actor)
    finalize(actor)
    return actor
end

function ActorEquilibrium(dd::IMAS.dd, par::FUSEparameters__ActorEquilibrium, act::ParametersAllActors; kw...)
    logging_actor_init(ActorEquilibrium)
    par = OverrideParameters(par; kw...)

    noop = ActorNoOperation(dd, act.ActorNoOperation)
    actor = ActorEquilibrium(dd, par, act, noop)

    if par.model == :FRESCO
        actor.eq_actor = ActorFRESCO(dd, act.ActorFRESCO, act)
    elseif par.model == :CHEASE
        actor.eq_actor = ActorCHEASE(dd, act.ActorCHEASE, act)
    elseif par.model == :TEQUILA
        actor.eq_actor = ActorTEQUILA(dd, act.ActorTEQUILA, act)
    elseif par.model == :EGGO
        actor.eq_actor = ActorEGGO(dd, act.ActorEGGO, act)
    elseif par.model == :replay
        actor.eq_actor = ActorReplay(dd, act.ActorReplay, actor)
    end

    return actor
end

"""
    _step(actor::ActorEquilibrium)

Clears and initializes data in eqt for equilibrium actors to run properly, then calls the `step()` function of the selected equilibrium actor
"""
function _step(actor::ActorEquilibrium)
    dd = actor.dd
    par = actor.par

    if par.do_plot
        if !isempty(dd.equilibrium.time_slice)
            plot(dd.equilibrium; label="before ActorEquilibrium")
        else
            plot()
        end
    end

    if par.model !== :none
        # initialize eqt for equilibrium actors
        prepare(actor)
    end

    # step selected equilibrium actor
    step(actor.eq_actor)

    return actor
end

"""
    _finalize(actor::ActorEquilibrium)

Calls the `finalize()` function of the selected equilibrium actor and populates flux surfaces information
"""
function _finalize(actor::ActorEquilibrium)
    dd = actor.dd
    par = actor.par

    # finalize selected equilibrium actor
    finalize(actor.eq_actor)

    if par.model ∉ (:none, :replay)
        eqt = dd.equilibrium.time_slice[]

        # symmetrize equilibrium if requested and number of X-points is even
        x_points = IMAS.x_points(dd.pulse_schedule.position_control.x_point)
        if par.symmetrize && mod(length(x_points), 2) != 1
            IMAS.symmetrize_equilibrium!(eqt)
        end

        # add flux surfaces information
        fw = IMAS.first_wall(dd.wall)
        try
            IMAS.flux_surfaces(eqt, fw.r, fw.z)
        catch e
            eqt2d = findfirst(:rectangular, eqt.profiles_2d)
            par.do_plot && display(current())
            contour(eqt2d.grid.dim1, eqt2d.grid.dim2, eqt2d.psi'; aspect_ratio=:equal)
            plot!(fw.r, fw.z; color=:gray)
            display(contour!(eqt2d.grid.dim1, eqt2d.grid.dim2, eqt2d.psi'; levels=[0], lw=3, color=:black, colorbar_entry=false))
            rethrow(e)
        end
    end

    if par.do_plot
        try
            display(plot!(dd.equilibrium; label="after ActorEquilibrium"))
        catch e
            if isa(e, BoundsError)
                display(plot(dd.equilibrium; label="after ActorEquilibrium"))
            else
                rethrow(e)
            end
        end
    end

    return actor
end

"""
    prepare(actor::ActorEquilibrium)

Prepare `dd.equilibrium` to run equilibrium actors

  - Clear equilibrium__time_slice
  - Set Ip, Bt from core_profiles, equilibrium, or pulse_schedule
  - Use position control from pulse_schedule
  - Use j_tor,pressure from core_profiles (for self-consistent iterations) or equilibrium (to re-solve equilibrium with different solver)
"""
function prepare(actor::ActorEquilibrium)
    dd = actor.dd
    par = actor.par

    ps = dd.pulse_schedule
    pc = ps.position_control

    # make sure j_tor and pressure on axis come in with zero gradient
    if par.j_p_from == :core_profiles
        @assert !isempty(dd.core_profiles.time)
        cp1d = dd.core_profiles.profiles_1d[]
        index = cp1d.grid.psi_norm .> 0.05
        psi0 = cp1d.grid.psi
        rho_tor_norm0 = cp1d.grid.rho_tor_norm
        rho_pol_norm_sqrt0 = vcat(-reverse(sqrt.(cp1d.grid.psi_norm[index])), sqrt.(cp1d.grid.psi_norm[index]))
        j_tor0 = vcat(reverse(cp1d.j_tor[index]), cp1d.j_tor[index])
        pressure0 = vcat(reverse(cp1d.pressure[index]), cp1d.pressure[index])
        j_itp = IMAS.interp1d(rho_pol_norm_sqrt0, j_tor0, :cubic)
        p_itp = IMAS.interp1d(rho_pol_norm_sqrt0, pressure0, :cubic)
    elseif par.j_p_from == :equilibrium
        @assert !isempty(dd.equilibrium.time)
        eqt1d = dd.equilibrium.time_slice[].profiles_1d
        psi0 = eqt1d.psi
        rho_tor_norm0 = eqt1d.rho_tor_norm
        rho_pol_norm_sqrt0 = sqrt.(eqt1d.psi_norm)
        j_tor0 = eqt1d.j_tor
        pressure0 = eqt1d.pressure
        j_itp = IMAS.interp1d(rho_pol_norm_sqrt0, j_tor0, :cubic)
        p_itp = IMAS.interp1d(rho_pol_norm_sqrt0, pressure0, :cubic)
    else
        @assert par.j_p_from in (:core_profiles, :equilibrium)
    end

    # get ip and b0 before wiping eqt in case ip_from=:equilibrium
    ip = IMAS.get_from(dd, Val{:ip}, actor.par.ip_from)
    r0, b0 = IMAS.get_from(dd, Val{:vacuum_r0_b0}, actor.par.vacuum_r0_b0_from)

    # geometric factors
    past_time_slice = false
    if !isempty(dd.equilibrium.time_slice)
        past_time_slice = true
        eqt = dd.equilibrium.time_slice[]
        psi = eqt.profiles_1d.psi
        gm1 = eqt.profiles_1d.gm1
        gm8 = eqt.profiles_1d.gm8
        gm9 = eqt.profiles_1d.gm9
    end

    # add/clear time-slice
    eqt = resize!(dd.equilibrium.time_slice)
    resize!(eqt.profiles_2d, 1)
    eqt1d = dd.equilibrium.time_slice[].profiles_1d

    # scalar quantities
    eqt.global_quantities.ip = ip
    eqt.global_quantities.vacuum_toroidal_field.b0 = b0
    eqt.global_quantities.vacuum_toroidal_field.r0 = r0

    # boundary from position control
    eqt.boundary.outline.r, eqt.boundary.outline.z = IMAS.boundary(pc)

    # boundary scalars from position control
    eqt.boundary.minor_radius = @ddtime(pc.minor_radius.reference)
    eqt.boundary.geometric_axis.r = @ddtime(pc.geometric_axis.r.reference)
    eqt.boundary.geometric_axis.z = @ddtime(pc.geometric_axis.z.reference)
    eqt.boundary.elongation = @ddtime(pc.elongation.reference)
    eqt.boundary.tilt = @ddtime(pc.tilt.reference)
    eqt.boundary.triangularity = @ddtime(pc.triangularity.reference)
    eqt.boundary.squareness = @ddtime(pc.squareness.reference)
    eqt.boundary.ovality = @ddtime(pc.ovality.reference)
    eqt.boundary.twist = @ddtime(pc.twist.reference)

    # x-points from position control
    # NOTE: we use get_time_array(...,:constant) instead of @ddtime because x-points can suddenly jump
    if !isempty(getproperty(pc, :x_point, []))
        n = 0
        for k in eachindex(pc.x_point)
            rx = IMAS.get_time_array(pc.x_point[k].r, :reference, :constant)
            zx = IMAS.get_time_array(pc.x_point[k].z, :reference, :constant)
            if rx > 0.0 && !isnan(rx) && !isnan(zx)
                n += 1
                resize!(eqt.boundary.x_point, n)
                eqt.boundary.x_point[n].r = rx
                eqt.boundary.x_point[n].z = zx
            end
        end
    end

    # stike-points from position control
    # NOTE: we use get_time_array(...,:constant) instead of @ddtime because strike-points can suddenly jump
    if !isempty(getproperty(pc, :strike_point, []))
        n = 0
        for k in eachindex(pc.strike_point)
            rs = IMAS.get_time_array(pc.strike_point[k].r, :reference, :constant)
            zs = IMAS.get_time_array(pc.strike_point[k].z, :reference, :constant)
            if rs > 0.0 && !isnan(rs) && !isnan(zs)
                n += 1
                resize!(eqt.boundary.strike_point, n)
                eqt.boundary.strike_point[n].r = rs
                eqt.boundary.strike_point[n].z = zs
            end
        end
    end

    # set j_tor and pressure, forcing zero derivative on axis
    eqt1d = dd.equilibrium.time_slice[].profiles_1d
    eqt1d.psi = psi0
    eqt1d.rho_tor_norm = rho_tor_norm0
    eqt1d.j_tor = j_itp.(sqrt.(eqt1d.psi_norm))
    eqt1d.pressure = p_itp.(sqrt.(eqt1d.psi_norm))

    # calculate pressure and j_tor using geometry from previous iteration
    # this is for equilibrium codes that cannot solve directly from pressure and current
    if past_time_slice
        pressure = IMAS.interp1d(psi0, eqt1d.pressure).(psi)
        j_tor = IMAS.interp1d(psi0, eqt1d.j_tor).(psi)
        tmp = IMAS.calc_pprime_ffprim_f(psi, gm8, gm9, gm1, r0, b0; pressure, j_tor)
        eqt1d.dpressure_dpsi = IMAS.interp1d(psi, tmp.dpressure_dpsi).(psi0)
        eqt1d.f_df_dpsi = IMAS.interp1d(psi, tmp.f_df_dpsi).(psi0)
        eqt1d.f = IMAS.interp1d(psi, tmp.f).(psi0)
    end

    # if sign(maximum(eqt1d.j_tor)) != sign(minimum(eqt1d.j_tor))
    #     j_tor = eqt1d.j_tor
    #     s = sign(sum(j_tor))
    #     j_tor = s .* j_tor
    #     min_j = sum(j_tor[j_tor.>0.0]) / length(j_tor) / 100.0
    #     j_tor[j_tor.<=min_j] .= min_j
    #     eqt1d.j_tor = s .* j_tor
    # end

    return dd
end

"""
    latest_equilibrium_grids!(dd)

Empties grids at the `dd.global_time` of `core_profiles`, `core_sources` and `core_transport`
so that the IMAS.jl expressions take the grid info from latest equilibrium. This is necessary
when iterating between equilibrium and other actors.

See: IMAS/src/expressions/onetime.jl
"""
function latest_equilibrium_grids!(dd::IMAS.dd)
    # core_profiles
    old_rho_tor_norm = dd.core_profiles.profiles_1d[].grid.rho_tor_norm
    empty!(dd.core_profiles.profiles_1d[].grid)
    dd.core_profiles.profiles_1d[].grid.rho_tor_norm = old_rho_tor_norm

    # core_sources
    for source in dd.core_sources.source
        old_rho_tor_norm = source.profiles_1d[].grid.rho_tor_norm
        empty!(source.profiles_1d[].grid)
        source.profiles_1d[].grid.rho_tor_norm = old_rho_tor_norm
    end

    # core_transport
    for model in dd.core_transport.model
        old_rho_tor_norm = model.profiles_1d[].grid_flux.rho_tor_norm
        empty!(model.profiles_1d[].grid_flux)
        model.profiles_1d[].grid_flux.rho_tor_norm = old_rho_tor_norm
    end
end

"""
    IMAS2Equilibrium(eqt::IMAS.equilibrium__time_slice)

Convert IMAS.equilibrium__time_slice to MXHEquilibrium.jl EFIT structure
"""
function IMAS2Equilibrium(eqt::IMAS.equilibrium__time_slice)
    eqt2d = findfirst(:rectangular, eqt.profiles_2d)
    dim1 = range(eqt2d.grid.dim1[1], eqt2d.grid.dim1[end], length(eqt2d.grid.dim1))
    @assert norm(dim1 .- eqt2d.grid.dim1) / norm(eqt2d.grid.dim1) < 1E-3
    dim2 = range(eqt2d.grid.dim2[1], eqt2d.grid.dim2[end], length(eqt2d.grid.dim2))
    @assert norm(dim2 .- eqt2d.grid.dim2) / norm(eqt2d.grid.dim2) < 1E-3
    psi = range(eqt.profiles_1d.psi[1], eqt.profiles_1d.psi[end], length(eqt.profiles_1d.psi))
    @assert norm(psi .- eqt.profiles_1d.psi) / norm(eqt.profiles_1d.psi) < 1E-3

    return MXHEquilibrium.efit(
        MXHEquilibrium.cocos(11), # COCOS
        dim1, # Radius/R range
        dim2, # Elevation/Z range
        psi, # Polodial Flux range (polodial flux from magnetic axis)
        eqt2d.psi, # Polodial Flux on RZ grid (polodial flux from magnetic axis)
        eqt.profiles_1d.f, # Polodial Current
        eqt.profiles_1d.pressure, # Plasma pressure
        eqt.profiles_1d.q, # Q profile
        eqt.profiles_1d.psi .* 0, # Electric Potential
        (eqt.global_quantities.magnetic_axis.r, eqt.global_quantities.magnetic_axis.z), # Magnetic Axis (raxis,zaxis)
        Int(sign(eqt.profiles_1d.f[end]) * sign(eqt.global_quantities.ip)) # sign(dot(J,B))
    )
end

function _step(replay_actor::ActorReplay, actor::ActorEquilibrium, replay_dd::IMAS.dd)
    IMAS.copy_timeslice!(actor.dd.equilibrium, replay_dd.equilibrium, actor.dd.global_time)
    return replay_actor
end