#= ================ =#
#  ActorEquilibrium  #
#= ================ =#
Base.@kwdef mutable struct FUSEparameters__ActorEquilibrium{T<:Real} <: ParametersActor{T}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    _time::Float64 = NaN
    #== actor parameters ==#
    model::Switch{Symbol} = Switch{Symbol}([:Solovev, :CHEASE, :TEQUILA, :none], "-", "Equilibrium actor to run"; default=:TEQUILA)
    symmetrize::Entry{Bool} = Entry{Bool}("-", "Force equilibrium up-down symmetry with respect to magnetic axis"; default=false)
    #== data flow parameters ==#
    j_p_from::Switch{Symbol} = Switch{Symbol}([:equilibrium, :core_profiles], "-", "Take j_tor and pressure profiles from this IDS"; default=:core_profiles)
    ip_from::Switch{Symbol} = switch_get_from(:ip)
    vacuum_r0_b0_from::Switch{Symbol} = switch_get_from(:vacuum_r0_b0)
    #== display and debugging parameters ==#
    do_plot::Entry{Bool} = act_common_parameters(; do_plot=false)
end

mutable struct ActorEquilibrium{D,P} <: CompoundAbstractActor{D,P}
    dd::IMAS.dd{D}
    par::FUSEparameters__ActorEquilibrium{P}
    act::ParametersAllActors
    eq_actor::Union{Nothing,ActorSolovev{D,P},ActorCHEASE{D,P},ActorTEQUILA{D,P},ActorNoOperation{D,P}}
end

"""
    ActorEquilibrium(dd::IMAS.dd, act::ParametersAllActors; kw...)

Provides a common interface to run multiple equilibrium actors
"""
function ActorEquilibrium(dd::IMAS.dd, act::ParametersAllActors; kw...)
    actor = ActorEquilibrium(dd, act.ActorEquilibrium, act; kw...)
    step(actor)
    finalize(actor)
    return actor
end

function ActorEquilibrium(dd::IMAS.dd, par::FUSEparameters__ActorEquilibrium, act::ParametersAllActors; kw...)
    logging_actor_init(ActorEquilibrium)
    par = par(kw...)
    if par.model == :Solovev
        eq_actor = ActorSolovev(dd, act.ActorSolovev; par.ip_from)
    elseif par.model == :CHEASE
        eq_actor = ActorCHEASE(dd, act.ActorCHEASE, act; par.ip_from)
    elseif par.model == :TEQUILA
        eq_actor = ActorTEQUILA(dd, act.ActorTEQUILA, act; par.ip_from)
    elseif par.model == :none
        eq_actor = ActorNoOperation(dd, act.ActorNoOperation)
    else
        error("ActorEquilibrium: model = `$(par.model)` can only be one of [:Solovev, :CHEASE, :TEQUILA, :none]")
    end
    return ActorEquilibrium(dd, par, act, eq_actor)
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

    if par.model !== :none
        eqt = dd.equilibrium.time_slice[]

        # symmetrize equilibrium if requested and number of X-points is even
        x_points = IMAS.x_points(dd.pulse_schedule.position_control.x_point)
        if par.symmetrize && mod(length(x_points), 2) != 1
            IMAS.symmetrize_equilibrium!(eqt)
        end

        # add flux surfaces information
        try
            fw = IMAS.first_wall(dd.wall)
            IMAS.flux_surfaces(eqt, fw.r, fw.z)
        catch e
            eqt2d = findfirst(:rectangular, eqt.profiles_2d)
            par.do_plot && display(current())
            contour(eqt2d.grid.dim1, eqt2d.grid.dim2, eqt2d.psi'; aspect_ratio=:equal)
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
  - Use j_tor,pressure from core_profiles or equilibrium
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
        rho_pol_norm0 = vcat(-reverse(sqrt.(cp1d.grid.psi_norm[index])), sqrt.(cp1d.grid.psi_norm[index]))
        j_tor0 = vcat(reverse(cp1d.j_tor[index]), cp1d.j_tor[index])
        pressure0 = vcat(reverse(cp1d.pressure[index]), cp1d.pressure[index])
        j_itp = IMAS.interp1d(rho_pol_norm0, j_tor0, :cubic)
        p_itp = IMAS.interp1d(rho_pol_norm0, pressure0, :cubic)
    elseif par.j_p_from == :equilibrium
        @assert !isempty(dd.equilibrium.time)
        eqt1d = dd.equilibrium.time_slice[].profiles_1d
        psi0 = eqt1d.psi
        rho_tor_norm0 = eqt1d.rho_tor_norm
        rho_pol_norm0 = sqrt.(eqt1d.psi_norm)
        j_tor0 = eqt1d.j_tor
        pressure0 = eqt1d.pressure
        j_itp = IMAS.interp1d(rho_pol_norm0, j_tor0, :cubic)
        p_itp = IMAS.interp1d(rho_pol_norm0, pressure0, :cubic)
    else
        @assert par.j_p_from in (:core_profiles, :equilibrium)
    end

    # get ip and b0 before wiping eqt in case ip_from=:equilibrium
    ip = IMAS.get_from(dd, Val{:ip}, actor.par.ip_from)
    r0, b0 = IMAS.get_from(dd, Val{:vacuum_r0_b0}, actor.par.vacuum_r0_b0_from)

    # add/clear time-slice
    eqt = resize!(dd.equilibrium.time_slice)
    resize!(eqt.profiles_2d, 1)
    eq1d = dd.equilibrium.time_slice[].profiles_1d

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
    if !isempty(getproperty(pc, :x_point, []))
        n = 0
        for k in eachindex(pc.x_point)
            rx = @ddtime(pc.x_point[k].r.reference)
            zx = @ddtime(pc.x_point[k].z.reference)
            if rx > 0.0 && !isnan(rx) && !isnan(zx)
                n += 1
                resize!(eqt.boundary.x_point, n)
                eqt.boundary.x_point[n].r = rx
                eqt.boundary.x_point[n].z = zx
            end
        end
    end

    # stike-points from position control
    if !isempty(getproperty(pc, :strike_point, []))
        n = 0
        for k in eachindex(pc.strike_point)
            rs = @ddtime(pc.strike_point[k].r.reference)
            zs = @ddtime(pc.strike_point[k].z.reference)
            if rs > 0.0 && !isnan(rs) && !isnan(zs)
                n += 1
                resize!(eqt.boundary.strike_point, n)
                eqt.boundary.strike_point[n].r = rs
                eqt.boundary.strike_point[n].z = zs
            end
        end
    end

    # set j_tor and pressure, forcing zero derivative on axis
    eq1d = dd.equilibrium.time_slice[].profiles_1d
    eq1d.psi = psi0
    eq1d.rho_tor_norm = rho_tor_norm0
    eq1d.j_tor = j_itp.(sqrt.(eq1d.psi_norm))
    eq1d.pressure = p_itp.(sqrt.(eq1d.psi_norm))

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
    @assert collect(dim1) ≈ eqt2d.grid.dim1
    dim2 = range(eqt2d.grid.dim2[1], eqt2d.grid.dim2[end], length(eqt2d.grid.dim2))
    @assert collect(dim2) ≈ eqt2d.grid.dim2
    psi = range(eqt.profiles_1d.psi[1], eqt.profiles_1d.psi[end], length(eqt.profiles_1d.psi))
    @assert collect(psi) ≈ eqt.profiles_1d.psi

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