#= ================ =#
#  ActorEquilibrium  #
#= ================ =#
Base.@kwdef mutable struct FUSEparameters__ActorEquilibrium{T} <: ParametersActor where {T<:Real}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    #== actor parameters ==#
    model::Switch{Symbol} = Switch{Symbol}([:Solovev, :CHEASE, :TEQUILA], "-", "Equilibrium actor to run"; default=:TEQUILA)
    symmetrize::Entry{Bool} = Entry{Bool}("-", "Force equilibrium up-down symmetry with respect to magnetic axis"; default=false)
    #== data flow parameters ==#
    ip_from::Switch{Symbol} = switch_get_from(:ip)
    #== display and debugging parameters ==#
    do_plot::Entry{Bool} = Entry{Bool}("-", "Plot before and after actor"; default=false)
end

mutable struct ActorEquilibrium{D,P} <: PlasmaAbstractActor{D,P}
    dd::IMAS.dd{D}
    par::FUSEparameters__ActorEquilibrium{P}
    act::ParametersAllActors
    eq_actor::Union{Nothing,ActorSolovev{D,P},ActorCHEASE{D,P},ActorTEQUILA{D,P}}
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
        eq_actor = ActorCHEASE(dd, act.ActorCHEASE; par.ip_from)
    elseif par.model == :TEQUILA
        eq_actor = ActorTEQUILA(dd, act.ActorTEQUILA; par.ip_from)
    else
        error("ActorEquilibrium: model = `$(par.model)` can only be `:Solovev` or `:CHEASE`")
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
            plot(dd.equilibrium; cx=true, label="before ActorEquilibrium")
        else
            plot()
        end
    end

    # initialize eqt for equilibrium actors
    prepare(actor)

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

    # symmetrize equilibrium if requested and number of X-points is even
    x_points = IMAS.x_points(dd.pulse_schedule.position_control.x_point)
    if par.symmetrize && mod(length(x_points), 2) != 1
        IMAS.symmetrize_equilibrium!(dd.equilibrium.time_slice[])
    end

    # add flux surfaces information
    IMAS.flux_surfaces(dd.equilibrium.time_slice[])

    if par.do_plot
        plot!(dd.equilibrium; label="after ActorEquilibrium")
    end

    return actor
end

"""
    prepare(actor::ActorEquilibrium)

Prepare `dd.equilibrium` to run equilibrium actors

  - clear equilibrium__time_slice
  - set Ip, Bt, position control from pulse_schedule
  - Copy pressure from core_profiles to equilibrium
  - Copy j_tor from core_profiles to equilibrium
"""
function prepare(actor::ActorEquilibrium)
    dd = actor.dd
    ps = dd.pulse_schedule
    pc = ps.position_control

    # freeze core_profiles before wiping eqt and get ip
    ip = IMAS.get_from(dd, Val{:ip}, actor.par.ip_from)
    cp1d = dd.core_profiles.profiles_1d[]
    psi = cp1d.grid.psi
    if true
        index = cp1d.grid.psi_norm .> 0.05
        rho_pol_norm0 = vcat(-reverse(sqrt.(cp1d.grid.psi_norm[index])), sqrt.(cp1d.grid.psi_norm[index]))
        j_tor0 = vcat(reverse(cp1d.j_tor[index]), cp1d.j_tor[index])
        pressure0 = vcat(reverse(cp1d.pressure[index]), cp1d.pressure[index])
    else
        rho_pol_norm0 = sqrt.(cp1d.grid.psi_norm)
        j_tor0 = cp1d.j_tor
        pressure0 = cp1d.pressure
    end

    # add/clear time-slice
    eqt = resize!(dd.equilibrium.time_slice)
    resize!(eqt.profiles_2d, 1)
    eq1d = dd.equilibrium.time_slice[].profiles_1d

    # scalar quantities
    eqt.global_quantities.ip = ip

    R0 = dd.equilibrium.vacuum_toroidal_field.r0
    B0 = @ddtime(ps.tf.b_field_tor_vacuum_r.reference.data) / R0
    @ddtime(dd.equilibrium.vacuum_toroidal_field.b0 = B0)

    # position control
    eqt.boundary.minor_radius = @ddtime(pc.minor_radius.reference.data)
    eqt.boundary.geometric_axis.r = @ddtime(pc.geometric_axis.r.reference.data)
    eqt.boundary.geometric_axis.z = @ddtime(pc.geometric_axis.z.reference.data)
    eqt.boundary.elongation = @ddtime(pc.elongation.reference.data)
    eqt.boundary.triangularity = @ddtime(pc.triangularity.reference.data)
    eqt.boundary.squareness = @ddtime(pc.squareness.reference.data)

    # boundary
    eqt.boundary.outline.r, eqt.boundary.outline.z = IMAS.boundary(pc)

    # x-points
    if length(getproperty(pc, :x_point, [])) >= 1
        n = 0
        for k in eachindex(pc.x_point)
            rx = @ddtime(pc.x_point[k].r.reference.data)
            zx = @ddtime(pc.x_point[k].z.reference.data)
            if rx > 0.0 && !isnan(rx) && !isnan(zx)
                n += 1
                resize!(eqt.boundary.x_point, n)
                eqt.boundary.x_point[n].r = rx
                eqt.boundary.x_point[n].z = zx
            end
        end
    end

    # set j_tor and pressure, forcing zero derivative on axis
    eq1d = dd.equilibrium.time_slice[].profiles_1d
    eq1d.psi = psi
    eq1d.j_tor = IMAS.interp1d(rho_pol_norm0, j_tor0, :cubic).(sqrt.(eq1d.psi_norm))
    eq1d.pressure = IMAS.interp1d(rho_pol_norm0, pressure0, :cubic).(sqrt.(eq1d.psi_norm))

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
    dim1 = LinRange(eqt2d.grid.dim1[1], eqt2d.grid.dim1[end], length(eqt2d.grid.dim1))
    @assert collect(dim1) ≈ eqt2d.grid.dim1
    dim2 = LinRange(eqt2d.grid.dim2[1], eqt2d.grid.dim2[end], length(eqt2d.grid.dim2))
    @assert collect(dim2) ≈ eqt2d.grid.dim2
    psi = LinRange(eqt.profiles_1d.psi[1], eqt.profiles_1d.psi[end], length(eqt.profiles_1d.psi))
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