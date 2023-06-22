#= ================ =#
#  ActorEquilibrium  #
#= ================ =#
Base.@kwdef mutable struct FUSEparameters__ActorEquilibrium{T} <: ParametersActor where {T<:Real}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    model::Switch{Symbol} = Switch{Symbol}([:Solovev, :CHEASE, :TEQUILA], "-", "Equilibrium actor to run"; default=:TEQUILA)
    symmetrize::Entry{Bool} = Entry{Bool}("-", "Force equilibrium up-down symmetry with respect to magnetic axis"; default=false)
end

mutable struct ActorEquilibrium{D,P} <: PlasmaAbstractActor
    dd::IMAS.dd{D}
    par::FUSEparameters__ActorEquilibrium{P}
    act::ParametersAllActors
    eq_actor::Union{Nothing,ActorSolovev{D,P},ActorCHEASE{D,P},ActorTEQUILA{D,P}}
    ip_from::Symbol
end

"""
    ActorEquilibrium(dd::IMAS.dd, act::ParametersAllActors; kw...)

Provides a common interface to run multiple equilibrium actors
"""
function ActorEquilibrium(dd::IMAS.dd, act::ParametersAllActors; ip_from::Symbol=:core_profiles, kw...)
    actor = ActorEquilibrium(dd, act.ActorEquilibrium, act, ip_from; kw...)
    step(actor)
    finalize(actor)
    return actor
end

function ActorEquilibrium(dd::IMAS.dd, par::FUSEparameters__ActorEquilibrium, act::ParametersAllActors, ip_from::Symbol; kw...)
    logging_actor_init(ActorEquilibrium)
    par = par(kw...)
    @show ip_from, "actEq", par.model
    if par.model == :Solovev
        eq_actor = ActorSolovev(dd, act.ActorSolovev; ip_from)
    elseif par.model == :CHEASE
        eq_actor = ActorCHEASE(dd, act.ActorCHEASE; ip_from)
    elseif par.model == :TEQUILA
        eq_actor = ActorTEQUILA(dd, act.ActorTEQUILA; ip_from)
    else
        error("ActorEquilibrium: model = `$(par.model)` can only be `:Solovev` or `:CHEASE`")
    end
    return ActorEquilibrium(dd, par, act, eq_actor, ip_from)
end

"""
    step(actor::ActorEquilibrium)

Runs through the selected equilibrium actor's step
"""
function _step(actor::ActorEquilibrium)
    step(actor.eq_actor)
    return actor
end

"""
    finalize(actor::ActorEquilibrium)

Finalizes the selected equilibrium actor
"""
function _finalize(actor::ActorEquilibrium)
    finalize(actor.eq_actor)
    if actor.par.symmetrize && mod(length(actor.dd.pulse_schedule.position_control.x_point), 2) != 1
        IMAS.symmetrize_equilibrium!(actor.dd.equilibrium.time_slice[])
        IMAS.flux_surfaces(actor.dd.equilibrium.time_slice[])
    end
    return actor
end

"""
    prepare_eq(dd::IMAS.dd)

Prepare `dd.equilbrium` to run equilibrium actors.
* clear equilibrium__time_slice
* set Ip, Bt, position control from pulse_schedule
* Copy pressure from core_profiles to equilibrium
* Copy j_tor from core_profiles to equilibrium

NOTE: prepare_eq(dd) must be called at the _step() stage
of all equilibrium actors to ensure that they work properly
when used in a transport-equilibrium loop.
"""
function prepare_eq(dd::IMAS.dd, ip_from::Symbol)
    ps = dd.pulse_schedule
    pc = ps.position_control

    # freeze cp1d before wiping eqt
    cp1d = IMAS.freeze(dd.core_profiles.profiles_1d[])

    # add/clear time-slice
    eqt = resize!(dd.equilibrium.time_slice)
    resize!(eqt.profiles_2d, 1)
    eq1d = dd.equilibrium.time_slice[].profiles_1d

    # scalar quantities
    eqt.global_quantities.ip = IMAS.get_from(dd, :ip, ip_from)
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
    resize!(eqt.boundary.x_point, length(pc.x_point))
    for k in eachindex(pc.x_point)
        eqt.boundary.x_point[k].r = @ddtime(pc.x_point[k].r.reference.data)
        eqt.boundary.x_point[k].z = @ddtime(pc.x_point[k].z.reference.data)
    end

    # set j_tor and pressure
    eq1d = dd.equilibrium.time_slice[].profiles_1d
    eq1d.psi = cp1d.grid.psi
    index = cp1d.grid.psi_norm .> 0.05 # force zero derivative current/pressure on axis
    rho_pol_norm0 = vcat(-reverse(sqrt.(cp1d.grid.psi_norm[index])), sqrt.(cp1d.grid.psi_norm[index]))
    j_tor0 = vcat(reverse(cp1d.j_tor[index]), cp1d.j_tor[index])
    pressure0 = vcat(reverse(cp1d.pressure[index]), cp1d.pressure[index])
    eq1d.j_tor = IMAS.interp1d(rho_pol_norm0, j_tor0, :cubic).(sqrt.(eq1d.psi_norm))
    eq1d.pressure = IMAS.interp1d(rho_pol_norm0, pressure0, :cubic).(sqrt.(eq1d.psi_norm))

    return dd
end

"""
    IMAS2Equilibrium(eqt::IMAS.equilibrium__time_slice)

Convert IMAS.equilibrium__time_slice to MXHEquilibrium.jl EFIT structure
"""
function IMAS2Equilibrium(eqt::IMAS.equilibrium__time_slice)
    dim1 = range(eqt.profiles_2d[1].grid.dim1[1], eqt.profiles_2d[1].grid.dim1[end], length=length(eqt.profiles_2d[1].grid.dim1))
    @assert collect(dim1) ≈ eqt.profiles_2d[1].grid.dim1
    dim2 = range(eqt.profiles_2d[1].grid.dim2[1], eqt.profiles_2d[1].grid.dim2[end], length=length(eqt.profiles_2d[1].grid.dim2))
    @assert collect(dim2) ≈ eqt.profiles_2d[1].grid.dim2
    psi = range(eqt.profiles_1d.psi[1], eqt.profiles_1d.psi[end], length=length(eqt.profiles_1d.psi))
    @assert collect(psi) ≈ eqt.profiles_1d.psi

    MXHEquilibrium.efit(
        MXHEquilibrium.cocos(11), # COCOS
        dim1, # Radius/R range
        dim2, # Elevation/Z range
        psi, # Polodial Flux range (polodial flux from magnetic axis)
        eqt.profiles_2d[1].psi, # Polodial Flux on RZ grid (polodial flux from magnetic axis)
        eqt.profiles_1d.f, # Polodial Current
        eqt.profiles_1d.pressure, # Plasma pressure
        eqt.profiles_1d.q, # Q profile
        eqt.profiles_1d.psi .* 0, # Electric Potential
        (eqt.global_quantities.magnetic_axis.r, eqt.global_quantities.magnetic_axis.z), # Magnetic Axis (raxis,zaxis)
        Int(sign(eqt.profiles_1d.f[end]) * sign(eqt.global_quantities.ip)) # sign(dot(J,B))
    )
end