#= ================ =#
#  ActorEquilibrium  #
#= ================ =#
Base.@kwdef mutable struct FUSEparameters__ActorEquilibrium{T} <: ParametersActor where {T<:Real}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    model::Switch{Symbol} = Switch(Symbol, [:Solovev, :CHEASE], "-", "Equilibrium actor to run"; default=:Solovev)
    symmetrize::Entry{Bool} = Entry(Bool, "-", "Force equilibrium up-down symmetry with respect to magnetic axis"; default=false)
end

mutable struct ActorEquilibrium <: PlasmaAbstractActor
    dd::IMAS.dd
    par::FUSEparameters__ActorEquilibrium
    act::ParametersAllActors
    eq_actor::Union{Nothing,ActorSolovev,ActorCHEASE}
end

"""
    ActorEquilibrium(dd::IMAS.dd, act::ParametersAllActors; kw...)

Provides a common interface to run multiple equilibrium actors
"""
function ActorEquilibrium(dd::IMAS.dd, act::ParametersAllActors; kw...)
    par = act.ActorEquilibrium
    actor = ActorEquilibrium(dd, par, act; kw...)
    step(actor)
    finalize(actor)
    return actor
end

function ActorEquilibrium(dd::IMAS.dd, par::FUSEparameters__ActorEquilibrium, act::ParametersAllActors; kw...)
    logging_actor_init(ActorEquilibrium)
    par = par(kw...)
    return ActorEquilibrium(dd, par, act, nothing)
end

"""
    step(actor::ActorEquilibrium)

Runs through the selected equilibrium actor's step
"""
function _step(actor::ActorEquilibrium)
    dd = actor.dd
    par = actor.par
    act = actor.act
    
    prepare(actor)

    if par.model == :Solovev
        actor.eq_actor = ActorSolovev(dd, act.ActorSolovev)
    elseif par.model == :CHEASE
        actor.eq_actor = ActorCHEASE(dd, act.ActorCHEASE)
    else
        error("ActorEquilibrium: model = `$(par.model)` can only be :Solovev or :CHEASE")
    end
    
    step(actor.eq_actor)
    
    return actor
end

"""
    finalize(actor::ActorEquilibrium)

Finalizes the selected equilibrium actor
"""
function _finalize(actor::ActorEquilibrium)
    finalize(actor.eq_actor)
    if actor.par.symmetrize
        IMAS.symmetrize_equilibrium!(actor.dd.equilibrium.time_slice[])
        IMAS.flux_surfaces(actor.dd.equilibrium.time_slice[])
    end
    return actor
end

"""
    prepare(actor::ActorEquilibrium)

Prepare `dd.equilbrium` to run ActorEquilibrium
* clear equilibrium__time_slice
* set Ip, Bt, position control from pulse_schedule
* Copy pressure from core_profiles to equilibrium
* Copy j_parallel from core_profiles to equilibrium
"""
function prepare(actor::ActorEquilibrium)
    dd = actor.dd
    ps = dd.pulse_schedule
    pc = ps.position_control

    # pressure and current profiles
    # NOTE: this is done this way because j_tor is an expression that depends on equilibrium (via j_bootstrap)
    #       so we cannot add/clear time-slice before getting j_tor (and pressure)
    if isempty(dd.equilibrium.time_slice)
        resize!(dd.equilibrium.time_slice)
    end
    eq1d = dd.equilibrium.time_slice[].profiles_1d
    cp1d = dd.core_profiles.profiles_1d[]
    if ismissing(eq1d, :psi)
        eq1d.psi = cp1d.grid.psi
    end
    psi = eq1d.psi
    j_tor = IMAS.interp1d(cp1d.grid.psi_norm, cp1d.j_tor).(eq1d.psi_norm)
    pressure = IMAS.interp1d(cp1d.grid.psi_norm, cp1d.pressure).(eq1d.psi_norm)

    # add/clear time-slice
    eqt = resize!(dd.equilibrium.time_slice)
    eq1d = dd.equilibrium.time_slice[].profiles_1d

    # scalar quantities
    eqt.global_quantities.ip = @ddtime(ps.flux_control.i_plasma.reference.data)
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
    eqt.boundary.outline.r = [@ddtime(pcb.r.reference.data) for pcb in pc.boundary_outline]
    eqt.boundary.outline.z = [@ddtime(pcb.z.reference.data) for pcb in pc.boundary_outline]

    # set j_tor and pressure
    eq1d.psi = psi
    eq1d.j_tor = j_tor
    eq1d.pressure = pressure

    return actor
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