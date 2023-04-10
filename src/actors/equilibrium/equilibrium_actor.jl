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
    eq_actor::Union{ActorSolovev,ActorCHEASE}
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
    if par.model == :Solovev
        eq_actor = ActorSolovev(dd, act.ActorSolovev)
    elseif par.model == :CHEASE
        eq_actor = ActorCHEASE(dd, act.ActorCHEASE)
    else
        error("ActorEquilibrium: model = $(par.model) is unknown")
    end
    return ActorEquilibrium(dd, par, eq_actor)
end

"""
    prepare(dd::IMAS.dd, :ActorEquilibrium, act::ParametersAllActors; kw...)

Prepare dd to run ActorEquilibrium
* call prapare function of the different equilibrium models
"""
function prepare(dd::IMAS.dd, ::Type{Val{:ActorEquilibrium}}, act::ParametersAllActors; kw...)
    par = act.ActorEquilibrium(kw...)
    if par.model == :Solovev
        return prepare(dd, :ActorSolovev, act; kw...)
    elseif par.model == :CHEASE
        return prepare(dd, :ActorCHEASE, act; kw...)
    else
        error("ActorEquilibrium: model = $(par.model) is unknown")
    end
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
    if actor.par.symmetrize
        IMAS.symmetrize_equilibrium!(actor.dd.equilibrium.time_slice[])
        IMAS.flux_surfaces(actor.dd.equilibrium.time_slice[])
    end
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