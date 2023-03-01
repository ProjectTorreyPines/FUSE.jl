import QED

#= =============== =#
#  ActorQEDcurrent  #
#= =============== =#
Base.@kwdef mutable struct FUSEparameters__ActorQEDcurrent{T} <: ParametersActor where {T<:Real}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
end

mutable struct ActorQEDcurrent <: PlasmaAbstractActor
    dd::IMAS.dd
    par::FUSEparameters__ActorQEDcurrent
    QI::QED.QED_state
    η#::Base.Callable
    QO::Union{QED.QED_state,Missing}
    initial_time
    tmax
end

"""
    ActorQEDcurrent(dd::IMAS.dd, act::ParametersAllActors; kw...)

Evolves the plasma current using the QED current diffusion solver

!!! note 
    Stores data in `dd.equilibrium`
"""
function ActorQEDcurrent(dd::IMAS.dd, act::ParametersAllActors; kw...)
    par = act.ActorQEDcurrent(kw...)
    actor = ActorQEDcurrent(dd, par)
    step(actor)
    finalize(actor)
    return actor
end

function ActorQEDcurrent(dd::IMAS.dd, par::FUSEparameters__ActorQEDcurrent; kw...)
    logging_actor_init(ActorQEDcurrent)
    par = par(kw...)
    return ActorQEDcurrent(dd, par, from_imas(dd), η_imas(dd), missing, @ddtime(dd.equilibrium.time), dd.global_time)
end

function _step(actor::ActorQEDcurrent, tmax, Nt, Vedge=nothing, Ip=nothing; resume=false)
    if resume
        if actor.QO !== missing
            actor.QI = actor.QO
            actor.QO = missing
        end
        actor.tmax += tmax
    else
        actor.tmax = tmax
    end
    actor.QO = QED.diffuse(actor.QI, actor.η, tmax, Nt; Vedge, Ip)
end

function _finalize(actor::ActorQEDcurrent)
    dd = actor.dd

    eqt = dd.equilibrium.time_slice[]
    newtime = actor.initial_time + actor.tmax
    resize!(dd.equilibrium.time_slice, Float64(newtime))
    dd.equilibrium.time_slice[newtime] = eqt_new = deepcopy(eqt)

    dΡ_dρ = eqt_new.profiles_1d.rho_tor[end]
    ρ = eqt_new.profiles_1d.rho_tor / dΡ_dρ

    eqt_new.profiles_1d.q = 1.0 ./ actor.QO.ι.(ρ)
    eqt_new.profiles_1d.j_tor = actor.QO.JtoR.(ρ) ./ eqt_new.profiles_1d.gm9
    eqt_new.time = newtime

    # if dd has core_profiles set cp1d.j_total from equilibrium (and set j_ohmic as expression)
    if !isempty(dd.core_profiles.profiles_1d)
        cp1d_new = cp1d = dd.core_profiles.profiles_1d[]
        if cp1d.time < newtime
            resize!(dd.core_profiles.profiles_1d, Float64(newtime))
            dd.core_profiles.profiles_1d[newtime] = cp1d_new = deepcopy(cp1d)
        end
        IMAS.j_total_from_equilibrium!(eqt_new, cp1d_new)
    end
    # update core_sources related to current
    IMAS.bootstrap_source!(dd)
    IMAS.ohmic_source!(dd)
    return dd
end

# utils
function from_imas(dd::IMAS.dd)
    eqt = dd.equilibrium.time_slice[]
    rho_tor = eqt.profiles_1d.rho_tor
    B0 = @ddtime(dd.equilibrium.vacuum_toroidal_field.b0)
    gm1 = eqt.profiles_1d.gm1
    f = eqt.profiles_1d.f
    dvolume_drho_tor = eqt.profiles_1d.dvolume_drho_tor
    q = eqt.profiles_1d.q
    j_tor = eqt.profiles_1d.j_tor
    gm9 = eqt.profiles_1d.gm9

    if !isempty(dd.core_profiles.profiles_1d) && !ismissing(dd.core_profiles.profiles_1d[], :j_non_inductive)
        prof1d = dd.core_profiles.profiles_1d[]
        ρ_j_non_inductive = (prof1d.grid.rho_tor_norm, prof1d.j_non_inductive)
    else
        ρ_j_non_inductive = nothing
    end

    return QED.initialize(rho_tor, B0, gm1, f, dvolume_drho_tor, q, j_tor, gm9; ρ_j_non_inductive)
end

function η_imas(dd::IMAS.dd; use_log=true)
    prof1d = dd.core_profiles.profiles_1d[]
    rho = prof1d.grid.rho_tor_norm
    η = 1.0 ./ prof1d.conductivity_parallel
    return QED.η_FE(rho, η; use_log)
end

#= ======================= =#
#  ActorSteadyStateCurrent  #
#= ======================= =#
Base.@kwdef mutable struct FUSEparameters__ActorSteadyStateCurrent{T} <: ParametersActor where {T<:Real}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
end

mutable struct ActorSteadyStateCurrent <: PlasmaAbstractActor
    dd::IMAS.dd
    par::FUSEparameters__ActorSteadyStateCurrent
    function ActorSteadyStateCurrent(dd::IMAS.dd, par::FUSEparameters__ActorSteadyStateCurrent; kw...)
        logging_actor_init(ActorSteadyStateCurrent)
        par = par(kw...)
        return new(dd, par)
    end
end

"""
    ActorSteadyStateCurrent(dd::IMAS.dd, act::ParametersAllActors; kw...)

Evolves the current to steady state using the conductivity from `dd.core_profiles` and current profile form `dd.equilibrium`.

Also sets the ohmic, bootstrap and non-inductive current profiles in `dd.core_profiles`

!!! note 
    Stores data in `dd.core_profiles, dd.equilbrium`
"""
function ActorSteadyStateCurrent(dd::IMAS.dd, act::ParametersAllActors; kw...)
    par = act.ActorSteadyStateCurrent(kw...)
    actor = ActorSteadyStateCurrent(dd, par)
    step(actor)
    finalize(actor)
    return actor
end

function _step(actor::ActorSteadyStateCurrent)
    dd = actor.dd
    IMAS.j_ohmic_steady_state!(dd.equilibrium.time_slice[], dd.core_profiles.profiles_1d[])
    # update core_sources related to current
    IMAS.bootstrap_source!(dd)
    IMAS.ohmic_source!(dd)
    return actor
end
