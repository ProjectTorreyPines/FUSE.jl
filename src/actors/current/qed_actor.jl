import QED
import FiniteElementHermite

#= ======== =#
#  ActorQED  #
#= ======== =#
Base.@kwdef mutable struct FUSEparameters__ActorQED{T} <: ParametersActor where {T<:Real}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    Δt::Entry{Float64} = Entry{Float64}("s", "Evolve for Δt")
    Nt::Entry{Int} = Entry{Int}("-", "Number of time steps during evolution"; default=100)
    solve_for::Switch{Symbol} = Switch{Symbol}([:ip, :vloop], "-", "Solve for specified Ip or Vloop"; default=:ip)
end

mutable struct ActorQED{D,P} <: PlasmaAbstractActor
    dd::IMAS.dd{D}
    par::FUSEparameters__ActorQED{P}
    QO::Union{Nothing,QED.QED_state}
end

"""
    ActorQED(dd::IMAS.dd, act::ParametersAllActors; kw...)

Evolves the plasma current using the QED current diffusion solver.

!!! note
    This actor wants the driver to do, before calling the actor:

       IMAS.new_timeslice!(dd, dd.global_time + Δt)
       dd.global_time += Δt
       ActorQED(dd, act)

!!! note 
    Stores data in `dd.equilibrium`, `dd.core_profiles`, `dd.core_sources`
"""
function ActorQED(dd::IMAS.dd, act::ParametersAllActors; kw...)
    actor = ActorQED(dd, act.ActorQED; kw...)
    step(actor)
    finalize(actor)
    return actor
end

function ActorQED(dd::IMAS.dd, par::FUSEparameters__ActorQED; kw...)
    logging_actor_init(ActorQED)
    par = par(kw...)
    return ActorQED(dd, par, nothing)
end

function _step(actor::ActorQED)
    dd = actor.dd
    par = actor.par

    eqt = dd.equilibrium.time_slice[]
    cp1d = dd.core_profiles.profiles_1d[]

    t0 = dd.global_time
    t1 = dd.global_time + par.Δt
    δt = par.Δt / par.Nt

    # initialization
    actor.QO = qed_init_from_imas(eqt, cp1d)

    if false
        # staircase approach to track current ramps: one QED diffuse call for each time step
        for tnow in LinRange(t0, t1, par.Nt + 1)[2:end]
            if par.solve_for == :ip
                Ip = IMAS.get_time_array(dd.pulse_schedule.flux_control.i_plasma.reference, :data, tnow, :linear)
                Vedge = nothing
            else
                error("Vloop advance not supported")
            end
            actor.QO = QED.diffuse(actor.QO, η_imas(dd.core_profiles.profiles_1d[tnow]), δt, 1; Vedge, Ip)
        end
    else
        if par.solve_for == :ip
            Ip = IMAS.get_time_array(dd.pulse_schedule.flux_control.i_plasma.reference, :data, t1, :linear)
            Vedge = nothing
        else
            error("Vloop advance not supported")
        end
        actor.QO = QED.diffuse(actor.QO, η_imas(dd.core_profiles.profiles_1d[t1]), par.Δt, par.Nt; Vedge, Ip)
    end

    return actor
end

function _finalize(actor::ActorQED)
    dd = actor.dd

    # set the total toroidal current in both equilibrium as well as core_profiles IDSs
    # NOTE: Here really we only care about core_profiles, since when the equilibrium actor is run,
    # then the new equilibrium time slice will be prepared based on the core_profiles current
    eqt = dd.equilibrium.time_slice[]
    ρ = eqt.profiles_1d.rho_tor_norm
    eqt.profiles_1d.q = 1.0 ./ actor.QO.ι.(ρ)
    eqt.profiles_1d.j_tor = actor.QO.JtoR.(ρ) ./ eqt.profiles_1d.gm9

    # update dd.core_profiles
    cp1d = dd.core_profiles.profiles_1d[]
    IMAS.j_total_from_equilibrium!(eqt, cp1d)

    # update dd.core_sources related to current
    IMAS.bootstrap_source!(dd)
    IMAS.ohmic_source!(dd)

    return actor
end

# utils
"""
    qed_init_from_imas(eqt::IMAS.equilibrium__time_slice, cp1d::IMAS.core_profiles__profiles_1d)

Setup QED from data in IMAS `dd`

NOTE: QED is initalized from equilibrium and not core_profiles because
it needs both `q` and `j_tor`, and equilibrium is the only place where
the two ought to be self-consistent
"""
function qed_init_from_imas(eqt::IMAS.equilibrium__time_slice, cp1d::IMAS.core_profiles__profiles_1d)
    R0, B0 = IMAS.vacuum_r0_b0(eqt)

    rho_tor = eqt.profiles_1d.rho_tor
    gm1 = eqt.profiles_1d.gm1
    f = eqt.profiles_1d.f
    dvolume_drho_tor = eqt.profiles_1d.dvolume_drho_tor
    q = eqt.profiles_1d.q
    j_tor = eqt.profiles_1d.j_tor
    gm9 = eqt.profiles_1d.gm9

    if ismissing(cp1d, :j_non_inductive)
        ρ_j_non_inductive = nothing
    else
        ρ_j_non_inductive = (cp1d.grid.rho_tor_norm, cp1d.j_non_inductive)
    end

    return QED.initialize(rho_tor, B0, gm1, f, dvolume_drho_tor, q, j_tor, gm9; ρ_j_non_inductive)
end

function η_imas(cp1d::IMAS.core_profiles__profiles_1d; use_log::Bool=true)
    rho = cp1d.grid.rho_tor_norm
    η = 1.0 ./ cp1d.conductivity_parallel
    return QED.η_FE(rho, η; use_log)
end
