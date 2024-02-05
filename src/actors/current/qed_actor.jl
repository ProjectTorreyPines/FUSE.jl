import QED
import FiniteElementHermite

#= ======== =#
#  ActorQED  #
#= ======== =#
Base.@kwdef mutable struct FUSEparameters__ActorQED{T} <: ParametersActor where {T<:Real}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    _time::Float64 = NaN
    Δt::Entry{Float64} = Entry{Float64}("s", "Evolve for Δt (Inf for steady state)"; default=Inf)
    Nt::Entry{Int} = Entry{Int}("-", "Number of time steps during evolution"; default=100)
    solve_for::Switch{Symbol} = Switch{Symbol}([:ip, :vloop], "-", "Solve for specified Ip or Vloop"; default=:ip)
    #== data flow parameters ==#
    ip_from::Switch{Symbol} = switch_get_from(:ip)
    vloop_from::Switch{Symbol} = switch_get_from(:vloop)
end

mutable struct ActorQED{D,P} <: PlasmaAbstractActor{D,P}
    dd::IMAS.dd{D}
    par::FUSEparameters__ActorQED{P}
    QO::Union{Nothing,QED.QED_state}
end

"""
    ActorQED(dd::IMAS.dd, act::ParametersAllActors; kw...)

Evolves the plasma current using the QED current diffusion solver.

!!! note

    This actor operates at "dd.global_time", any time advance must be done outside of the actor

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

    if par.Nt == 0
        # @info("QED: initialize")
        # initialize
        actor.QO = qed_init_from_imas(eqt, cp1d; uniform_rho=true)

    elseif par.Δt == Inf
        # @info("QED: steady state")
        # steady state solution
        actor.QO = qed_init_from_imas(eqt, cp1d; uniform_rho=false)

        if par.solve_for == :ip
            Ip = IMAS.get_from(dd, Val{:ip}, par.ip_from)
            Vedge = nothing
        else
            Ip = nothing
            Vedge = IMAS.get_from(dd, Val{:vloop}, par.vloop_from)
        end

        actor.QO = QED.steady_state(actor.QO, η_imas(dd.core_profiles.profiles_1d[]); Vedge, Ip)

    elseif par.Δt < 0.0
        # @info("QED: fake rampup for $(par.Δt) [s]")
        # scales the conductivity and non-inductive sources from zero to current value with some fudge law
        α_ni_evo_ip = 1.0 / 5.0 # early application of non-inductive current
        α_η_evo_ip = 5.0 # later transition to H-mode

        actor.QO = qed_init_from_imas(eqt, cp1d; uniform_rho=true)

        t0 = dd.global_time + par.Δt
        t1 = dd.global_time

        δt = abs(par.Δt / par.Nt)
        No = par.Nt
        Ni = 100

        rho = actor.QO.ι.x
        JBni1 = actor.QO.JBni.(rho)
        Ip1 = IMAS.get_from(dd, Val{:ip}, :pulse_schedule; time0=t1)
        η1 = η_imas(dd.core_profiles.profiles_1d[dd.global_time]).(rho)
        Vedge = nothing

        # start with flat resistivity and zero non-inductive sources
        η = QED.η_FE(rho, rho .* 0.0 .+ 1.0; use_log=true)
        actor.QO.JBni = QED.FE(rho, JBni1 .* 0.0)
        actor.QO = QED.steady_state(actor.QO, η; Vedge, Ip=1E3)

        #p = hline([1.0]; ls=:dash, color=:black, label="")
        for (k, tnow) in enumerate(range(t0 + δt / 2.0, t1 + δt / 2.0, No + 1)[1:end-1])
            Ip = IMAS.get_from(dd, Val{:ip}, :pulse_schedule; time0=tnow)
            evo_t = k / No
            evo_ip = Ip / Ip1
            actor.QO.JBni = QED.FE(rho, JBni1 .* evo_ip .^ α_ni_evo_ip) #
            η = QED.η_FE(rho, η1 ./ evo_ip .^ α_η_evo_ip; use_log=true)
            actor.QO = QED.diffuse(actor.QO, η, δt, Ni; Vedge, Ip)
        #    plot!(p, actor.QO; what=:q, label="", xlim=(0, 1), ylim=(0, 10))
        end
        #display(p)

    else
        # @info("QED: current diffusion for $(par.Δt) [s]")
        actor.QO = qed_init_from_imas(eqt, cp1d; uniform_rho=false)

        t0 = dd.global_time
        t1 = dd.global_time + par.Δt

        if false
            # staircase approach to track current ramps: one QED diffuse call for each time step
            δt = par.Δt / par.Nt
            No = par.Nt
            Ni = 1
        else
            δt = t1 - t0
            No = 1
            Ni = par.Nt
        end

        for tnow in range(t0 + δt / 2.0, t1 + δt / 2.0, No + 1)[1:end-1]
            if par.solve_for == :ip
                Ip = IMAS.get_from(dd, Val{:ip}, par.ip_from; time0=tnow)
                Vedge = nothing
            else
                Ip = nothing
                Vedge = IMAS.get_from(dd, Val{:vloop}, par.vloop_from; time0=tnow)
            end
            actor.QO = QED.diffuse(actor.QO, η_imas(dd.core_profiles.profiles_1d[tnow]), δt, Ni; Vedge, Ip, debug=false)
        end
    end

    return actor
end

function _finalize(actor::ActorQED)
    dd = actor.dd
    par = actor.par

    # set the total toroidal current in both equilibrium as well as core_profiles IDSs
    # NOTE: Here really we only care about core_profiles, since when the equilibrium actor is run,
    # then the new equilibrium time slice will be prepared based on the core_profiles current
    eqt = dd.equilibrium.time_slice[]
    ρ = eqt.profiles_1d.rho_tor_norm
    eqt.profiles_1d.q = 1.0 ./ actor.QO.ι.(ρ)
    eqt.profiles_1d.j_tor = actor.QO.JtoR.(ρ) ./ eqt.profiles_1d.gm9
    _, B0 = IMAS.vacuum_r0_b0(eqt)
    if true # we prefer using an expression to ensure consistency
        empty!(eqt.profiles_1d, :j_parallel) # restore expression
    else
        # more accurate
        eqt.profiles_1d.j_parallel = QED.JB(actor.QO; ρ) ./ B0
    end

    # update dd.core_profiles
    cp1d = dd.core_profiles.profiles_1d[]
    cp1d.j_total = QED.JB(actor.QO; ρ=cp1d.grid.rho_tor_norm) ./ B0

    if ismissing(cp1d, :j_non_inductive)
        cp1d.j_ohmic = cp1d.j_total
    else
        cp1d.j_ohmic = cp1d.j_total .- cp1d.j_non_inductive
    end
    empty!(cp1d, :j_tor) # restore expression

    # update dd.core_sources related to current
    IMAS.bootstrap_source!(dd)
    IMAS.ohmic_source!(dd)

    # add vloop info to pulse_schedule
    if par.solve_for == :ip
        vloop = IMAS.get_from(dd, Val{:vloop}, :core_profiles)
        @ddtime(dd.pulse_schedule.flux_control.loop_voltage.reference = vloop)
    end

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
function qed_init_from_imas(eqt::IMAS.equilibrium__time_slice, cp1d::IMAS.core_profiles__profiles_1d; uniform_rho::Bool=false)
    R0, B0 = IMAS.vacuum_r0_b0(eqt)

    rho_tor = eqt.profiles_1d.rho_tor
    gm1 = eqt.profiles_1d.gm1
    f = eqt.profiles_1d.f
    dvolume_drho_tor = eqt.profiles_1d.dvolume_drho_tor
    q = eqt.profiles_1d.q
    gm9 = eqt.profiles_1d.gm9

    # DO NOT use the equilibrium j_tor, since it's quality depends on the quality/resolution of the equilibrium solver
    # better to use the j_tor from core_profiles, which is the same quantity that is input in the equilibrium solver
    if false
        j_tor = eqt.profiles_1d.j_tor
    else
        j_tor = IMAS.interp1d(cp1d.grid.rho_tor_norm, cp1d.j_tor, :cubic).(IMAS.norm01(rho_tor))
    end

    y = log10.(1.0 ./ cp1d.conductivity_parallel) # `y` is used for packing points
    if ismissing(cp1d, :j_non_inductive)
        ρ_j_non_inductive = nothing
    else
        ρ_j_non_inductive = (cp1d.grid.rho_tor_norm, cp1d.j_non_inductive)
        y .*= ρ_j_non_inductive[2]
    end
    if uniform_rho
        ρ_grid = collect(range(0.0, 1.0, 2001))
    else
        ρ_grid = IMAS.pack_grid_gradients(cp1d.grid.rho_tor_norm, y; l=1E-2)
    end

    return QED.initialize(rho_tor, B0, gm1, f, dvolume_drho_tor, q, j_tor, gm9; ρ_j_non_inductive, ρ_grid)
end

function η_imas(cp1d::IMAS.core_profiles__profiles_1d; use_log::Bool=true)
    rho = cp1d.grid.rho_tor_norm
    η = 1.0 ./ cp1d.conductivity_parallel
    return QED.η_FE(rho, η; use_log)
end
