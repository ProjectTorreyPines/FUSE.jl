import QED

#= ======== =#
#  ActorQED  #
#= ======== =#
Base.@kwdef mutable struct FUSEparameters__ActorQED{T<:Real} <: ParametersActor{T}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    _time::Float64 = NaN
    Δt::Entry{Float64} = Entry{Float64}("s", "Evolve for Δt (Inf for steady state)"; default=Inf, check=x->@assert x>=0 "Δt must be >= 0.0")
    Nt::Entry{Int} = Entry{Int}("-", "Number of time steps during evolution"; default=100)
    solve_for::Switch{Symbol} = Switch{Symbol}([:ip, :vloop], "-", "Solve for specified Ip or Vloop"; default=:ip)
    allow_floating_plasma_current::Entry{Bool} = Entry{Bool}("-", "Zero loop voltage if non-inductive fraction exceeds 100% of the target Ip")
    #== data flow parameters ==#
    ip_from::Switch{Symbol} = switch_get_from(:ip)
    vloop_from::Switch{Symbol} = switch_get_from(:vloop)
end

mutable struct ActorQED{D,P} <: SingleAbstractActor{D,P}
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

    # non_inductive contribution
    Ip_non_inductive = trapz(cp1d.grid.area, cp1d.j_non_inductive)
    B0 = eqt.global_quantities.vacuum_toroidal_field.b0
    JBni = QED.FE(cp1d.grid.rho_tor_norm, cp1d.j_non_inductive .* B0)

    # initialize QED
    if actor.QO === nothing || par.Δt == Inf
        actor.QO = qed_init_from_imas(eqt, cp1d; uniform_rho = 501)
    else
        actor.QO.JBni = JBni
    end

    if par.Nt == 0
        # only initialize, nothing to do

    elseif par.Δt != Inf
        # current diffusion
        t0 = dd.global_time
        t1 = t0 + par.Δt

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

        for time0 in range(t0 + δt / 2.0, t1 + δt / 2.0, No + 1)[1:end-1]
            if par.solve_for == :ip
                Ip = IMAS.get_from(dd, Val{:ip}, par.ip_from; time0)
                Vedge = nothing
                if par.allow_floating_plasma_current && abs(Ip) < abs(Ip_non_inductive)
                    Ip = nothing
                    Vedge = 0.0
                end
            else
                Ip = nothing
                Vedge = IMAS.get_from(dd, Val{:vloop}, par.vloop_from; time0)
            end
            actor.QO = QED.diffuse(actor.QO, η_imas(dd.core_profiles.profiles_1d[time0]), δt, Ni; Vedge, Ip, debug=false)
        end

    elseif par.Δt == Inf
        # steady state solution
        if par.solve_for == :ip
            Ip = IMAS.get_from(dd, Val{:ip}, par.ip_from)
            Vedge = nothing
            if par.allow_floating_plasma_current && abs(Ip) < abs(Ip_non_inductive)
                Ip = nothing
                Vedge = 0.0
            end
        else
            Ip = nothing
            Vedge = IMAS.get_from(dd, Val{:vloop}, par.vloop_from)
        end

        actor.QO = QED.steady_state(actor.QO, η_imas(dd.core_profiles.profiles_1d[]); Vedge, Ip)
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
    _, B0 = eqt.global_quantities.vacuum_toroidal_field.r0, eqt.global_quantities.vacuum_toroidal_field.b0
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

    return actor
end

# utils
"""
    qed_init_from_imas(eqt::IMAS.equilibrium__time_slice, cp1d::IMAS.core_profiles__profiles_1d; uniform_rho::Int)

Setup QED from data in IMAS `dd`

NOTE: QED is initalized from equilibrium and not core_profiles because
it needs both `q` and `j_tor`, and equilibrium is the only place where
the two ought to be self-consistent
"""
function qed_init_from_imas(eqt::IMAS.equilibrium__time_slice, cp1d::IMAS.core_profiles__profiles_1d; uniform_rho::Int)
    R0, B0 = eqt.global_quantities.vacuum_toroidal_field.r0, eqt.global_quantities.vacuum_toroidal_field.b0

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

    # uniform_rho just works better, and QED is fast enough that it can handle many radial points
    if uniform_rho > 0
        ρ_grid = collect(range(0.0, 1.0, uniform_rho))
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
