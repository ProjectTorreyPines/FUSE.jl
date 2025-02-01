import QED

#= ======== =#
#  ActorQED  #
#= ======== =#
Base.@kwdef mutable struct FUSEparameters__ActorQED{T<:Real} <: ParametersActor{T}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    _time::Float64 = NaN
    Δt::Entry{Float64} = Entry{Float64}("s", "Evolve for Δt (Inf for steady state)"; default=Inf, check=x -> @assert x >= 0 "Δt must be >= 0.0")
    Nt::Entry{Int} = Entry{Int}("-", "Number of time steps during evolution"; default=100, check=x -> @assert x > 0 "Nt must be > 0")
    solve_for::Switch{Symbol} = Switch{Symbol}([:ip, :vloop], "-", "Solve for specified Ip or Vloop"; default=:ip)
    allow_floating_plasma_current::Entry{Bool} = Entry{Bool}("-", "Zero loop voltage if non-inductive fraction exceeds 100% of the target Ip")
    qmin_desired::Entry{Float64} = Entry{Float64}("-", "Keep the minimum magnitude of the q-profile above this value"; default=1.0, check=x -> @assert x >= 0 "qmin_desired >= 0")
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

    Stores data in `dd.core_profiles.profiles_1d[].j_ohmic`
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
    B0 = eqt.global_quantities.vacuum_toroidal_field.b0
    JBni = QED.FE(cp1d.grid.rho_tor_norm, cp1d.j_non_inductive .* B0)
    if par.solve_for == :ip && par.allow_floating_plasma_current
        ip_non_inductive = IMAS.Ip_non_inductive(cp1d, eqt)
    end

    # initialize QED
    if actor.QO === nothing || par.Δt == Inf
        actor.QO = qed_init_from_imas(eqt, cp1d; uniform_rho=501)
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
                if par.allow_floating_plasma_current && abs(Ip) < abs(ip_non_inductive)
                    Ip = nothing
                    Vedge = 0.0
                end
            else
                Ip = nothing
                Vedge = IMAS.get_from(dd, Val{:vloop}, par.vloop_from; time0)
            end

            # check where q<1 based on the the q-profile at the previous
            # time-step and change the resisitivity to keep q>1
            qval = 1.0 ./ abs.(actor.QO.ι.(cp1d.grid.rho_tor_norm))
            i_qdes = findlast(qval .< par.qmin_desired)
            if i_qdes === nothing
                i_qdes = 0
            end

            actor.QO = QED.diffuse(actor.QO, η_imas(dd.core_profiles.profiles_1d[time0], i_qdes), δt, Ni; Vedge, Ip, debug=false)
        end

    elseif par.Δt == Inf
        # steady state solution
        if par.solve_for == :ip
            Ip = IMAS.get_from(dd, Val{:ip}, par.ip_from)
            Vedge = nothing
            if par.allow_floating_plasma_current && abs(Ip) < abs(ip_non_inductive)
                Ip = nothing
                Vedge = 0.0
            end
        else
            Ip = nothing
            Vedge = IMAS.get_from(dd, Val{:vloop}, par.vloop_from)
        end

        # we need to run steady state twice, the first time to find the q-profile when the
        # current fully relaxes, and the second time we change the resisitivity to keep q>1
        i_qdes = 0
        for _ in (1, 2)
            actor.QO = QED.steady_state(actor.QO, η_imas(dd.core_profiles.profiles_1d[], i_qdes); Vedge, Ip)
            # check where q<1
            qval = 1.0 ./ abs.(actor.QO.ι.(cp1d.grid.rho_tor_norm))
            i_qdes = findlast(qval .< par.qmin_desired)
            if i_qdes === nothing
                break
            end
        end
    end

    return actor
end

function _finalize(actor::ActorQED)
    dd = actor.dd

    eqt = dd.equilibrium.time_slice[]
    B0 = eqt.global_quantities.vacuum_toroidal_field.b0

    cp1d = dd.core_profiles.profiles_1d[]
    j_total = QED.JB(actor.QO; ρ=cp1d.grid.rho_tor_norm) ./ B0

    if ismissing(cp1d, :j_non_inductive)
        cp1d.j_ohmic = j_total
    else
        cp1d.j_ohmic = j_total .- cp1d.j_non_inductive
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
    B0 = eqt.global_quantities.vacuum_toroidal_field.b0

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

"""
    η_imas(cp1d::IMAS.core_profiles__profiles_1d, i_qdes::Int; use_log::Bool=true)

Jardin's model for stationary sawteeth changes the plasma resistivity to raise q>1
"""
function η_imas(cp1d::IMAS.core_profiles__profiles_1d, i_qdes::Int; use_log::Bool=true)
    η = 1.0 ./ cp1d.conductivity_parallel
    if i_qdes != 0
        η[1:i_qdes] .= η[i_qdes]
    end
    return QED.η_FE(cp1d.grid.rho_tor_norm, η; use_log)
end
