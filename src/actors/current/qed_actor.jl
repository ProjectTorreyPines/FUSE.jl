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
    qmin_desired::Entry{Float64} = Entry{Float64}("-", "Keep the minimum magnitude of the q-profile above this value"; default=1.0, check=x -> @assert x >= 0 "qmin_desired >= 0")
    #== data flow parameters ==#
    ip_from::Switch{Symbol} = switch_get_from(:ip)
    vloop_from::Switch{Symbol} = switch_get_from(:vloop)
end

mutable struct ActorQED{D,P} <: CompoundAbstractActor{D,P}
    dd::IMAS.dd{D}
    par::OverrideParameters{P,FUSEparameters__ActorQED{P}}
    act::ParametersAllActors{P}
    ip_controller::ActorControllerIp{D,P}
    QO::Union{Nothing,QED.QED_state}
end

"""
    ActorQED(dd::IMAS.dd, act::ParametersAllActors; kw...)

Evolves the plasma current using the QED current diffusion solver.

The fundamental quantitiy being solved is `j_total` in `dd.core_profiles.profiles_1d[]`

!!! note

    This actor operates at "dd.global_time", any time advance must be done outside of the actor

        IMAS.new_timeslice!(dd, dd.global_time + Δt)
        dd.global_time += Δt
        ActorQED(dd, act)
"""
function ActorQED(dd::IMAS.dd, act::ParametersAllActors; kw...)
    actor = ActorQED(dd, act.ActorQED, act; kw...)
    step(actor)
    finalize(actor)
    return actor
end

function ActorQED(dd::IMAS.dd, par::FUSEparameters__ActorQED, act::ParametersAllActors; kw...)
    logging_actor_init(ActorQED)
    par = OverrideParameters(par; kw...)
    ip_controller = ActorControllerIp(dd, act.ActorControllerIp)
    return ActorQED(dd, par, act, ip_controller, nothing)
end

function _step(actor::ActorQED)
    dd = actor.dd
    par = actor.par

    eqt = dd.equilibrium.time_slice[]
    cp1d = dd.core_profiles.profiles_1d[]

    B0 = eqt.global_quantities.vacuum_toroidal_field.b0
    # no ohmic, no sawteeth, no time dependent
    j_non_inductive = IMAS.total_sources(dd.core_sources, cp1d; time0=dd.global_time, exclude_indexes=[7, 409, 701], fields=[:j_parallel]).j_parallel
    conductivity_parallel = IMAS.neo_conductivity(eqt, cp1d)

    # initialize QED
    # we must reinitialize to update the equilibrium metrics
    actor.QO = QED.initialize(dd, par.qmin_desired; uniform_rho=501)

    if par.Δt > 0.0 && par.Δt < Inf
        # current diffusion
        t0 = dd.global_time - par.Δt
        t1 = dd.global_time

        if par.solve_for == :vloop && par.vloop_from == :controllers__ip
            # staircase approach to call Ip control at each step of the current ramp: one QED diffuse call for each time step
            # NOTE: `QED.diffuse` is inefficient, would be beneficial to have a `QED.diffuse!` function without allocations
            δt = par.Δt / par.Nt
            No = par.Nt
            Ni = 1
        else
            δt = t1 - t0
            No = 1
            Ni = par.Nt
        end

        i_qdes = nothing
        flattened_j_non_inductive = j_non_inductive
        for tt in range(t0, t1, No + 1)[1:end-1]
            if par.solve_for == :ip
                Ip = IMAS.get_from(dd, Val{:ip}, par.ip_from; time0=tt + δt)
                Vedge = nothing
            else
                # run Ip controller if vloop_from == :controllers__ip
                if par.vloop_from == :controllers__ip
                    finalize(step(actor.ip_controller; time0=tt + δt))
                end
                Ip = nothing
                Vedge = IMAS.get_from(dd, Val{:vloop}, par.vloop_from; time0=tt + δt)
            end

            # check where q<1 based on the the q-profile at the previous
            # time-step and change the resisitivity to keep q>1
            qval = 1.0 ./ abs.(actor.QO.ι.(cp1d.grid.rho_tor_norm))
            i_qdes = findlast(qval .< par.qmin_desired)
            if i_qdes === nothing
                rho_qdes = -1.0
            else
                rho_qdes = cp1d.grid.rho_tor_norm[i_qdes]
            end

            η_jardin, flattened_j_non_inductive = QED.η_JBni_sawteeth(cp1d, j_non_inductive, rho_qdes; conductivity_parallel)
            actor.QO.JBni = QED.FE(cp1d.grid.rho_tor_norm, flattened_j_non_inductive .* B0)

            actor.QO = QED.diffuse(actor.QO, η_jardin, δt, Ni; Vedge, Ip, debug=false)
        end

    elseif par.Δt == Inf
        # steady state solution
        if par.solve_for == :ip
            Ip = IMAS.get_from(dd, Val{:ip}, par.ip_from)
            Vedge = nothing
        else
            Ip = nothing
            Vedge = IMAS.get_from(dd, Val{:vloop}, par.vloop_from)
        end

        # fist try full relaxation
        rho_qdes = 0.0
        η_jardin, flattened_j_non_inductive = QED.η_JBni_sawteeth(cp1d, j_non_inductive, rho_qdes; conductivity_parallel)
        actor.QO.JBni = QED.FE(cp1d.grid.rho_tor_norm, flattened_j_non_inductive .* B0)
        actor.QO = QED.steady_state(actor.QO, η_jardin; Vedge, Ip)
        qval = 1.0 ./ abs.(actor.QO.ι.(cp1d.grid.rho_tor_norm))
        i_qdes = findlast(qval .< par.qmin_desired)

        # then identify inversion radius
        if i_qdes !== nothing
            for i_qdes in 1:i_qdes
                rho_qdes = cp1d.grid.rho_tor_norm[i_qdes]
                η_jardin, flattened_j_non_inductive = QED.η_JBni_sawteeth(cp1d, j_non_inductive, rho_qdes; conductivity_parallel)
                actor.QO.JBni = QED.FE(cp1d.grid.rho_tor_norm, flattened_j_non_inductive .* B0)
                actor.QO = QED.steady_state(actor.QO, η_jardin; Vedge, Ip)
                qval = 1.0 ./ abs.(actor.QO.ι.(cp1d.grid.rho_tor_norm))
                if findlast(qval .< par.qmin_desired) === nothing
                    break
                end
            end
        end
    else

        error("act.ActorQED.Δt = $(par.Δt) is not valid")
    end

    cp1d.j_total = QED.JB(actor.QO; ρ=cp1d.grid.rho_tor_norm) ./ B0
    cp1d.j_non_inductive = flattened_j_non_inductive
    cp1d.q = 1.0 ./ abs.(actor.QO.ι.(cp1d.grid.rho_tor_norm))

    return actor
end
