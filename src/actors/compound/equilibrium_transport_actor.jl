#= ========================= =#
#  ActorEquilibriumTransport  #
#= ========================= =#
Base.@kwdef mutable struct FUSEparameters__ActorEquilibriumTransport{T} <: ParametersActor where {T<:Real}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    do_plot::Entry{Bool} = Entry(Bool, "-", "plot"; default=false)
    max_iter::Entry{Int} = Entry(Int, "-", "max number of transport-equilibrium iterations"; default=5)
    convergence_error::Entry{T} = Entry(T, "-", "Convergence error threshold (relative change in current and pressure profiles)"; default=1E-2)
end

mutable struct ActorEquilibriumTransport <: PlasmaAbstractActor
    dd::IMAS.dd
    par::FUSEparameters__ActorEquilibriumTransport
    act::ParametersAllActors
    actor_tr::ActorCoreTransport
    actor_hc::ActorHCD
    actor_jt::ActorSteadyStateCurrent
    actor_eq::ActorEquilibrium
end

"""
    ActorEquilibriumTransport(dd::IMAS.dd, act::ParametersAllActors; kw...)

Compound actor that runs the following actors in succesion:
* ActorSteadyStateCurrent
* ActorHCD
* ActorCoreTransport
* ActorEquilibrium

!!! note
    Stores data in `dd.equilibrium, dd.core_profiles, dd.core_sources`
"""
function ActorEquilibriumTransport(dd::IMAS.dd, act::ParametersAllActors; kw...)
    par = act.ActorEquilibriumTransport
    actor = ActorEquilibriumTransport(dd, par, act; kw...)
    step(actor)
    finalize(actor)
    return actor
end

function ActorEquilibriumTransport(dd::IMAS.dd, par::FUSEparameters__ActorEquilibriumTransport, act::ParametersAllActors; kw...)
    logging_actor_init(ActorEquilibriumTransport)
    par = par(kw...)
    actor_tr = ActorCoreTransport(dd, act.ActorCoreTransport, act)
    actor_hc = ActorHCD(dd, act.ActorHCD, act)
    actor_jt = ActorSteadyStateCurrent(dd, act.ActorSteadyStateCurrent)
    actor_eq = ActorEquilibrium(dd, act.ActorEquilibrium, act)
    return ActorEquilibriumTransport(dd, par, act, actor_tr, actor_hc, actor_jt, actor_eq)
end

function _step(actor::ActorEquilibriumTransport)
    dd = actor.dd
    par = actor.par

    if par.do_plot
        pe = plot(dd.equilibrium; color=:gray, label=" (before)", coordinate=:rho_tor_norm)
        pp = plot(dd.core_profiles; color=:gray, label=" (before)")
        ps = plot(dd.core_sources; color=:gray, label=" (before)")
    end

    # run HCD to get updated current drive
    finalize(step(actor.actor_hc))

    # set j_ohmic to steady state
    finalize(step(actor.actor_jt))

    # set actors switches specific to this workflow
    if actor.actor_eq.par.model == :CHEASE
        chease_par = actor.actor_eq.eq_actor.par
        orig_par_chease = deepcopy(chease_par)
        chease_par.free_boundary = false
        chease_par.rescale_eq_to_ip = true
    end

    try
        total_error = Float64[]
        while isempty(total_error) || (total_error[end] > par.convergence_error)
            # get current and pressure profiles before updating them
            j_tor_before = dd.core_profiles.profiles_1d[].j_tor
            pressure_before = dd.core_profiles.profiles_1d[].pressure

            # run transport actor
            finalize(step(actor.actor_tr))

            # run HCD to get updated current drive
            finalize(step(actor.actor_hc))

            # set j_ohmic to steady state
            finalize(step(actor.actor_jt))

            # run equilibrium actor with the updated beta
            finalize(step(actor.actor_eq))

            # evaluate change in current and pressure profiles after the update
            j_tor_after = dd.core_profiles.profiles_1d[].j_tor
            pressure_after = dd.core_profiles.profiles_1d[].pressure
            error_jtor = sum((j_tor_after .- j_tor_before) .^ 2) / sum(j_tor_before .^ 2)
            error_pressure = sum((pressure_after .- pressure_before) .^ 2) / sum(pressure_before .^ 2)
            push!(total_error, sqrt(error_jtor + error_pressure) / 2.0)

            if par.do_plot
                @info("Iteration = $(length(total_error)) , convergence error = $(round(total_error[end],digits = 3)), threshold = $(par.convergence_error)")
            end

            if (total_error[end] > par.convergence_error) && (length(total_error) == par.max_iter)
                @warn "Max number of iterations ($(par.max_iter)) has been reached with convergence error of $(round(total_error[end],digits = 3)) compared to threshold of $(par.convergence_error)"
                break
            end
        end

    finally
        if actor.actor_eq.par.model == :CHEASE
            actor.actor_eq.eq_actor.par = orig_par_chease
            if orig_par_chease.free_boundary
                finalize(step(actor.actor_eq))
            end
        end
    end

    if par.do_plot
        display(plot!(pe, dd.equilibrium, coordinate=:rho_tor_norm, label=" (after)"))
        display(plot!(pp, dd.core_profiles, label=" (after)"))
        display(plot!(ps, dd.core_sources, label=" (after)"))
    end

    return actor
end
