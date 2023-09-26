#= ===================== =#
#  ActorStationaryPlasma  #
#= ===================== =#
Base.@kwdef mutable struct FUSEparameters__ActorStationaryPlasma{T} <: ParametersActor where {T<:Real}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    do_plot::Entry{Bool} = Entry{Bool}("-", "Plot"; default=false)
    update_pedestal::Entry{Bool} = Entry{Bool}("-", "Updates the pedestal before and after the transport step"; default=false)
    max_iter::Entry{Int} = Entry{Int}("-", "max number of transport-equilibrium iterations"; default=5)
    convergence_error::Entry{T} = Entry{T}("-", "Convergence error threshold (relative change in current and pressure profiles)"; default=5E-2)
end

mutable struct ActorStationaryPlasma{D,P} <: PlasmaAbstractActor
    dd::IMAS.dd{D}
    par::FUSEparameters__ActorStationaryPlasma{P}
    act::ParametersAllActors
    actor_tr::ActorCoreTransport{D,P}
    actor_hc::ActorHCD{D,P}
    actor_jt::ActorCurrent{D,P}
    actor_eq::ActorEquilibrium{D,P}
    actor_ped::ActorPedestal{D,P}
end

"""
    ActorStationaryPlasma(dd::IMAS.dd, act::ParametersAllActors; kw...)

Compound actor that runs the following actors in succesion:
* ActorCurrent
* ActorHCD
* ActorCoreTransport
* ActorEquilibrium

!!! note
    Stores data in `dd.equilibrium`, `dd.core_profiles`, `dd.core_sources`, `dd.core_transport`
"""
function ActorStationaryPlasma(dd::IMAS.dd, act::ParametersAllActors; kw...)
    actor = ActorStationaryPlasma(dd, act.ActorStationaryPlasma, act; kw...)
    step(actor)
    finalize(actor)
    return actor
end

function ActorStationaryPlasma(dd::IMAS.dd, par::FUSEparameters__ActorStationaryPlasma, act::ParametersAllActors; kw...)
    logging_actor_init(ActorStationaryPlasma)
    par = par(kw...)
    actor_tr = ActorCoreTransport(dd, act.ActorCoreTransport, act)
    actor_hc = ActorHCD(dd, act.ActorHCD, act)
    actor_jt = ActorCurrent(dd, act.ActorCurrent, act; ip_from=:pulse_schedule, vloop_from=:pulse_schedule)
    actor_eq = ActorEquilibrium(dd, act.ActorEquilibrium, act; ip_from=:pulse_schedule)
    actor_ped = ActorPedestal(dd, act.ActorPedestal; ip_from=:equilibrium, βn_from=:equilibrium)
    return ActorStationaryPlasma(dd, par, act, actor_tr, actor_hc, actor_jt, actor_eq, actor_ped)
end

function _step(actor::ActorStationaryPlasma)
    dd = actor.dd
    par = actor.par
    act = actor.act

    if par.do_plot
        pe = plot(dd.equilibrium; color=:gray, label=" (before)", coordinate=:rho_tor_norm)
        pp = plot(dd.core_profiles; color=:gray, label=" (before)")
        ps = plot(dd.core_sources; color=:gray, label=" (before)")

        @printf("Jtor0_before = %.2f MA/m²\n", getproperty(dd.core_profiles.profiles_1d[],:j_tor,[0.0])[1]/1e6)
        @printf("P0_before = %.2f kPa\n", getproperty(dd.core_profiles.profiles_1d[], :pressure, [0.0])[1]/1e3)
        @printf("βn_MHD = %.2f\n", dd.equilibrium.time_slice[].global_quantities.beta_normal)
        @printf("βn_tot = %.2f\n", @ddtime(dd.summary.global_quantities.beta_tor_norm.value))
        @printf("Te_ped = %.2e eV\n", @ddtime(dd.summary.local.pedestal.t_e.value))
        @printf("rho_ped = %.4f\n", @ddtime(dd.summary.local.pedestal.position.rho_tor_norm))
    end

    # set actors switches specific to this workflow
    if typeof(actor.actor_eq) <: ActorCHEASE
        chease_par = actor.actor_eq.eq_actor.par
        orig_par_chease = deepcopy(chease_par)
        chease_par.rescale_eq_to_ip = true
    end

    try
        # run HCD to get updated current drive
        finalize(step(actor.actor_hc))

        # evolve j_ohmic
        finalize(step(actor.actor_jt))

        total_error = Float64[]
        while isempty(total_error) || (total_error[end] > par.convergence_error)
            # get current and pressure profiles before updating them
            j_tor_before = dd.core_profiles.profiles_1d[].j_tor
            pressure_before = dd.core_profiles.profiles_1d[].pressure

            # core_profiles, core_sources, core_transport grids from latest equilibrium
            latest_equilibrium_grids!(dd)

            # run transport actor
            finalize(step(actor.actor_tr))

            # run pedestal actor
            if par.update_pedestal
                finalize(step(actor.actor_ped))
            end

            # run HCD to get updated current drive
            finalize(step(actor.actor_hc))

            # evolve j_ohmic
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
                plot!(pe, dd.equilibrium, coordinate=:rho_tor_norm, label="i=$(length(total_error))")
                plot!(pp, dd.core_profiles, label="i=$(length(total_error))")
                plot!(ps, dd.core_sources, label="i=$(length(total_error))")

                @printf("\n")
                @printf("Jtor0_after = %.2f MA\n", j_tor_after[1]/1e6)
                @printf("P0_after = %.2f kPa\n", pressure_after[1]/1e3)
                @printf("βn_MHD = %.2f\n", dd.equilibrium.time_slice[].global_quantities.beta_normal)
                @printf("βn_tot = %.2f\n", @ddtime(dd.summary.global_quantities.beta_tor_norm.value))
                @printf("Te_ped = %.2e eV\n", @ddtime(dd.summary.local.pedestal.t_e.value))
                @printf("rho_ped = %.4f\n", @ddtime(dd.summary.local.pedestal.position.rho_tor_norm))
                @info("Iteration = $(length(total_error)) , convergence error = $(round(total_error[end],digits = 5)), threshold = $(par.convergence_error)")
            end

            if (total_error[end] > par.convergence_error) && (length(total_error) == par.max_iter)
                @warn "Max number of iterations ($(par.max_iter)) has been reached with convergence error of $(collect(map(x->round(x,digits = 3),total_error))) compared to threshold of $(par.convergence_error)"
                break
            end
        end

        # run HCD to get updated current drive
        finalize(step(actor.actor_hc))

        # evolve j_ohmic
        finalize(step(actor.actor_jt))

    finally
        if typeof(actor.actor_eq) <: ActorCHEASE
            actor.actor_eq.eq_actor.par = orig_par_chease
        end
    end

    if par.do_plot
        display(pe)
        display(pp)
        display(ps)
    end

    return actor
end
