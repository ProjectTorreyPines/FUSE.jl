#= ===================== =#
#  ActorStationaryPlasma  #
#= ===================== =#
Base.@kwdef mutable struct FUSEparameters__ActorStationaryPlasma{T<:Real} <: ParametersActor{T}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    _time::Float64 = NaN
    max_iter::Entry{Int} = Entry{Int}("-", "max number of transport-equilibrium iterations"; default=5)
    convergence_error::Entry{T} = Entry{T}("-", "Convergence error threshold (relative change in current and pressure profiles)"; default=5E-2)
    #== display and debugging parameters ==#
    do_plot::Entry{Bool} = act_common_parameters(; do_plot=false)
    verbose::Entry{Bool} = act_common_parameters(; verbose=false)
end

mutable struct ActorStationaryPlasma{D,P} <: CompoundAbstractActor{D,P}
    dd::IMAS.dd{D}
    par::FUSEparameters__ActorStationaryPlasma{P}
    act::ParametersAllActors
    actor_tr::ActorCoreTransport{D,P}
    actor_ped::Union{ActorPedestal{D,P},ActorNoOperation{D,P}}
    actor_hc::ActorHCD{D,P}
    actor_jt::ActorCurrent{D,P}
    actor_eq::ActorEquilibrium{D,P}
end

"""
    ActorStationaryPlasma(dd::IMAS.dd, act::ParametersAllActors; kw...)

Compound actor that runs the following actors in succesion:

  - ActorCurrent
  - ActorHCD
  - ActorCoreTransport
  - ActorEquilibrium

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

    if act.ActorCoreTransport.model == :FluxMatcher
        actor_ped = ActorPedestal(
            dd,
            act.ActorPedestal,
            act;
            ip_from=:core_profiles,
            βn_from=:core_profiles,
            ne_from=:pulse_schedule,
            zeff_ped_from=:pulse_schedule,
            rho_nml=actor_tr.tr_actor.par.rho_transport[end-1],
            rho_ped=actor_tr.tr_actor.par.rho_transport[end])
    else
        actor_ped = ActorNoOperation(dd, act.ActorNoOperation)
    end

    actor_hc = ActorHCD(dd, act.ActorHCD, act)

    actor_jt = ActorCurrent(dd, act.ActorCurrent, act; model=:QED, ip_from=:pulse_schedule, vloop_from=:pulse_schedule)

    actor_eq = ActorEquilibrium(dd, act.ActorEquilibrium, act; ip_from=:core_profiles)

    return ActorStationaryPlasma(dd, par, act, actor_tr, actor_ped, actor_hc, actor_jt, actor_eq)
end

function _step(actor::ActorStationaryPlasma)
    dd = actor.dd
    par = actor.par

    if par.do_plot
        pe = plot(dd.equilibrium; color=:gray, label=" (before)", coordinate=:rho_tor_norm)
        pp = plot(dd.core_profiles; color=:gray, label=" (before)")
        ps = plot(dd.core_sources; color=:gray, label=" (before)")

        println("initial")
        @printf("Jtor0  = %.2f MA/m²\n", getproperty(dd.core_profiles.profiles_1d[], :j_tor, [0.0])[1] / 1e6)
        @printf("P0     = %.2f kPa\n", getproperty(dd.core_profiles.profiles_1d[], :pressure, [0.0])[1] / 1e3)
        @printf("βn_MHD = %.2f\n", dd.equilibrium.time_slice[].global_quantities.beta_normal)
        @printf("βn_tot = %.2f\n", @ddtime(dd.summary.global_quantities.beta_tor_norm.value))
        @printf("ne_ped = %.2e m⁻³\n", @ddtime(dd.summary.local.pedestal.n_e.value))
        @printf("Te_ped = %.2e eV\n", @ddtime(dd.summary.local.pedestal.t_e.value))
        @printf(" ρ_ped = %.4f\n", @ddtime(dd.summary.local.pedestal.position.rho_tor_norm))
    end

    # set actors switches specific to this workflow
    if typeof(actor.actor_eq) <: ActorCHEASE
        chease_par = actor.actor_eq.eq_actor.par
        orig_par_chease = deepcopy(chease_par)
        chease_par.rescale_eq_to_ip = true
    end

    # set Δt of the time-dependent actors
    if actor.actor_jt.par.model == :QED
        actor.actor_jt.jt_actor.par.Δt = Inf
    end
    if actor.actor_tr.par.model == :FluxMatcher
        actor.actor_tr.tr_actor.par.Δt = Inf
    end

    ProgressMeter.ijulia_behavior(:clear)
    prog = ProgressMeter.Progress((par.max_iter + 1) * 5 + 2; dt=0.0, showspeed=true, enabled=par.verbose && !par.do_plot)
    old_logging = actor_logging(dd, !(par.verbose && !par.do_plot))
    total_error = Float64[]
    cp1d = dd.core_profiles.profiles_1d[]
    try

        if !(par.verbose && !par.do_plot)
            logging(Logging.Info, :actors, " "^workflow_depth(actor.dd) * "--------------- 1/$(par.max_iter)")
        end

        # unless `par.max_iter==1` we want to iterate at least twice to ensure consistency between equilibrium and profiles
        while length(total_error) < 2 || (total_error[end] > par.convergence_error)

            # get current and pressure profiles before updating them
            j_tor_before = cp1d.j_tor
            pressure_before = cp1d.pressure

            # core_profiles, core_sources, core_transport grids from latest equilibrium
            latest_equilibrium_grids!(dd)

            # run HCD to get updated current drive
            ProgressMeter.next!(prog; showvalues=progress_ActorStationaryPlasma(total_error, actor, actor.actor_hc))
            finalize(step(actor.actor_hc))

            # evolve j_ohmic (because hcd has changed non-inductive current drive)
            ProgressMeter.next!(prog; showvalues=progress_ActorStationaryPlasma(total_error, actor, actor.actor_jt))
            finalize(step(actor.actor_jt))

            # run pedestal actor
            ProgressMeter.next!(prog; showvalues=progress_ActorStationaryPlasma(total_error, actor, actor.actor_ped))
            finalize(step(actor.actor_ped))

            # run transport actor
            ProgressMeter.next!(prog; showvalues=progress_ActorStationaryPlasma(total_error, actor, actor.actor_tr))
            finalize(step(actor.actor_tr))

            # evolve j_ohmic (because transport and pedestal have updated my bootstrap)
            ProgressMeter.next!(prog; showvalues=progress_ActorStationaryPlasma(total_error, actor, actor.actor_jt))
            finalize(step(actor.actor_jt))

            # run equilibrium actor with the updated beta
            ProgressMeter.next!(prog; showvalues=progress_ActorStationaryPlasma(total_error, actor, actor.actor_eq))
            finalize(step(actor.actor_eq))

            # evaluate change in current and pressure profiles after the update
            j_tor = cp1d.j_tor
            dj2 = (k, x) -> (j_tor[k] .- j_tor_before[k])^2
            j2 = (k, x) -> j_tor_before[k]^2
            error_jtor = trapz(cp1d.grid.area, dj2) / trapz(cp1d.grid.area, j2)

            pressure = cp1d.pressure
            dp2 = (k, x) -> (pressure[k] .- pressure_before[k])^2
            p2 = (k, x) -> pressure_before[k]^2
            error_pressure = trapz(cp1d.grid.volume, dp2) / trapz(cp1d.grid.volume, p2)
            push!(total_error, sqrt(error_jtor + error_pressure) / 2.0)

            if par.do_plot
                plot!(pe, dd.equilibrium; coordinate=:rho_tor_norm, label="i=$(length(total_error))")
                plot!(pp, dd.core_profiles; label="i=$(length(total_error))")
                plot!(ps, dd.core_sources; label="i=$(length(total_error))")

                @printf("\n")
                @printf(" Jtor0 = %.2f MA m²\n", cp1d.j_tor[1] / 1e6)
                @printf("    P0 = %.2f kPa\n", cp1d.pressure[1] / 1e3)
                @printf("βn_MHD = %.2f\n", dd.equilibrium.time_slice[].global_quantities.beta_normal)
                @printf("βn_tot = %.2f\n", @ddtime(dd.summary.global_quantities.beta_tor_norm.value))
                @printf("ne_ped = %.2e m⁻³\n", @ddtime(dd.summary.local.pedestal.n_e.value))
                @printf("Te_ped = %.2e eV\n", @ddtime(dd.summary.local.pedestal.t_e.value))
                @printf(" ρ_ped = %.4f\n", @ddtime(dd.summary.local.pedestal.position.rho_tor_norm))
                @info("Iteration = $(length(total_error)) , convergence error = $(round(total_error[end],digits = 5)), threshold = $(par.convergence_error)")
            end

            if !(par.verbose && !par.do_plot)
                logging(
                    Logging.Info,
                    :actors,
                    " "^workflow_depth(actor.dd) * "--------------- $(length(total_error))/$(par.max_iter) @ $(@sprintf("%3.2f",100*total_error[end]/par.convergence_error))%"
                )
            end

            if (total_error[end] > par.convergence_error) && (length(total_error) == par.max_iter)
                @warn "Max number of iterations ($(par.max_iter)) has been reached with convergence error of (1)$(collect(map(x->round(x,digits = 3),total_error)))($(length(total_error))) compared to threshold of $(par.convergence_error)"
                break
            elseif par.max_iter == 1
                break
            end
        end

    finally
        if typeof(actor.actor_eq) <: ActorCHEASE
            actor.actor_eq.eq_actor.par = orig_par_chease
        end
        actor_logging(dd, old_logging)
    end
    ProgressMeter.finish!(prog; showvalues=progress_ActorStationaryPlasma(total_error, actor))

    if par.do_plot
        display(pe)
        display(pp)
        display(ps)
    end

    return actor
end

function progress_ActorStationaryPlasma(total_error::Vector{Float64}, actor::ActorStationaryPlasma, step_actor::Union{Nothing,AbstractActor}=nothing)
    dd = actor.dd
    par = actor.par
    cp1d = dd.core_profiles.profiles_1d[]
    tmp = [
        (par.max_iter == 1 ? "                 iteration" : "         iteration (min 2)", "$(length(total_error))/$(par.max_iter)"),
        ("required convergence error", par.convergence_error),
        ("       convergence history", isempty(total_error) ? "N/A" : reverse(total_error)),
        ("                     stage", step_actor === nothing ? "N/A" : "$(name(step_actor))"),
        ("                   Ip [MA]", IMAS.get_from(dd, Val{:ip}, :equilibrium) / 1E6),
        ("                 Ti0 [keV]", cp1d.t_i_average[1] / 1E3),
        ("                 Te0 [keV]", cp1d.electrons.temperature[1] / 1E3),
        ("            ne0 [10²⁰ m⁻³]", cp1d.electrons.density_thermal[1] / 1E20)
    ]
    return tuple(tmp...)
end
