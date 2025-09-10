#= ===================== =#
#  ActorStationaryPlasma  #
#= ===================== =#
@actor_parameters_struct ActorStationaryPlasma{T} begin
    max_iterations::Entry{Int} = Entry{Int}("-", "max number of transport-equilibrium iterations"; default=5)
    convergence_error::Entry{T} = Entry{T}("-", "Convergence error threshold (relative change in current and pressure profiles)"; default=5E-2)
    #== display and debugging parameters ==#
    do_plot::Entry{Bool} = act_common_parameters(; do_plot=false)
    verbose::Entry{Bool} = act_common_parameters(; verbose=false)
end

mutable struct ActorStationaryPlasma{D,P} <: CompoundAbstractActor{D,P}
    dd::IMAS.dd{D}
    par::OverrideParameters{P,FUSEparameters__ActorStationaryPlasma{P}}
    act::ParametersAllActors{P}
    actor_tr::ActorCoreTransport{D,P}
    actor_ped::Union{ActorPedestal{D,P},ActorNoOperation{D,P}}
    actor_hc::ActorHCD{D,P}
    actor_jt::ActorCurrent{D,P}
    actor_eq::ActorEquilibrium{D,P}
    actor_saw::ActorSawteeth{D,P}
end

"""
    ActorStationaryPlasma(dd::IMAS.dd, act::ParametersAllActors; kw...)

Compound actor that runs the following actors in succesion to find a self-consistent stationary plasma solution

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
    par = OverrideParameters(par; kw...)

    actor_tr = ActorCoreTransport(dd, act.ActorCoreTransport, act)

    # allows users to hardwire `rho_nml` and `rho_ped` (same logic here as in ActorDynamicPlasma)
    if act.ActorCoreTransport.model == :FluxMatcher && ismissing(act.ActorPedestal, :rho_nml)
        rho_nml = actor_tr.tr_actor.par.rho_transport[end-1]
    elseif ismissing(act.ActorPedestal, :rho_nml)
        rho_nml = act.ActorFluxMatcher.rho_transport[end-1]
    else
        rho_nml = act.ActorPedestal.rho_nml
    end
    if act.ActorCoreTransport.model == :FluxMatcher && ismissing(act.ActorPedestal, :rho_ped)
        rho_ped = actor_tr.tr_actor.par.rho_transport[end]
    elseif ismissing(act.ActorPedestal, :rho_ped)
        rho_ped = act.ActorFluxMatcher.rho_transport[end]
    else
        rho_ped = act.ActorPedestal.rho_ped
    end
    actor_ped = ActorPedestal(
        dd,
        act.ActorPedestal,
        act;
        ip_from=:core_profiles,
        βn_from=:core_profiles,
        ne_from=:pulse_schedule,
        zeff_from=:pulse_schedule,
        rho_nml,
        rho_ped)

    actor_hc = ActorHCD(dd, act.ActorHCD, act)

    actor_jt = ActorCurrent(dd, act.ActorCurrent, act; ip_from=:pulse_schedule, vloop_from=:pulse_schedule)

    actor_eq = ActorEquilibrium(dd, act.ActorEquilibrium, act; ip_from=:core_profiles)

    actor_saw = ActorSawteeth(dd, act.ActorSawteeth)

    return ActorStationaryPlasma(dd, par, act, actor_tr, actor_ped, actor_hc, actor_jt, actor_eq, actor_saw)
end

function _step(actor::ActorStationaryPlasma)
    dd = actor.dd
    par = actor.par

    if par.do_plot
        pe = plot(dd.equilibrium; color=:gray, label="", coordinate=:rho_tor_norm)
        pp = plot(dd.core_profiles; color=:gray, label="")
        ps = plot(dd.core_sources; color=:gray, label="")

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
    was_logging = actor_logging(dd)
    is_logging = was_logging && !(par.verbose && !par.do_plot)
    actor_logging(dd, is_logging)
    prog = ProgressMeter.Progress((par.max_iterations + 1) * 6 + 2; dt=0.0, showspeed=true, enabled=par.verbose && !par.do_plot)
    total_error = Float64[]
    cp1d = dd.core_profiles.profiles_1d[]
    try

        if is_logging
            logging(Logging.Info, :actors, " "^workflow_depth(actor.dd) * "--------------- 1/$(par.max_iterations)")
        end

        # unless `par.max_iterations==1` we want to iterate at least twice to ensure consistency between equilibrium and profiles
        while length(total_error) < 2 || (total_error[end] > par.convergence_error)

            # get current and pressure profiles before updating them
            j_tor_before = cp1d.j_tor
            pressure_before = cp1d.pressure

            # core_profiles, core_sources, core_transport grids from latest equilibrium
            latest_equilibrium_grids!(dd)

            # run HCD to get updated current drive
            ProgressMeter.next!(prog; showvalues=progress_ActorStationaryPlasma(total_error, actor, actor.actor_hc))
            finalize(step(actor.actor_hc))

            # run pedestal actor
            ProgressMeter.next!(prog; showvalues=progress_ActorStationaryPlasma(total_error, actor, actor.actor_ped))
            finalize(step(actor.actor_ped))

            # run transport actor
            ProgressMeter.next!(prog; showvalues=progress_ActorStationaryPlasma(total_error, actor, actor.actor_tr))
            finalize(step(actor.actor_tr))

            # evolve j_ohmic
            ProgressMeter.next!(prog; showvalues=progress_ActorStationaryPlasma(total_error, actor, actor.actor_jt))
            finalize(step(actor.actor_jt))

            # run sawteeth actor
            ProgressMeter.next!(prog; showvalues=progress_ActorStationaryPlasma(total_error, actor, actor.actor_saw))
            finalize(step(actor.actor_saw))

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
                plot!(pe, dd.equilibrium; label="", coordinate=:rho_tor_norm)
                plot!(pp, dd.core_profiles; label="", legend=nothing)
                plot!(ps, dd.core_sources; label="", legend=nothing)

                @printf("\n")
                @printf(" Jtor0 = %.2f MA m²\n", cp1d.j_tor[1] / 1e6)
                @printf("    P0 = %.2f kPa\n", cp1d.pressure[1] / 1e3)
                @printf("βn_MHD = %.2f\n", dd.equilibrium.time_slice[].global_quantities.beta_normal)
                @printf("βn_tot = %.2f\n", @ddtime(dd.summary.global_quantities.beta_tor_norm.value))
                @printf("ne_ped = %.2e m⁻³\n", @ddtime(dd.summary.local.pedestal.n_e.value))
                @printf("Te_ped = %.2e eV\n", @ddtime(dd.summary.local.pedestal.t_e.value))
                @printf(" ρ_ped = %.4f\n", @ddtime(dd.summary.local.pedestal.position.rho_tor_norm))
                @printf(" ϵ jtor = %.4f\n", error_jtor)
                @printf(" ϵ pres = %.4f\n", error_pressure)
                @info("Iteration = $(length(total_error)) , convergence error = $(round(total_error[end],digits = 5)), threshold = $(par.convergence_error)")
            end

            if is_logging
                logging(
                    Logging.Info,
                    :actors,
                    " "^workflow_depth(actor.dd) *
                    "--------------- $(length(total_error))/$(par.max_iterations) @ $(@sprintf("%3.2f",100*total_error[end]/par.convergence_error))%"
                )
            end

            callback(actor, :iteration_end; total_error)

            if (total_error[end] > par.convergence_error) && (length(total_error) == par.max_iterations)
                if is_logging
                    @warn "Max number of iterations ($(par.max_iterations)) has been reached with convergence error of (1)$(collect(map(x->round(x,digits = 3),total_error)))($(length(total_error))) compared to threshold of $(par.convergence_error)"
                end
                break
            elseif par.max_iterations == 1
                break
            end
        end

    finally
        if typeof(actor.actor_eq) <: ActorCHEASE
            actor.actor_eq.eq_actor.par = orig_par_chease
        end
        actor_logging(dd, was_logging)
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
        (par.max_iterations == 1 ? "                 iteration" : "         iteration (min 2)", "$(length(total_error))/$(par.max_iterations)"),
        ("required convergence error", par.convergence_error),
        ("       convergence history", isempty(total_error) ? "N/A" : reverse(total_error)),
        ("                     stage", step_actor === nothing ? "N/A" : "$(name(step_actor))"),
        ("                   Ip [MA]", IMAS.get_from(dd, Val(:ip), :equilibrium) / 1E6),
        ("                 Ti0 [keV]", cp1d.t_i_average[1] / 1E3),
        ("                 Te0 [keV]", cp1d.electrons.temperature[1] / 1E3),
        ("            ne0 [10²⁰ m⁻³]", cp1d.electrons.density_thermal[1] / 1E20),
        ("                 max(zeff)", maximum(cp1d.zeff)),
        ("               ω0 [krad/s]", cp1d.rotation_frequency_tor_sonic[1] / 1E3)
    ]
    return tuple(tmp...)
end
