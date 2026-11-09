#= ===================== =#
#  ActorStationaryPlasma  #
#= ===================== =#
@actor_parameters_struct ActorStationaryPlasma{T} begin
    max_iterations::Entry{Int} = Entry{Int}("-", "max number of transport-equilibrium iterations"; default=5)
    convergence_error::Entry{T} = Entry{T}("-", "Convergence error threshold (relative change in current and pressure profiles)"; default=5E-2)
    couple_sol::Entry{Bool} = Entry{Bool}("-", "Run ActorSOL at each iteration and use its separatrix Te as the pedestal boundary condition"; default=false)
    couple_sol_ti::Entry{Bool} = Entry{Bool}(
        "-",
        "With `couple_sol=true`, also use the SOL separatrix Ti as ActorPedestal `Ti_sep` (only the separatrix point: the pedestal Ti height still follows `T_ratio_pedestal`)";
        default=false)
    sol_relaxation::Entry{T} = Entry{T}(
        "-",
        "Under-relaxation of the SOL-coupled Te_sep (and Ti_sep): Te_sep ← Te_sep + α (Te_sep_SOL - Te_sep). The first iteration always takes the SOL value (α=1)";
        default=1.0,
        check=x -> @assert 0.0 < x <= 1.0 "sol_relaxation must be in (0, 1]")
    #== display and debugging parameters ==#
    do_plot::Entry{Bool} = act_common_parameters(; do_plot=false)
    verbose::Entry{Bool} = act_common_parameters(; verbose=false)
end

mutable struct ActorStationaryPlasma{D,P} <: CompoundAbstractActor{D,P}
    dd::IMAS.DD{D}
    par::OverrideParameters{P,FUSEparameters__ActorStationaryPlasma{P}}
    act::ParametersAllActors{P}
    actor_tr::ActorCoreTransport{D,P}
    actor_ped::Union{ActorPedestal{D,P},ActorNoOperation{D,P}}
    actor_src::ActorSources{D,P}
    actor_jt::ActorCurrent{D,P}
    actor_eq::ActorEquilibrium{D,P}
    actor_saw::ActorSawteethSource{D,P}
    actor_sol::Union{Nothing,ActorSOL{D,P}}
end

"""
    ActorStationaryPlasma(dd::IMAS.DD, act::ParametersAllActors; kw...)

Iteratively solves for self-consistent stationary plasma equilibrium and transport.

This compound actor finds steady-state plasma solutions by iterating between plasma
physics models until convergence. It balances current drive, heating, transport,
and equilibrium to achieve a consistent stationary state.

Core iteration workflow:
1. **Current Evolution**: Updates current density from current drive and resistive diffusion
2. **Heating & Current Drive**: Calculates power deposition and driven current from external sources
3. **Pedestal**: Determines edge boundary conditions and pedestal profiles
4. **Core Transport**: Evolves temperature and density profiles based on transport fluxes
5. **Sawteeth**: Applies sawtooth mixing to maintain realistic q-profiles
6. **Equilibrium**: Solves MHD equilibrium with updated pressure and current profiles

The iteration continues until the relative changes in current density and pressure
profiles fall below the convergence threshold. Each iteration updates the flux surface
geometry based on the latest equilibrium solution.

Key convergence criteria:
- Relative change in current density profile (volume-averaged)
- Relative change in pressure profile (volume-averaged)
- Combined error must be below `convergence_error` threshold
- Maximum iterations limit prevents infinite loops

Actors coordinated:
- **ActorCurrent**: Current density evolution (Ohmic + driven + bootstrap)
- **ActorSources**: Heating and current drive power deposition
- **ActorPedestal**: Edge pressure and temperature boundary conditions
- **ActorCoreTransport**: Core transport flux calculations
- **ActorSawteethSource**: Sawtooth crash mixing and q-profile flattening
- **ActorEquilibrium**: MHD equilibrium solution with updated profiles
- **ActorSOL** (only if `couple_sol=true`): separatrix Te from a SOL model (e.g. `:solwidth` or `:box`),
  run after the sources and fed to ActorPedestal as `Te_sep` (under-relaxed by `sol_relaxation`).
  With `couple_sol_ti=true` the SOL separatrix Ti is also passed as `Ti_sep` (separatrix point only).
  The relative change of `Te_sep` (and `Ti_sep`) is included in the convergence error.
  Only pedestal models that use `Te_sep` (`:EPED`, `:WPED`, `:analytic`, `:dynamic`) are affected by it.

!!! note

    Stores converged plasma state in `dd.equilibrium`, `dd.core_profiles`,
    `dd.core_sources`, `dd.core_transport`
"""
function ActorStationaryPlasma(dd::IMAS.DD, act::ParametersAllActors; kw...)
    actor = ActorStationaryPlasma(dd, act.ActorStationaryPlasma, act; kw...)
    step(actor)
    finalize(actor)
    return actor
end

function ActorStationaryPlasma(dd::IMAS.DD, par::FUSEparameters__ActorStationaryPlasma, act::ParametersAllActors; kw...)
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

    actor_src = ActorSources(dd, act.ActorSources, act)

    actor_jt = ActorCurrent(dd, act.ActorCurrent, act; ip_from=:pulse_schedule, vloop_from=:pulse_schedule)

    actor_eq = ActorEquilibrium(dd, act.ActorEquilibrium, act; ip_from=:core_profiles)

    actor_saw = ActorSawteethSource(dd, act.ActorSawteethSource)

    if par.couple_sol
        @assert actor_ped isa ActorPedestal "ActorStationaryPlasma: couple_sol=true requires ActorPedestal"
        @assert act.ActorSOL.model ∉ (:none, :solps_nn) "ActorStationaryPlasma: couple_sol=true requires an ActorSOL model that sets the separatrix Te in dd.edge_profiles.profiles_1d (e.g. :solwidth or :box), not :$(act.ActorSOL.model)"
        if actor_ped.par.model ∈ (:replay, :none)
            @warn "ActorStationaryPlasma: couple_sol=true has no effect with ActorPedestal.model=:$(actor_ped.par.model), which does not use Te_sep"
        end
        actor_sol = ActorSOL(dd, act.ActorSOL, act)
    else
        actor_sol = nothing
    end

    return ActorStationaryPlasma(dd, par, act, actor_tr, actor_ped, actor_src, actor_jt, actor_eq, actor_saw, actor_sol)
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

    (par.verbose && !par.do_plot) && ProgressMeter.ijulia_behavior(:clear)
    was_logging = actor_logging(dd)
    is_logging = was_logging && !(par.verbose && !par.do_plot)
    actor_logging(dd, is_logging)
    prog = ProgressMeter.Progress((par.max_iterations + 1) * (6 + par.couple_sol) + 2; dt=0.0, showspeed=true, enabled=par.verbose && !par.do_plot)
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

            # run sources to get updated current drive
            ProgressMeter.next!(prog; showvalues=progress_ActorStationaryPlasma(total_error, actor, actor.actor_src))
            finalize(step(actor.actor_src))

            # run SOL actor and use its separatrix Te as pedestal boundary condition
            error_te_sep = 0.0
            if actor.actor_sol !== nothing
                ProgressMeter.next!(prog; showvalues=progress_ActorStationaryPlasma(total_error, actor, actor.actor_sol))
                finalize(step(actor.actor_sol))
                Te_sep_before = actor.actor_ped.par.Te_sep
                Ti_sep_before = getproperty(actor.actor_ped.par, :Ti_sep, missing)
                α = isempty(total_error) ? 1.0 : par.sol_relaxation
                sep = couple_separatrix!(actor.actor_ped, dd; relaxation=α, couple_ti=par.couple_sol_ti)
                if !isempty(total_error)
                    error_te_sep = ((sep.Te_sep - Te_sep_before) / Te_sep_before)^2
                    if par.couple_sol_ti && !ismissing(Ti_sep_before)
                        error_te_sep += ((sep.Ti_sep - Ti_sep_before) / Ti_sep_before)^2
                    end
                end
            end

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
            push!(total_error, sqrt(error_jtor + error_pressure + error_te_sep) / 2.0)

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
                if actor.actor_sol !== nothing
                    @printf("Te_sep = %.1f eV\n", actor.actor_ped.par.Te_sep)
                    par.couple_sol_ti && @printf("Ti_sep = %.1f eV\n", actor.actor_ped.par.Ti_sep)
                    @printf("ϵ Tsep = %.4f\n", error_te_sep)
                end
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

"""
    couple_separatrix!(actor_ped::ActorPedestal, dd::IMAS.DD; relaxation::Real=1.0, couple_ti::Bool=false)

Sets `Te_sep` (and `Ti_sep` if `couple_ti`) of `actor_ped` from the separatrix point (`rho_pol_norm=1.0`)
of `dd.edge_profiles.profiles_1d[]`, under-relaxed as `X ← X + relaxation * (X_SOL - X)`.
A `Ti_sep` that was not set before is taken directly from the SOL model.

Returns `(Te_sep=..., Ti_sep=...)` with the new values (`Ti_sep` is `missing` if not set).
"""
function couple_separatrix!(actor_ped::ActorPedestal, dd::IMAS.DD; relaxation::Real=1.0, couple_ti::Bool=false)
    @assert 0.0 < relaxation <= 1.0 "relaxation must be in (0, 1]"
    ep1d = dd.edge_profiles.profiles_1d[]
    @assert !ismissing(ep1d.grid, :rho_pol_norm) && ep1d.grid.rho_pol_norm[1] == 1.0 "SOL model did not populate the separatrix in dd.edge_profiles"
    par = actor_ped.par

    Te_sep_sol = ep1d.electrons.temperature[1]
    @assert isfinite(Te_sep_sol) && Te_sep_sol > 0.0 "SOL model returned invalid Te_sep = $Te_sep_sol"
    par.Te_sep = par.Te_sep + relaxation * (Te_sep_sol - par.Te_sep)

    if couple_ti
        @assert !ismissing(ep1d, :t_i_average) "couple_sol_ti=true but the SOL model did not set the separatrix Ti in dd.edge_profiles (use e.g. ActorSOL.model=:solwidth)"
        Ti_sep_sol = ep1d.t_i_average[1]
        @assert isfinite(Ti_sep_sol) && Ti_sep_sol > 0.0 "SOL model returned invalid Ti_sep = $Ti_sep_sol"
        Ti_sep_old = getproperty(par, :Ti_sep, missing)
        par.Ti_sep = ismissing(Ti_sep_old) ? Ti_sep_sol : Ti_sep_old + relaxation * (Ti_sep_sol - Ti_sep_old)
    end

    return (Te_sep=par.Te_sep, Ti_sep=getproperty(par, :Ti_sep, missing))
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
