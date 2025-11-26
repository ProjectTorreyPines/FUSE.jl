import NLsolve
using LinearAlgebra
import TJLF: InputTJLF

import NonlinearSolve, FixedPointAcceleration

#= ================ =#
#  ActorFluxMatcher  #
#= ================ =#
@actor_parameters_struct ActorFluxMatcher{T} begin
    rho_transport::Entry{AbstractVector{T}} = Entry{AbstractVector{T}}("-", "ρ transport grid"; default=0.25:0.1:0.85)
    evolve_Ti::Switch{Symbol} = Switch{Symbol}([:flux_match, :fixed, :replay], "-", "Ion temperature `:flux_match`, keep `:fixed`, or `:replay` from replay_dd"; default=:flux_match)
    evolve_Te::Switch{Symbol} = Switch{Symbol}([:flux_match, :fixed, :replay], "-", "Electron temperature `:flux_match`, keep `:fixed`, or `:replay` from replay_dd"; default=:flux_match)
    evolve_densities::Entry{Union{AbstractDict,Symbol}} =
        Entry{Union{AbstractDict,Symbol}}(
            "-",
            """
            * `:fixed`: no density evolution
            * `:flux_match:` electron flux-match and ions match ne scale while keeping quasi-neutrality, giving a flat zeff
            * `:zeff:` electron flux-match and ions scaled to keep zeff constant while keeping quasi-neutrality
            * `:replay:` replay densities from replay_dd
            * Dict to specify which species are `:flux_match`, kept `:fixed`, used to enforce `:quasi_neutrality`, scaled to `:match_ne_scale`, or `:replay`""";
            default=:flux_match
        )
    evolve_rotation::Switch{Symbol} = Switch{Symbol}([:flux_match, :fixed, :replay], "-", "Rotation `:flux_match`, keep `:fixed`, or `:replay` from replay_dd"; default=:fixed)
    evolve_pedestal::Entry{Bool} = Entry{Bool}("-", "Evolve the pedestal at each iteration"; default=true)
    evolve_plasma_sources::Entry{Bool} = Entry{Bool}("-", "Update the plasma sources at each iteration"; default=true)
    find_widths::Entry{Bool} = Entry{Bool}("-", "Runs turbulent transport actor TJLF finding widths after first iteration"; default=true)
    max_iterations::Entry{Int} = Entry{Int}("-", "Maximum optimizer iterations"; default=0)
    xtol::Entry{T} = Entry{T}("-", "Tolerance on the solution vector"; default=1E-3, check=x -> @assert x > 0.0 "must be: xtol > 0.0")
    algorithm::Switch{Symbol} =
        Switch{Symbol}(
            [:default, :basic_polyalg, :polyalg, :broyden, :anderson, :simple_trust, :simple_dfsane, :trust, :simple, :old_anderson, :custom, :none],
            "-",
            "Optimizing algorithm used for the flux matching";
            default=:default
        )
    custom_algorithm::Entry{NonlinearSolve.AbstractNonlinearSolveAlgorithm} =
        Entry{NonlinearSolve.AbstractNonlinearSolveAlgorithm}("-", "User-defined custom solver from NonlinearSolve")
    step_size::Entry{T} = Entry{T}(
        "-",
        "Step size for each algorithm iteration (note this has a different meaning for each algorithm)";
        default=1.0,
        check=x -> @assert x > 0.0 "must be: step_size > 0.0"
    )
    Δt::Entry{Float64} = Entry{Float64}("s", "Evolve for Δt (Inf for steady state)"; default=Inf)
    relax::Entry{Float64} = Entry{Float64}("-", "Relaxation on the final solution"; default=1.0, check=x -> @assert 0.0 <= x <= 1.0 "must be: 0.0 <= relax <= 1.0")
    scale_turbulence_law::Switch{Symbol} =
        Switch{Symbol}([:h98, :ds03], "-", "Scale turbulent transport to achieve a desired confinement law: NOTE: scale law evaluated with include_radiation=false")
    scale_turbulence_value::Entry{T} = Entry{T}(
        "-",
        "Scale turbulent transport to achieve a desired confinement value for the `scale_turbulence_law`";
        default=1.0,
        check=x -> @assert x > 0.0 "must be: turbulence_scale_value > 0.0"
    )
    do_plot::Entry{Bool} = act_common_parameters(; do_plot=false)
    verbose::Entry{Bool} = act_common_parameters(; verbose=false)
    show_trace::Entry{Bool} = Entry{Bool}("-", "Show convergence trace of nonlinear solver"; default=false)
end

mutable struct ActorFluxMatcher{D,P} <: CompoundAbstractActor{D,P}
    dd::IMAS.dd{D}
    par::OverrideParameters{P,FUSEparameters__ActorFluxMatcher{P}}
    act::ParametersAllActors{P}
    actor_ct::ActorFluxCalculator{D,P}
    actor_replay::ActorReplay{D,P}
    actor_ped::ActorPedestal{D,P}
    norms::Vector{D}
    error::D
end

"""
    ActorFluxMatcher(dd::IMAS.dd, act::ParametersAllActors; kw...)

Performs self-consistent transport evolution by matching turbulent/neoclassical transport fluxes to source fluxes.

This actor solves the core transport equation ∇·Γ = S by iteratively adjusting plasma profile 
gradients until transport fluxes balance particle/energy sources. The process:

1. **Profile Evolution**: Adjusts temperature/density/rotation profile gradients based on user configuration
2. **Transport Calculation**: Evaluates turbulent and neoclassical transport fluxes using ActorFluxCalculator  
3. **Source Calculation**: Updates plasma heating, particle, and momentum sources
4. **Flux Matching**: Minimizes residual (transport_flux - source_flux) using nonlinear solvers
5. **Pedestal Coupling**: Optionally evolves pedestal conditions during the iteration process

Evolution options per channel:
- `:flux_match`: Evolve profile gradients to match transport and source fluxes  
- `:fixed`: Keep profiles fixed
- `:replay`: Use profiles from experimental data

Advanced features include turbulence scaling to target confinement laws (H98, DS03),
time-dependent evolution with ∂/∂t terms, and various nonlinear solver algorithms
optimized for different transport model characteristics.
"""
function ActorFluxMatcher(dd::IMAS.dd, act::ParametersAllActors; kw...)
    actor = ActorFluxMatcher(dd, act.ActorFluxMatcher, act; kw...)
    step(actor)
    finalize(actor)
    return actor
end

function ActorFluxMatcher(dd::IMAS.dd{D}, par::FUSEparameters__ActorFluxMatcher{P}, act::ParametersAllActors{P}; kw...) where {D<:Real,P<:Real}
    logging_actor_init(ActorFluxMatcher)
    par = OverrideParameters(par; kw...)
    actor_ct = ActorFluxCalculator(dd, act.ActorFluxCalculator, act; par.rho_transport)
    actor_replay = ActorReplay(dd, act.ActorReplay, ActorNoOperation(dd, act.ActorNoOperation))
    actor_ped = ActorPedestal(
        dd,
        act.ActorPedestal,
        act;
        ip_from=:core_profiles,
        βn_from=:core_profiles,
        ne_from=:pulse_schedule,
        zeff_from=:pulse_schedule,
        rho_nml=par.rho_transport[end-1],
        rho_ped=par.rho_transport[end])
    actor = ActorFluxMatcher(dd, par, act, actor_ct, actor_replay, actor_ped, D[], D(Inf))
    actor.actor_replay = ActorReplay(dd, act.ActorReplay, actor)
    return actor
end

"""
    _step(actor::ActorFluxMatcher)

ActorFluxMatcher step
"""
function _step(actor::ActorFluxMatcher{D,P}) where {D<:Real,P<:Real}
    dd = actor.dd
    par = actor.par
    cp1d = dd.core_profiles.profiles_1d[]

    # Apply replay profiles if any evolve options are set to :replay
    evolve_densities = evolve_densities_dictionary(cp1d, par)
    if par.evolve_Te == :replay || par.evolve_Ti == :replay || par.evolve_rotation == :replay ||
                 (!isempty(evolve_densities) && any(evolve == :replay for (_, evolve) in evolve_densities))
        finalize(step(actor.actor_replay))
    end

    IMAS.sources!(dd)
    # freeze current expressions for speed
    IMAS.refreeze!(cp1d, :j_ohmic)
    IMAS.refreeze!(cp1d, :j_non_inductive)
    IMAS.freeze!(cp1d, :rotation_frequency_tor_sonic)

    if !isinf(par.Δt)
        # "∂/∂t" is to account to changes in the profiles that
        # have happened between t-1 and tnow outside of the flux matcher.
        # The flux_macher then temporarily adds a "∂/∂t implicit" term to this to
        # account for changes that occur because of transport predictions.
        IMAS.time_derivative_source!(dd)
    end

    initial_cp1d = cp1d_copy_primary_quantities(cp1d)

    @assert nand(typeof(actor.actor_ct.actor_neoc) <: ActorNoOperation, typeof(actor.actor_ct.actor_turb) <: ActorNoOperation) "Unable to fluxmatch when all transport actors are turned off"

    z_init, profiles_paths, fluxes_paths = pack_z_profiles(cp1d, par)
    z_init_scaled = scale_z_profiles(z_init) # scale z_profiles to get smaller stepping using NLsolve
    N_radii = length(par.rho_transport)
    N_channels = round(Int, length(z_init) / N_radii, RoundDown)

    z_scaled_history = []
    err_history = Vector{Vector{D}}()

    actor.norms = fill(D(NaN), N_channels)

    ProgressMeter.ijulia_behavior(:clear)
    prog = ProgressMeter.ProgressUnknown(; dt=0.1, desc="Calls:", enabled=par.verbose)
    old_logging = actor_logging(dd, false)

    if ismissing(par, :scale_turbulence_law)
        opt_parameters = z_init_scaled
    else
        opt_parameters = [1.0; z_init_scaled]
    end

    autodiff = NonlinearSolve.ADTypes.AutoFiniteDiff()
    if D <: ForwardDiff.Dual
        autodiff = NonlinearSolve.ADTypes.AutoForwardDiff()
    end

    algorithm = if par.algorithm === :default
        if actor.actor_ct.actor_turb.par.model === :TGLFNN
            # combines speed and robustness, but needs smooth derivatives
            :basic_polyalg
        else
            # derivative-free method
            :simple_dfsane
        end
    else
        par.algorithm
    end

    # Different defaults for gradient-based methods
    if par.max_iterations != 0
        max_iterations = abs(par.max_iterations)
    elseif algorithm in (:trust, :simple_trust, :broyden, :polyalg)
        max_iterations = 50
    else
        max_iterations = 500
    end

    res = try
        if algorithm == :none
            res = (zero=opt_parameters,)

        elseif algorithm == :simple
            ftol = Inf # always use xtol condition as in NonlinearSolve.jl
            res = flux_match_simple(actor, opt_parameters, initial_cp1d, z_scaled_history, err_history, max_iterations, ftol, par.xtol, prog)

        else
            # 1. In-place residual
            function f!(F, u, initial_cp1d)
                try
                    F .= flux_match_errors(actor, u, initial_cp1d; z_scaled_history, err_history, prog).errors
                catch e
                    if isa(e, InterruptException)
                        rethrow(e)
                    end
                    F .= Inf
                end
            end

            # 2. Problem definition
            problem = NonlinearSolve.NonlinearProblem(f!, opt_parameters, initial_cp1d)

            # 3. Algorithm selection
            alg = if algorithm == :old_anderson
                NonlinearSolve.NLsolveJL(; method=:anderson, m=4, beta=-par.step_size * 5)

            elseif algorithm == :anderson
                NonlinearSolve.FixedPointAccelerationJL(;
                    algorithm=:Anderson,
                    m=5, # history length
                    dampening=-par.step_size * 0.5,
                    condition_number_threshold=1e8,
                    replace_invalids=:ReplaceVector
                )

            elseif algorithm == :broyden
                NonlinearSolve.Broyden(; autodiff)

            elseif algorithm == :trust
                NonlinearSolve.TrustRegion(; autodiff)

            elseif algorithm == :simple_trust
                NonlinearSolve.SimpleTrustRegion(; autodiff)

            elseif algorithm == :simple_dfsane
                NonlinearSolve.SimpleDFSane()

            elseif algorithm === :basic_polyalg
                NonlinearSolve.NonlinearSolvePolyAlgorithm((NonlinearSolve.Broyden(; autodiff),
                    NonlinearSolve.SimpleTrustRegion(; autodiff)))

            elseif algorithm == :polyalg
                # Default NonlinearSolve algorithm
                NonlinearSolve.FastShortcutNonlinearPolyalg(; autodiff)

            elseif algorithm == :custom
                ismissing(par.custom_algorithm) && error("custom_algorithm must be set to a NonlinearSolve algorithm for algorithm=:custom")
                par.custom_algorithm

            else
                error("Unsupported algorithm: $(algorithm)")
            end

            # 4. Solve with matching tolerances and iteration limits
            # NonlinearSolve abstol is meant to be on u, but actually gets
            #   passed to ftol in NLsolve which is an error on the residual
            # See https://github.com/SciML/NonlinearSolve.jl/issues/593
            abstol = par.xtol

            NonlinearSolve.solve(
                problem, alg;
                abstol,
                maxiters=max_iterations,
                show_trace=Val(par.show_trace),
                store_trace=Val(false),
                verbose=false
            )

            # NonlinearSolve returns the first value if the optimization was not successful
            # but we want it to return the best solution, even if the optimization did not
            # reach the desired tolerance
            if !isempty(err_history) && length(err_history) > 1 && norm(err_history[end]) == norm(err_history[1])
                pop!(err_history)
                pop!(z_scaled_history)
            end

            # Extract the solution vector
            if !isempty(err_history)
                k = argmin(map(norm, err_history))
                res = (zero=z_scaled_history[k],)
            else
                res = (zero=opt_parameters,)
            end
        end

    finally
        actor_logging(dd, old_logging)
    end

    # detect cases where all optimization calls failed
    if algorithm != :none && isempty(err_history)
        flux_match_errors(actor, opt_parameters, initial_cp1d)
        error("FluxMatcher failed")
    end

    # evaluate profiles at the best-matching gradients
    out = flux_match_errors(actor, collect(res.zero), initial_cp1d) # z_profiles for the smallest error iteration

    # statistics
    actor.error = norm(out.errors)
    @ddtime(dd.transport_solver_numerics.convergence.time_step.time = dd.global_time)
    @ddtime(dd.transport_solver_numerics.convergence.time_step.data = actor.error)
    dd.transport_solver_numerics.ids_properties.name = "FluxMatcher"
    ProgressMeter.finish!(prog; showvalues=progress_ActorFluxMatcher(dd, norm(out.errors)))

    # plotting of the channels that have been evolved
    if par.do_plot
        if ismissing(par, :scale_turbulence_law)
            specials = String[]
        else
            specials = String["scale_turbulence_law"]
        end
        N_specials = length(specials)

        p = plot()
        normed_errors = vcat(map(norm, err_history)...)
        plot!(
            normed_errors;
            color=:black,
            lw=2.0,
            yscale=:log10,
            ylabel="Convergence errror",
            xlabel="Iterations",
            label=@sprintf("total [error = %.3e]", (minimum(normed_errors)))
        )
        for sp in 1:N_specials
            plot!(vcat(map(x -> errors_by_channel(abs.(x), N_radii, N_channels, N_specials).specials, err_history)...)[:, sp]; label=specials[sp])
        end
        for (ch, profiles_path) in enumerate(profiles_paths)
            label = profiles_title(cp1d, profiles_path)
            plot!(vcat(map(x -> errors_by_channel(x, N_radii, N_channels, N_specials).radii_channels, err_history)...)[:, ch]; label)
        end
        display(p)

        channels_evolution = transpose(hcat(map(z -> collect(unscale_z_profiles(z)), z_scaled_history)...))
        data = reshape(channels_evolution, (length(err_history), N_radii, N_channels))
        p = plot()
        for (ch, profiles_path) in enumerate(profiles_paths)
            label = profiles_title(cp1d, profiles_path)
            for kr in 1:length(par.rho_transport)
                plot!(data[:, kr, ch]; ylabel="Inverse scale length [m⁻¹]", xlabel="Iterations", primary=kr == 1, lw=kr, label)
            end
        end
        display(p)

        p = plot(; layout=(N_channels, 2), size=(1000, 300 * N_channels))
        total_flux1d = IMAS.total_fluxes(dd.core_transport, cp1d, par.rho_transport; time0=dd.global_time)
        total_source1d = IMAS.total_sources(dd.core_sources, cp1d; time0=dd.global_time)
        model_type = IMAS.name_2_index(dd.core_transport.model)
        for (ch, (profiles_path, fluxes_path)) in enumerate(zip(profiles_paths, fluxes_paths))
            plot!(
                IMAS.goto(total_source1d, fluxes_path[1:end-1]),
                Val(fluxes_path[end]);
                flux=true,
                subplot=2 * ch - 1,
                name="total sources",
                color=:blue,
                linewidth=2,
                show_zeros=true
            )
            for model in dd.core_transport.model
                if !ismissing(IMAS.goto(model.profiles_1d[], fluxes_path), :flux)
                    if model.identifier.index == model_type[:anomalous]
                        plot_opts = Dict(:markershape => :diamond, :markerstrokewidth => 0.5, :linewidth => 0, :color => :orange)
                    elseif model.identifier.index == model_type[:neoclassical]
                        plot_opts = Dict(:markershape => :cross, :linewidth => 0, :color => :purple)
                    else
                        continue
                    end
                    plot!(IMAS.goto(model.profiles_1d[], fluxes_path), Val(:flux); subplot=2 * ch - 1, plot_opts...)
                end
            end
            plot!(IMAS.goto(total_flux1d, fluxes_path), Val(:flux); subplot=2 * ch - 1, color=:red, label="total transport", linewidth=2)

            title = profiles_title(cp1d, profiles_path)
            plot!(IMAS.goto(initial_cp1d, profiles_path[1:end-1]), profiles_path[end]; subplot=2 * ch, label="before", linestyle=:dash, color=:black)
            plot!(IMAS.goto(cp1d, profiles_path[1:end-1]), profiles_path[end]; subplot=2 * ch, label="after", title)
            if profiles_path[end] != :momentum_tor
                plot!(; subplot=2 * ch, ylim=[0, Inf])
            end
        end
        display(p)
    end

    IMAS.unfreeze!(cp1d, :j_ohmic)
    IMAS.unfreeze!(cp1d, :j_non_inductive)

    # final relaxation of profiles
    if par.relax < 1.0
        for profiles_path in profiles_paths
            field = profiles_path[end]
            ids1 = IMAS.goto(cp1d, profiles_path[1:end-1])
            ids2 = IMAS.goto(initial_cp1d, profiles_path[1:end-1])
            if !ismissing(ids1, field) && !ismissing(ids2, field)
                value1 = getproperty(ids1, field)
                value2 = getproperty(ids2, field)
                value = par.relax * value1 + (1.0 - par.relax) * value2
                setproperty!(ids1, field, value)
            end
        end

        # transfer t_i_average to individual ions temperature and restore t_i_average as expression
        for ion in cp1d.ion
            ion.temperature = cp1d.t_i_average
        end
        IMAS.unfreeze!(cp1d, :t_i_average)

        # refresh sources with relatex profiles
        IMAS.sources!(dd)

        # Ensure quasi neutrality if densities are evolved
        # NOTE: check_evolve_densities() takes care of doing proper error handling for user inputs
        for (species, evolve) in evolve_densities
            if evolve == :quasi_neutrality
                IMAS.enforce_quasi_neutrality!(cp1d, species)
                break
            end
        end
    end

    # remove the temporary `∂/∂t implicit` term and update the full `∂/∂t` term
    if !isinf(par.Δt)
        deleteat!(dd.core_sources.source, :time_derivative, "identifier.name" => "∂/∂t implicit")
        IMAS.time_derivative_source!(dd)
    end

    # free pressures expressions
    IMAS.unfreeze!(cp1d.electrons, :pressure_thermal)
    IMAS.unfreeze!(cp1d.electrons, :pressure)
    for ion in cp1d.ion
        IMAS.unfreeze!(ion, :rotation_frequency_tor)
        IMAS.unfreeze!(ion, :pressure_thermal)
        IMAS.unfreeze!(ion, :pressure)
    end
    for field in [:pressure_ion_total, :pressure_thermal, :pressure]
        IMAS.unfreeze!(cp1d, field)
    end

    return actor
end

function profiles_title(cp1d, profiles_path)
    if profiles_path[1] == :electrons
        if profiles_path[2] == :temperature
            return "Electron temperature"
        else
            return "Electron density"
        end
    elseif profiles_path[1] == :rotation_frequency_tor_sonic
        return "Rotation"
    elseif profiles_path[1] == :t_i_average
        return "Ions temperature"
    elseif profiles_path[1] == :ion
        ion = cp1d.ion[profiles_path[2]]
        return "$(ion.label) density"
    else
        error("$(profiles_path) is not a recognized core_profiles.profiles_1d path")
    end
end

function errors_by_channel(errors::Vector{T}, N_radii::Int, N_channels::Int, N_specials::Int) where {T<:Real}
    return (specials=errors[1+N_specials], radii_channels=mapslices(norm, reshape(errors[1+N_specials:end], (N_radii, N_channels)); dims=1))
end

"""
    scale_z_profiles(z_profiles)

scale z_profiles to get smaller stepping in nlsolve
"""
function scale_z_profiles(z_profiles)
    return z_profiles .* 100
end

function unscale_z_profiles(z_profiles_scaled)
    return z_profiles_scaled ./ 100
end

"""
    flux_match_errors(
        actor::ActorFluxMatcher{D,P},
        opt_parameters::Vector{<:Real},
        initial_cp1d::IMAS.core_profiles__profiles_1d;
        z_scaled_history::Vector=[],
        err_history::Vector{Vector{D}}=Vector{Vector{D}}(),
        prog::Any=nothing) where {D<:Real,P<:Real}

Update the profiles, evaluates neoclassical and turbulent fluxes, sources (ie target fluxes), and returns named tuple with (targets, fluxes, errors)

NOTE: flux matching is done in physical units
"""
function flux_match_errors(
    actor::ActorFluxMatcher{D,P},
    opt_parameters::Vector{<:Real},
    initial_cp1d::IMAS.core_profiles__profiles_1d;
    z_scaled_history::Vector=[],
    err_history::Vector{Vector{D}}=Vector{Vector{D}}(),
    prog::Any=nothing) where {D<:Real,P<:Real}

    dd = actor.dd
    par = actor.par
    cp1d = dd.core_profiles.profiles_1d[]

    if ismissing(par, :scale_turbulence_law)
        z_profiles_scaled = deepcopy(opt_parameters)
    else
        turbulence_scale = opt_parameters[1]
        z_profiles_scaled = deepcopy(opt_parameters[2:end])
    end

    # unscale z_profiles
    z_profiles = unscale_z_profiles(z_profiles_scaled)

    # restore profiles at initial conditions
    cp1d_copy_primary_quantities!(cp1d, initial_cp1d)

    # evolve pedestal
    if par.evolve_pedestal
        # modify cp1d with new z_profiles
        unpack_z_profiles(cp1d, par, z_profiles)
        # run pedestal
        actor.actor_ped.par.βn_from = :core_profiles
        finalize(step(actor.actor_ped))
    end

    # modify cp1d with new z_profiles
    unpack_z_profiles(cp1d, par, z_profiles)

    # evaluate sources (ie. target fluxes)
    if par.evolve_plasma_sources
        IMAS.sources!(dd; bootstrap=false)
    end
    if par.Δt < Inf
        IMAS.time_derivative_source!(dd, initial_cp1d, par.Δt; name="∂/∂t implicit")
    end

    # handle fixed widths in TJLF
    if actor.actor_ct.par.turbulence_model == :TGLF && actor.actor_ct.actor_turb.par.model == :TJLF && !par.find_widths && !isempty(err_history)
        for input_tglf in actor.actor_ct.actor_turb.input_tglfs
            input_tglf.FIND_WIDTH = false
        end
    end

    # evaluate neoclassical + turbulent fluxes
    finalize(step(actor.actor_ct))

    # scale turbulent fluxes, if par.scale_turbulence_law is set
    if !ismissing(par, :scale_turbulence_law)
        m1d = dd.core_transport.model[:anomalous].profiles_1d[]
        for leaf in AbstractTrees.Leaves(m1d)
            if leaf.field == :flux
                getproperty(leaf.ids, leaf.field) .*= turbulence_scale
            end
        end
    end

    # get transport fluxes and sources
    fluxes = flux_match_fluxes(dd, par)
    targets = flux_match_targets(dd, par)

    cp_gridpoints = [argmin_abs(cp1d.grid.rho_tor_norm, rho_x) for rho_x in par.rho_transport]
    surface0 = cp1d.grid.surface[cp_gridpoints] ./ cp1d.grid.surface[end]

    # Evaluate the flux_matching errors
    nrho = length(par.rho_transport)
    errors = similar(fluxes)
    for (inorm, norm0) in enumerate(actor.norms)
        index = (inorm-1)*nrho+1:inorm*nrho
        if isnan(norm0)
            # this sets the norm to be the based on the average of the fluxes and targets
            # actor.norm is only used if the targets are zero
            actor.norms[inorm] = norm0 = (norm(fluxes[index] .* surface0) + norm(targets[index] .* surface0)) / 2.0
        end
        errors[index] .= @views (targets[index] .- fluxes[index]) ./ norm0 .* surface0
    end

    # update progress meter
    if prog !== nothing
        ProgressMeter.next!(prog; showvalues=progress_ActorFluxMatcher(dd, norm(errors)))
    end

    # add error toward achieving desired scaling law value
    if !ismissing(par, :scale_turbulence_law)
        tau_th = IMAS.tau_e_thermal(dd; include_radiation=false)
        if par.scale_turbulence_law == :h98
            tauH = IMAS.tau_e_h98(dd; include_radiation=false)
        elseif par.scale_turbulence_law == :ds03
            tauH = IMAS.tau_e_ds03(dd; include_radiation=false)
        else
            error("act.ActorFluxMatcher.scale_turbulence_law=$(par.scale_turbulence_law) not recognized: valid options are :h98 or :ds03")
        end
        H_value = tau_th / tauH
        H_target = par.scale_turbulence_value
        errors = [(H_value - H_target) / H_target * nrho; errors]
    end

    # update history
    push!(z_scaled_history, z_profiles_scaled)
    push!(err_history, errors)

    return (targets=targets, fluxes=fluxes, errors=errors)
end

function norm_transformation(norm_source::Vector{T}, norm_transp::Vector{T}) where {T<:Real}
    return (sum(abs.(norm_source)) .+ sum(abs.(norm_transp))) / length(norm_source) / 2.0
end

"""
    flux_match_targets(dd::IMAS.dd, par::OverrideParameters{P,FUSEparameters__ActorFluxMatcher{P}}) where {P<:Real}

Evaluates the flux_matching targets for the :flux_match species and channels

NOTE: flux matching is done in physical units
"""
function flux_match_targets(dd::IMAS.dd{D}, par::OverrideParameters{P,FUSEparameters__ActorFluxMatcher{P}}) where {D<:Real,P<:Real}
    cp1d = dd.core_profiles.profiles_1d[]

    total_source = resize!(dd.core_sources.source, :total; wipe=false)
    total_source1d = resize!(total_source.profiles_1d; wipe=false)
    IMAS.total_sources!(total_source1d, dd.core_sources, cp1d; time0=dd.global_time, fields=[:total_ion_power_inside, :power_inside, :particles_inside, :torque_tor_inside])
    cs_gridpoints = [argmin_abs(total_source1d.grid.rho_tor_norm, rho_x) for rho_x in par.rho_transport]

    targets = D[]

    if par.evolve_Te == :flux_match
        target = total_source1d.electrons.power_inside[cs_gridpoints] ./ total_source1d.grid.surface[cs_gridpoints]
        append!(targets, target)
    end

    if par.evolve_Ti == :flux_match
        target = total_source1d.total_ion_power_inside[cs_gridpoints] ./ total_source1d.grid.surface[cs_gridpoints]
        append!(targets, target)
    end

    if par.evolve_rotation == :flux_match
        target = total_source1d.torque_tor_inside[cs_gridpoints] ./ total_source1d.grid.surface[cs_gridpoints]
        append!(targets, target)
    end

    evolve_densities = evolve_densities_dictionary(cp1d, par)
    if !isempty(evolve_densities)
        if evolve_densities[:electrons] == :flux_match
            target = total_source1d.electrons.particles_inside[cs_gridpoints] ./ total_source1d.grid.surface[cs_gridpoints]
            append!(targets, target)
        end
        for (k, ion) in enumerate(cp1d.ion)
            if evolve_densities[Symbol(ion.label)] == :flux_match
                target = total_source1d.ion[k].particles_inside[cs_gridpoints] ./ total_source1d.grid.surface[cs_gridpoints]
                append!(targets, target)
            end
        end
    end

    return targets
end

"""
    flux_match_fluxes(dd::IMAS.dd{D}, par::OverrideParameters{P,FUSEparameters__ActorFluxMatcher{P}}) where {D<:Real, P<:Real}

Evaluates the flux_matching fluxes for the :flux_match species and channels

NOTE: flux matching is done in physical units
"""
function flux_match_fluxes(dd::IMAS.dd{D}, par::OverrideParameters{P,FUSEparameters__ActorFluxMatcher{P}}) where {D<:Real,P<:Real}
    cp1d = dd.core_profiles.profiles_1d[]

    total_flux = resize!(dd.core_transport.model, :combined; wipe=false)
    total_flux1d = resize!(total_flux.profiles_1d; wipe=false)
    IMAS.total_fluxes!(total_flux1d, dd.core_transport, cp1d, par.rho_transport; time0=dd.global_time)

    fluxes = D[]

    if par.evolve_Te == :flux_match
        flux = total_flux1d.electrons.energy.flux
        check_output_fluxes(flux, "electrons.energy")
        append!(fluxes, flux)
    end

    if par.evolve_Ti == :flux_match
        flux = total_flux1d.total_ion_energy.flux
        check_output_fluxes(flux, "total_ion_energy")
        append!(fluxes, flux)
    end

    if par.evolve_rotation == :flux_match
        flux = total_flux1d.momentum_tor.flux
        check_output_fluxes(flux, "momentum_tor")
        append!(fluxes, flux)
    end

    evolve_densities = evolve_densities_dictionary(cp1d, par)
    if !isempty(evolve_densities)
        if evolve_densities[:electrons] == :flux_match
            flux = total_flux1d.electrons.particles.flux
            check_output_fluxes(flux, "electrons.particles")
            append!(fluxes, flux)
        end
        for (k, ion) in enumerate(cp1d.ion)
            if evolve_densities[Symbol(ion.label)] == :flux_match
                flux = total_flux1d.ion[k].particles.flux
                check_output_fluxes(flux, "ion[$k].particles")
                append!(fluxes, flux)
            end
        end
    end

    return fluxes
end

"""
    flux_match_simple(
        actor::ActorFluxMatcher{D,P},
        opt_parameters::Vector{<:Real},
        initial_cp1d::IMAS.core_profiles__profiles_1d,
        z_scaled_history::Vector,
        err_history::Vector{Vector{D}},
        max_iterations::Int,
        ftol::Float64,
        xtol::Float64,
        prog::Any) where {D<:Real,P<:Real}

Updates zprofiles based on TGYRO simple algorithm
"""
function flux_match_simple(
    actor::ActorFluxMatcher{D,P},
    opt_parameters::Vector{<:Real},
    initial_cp1d::IMAS.core_profiles__profiles_1d,
    z_scaled_history::Vector,
    err_history::Vector{Vector{D}},
    max_iterations::Int,
    ftol::Float64,
    xtol::Float64,
    prog::Any) where {D<:Real,P<:Real}

    par = actor.par

    i = 0

    if ismissing(par, :scale_turbulence_law)
        z_init_scaled = opt_parameters
    else
        turbulence_scale = opt_parameters[1]
        @assert turbulence_scale == 1.0
        z_init_scaled = @views opt_parameters[2:end]
    end

    zprofiles_old = unscale_z_profiles(z_init_scaled)
    targets, fluxes, errors = flux_match_errors(actor, opt_parameters, initial_cp1d; z_scaled_history, err_history, prog)
    ferror = norm(errors)
    xerror = Inf
    step_size = par.step_size
    while (ferror > ftol) || (xerror .> xtol)
        i += 1
        if (i > max_iterations)
            # if max_iterations > 0
            #     @info "Unable to flux-match within $(max_iterations) iterations (aerr) = $(ferror) (ftol=$ftol) (xerr) = $(xerror) (xtol = $xtol)"
            # end
            break
        end

        zprofiles = zprofiles_old .* (1.0 .+ step_size * 0.1 .* (targets .- fluxes) ./ sqrt.(1.0 .+ fluxes .^ 2 + targets .^ 2))
        if ismissing(par, :scale_turbulence_law)
            targets, fluxes, errors = flux_match_errors(actor, scale_z_profiles(zprofiles), initial_cp1d; z_scaled_history, err_history, prog)
        else
            turbulence_scale += errors[1] * step_size
            targets, fluxes, errors =
                flux_match_errors(actor, [turbulence_scale; scale_z_profiles(zprofiles)], initial_cp1d; z_scaled_history, err_history, prog)
        end
        xerror = maximum(abs.(zprofiles .- zprofiles_old)) / step_size
        ferror = norm(errors)
        zprofiles_old = zprofiles
    end

    if ismissing(par, :scale_turbulence_law)
        return (zero=z_scaled_history[argmin(map(norm, err_history))],)
    else
        return (zero=[turbulence_scale; z_scaled_history[argmin(map(norm, err_history))]],)
    end
end

function progress_ActorFluxMatcher(dd::IMAS.dd, error::Real)
    cp1d = dd.core_profiles.profiles_1d[]
    out = Tuple{String,Float64}[]
    push!(out, ("         error", IMAS.force_float64(error)))
    pfus = IMAS.fusion_power(cp1d)
    if pfus > 1E3
        push!(out, ("  Pfusion [MW]", pfus / 1E6))
    end
    push!(out, ("     Ti0 [keV]", IMAS.force_float64(cp1d.t_i_average[1] / 1E3)))
    push!(out, ("     Te0 [keV]", IMAS.force_float64(cp1d.electrons.temperature[1] / 1E3)))
    push!(out, ("ne0 [10²⁰ m⁻³]", IMAS.force_float64(cp1d.electrons.density_thermal[1] / 1E20)))
    return out
end

function evolve_densities_dictionary(cp1d::IMAS.core_profiles__profiles_1d, par::OverrideParameters{P,FUSEparameters__ActorFluxMatcher{P}}) where {P<:Real}
    if par.evolve_densities == :fixed
        return setup_density_evolution_fixed(cp1d)
    elseif par.evolve_densities == :flux_match
        return setup_density_evolution_electron_flux_match_rest_ne_scale(cp1d)
    elseif par.evolve_densities == :zeff
        return setup_density_evolution_electron_flux_match_fixed_zeff(cp1d)
    elseif par.evolve_densities == :replay
        return setup_density_evolution_replay(cp1d)
    elseif typeof(par.evolve_densities) <: AbstractDict{Symbol,Symbol}
        for (k, v) in par.evolve_densities
            @assert v in (:quasi_neutrality, :match_ne_scale, :fixed, :zeff, :flux_match, :replay) "evolve_density `:$(k) => :$(v)` is not allowed. Choose one of [:quasi_neutrality, :match_ne_scale, :fixed, :zeff, :flux_match, :replay]"
        end
        return par.evolve_densities
    else
        error(
            "act.ActorFluxMatcher.evolve_densities cannot be `$(repr(par.evolve_densities))`. Use `:fixed`, `:flux_match`, `:replay` or a dictionary specifiying [:quasi_neutrality, :match_ne_scale, :fixed, :zeff, :flux_match, :replay] for each specie"
        )
    end
end

"""
    pack_z_profiles(cp1d::IMAS.core_profiles__profiles_1d{D}, par::OverrideParameters{P,FUSEparameters__ActorFluxMatcher{P}}) where {D<:Real, P<:Real}

Packs the z_profiles based on evolution parameters

NOTE: the order for packing and unpacking is always: [Te, Ti, Rotation, ne, nis...]
"""
function pack_z_profiles(cp1d::IMAS.core_profiles__profiles_1d{D}, par::OverrideParameters{P,FUSEparameters__ActorFluxMatcher{P}}) where {D<:Real,P<:Real}
    cp_gridpoints = [argmin_abs(cp1d.grid.rho_tor_norm, rho_x) for rho_x in par.rho_transport]

    z_profiles = D[]
    profiles_paths = []
    fluxes_paths = []

    if par.evolve_Te == :flux_match
        z_Te = IMAS.calc_z(cp1d.grid.rho_tor_norm, cp1d.electrons.temperature, :backward)[cp_gridpoints]
        append!(z_profiles, z_Te)
        push!(profiles_paths, (:electrons, :temperature))
        push!(fluxes_paths, (:electrons, :energy))
    end

    if par.evolve_Ti == :flux_match
        z_Ti = IMAS.calc_z(cp1d.grid.rho_tor_norm, cp1d.t_i_average, :backward)[cp_gridpoints]
        append!(z_profiles, z_Ti)
        push!(profiles_paths, (:t_i_average,))
        push!(fluxes_paths, (:total_ion_energy,))
    end

    if par.evolve_rotation == :flux_match
        max_rho_idx = IMAS.argmin_abs(cp1d.grid.rho_tor_norm, maximum(par.rho_transport))
        rot_check = cp1d.rotation_frequency_tor_sonic[1:max_rho_idx]
        if !(any(x -> x >= 0, rot_check) && any(x -> x < 0, rot_check))
            z_rot = IMAS.calc_z(cp1d.grid.rho_tor_norm, cp1d.rotation_frequency_tor_sonic, :backward)[cp_gridpoints]
            append!(z_profiles, z_rot)
        else
            # Use TGYRO approach: evolve normalized rotation shear f_rot = (dω/dr) / w0_norm
            w0_norm = calculate_w0_norm(cp1d.electrons.temperature[1])

            # Calculate rotation shear dω/dr
            dw_dr = IMAS.gradient(cp1d.grid.rho_tor_norm, cp1d.rotation_frequency_tor_sonic)[cp_gridpoints]

            # Convert to normalized rotation shear f_rot
            f_rot = dw_dr ./ w0_norm

            append!(z_profiles, f_rot)
        end
        push!(profiles_paths, (:rotation_frequency_tor_sonic,))
        push!(fluxes_paths, (:momentum_tor,))
    end

    evolve_densities = evolve_densities_dictionary(cp1d, par)
    if !isempty(evolve_densities)
        check_evolve_densities(cp1d, evolve_densities)
        if evolve_densities[:electrons] == :flux_match
            z_ne = IMAS.calc_z(cp1d.grid.rho_tor_norm, cp1d.electrons.density_thermal, :backward)[cp_gridpoints]
            append!(z_profiles, z_ne)
            push!(profiles_paths, (:electrons, :density_thermal))
            push!(fluxes_paths, (:electrons, :particles))
        end
        for (k, ion) in enumerate(cp1d.ion)
            if evolve_densities[Symbol(ion.label)] == :flux_match
                z_ni = IMAS.calc_z(cp1d.grid.rho_tor_norm, ion.density_thermal, :backward)[cp_gridpoints]
                append!(z_profiles, z_ni)
                push!(profiles_paths, (:ion, k, :density_thermal))
                push!(fluxes_paths, [:ion, k, :particles])
            end
        end
    end

    # note: we update the local actor.par.evolve_densities
    #       to make it easier to show what's going on
    par.evolve_densities = evolve_densities

    return (z_profiles=z_profiles, profiles_paths=profiles_paths, fluxes_paths=fluxes_paths)
end

"""
    unpack_z_profiles(
        cp1d::IMAS.core_profiles__profiles_1d,
        par::OverrideParameters{P,FUSEparameters__ActorFluxMatcher{P}},
        z_profiles::AbstractVector{<:Real}) where {P<:Real}

Unpacks z_profiles based on evolution parameters

NOTE: The order for packing and unpacking is always: [Ti, Te, Rotation, ne, nis...]
"""
function unpack_z_profiles(
    cp1d::IMAS.core_profiles__profiles_1d,
    par::OverrideParameters{P,FUSEparameters__ActorFluxMatcher{P}},
    z_profiles::AbstractVector{<:Real}) where {P<:Real}

    # bound range of accepted z_profiles to avoid issues during optimization
    z_max = 10.0
    z_profiles .= min.(max.(z_profiles, -z_max), z_max)

    cp_gridpoints = [argmin_abs(cp1d.grid.rho_tor_norm, rho_x) for rho_x in par.rho_transport]
    cp_rho_transport = cp1d.grid.rho_tor_norm[cp_gridpoints]

    N = length(par.rho_transport)
    counter = 0

    evolve_densities = evolve_densities_dictionary(cp1d, par)
    if any(evolve_densities[Symbol(ion.label)] == :zeff for ion in cp1d.ion)
        old_zeff = cp1d.zeff
    end

    if par.evolve_Te == :flux_match
        cp1d.electrons.temperature = IMAS.profile_from_z_transport(cp1d.electrons.temperature, cp1d.grid.rho_tor_norm, cp_rho_transport, z_profiles[counter+1:counter+N])
        counter += N
    end

    if par.evolve_Ti == :flux_match
        Ti_new = IMAS.profile_from_z_transport(cp1d.t_i_average, cp1d.grid.rho_tor_norm, cp_rho_transport, z_profiles[counter+1:counter+N])
        counter += N
        for ion in cp1d.ion
            ion.temperature = Ti_new
        end
    end

    if par.evolve_rotation == :flux_match
        max_rho_idx = IMAS.argmin_abs(cp1d.grid.rho_tor_norm, maximum(par.rho_transport))
        rot_check = cp1d.rotation_frequency_tor_sonic[1:max_rho_idx]
        if !(any(x -> x >= 0, rot_check) && any(x -> x < 0, rot_check))
            cp1d.rotation_frequency_tor_sonic = IMAS.profile_from_z_transport(cp1d.rotation_frequency_tor_sonic .+ 1, cp1d.grid.rho_tor_norm, cp_rho_transport, z_profiles[counter+1:counter+N])
        else
            # Use TGYRO approach: convert normalized rotation shear back to rotation frequency
            w0_norm = calculate_w0_norm(cp1d.electrons.temperature[1])

            # Get normalized rotation shear f_rot from z_profiles
            f_rot_evolved = z_profiles[counter+1:counter+N]

            # Convert to rotation shear: dω/dr = f_rot * w0_norm
            dw_dr_evolved = f_rot_evolved .* w0_norm

            # Use the new profile_from_rotation_shear_transport function
            cp1d.rotation_frequency_tor_sonic = IMAS.profile_from_rotation_shear_transport(
                cp1d.rotation_frequency_tor_sonic,
                cp1d.grid.rho_tor_norm,
                cp_rho_transport,
                dw_dr_evolved
            )
        end
        counter += N
    end

    if !isempty(evolve_densities)
        if evolve_densities[:electrons] == :flux_match
            z_ne = z_profiles[counter+1:counter+N]
            cp1d.electrons.density_thermal = IMAS.profile_from_z_transport(cp1d.electrons.density_thermal, cp1d.grid.rho_tor_norm, cp_rho_transport, z_ne)
            counter += N
        else
            z_ne = IMAS.calc_z(cp1d.grid.rho_tor_norm, cp1d.electrons.density_thermal, :backward)[cp_gridpoints]
        end
        for ion in cp1d.ion
            if evolve_densities[Symbol(ion.label)] == :zeff
                IMAS.scale_ion_densities_to_target_zeff!(cp1d, old_zeff)
                break
            elseif evolve_densities[Symbol(ion.label)] == :flux_match
                ion.density_thermal = IMAS.profile_from_z_transport(ion.density_thermal, cp1d.grid.rho_tor_norm, cp_rho_transport, z_profiles[counter+1:counter+N])
                counter += N
            elseif evolve_densities[Symbol(ion.label)] == :match_ne_scale
                ion.density_thermal = IMAS.profile_from_z_transport(ion.density_thermal, cp1d.grid.rho_tor_norm, cp_rho_transport, z_ne)
            end
        end
    end

    # Ensure quasi neutrality if densities are evolved
    # NOTE: check_evolve_densities() takes care of doing proper error handling for user inputs
    for (species, evolve) in evolve_densities
        if evolve == :quasi_neutrality
            IMAS.enforce_quasi_neutrality!(cp1d, species)
            break
        end
    end

    # re-freeze pressures expressions
    zero_value = zero(cp1d.grid.rho_tor_norm)
    IMAS.refreeze!(cp1d.electrons, :pressure_thermal, zero_value)
    IMAS.refreeze!(cp1d.electrons, :pressure, zero_value)
    for ion in cp1d.ion
        IMAS.refreeze!(ion, :pressure_thermal, zero_value)
        IMAS.refreeze!(ion, :pressure, zero_value)
    end
    for field in [:pressure_ion_total, :pressure_thermal, :pressure]
        IMAS.refreeze!(cp1d, field, zero_value)
    end

    return cp1d
end

"""
    check_evolve_densities(cp1d::IMAS.core_profiles__profiles_1d, evolve_densities::Dict)

Checks if the evolve_densities dictionary makes sense and return sensible errors if this is not the case
"""
function check_evolve_densities(cp1d::IMAS.core_profiles__profiles_1d, evolve_densities::AbstractDict)
    dd_species =
        Symbol[specie.name for specie in IMAS.species(cp1d; only_electrons_ions=:all, only_thermal_fast=:all, return_zero_densities=true)]

    # Check if evolve_densities contains all of dd thermal species
    @assert sort!([specie for (specie, evolve) in evolve_densities]) == sort!(dd_species) "Mismatch: dd species $(sort!(dd_species)) VS evolve_densities species : $(sort!(collect(keys(evolve_densities))))"

    # Check that either all species are fixed, or there is 1 quasi_neutrality specie when evolving densities
    if any(evolve == :zeff for (specie, evolve) in evolve_densities if specie != :electrons)
        txt = "When flux_matching densities, either none or all ion species must be :zeff"
        @assert all(evolve == :zeff for (specie, evolve) in evolve_densities if specie != :electrons)
    elseif all(evolve in (:fixed, :replay) for (specie, evolve) in evolve_densities if evolve != :quasi_neutrality)
        txt = "When using fixed/replay densities, no more than one species can be set to :quasi_neutrality"
        @assert length([specie for (specie, evolve) in evolve_densities if evolve == :quasi_neutrality]) <= 1 txt
    elseif any(evolve == :flux_match for (specie, evolve) in evolve_densities)
        txt = "When flux_matching densities, one an only one species must be set to :quasi_neutrality"
        @assert length([specie for (specie, evolve) in evolve_densities if evolve == :quasi_neutrality]) == 1 txt
    else
        error("Invalid evolve_densities = $(evolve_densities)")
    end
end

"""
    evolve_densities_dict_creation(flux_match_species::Vector, fixed_species::Vector, match_ne_scale_species::Vector, replay_species::Vector; quasi_neutrality_specie::Union{Symbol,Bool}=false)

Create the density_evolution dict based on input vectors: flux_match_species, fixed_species, match_ne_scale_species, replay_species, quasi_neutrality_specie
"""
function evolve_densities_dict_creation(
    flux_match_species::Vector;
    fixed_species::Vector{Symbol}=Symbol[],
    match_ne_scale_species::Vector{Symbol}=Symbol[],
    replay_species::Vector{Symbol}=Symbol[],
    quasi_neutrality_specie::Union{Symbol,Bool}=false)

    parse_list = vcat(
        [[sym, :flux_match] for sym in flux_match_species],
        [[sym, :match_ne_scale] for sym in match_ne_scale_species],
        [[sym, :fixed] for sym in fixed_species],
        [[sym, :replay] for sym in replay_species]
    )
    if typeof(quasi_neutrality_specie) <: Symbol
        parse_list = vcat(parse_list, [[quasi_neutrality_specie, :quasi_neutrality]])
    end

    return Dict(sym => evolve for (sym, evolve) in parse_list)
end


"""
    setup_density_evolution_electron_flux_match_rest_ne_scale(cp1d::IMAS.core_profiles__profiles_1d)

Sets up the evolve_density dict to evolve only ne and keep the rest matching the ne_scale lengths
"""
function setup_density_evolution_electron_flux_match_rest_ne_scale(cp1d::IMAS.core_profiles__profiles_1d)
    dd_thermal = Symbol[specie.name for specie in IMAS.species(cp1d; only_electrons_ions=:ions, only_thermal_fast=:thermal, return_zero_densities=true)]
    dd_fast = Symbol[specie.name for specie in IMAS.species(cp1d; only_electrons_ions=:all, only_thermal_fast=:fast, return_zero_densities=true)]
    quasi_neutrality_specie = :D
    if :DT ∈ dd_thermal
        quasi_neutrality_specie = :DT
    end

    return evolve_densities_dict_creation([:electrons]; fixed_species=dd_fast, match_ne_scale_species=dd_thermal, quasi_neutrality_specie)
end

"""
    setup_density_evolution_electron_flux_match_impurities_fixed(cp1d::IMAS.core_profiles__profiles_1d)

Sets up the evolve_density dict to evolve only ne, quasi_neutrality on main ion (:D or :DT) and fix the rest
"""
function setup_density_evolution_electron_flux_match_impurities_fixed(cp1d::IMAS.core_profiles__profiles_1d)
    dd_thermal = Symbol[specie.name for specie in IMAS.species(cp1d; only_electrons_ions=:ions, only_thermal_fast=:thermal, return_zero_densities=true)]
    dd_fast = Symbol[specie.name for specie in IMAS.species(cp1d; only_electrons_ions=:all, only_thermal_fast=:fast, return_zero_densities=true)]
    if :DT ∈ dd_thermal
        quasi_neutrality_specie = :DT
    else
        quasi_neutrality_specie = :D
    end
    return evolve_densities_dict_creation([:electrons]; fixed_species=[dd_thermal..., dd_fast...], quasi_neutrality_specie)
end

"""
    setup_density_evolution_fixed(cp1d::IMAS.core_profiles__profiles_1d)

Sets up the evolve_density dict to keep all species fixed
"""
function setup_density_evolution_fixed(cp1d::IMAS.core_profiles__profiles_1d)
    all_species = [specie.name for specie in IMAS.species(cp1d; return_zero_densities=true)]
    return evolve_densities_dict_creation(Symbol[]; fixed_species=all_species, quasi_neutrality_specie=false)
end

"""
    setup_density_evolution_electron_flux_match_fixed_zeff(cp1d::IMAS.core_profiles__profiles_1d)

Sets up the evolve_density dict to evolve only ne, and scale ion species in such a way to keep Zeff fixed
"""
function setup_density_evolution_electron_flux_match_fixed_zeff(cp1d::IMAS.core_profiles__profiles_1d)
    evolve_density = Dict{Symbol,Symbol}()
    evolve_density[:electrons] = :flux_match
    for ion_name in IMAS.species(cp1d; only_electrons_ions=:ions, return_zero_densities=true)
        evolve_density[ion_name.name] = :zeff
    end
    return evolve_density
end

"""
    setup_density_evolution_replay(cp1d::IMAS.core_profiles__profiles_1d)

Sets up the evolve_density dict to replay all species from replay_dd
"""
function setup_density_evolution_replay(cp1d::IMAS.core_profiles__profiles_1d)
    all_species = [specie.name for specie in IMAS.species(cp1d; return_zero_densities=true)]
    return evolve_densities_dict_creation(Symbol[]; replay_species=all_species, quasi_neutrality_specie=false)
end

"""
    check_output_fluxes(output::Vector{<:Real}, what::String)

Checks if there are any NaNs in the output
"""
function check_output_fluxes(output::Vector{<:Real}, what::String)
    @assert !any(isnan, output) "The transport flux is NaN check your transport model fluxes in dd.core_transport ($(what)): $(output)"
end

function cp1d_copy_primary_quantities(cp1d::IMAS.core_profiles__profiles_1d{T}) where {T<:Real}
    to_cp1d = IMAS.core_profiles__profiles_1d{T}()
    to_cp1d.grid.rho_tor_norm = deepcopy(cp1d.grid.rho_tor_norm)
    to_cp1d.electrons.density_thermal = deepcopy(cp1d.electrons.density_thermal)
    to_cp1d.electrons.temperature = deepcopy(cp1d.electrons.temperature)
    if !ismissing(cp1d.electrons, :density_fast)
        to_cp1d.electrons.density_fast = deepcopy(cp1d.electrons.density_fast)
    end
    resize!(to_cp1d.ion, length(cp1d.ion))
    for (initial_ion, ion) in zip(to_cp1d.ion, cp1d.ion)
        initial_ion.element = ion.element
        initial_ion.density_thermal = deepcopy(ion.density_thermal)
        if !ismissing(ion, :density_fast)
            initial_ion.density_fast = deepcopy(ion.density_fast)
        end
        initial_ion.temperature = deepcopy(ion.temperature)
    end
    to_cp1d.rotation_frequency_tor_sonic = deepcopy(cp1d.rotation_frequency_tor_sonic)
    return to_cp1d
end

function cp1d_copy_primary_quantities!(to_cp1d::T, cp1d::T) where {T<:IMAS.core_profiles__profiles_1d{<:Real}}
    @assert length(to_cp1d.ion) == length(cp1d.ion)
    to_cp1d.grid.rho_tor_norm .= cp1d.grid.rho_tor_norm
    to_cp1d.electrons.density_thermal .= cp1d.electrons.density_thermal
    to_cp1d.electrons.temperature .= cp1d.electrons.temperature
    if !ismissing(cp1d.electrons, :density_fast)
        to_cp1d.electrons.density_fast .= cp1d.electrons.density_fast
    end
    for (initial_ion, ion) in zip(to_cp1d.ion, cp1d.ion)
        initial_ion.density_thermal .= ion.density_thermal
        initial_ion.temperature .= ion.temperature
        if !ismissing(ion, :density_fast)
            initial_ion.density_fast .= ion.density_fast
        end
    end
    to_cp1d.rotation_frequency_tor_sonic .= cp1d.rotation_frequency_tor_sonic
    return to_cp1d
end

"""
    calculate_w0_norm(cp1d::IMAS.core_profiles__profiles_1d)

Calculate TGYRO normalization frequency w0_norm = cs/R₀ where:

    cs = sound speed at axis = √(kB*Te(0)/md)
"""
function calculate_w0_norm(Te_axis)
    cs_axis = sqrt(IMAS.mks.k_B * Te_axis / IMAS.mks.m_d)  # sound speed at axis
    R0 = 1.0 # 1 [m]
    return cs_axis / R0
end

"""
    _step(replay_actor::ActorReplay, actor::ActorFluxMatcher, replay_dd::IMAS.dd)

Replay profiles from replay_dd to current dd for channels set to :replay
"""
function _step(replay_actor::ActorReplay, actor::ActorFluxMatcher, replay_dd::IMAS.dd)
    dd = actor.dd
    par = actor.par

    time0 = dd.global_time
    cp1d = dd.core_profiles.profiles_1d[time0]
    replay_cp1d = replay_dd.core_profiles.profiles_1d[time0]
    rho = cp1d.grid.rho_tor_norm

    # Get the transport grid boundary for blending
    rho_nml = par.rho_transport[end]
    rho_ped = par.rho_transport[end]

    # Replay electron temperature if set to :replay
    if par.evolve_Te == :replay
        cp1d.electrons.temperature = IMAS.blend_core_edge(replay_cp1d.electrons.temperature, cp1d.electrons.temperature, rho, rho_nml, rho_ped)
    end

    # Replay ion temperature if set to :replay
    if par.evolve_Ti == :replay
        Ti_replay = IMAS.blend_core_edge(replay_cp1d.t_i_average, cp1d.t_i_average, rho, rho_nml, rho_ped)
        for ion in cp1d.ion
            ion.temperature = Ti_replay
        end
    end

    # Replay rotation if set to :replay
    if par.evolve_rotation == :replay
        i_nml = IMAS.argmin_abs(rho, rho_nml)
        ω_core = deepcopy(replay_cp1d.rotation_frequency_tor_sonic)
        ω_edge = IMAS.freeze!(cp1d, :rotation_frequency_tor_sonic)
        ω_core[i_nml+1:end] = ω_edge[i_nml+1:end]
        ω_core[1:i_nml] = ω_core[1:i_nml] .- ω_core[i_nml] .+ ω_edge[i_nml]
        cp1d.rotation_frequency_tor_sonic = ω_core
    end

    # Handle density replays
    evolve_densities = evolve_densities_dictionary(cp1d, par)
    if !isempty(evolve_densities)
        # Replay electron density if set to :replay
        if evolve_densities[:electrons] == :replay
            cp1d.electrons.density_thermal = IMAS.blend_core_edge(replay_cp1d.electrons.density_thermal, cp1d.electrons.density_thermal, rho, rho_nml, rho_ped; method=:scale)
        end

        # Replay ion densities if set to :replay
        for (ion, replay_ion) in zip(cp1d.ion, replay_cp1d.ion)
            ion_symbol = Symbol(ion.label)
            if haskey(evolve_densities, ion_symbol) && evolve_densities[ion_symbol] == :replay
                if !ismissing(ion, :density_thermal)
                    ion.density_thermal = IMAS.blend_core_edge(replay_ion.density_thermal, ion.density_thermal, rho, rho_nml, rho_ped; method=:scale)
                end
            end
        end

        # Ensure quasi neutrality if needed
        for (species, evolve) in evolve_densities
            if evolve == :quasi_neutrality
                IMAS.enforce_quasi_neutrality!(cp1d, species)
                break
            end
        end
    end

    return replay_actor
end
