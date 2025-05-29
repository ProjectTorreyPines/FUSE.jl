import NLsolve
using LinearAlgebra
import TJLF: InputTJLF

import NonlinearSolve, FixedPointAcceleration

#= ================ =#
#  ActorFluxMatcher  #
#= ================ =#
Base.@kwdef mutable struct FUSEparameters__ActorFluxMatcher{T<:Real} <: ParametersActor{T}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    _time::Float64 = NaN
    rho_transport::Entry{AbstractVector{T}} = Entry{AbstractVector{T}}("-", "ρ transport grid"; default=0.25:0.1:0.85)
    evolve_Ti::Switch{Symbol} = Switch{Symbol}([:flux_match, :fixed], "-", "Ion temperature `:flux_match` or keep `:fixed`"; default=:flux_match)
    evolve_Te::Switch{Symbol} = Switch{Symbol}([:flux_match, :fixed], "-", "Electron temperature `:flux_match` or keep `:fixed`"; default=:flux_match)
    evolve_densities::Entry{Union{AbstractDict,Symbol}} =
        Entry{Union{AbstractDict,Symbol}}(
            "-",
            """
            * `:fixed`: no density evolution
            * `:flux_match:` electron flux-match and ions match ne scale while keeping quasi-neutrality, giving a flat zeff
            * `:zeff:` electron flux-match and ions scaled to keep zeff constant while keeping quasi-neutrality
            * Dict to specify which species are `:flux_match`, kept `:fixed`, used to enforce `:quasi_neutrality`, or scaled to `:match_ne_scale`""";
            default=:flux_match
        )
    evolve_rotation::Switch{Symbol} = Switch{Symbol}([:flux_match, :fixed], "-", "Rotation `:flux_match` or keep `:fixed`"; default=:fixed)
    evolve_pedestal::Entry{Bool} = Entry{Bool}("-", "Evolve the pedestal at each iteration"; default=true)
    evolve_plasma_sources::Entry{Bool} = Entry{Bool}("-", "Update the plasma sources at each iteration"; default=true)
    find_widths::Entry{Bool} = Entry{Bool}("-", "Runs turbulent transport actor TJLF finding widths after first iteration"; default=true)
    max_iterations::Entry{Int} = Entry{Int}("-", "Maximum optimizer iterations"; default=-1)
    algorithm::Switch{Symbol} =
        Switch{Symbol}([:polyalg, :broyden, :anderson, :simple, :old_anderson, :none], "-", "Optimizing algorithm used for the flux matching"; default=:broyden)
    step_size::Entry{T} = Entry{T}(
        "-",
        "Step size for each algorithm iteration (note this has a different meaning for each algorithm)";
        default=1.0,
        check=x -> @assert x > 0.0 "must be: step_size > 0.0"
    )
    Δt::Entry{Float64} = Entry{Float64}("s", "Evolve for Δt (Inf for steady state)"; default=Inf)
    relax::Entry{Float64} = Entry{Float64}("-", "Relaxation on the final solution"; default=1.0, check=x -> @assert 0.0 <= x <= 1.0 "must be: 0.0 <= relax <= 1.0")
    scale_turbulence_law::Switch{Symbol} = Switch{Symbol}([:h98, :ds03], "-", "Scale turbulent transport to achieve a desired confinement law")
    scale_turbulence_value::Entry{Float64} = Entry{Float64}(
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
    actor_ped::ActorPedestal{D,P}
    norms::Vector{Float64}
    error::Float64
end

"""
    ActorFluxMatcher(dd::IMAS.dd, act::ParametersAllActors; kw...)

Evalutes the transport fluxes and source fluxes and minimizes the flux_match error
"""
function ActorFluxMatcher(dd::IMAS.dd, act::ParametersAllActors; kw...)
    actor = ActorFluxMatcher(dd, act.ActorFluxMatcher, act; kw...)
    step(actor)
    finalize(actor)
    return actor
end

function ActorFluxMatcher(dd::IMAS.dd, par::FUSEparameters__ActorFluxMatcher, act::ParametersAllActors; kw...)
    logging_actor_init(ActorFluxMatcher)
    par = OverrideParameters(par; kw...)
    actor_ct = ActorFluxCalculator(dd, act.ActorFluxCalculator, act; par.rho_transport)
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
    return ActorFluxMatcher(dd, par, act, actor_ct, actor_ped, Float64[], Inf)
end

"""
    _step(actor::ActorFluxMatcher)

ActorFluxMatcher step
"""
function _step(actor::ActorFluxMatcher{D,P}) where {D<:Real,P<:Real}
    dd = actor.dd
    par = actor.par
    cp1d = dd.core_profiles.profiles_1d[]

    IMAS.sources!(dd)

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
    N_channels = Int(floor(length(z_init) / N_radii))

    z_scaled_history = Vector{NTuple{length(z_init),Float64}}()
    err_history = Vector{Vector{Float64}}()

    actor.norms = fill(NaN, N_channels)

    ProgressMeter.ijulia_behavior(:clear)
    prog = ProgressMeter.ProgressUnknown(; dt=0.1, desc="Calls:", enabled=par.verbose)
    old_logging = actor_logging(dd, false)

    if ismissing(par, :scale_turbulence_law)
        opt_parameters = z_init_scaled
    else
        opt_parameters = [1.0; z_init_scaled]
    end

    out = try
        if par.algorithm == :none
            res = (zero=opt_parameters,)

        elseif par.algorithm == :simple
            ftol = 1E-2 # relative error
            xtol = 1E-3 # difference in input array
            res = flux_match_simple(actor, opt_parameters, initial_cp1d, z_scaled_history, err_history, ftol, xtol, prog)

        else
            # 1. In-place residual
            function f!(F, u, initial_cp1d)
                try
                    F .= flux_match_errors(actor, u, initial_cp1d; z_scaled_history, err_history, prog).errors
                catch error
                    if isa(error, InterruptException)
                        rethrow(error)
                    end
                    F .= Inf
                end
            end

            # 2. Problem definition
            problem = NonlinearSolve.NonlinearProblem(f!, opt_parameters, initial_cp1d)

            # 3. Algorithm selection
            alg = if par.algorithm == :old_anderson
                NonlinearSolve.NLsolveJL(; method=:anderson, m=4, beta=-par.step_size * 0.5)

            elseif par.algorithm == :anderson
                NonlinearSolve.FixedPointAccelerationJL(;
                    algorithm=:Anderson,
                    m=5,                        # history length
                    dampening=-par.step_size * 0.5,
                    condition_number_threshold=1e8,
                    replace_invalids=:ReplaceVector
                )

            elseif par.algorithm == :broyden
                alg = NonlinearSolve.Broyden(;
                    linesearch=NonlinearSolve.LineSearch.BackTracking(;
                        c_1=1e-4,   # Armijo parameter
                        ρ_hi=0.5,    # shrink factor
                        order=2,      # model order
                        autodiff=NonlinearSolve.ADTypes.AutoFiniteDiff()
                    ),
                    update_rule=Val(:good_broyden),
                    init_jacobian=Val(:true_jacobian),
                    autodiff=NonlinearSolve.ADTypes.AutoFiniteDiff()
                )
            elseif par.algorithm == :polyalg
                # Default NonlinearSolve algorithm
                NonlinearSolve.FastShortcutNonlinearPolyalg(; autodiff=NonlinearSolve.ADTypes.AutoFiniteDiff())

            else
                error("Unsupported algorithm: $(par.algorithm)")
            end

            if par.max_iterations > 0
                max_iterations = par.max_iterations
            else
                # Default for gradient methods 20, otherwise 500
                max_iterations = (par.algorithm in (:broyden, :polyalg)) ? 15 : 500
            end

            # 4. Solve with matching tolerances and iteration limits
            # NonlinearSolve abstol is meant to be on u, but actually gets
            #   passed to ftol in NLsolve which is an error on the residual
            # See https://github.com/SciML/NonlinearSolve.jl/issues/593
            abstol = 1E-2
            sol = NonlinearSolve.solve(problem, alg;
                abstol,
                maxiters=max_iterations,
                show_trace=Val(par.show_trace),
                store_trace=Val(false),
                verbose=false
            )

            # Extract the solution vector
            res = (zero=sol.u,)
        end

        flux_match_errors(actor, collect(res.zero), initial_cp1d) # z_profiles for the smallest error iteration

    finally

        actor_logging(dd, old_logging)
    end

    # evaluate profiles at the best-matching gradients
    actor.error = norm(out.errors)
    @ddtime(dd.transport_solver_numerics.convergence.time_step.time = dd.global_time)
    @ddtime(dd.transport_solver_numerics.convergence.time_step.data = actor.error)
    dd.transport_solver_numerics.ids_properties.name = "FluxMatcher"
    ProgressMeter.finish!(prog; showvalues=progress_ActorFluxMatcher(dd, norm(out.errors)))

    # plotting of the channels that have been evolved
    if par.do_plot
        plot()
        for (ch, profiles_path) in enumerate(profiles_paths)
            title = profiles_title(cp1d, profiles_path)
            plot!(vcat(map(x -> errors_by_channel(x, N_radii, N_channels), err_history)...)[:, ch]; label=title)
        end
        display(
            plot!(vcat(map(norm, err_history)...);
                color=:black,
                lw=0.5,
                yscale=:log10,
                ylabel="Log₁₀ of convergence errror",
                xlabel="Iterations",
                label=@sprintf("total [error = %.3e]", (minimum(map(norm, err_history))))))

        channels_evolution = transpose(hcat(map(z -> collect(unscale_z_profiles(z)), z_scaled_history)...))
        data = reshape(channels_evolution, (length(err_history), length(par.rho_transport), N_channels))
        p = plot()
        for (ch, profiles_path) in enumerate(profiles_paths)
            title = profiles_title(cp1d, profiles_path)
            for kr in 1:length(par.rho_transport)
                plot!(data[:, kr, ch]; ylabel="Inverse scale length [m⁻¹]", xlabel="Iterations", primary=kr == 1, lw=kr, label=title)
            end
        end
        display(p)

        p = plot(; layout=(N_channels, 2), size=(1000, 300 * N_channels))
        total_flux1d = IMAS.total_fluxes(dd.core_transport, cp1d, par.rho_transport; time0=dd.global_time)
        total_source1d = IMAS.total_sources(dd.core_sources, cp1d; time0=dd.global_time)
        model_type = IMAS.name_2_index(dd.core_transport.model)
        for (ch, (profiles_path, fluxes_path)) in enumerate(zip(profiles_paths, fluxes_paths))
            title = IMAS.p2i(collect(map(string, fluxes_path)))
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
    end

    # remove the temporary `∂/∂t implicit` term and update the full `∂/∂t` term
    if !isinf(par.Δt)
        deleteat!(dd.core_sources.source, :time_derivative, "identifier.name" => "∂/∂t implicit")
        IMAS.time_derivative_source!(dd)
    end

    # free total densities expressions
    IMAS.unfreeze!(cp1d.electrons, :density)
    for ion in cp1d.ion
        IMAS.unfreeze!(ion, :density)
    end

    # free pressures expressions
    IMAS.unfreeze!(cp1d.electrons, :pressure_thermal)
    IMAS.unfreeze!(cp1d.electrons, :pressure)
    for ion in cp1d.ion
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

function errors_by_channel(errors::Vector{Float64}, N_radii::Int, N_channels::Int)
    return mapslices(norm, reshape(errors, (N_radii, N_channels)); dims=1)
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
        actor::ActorFluxMatcher,
        z_profiles_scaled::Vector{<:Real},
        initial_cp1d::IMAS.core_profiles__profiles_1d,
        z_scaled_history::Vector=[],
        err_history::Vector{Vector{Float64}}=Vector{Vector{Float64}}(),
        prog::Any=nothing)

Update the profiles, evaluates neoclassical and turbulent fluxes, sources (ie target fluxes), and returns named tuple with (targets, fluxes, errors)

NOTE: flux matching is done in physical units
"""
function flux_match_errors(
    actor::ActorFluxMatcher,
    opt_parameters::Vector{<:Real},
    initial_cp1d::IMAS.core_profiles__profiles_1d;
    z_scaled_history::Vector=[],
    err_history::Vector{Vector{Float64}}=Vector{Vector{Float64}}(),
    prog::Any=nothing)

    dd = actor.dd
    par = actor.par
    cp1d = dd.core_profiles.profiles_1d[]

    if ismissing(par, :scale_turbulence_law)
        z_profiles_scaled = opt_parameters
    else
        turbulence_scale = opt_parameters[1]
        z_profiles_scaled = @views opt_parameters[2:end]
    end

    # unscale z_profiles
    push!(z_scaled_history, Tuple(z_profiles_scaled))
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
        m1d.total_ion_energy.flux .*= turbulence_scale
        m1d.electrons.energy.flux .*= turbulence_scale
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
        if norm0 === NaN
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
        tau_th = IMAS.tau_e_thermal(dd; include_radiation=true)
        if par.scale_turbulence_law == :h98
            tauH = IMAS.tau_e_h98(dd; include_radiation=true)
        elseif par.scale_turbulence_law == :ds03
            tauH = IMAS.tau_e_ds03(dd; include_radiation=true)
        else
            error("act.ActorFluxMatcher.scale_turbulence_law=$(par.scale_turbulence_law) not recognized: valid options are :h98 or :ds03")
        end
        H_value = tau_th / tauH
        H_target = par.scale_turbulence_value
        errors = [(H_value - H_target) / H_target; errors]
    end

    # update error history
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
function flux_match_targets(dd::IMAS.dd, par::OverrideParameters{P,FUSEparameters__ActorFluxMatcher{P}}) where {P<:Real}
    cp1d = dd.core_profiles.profiles_1d[]

    total_source = resize!(dd.core_sources.source, :total; wipe=false)
    total_source1d = resize!(total_source.profiles_1d; wipe=false)
    IMAS.total_sources!(total_source1d, dd.core_sources, cp1d; time0=dd.global_time, fields=[:total_ion_power_inside, :power_inside, :particles_inside, :torque_tor_inside])
    cs_gridpoints = [argmin_abs(total_source1d.grid.rho_tor_norm, rho_x) for rho_x in par.rho_transport]

    targets = Float64[]

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
    flux_match_fluxes(dd::IMAS.dd{T}, par::OverrideParameters{P,FUSEparameters__ActorFluxMatcher{P}}) where {T<:Real, P<:Real}

Evaluates the flux_matching fluxes for the :flux_match species and channels

NOTE: flux matching is done in physical units
"""
function flux_match_fluxes(dd::IMAS.dd{T}, par::OverrideParameters{P,FUSEparameters__ActorFluxMatcher{P}}) where {T<:Real,P<:Real}
    cp1d = dd.core_profiles.profiles_1d[]

    total_flux = resize!(dd.core_transport.model, :combined; wipe=false)
    total_flux1d = resize!(total_flux.profiles_1d; wipe=false)
    IMAS.total_fluxes!(total_flux1d, dd.core_transport, cp1d, par.rho_transport; time0=dd.global_time)

    fluxes = T[]

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
        actor::ActorFluxMatcher,
        z_init::Vector{<:Real},
        initial_cp1d::IMAS.core_profiles__profiles_1d,
        z_scaled_history::Vector,
        err_history::Vector{Vector{Float64}},
        ftol::Float64,
        xtol::Float64,
        prog::Any)

Updates zprofiles based on TGYRO simple algorithm
"""
function flux_match_simple(
    actor::ActorFluxMatcher,
    opt_parameters::Vector{<:Real},
    initial_cp1d::IMAS.core_profiles__profiles_1d,
    z_scaled_history::Vector,
    err_history::Vector{Vector{Float64}},
    ftol::Float64,
    xtol::Float64,
    prog::Any)

    par = actor.par

    i = 0

    if ismissing(par, :scale_turbulence_law)
        z_init_scaled = opt_parameters
    else
        turbulence_scale = opt_parameters[1]
        z_init_scaled = @views opt_parameters[2:end]
    end

    zprofiles_old = unscale_z_profiles(z_init_scaled)
    targets, fluxes, errors = flux_match_errors(actor, opt_parameters, initial_cp1d; z_scaled_history, err_history, prog)
    ferror = norm(errors)
    xerror = Inf
    step_size = par.step_size
    max_iterations = par.max_iterations
    while (ferror > ftol) || (xerror .> xtol)
        i += 1
        if (i > abs(max_iterations))
            if max_iterations > 0
                @info "Unable to flux-match within $(max_iterations) iterations (aerr) = $(ferror) (ftol=$ftol) (xerr) = $(xerror) (xtol = $xtol)"
            end
            break
        end

        zprofiles = zprofiles_old .* (1.0 .+ step_size * 0.1 .* (targets .- fluxes) ./ sqrt.(1.0 .+ fluxes .^ 2 + targets .^ 2))
        if ismissing(par, :scale_turbulence_law)
            targets, fluxes, errors = flux_match_errors(actor, scale_z_profiles(zprofiles), initial_cp1d; z_scaled_history, err_history, prog)
        else
            turbulence_scale += errors[1] / 10.0
            targets, fluxes, errors =
                flux_match_errors(actor, [turbulence_scale; scale_z_profiles(zprofiles)], initial_cp1d; z_scaled_history, err_history, prog)
        end
        xerror = maximum(abs.(zprofiles .- zprofiles_old)) / step_size
        ferror = norm(errors)
        zprofiles_old = zprofiles
    end

    if ismissing(par, :scale_turbulence_law)
        return (zero=collect(z_scaled_history[argmin(map(norm, err_history))]),)
    else
        return (zero=[turbulence_scale; collect(z_scaled_history[argmin(map(norm, err_history))])],)
    end
end

function progress_ActorFluxMatcher(dd::IMAS.dd, error::Float64)
    cp1d = dd.core_profiles.profiles_1d[]
    out = Tuple{String,Float64}[]
    push!(out, ("         error", error))
    pfus = IMAS.fusion_power(cp1d)
    if pfus > 1E3
        push!(out, ("  Pfusion [MW]", pfus / 1E6))
    end
    push!(out, ("     Ti0 [keV]", cp1d.t_i_average[1] / 1E3))
    push!(out, ("     Te0 [keV]", cp1d.electrons.temperature[1] / 1E3))
    push!(out, ("ne0 [10²⁰ m⁻³]", cp1d.electrons.density_thermal[1] / 1E20))
    return out
end

function evolve_densities_dictionary(cp1d::IMAS.core_profiles__profiles_1d, par::OverrideParameters{P,FUSEparameters__ActorFluxMatcher{P}}) where {P<:Real}
    if par.evolve_densities == :fixed
        return setup_density_evolution_fixed(cp1d)
    elseif par.evolve_densities == :flux_match
        return setup_density_evolution_electron_flux_match_rest_ne_scale(cp1d)
    elseif par.evolve_densities == :zeff
        return setup_density_evolution_electron_flux_match_fixed_zeff(cp1d)
    elseif typeof(par.evolve_densities) <: AbstractDict{Symbol,Symbol}
        for (k, v) in par.evolve_densities
            @assert v in (:quasi_neutrality, :match_ne_scale, :fixed, :zeff, :flux_match) "evolve_density `:$(k) => :$(v)` is not allowed. Choose one of [:quasi_neutrality, :match_ne_scale, :fixed, :zeff, :flux_match]"
        end
        return par.evolve_densities
    else
        error(
            "act.ActorFluxMatcher.evolve_densities cannot be `$(repr(par.evolve_densities))`. Use `:fixed`, `:flux_match` or a dictionary specifiying [:quasi_neutrality, :match_ne_scale, :fixed, :zeff, :flux_match] for each specie"
        )
    end
end

"""
    pack_z_profiles(cp1d::IMAS.core_profiles__profiles_1d, par::OverrideParameters{P,FUSEparameters__ActorFluxMatcher{P}}) where {P<:Real}

Packs the z_profiles based on evolution parameters

NOTE: the order for packing and unpacking is always: [Te, Ti, Rotation, ne, nis...]
"""
function pack_z_profiles(cp1d::IMAS.core_profiles__profiles_1d, par::OverrideParameters{P,FUSEparameters__ActorFluxMatcher{P}}) where {P<:Real}
    cp_gridpoints = [argmin_abs(cp1d.grid.rho_tor_norm, rho_x) for rho_x in par.rho_transport]

    z_profiles = Float64[]
    profiles_paths = []
    fluxes_paths = []

    if par.evolve_Te == :flux_match
        z_Te = IMAS.calc_z(cp1d.grid.rho_tor_norm, cp1d.electrons.temperature, :third_order)[cp_gridpoints]
        append!(z_profiles, z_Te)
        push!(profiles_paths, (:electrons, :temperature))
        push!(fluxes_paths, (:electrons, :energy))
    end

    if par.evolve_Ti == :flux_match
        z_Ti = IMAS.calc_z(cp1d.grid.rho_tor_norm, cp1d.t_i_average, :third_order)[cp_gridpoints]
        append!(z_profiles, z_Ti)
        push!(profiles_paths, (:t_i_average,))
        push!(fluxes_paths, (:total_ion_energy,))
    end

    if par.evolve_rotation == :flux_match
        z_rot = IMAS.calc_z(cp1d.grid.rho_tor_norm, cp1d.rotation_frequency_tor_sonic, :third_order)[cp_gridpoints]
        append!(z_profiles, z_rot)
        push!(profiles_paths, (:rotation_frequency_tor_sonic,))
        push!(fluxes_paths, (:momentum_tor,))
    end

    evolve_densities = evolve_densities_dictionary(cp1d, par)
    if !isempty(evolve_densities)
        check_evolve_densities(cp1d, evolve_densities)
        if evolve_densities[:electrons] == :flux_match
            z_ne = IMAS.calc_z(cp1d.grid.rho_tor_norm, cp1d.electrons.density_thermal, :third_order)[cp_gridpoints]
            append!(z_profiles, z_ne)
            push!(profiles_paths, (:electrons, :density_thermal))
            push!(fluxes_paths, (:electrons, :particles))
        end
        for (k, ion) in enumerate(cp1d.ion)
            if evolve_densities[Symbol(ion.label)] == :flux_match
                z_ni = IMAS.calc_z(cp1d.grid.rho_tor_norm, ion.density_thermal, :third_order)[cp_gridpoints]
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
        cp1d.rotation_frequency_tor_sonic =
            IMAS.profile_from_z_transport(cp1d.rotation_frequency_tor_sonic .+ 1, cp1d.grid.rho_tor_norm, cp_rho_transport, z_profiles[counter+1:counter+N])
        counter += N
    end

    if !isempty(evolve_densities)
        if evolve_densities[:electrons] == :flux_match
            z_ne = z_profiles[counter+1:counter+N]
            cp1d.electrons.density_thermal = IMAS.profile_from_z_transport(cp1d.electrons.density_thermal, cp1d.grid.rho_tor_norm, cp_rho_transport, z_ne)
            counter += N
        else
            z_ne = IMAS.calc_z(cp1d.grid.rho_tor_norm, cp1d.electrons.density_thermal, :third_order)[cp_gridpoints]
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
    dd_species = [
        :electrons;
        Symbol[specie.name for specie in IMAS.species(cp1d; only_electrons_ions=:ions, only_thermal_fast=:thermal)];
        Symbol[specie.name for specie in IMAS.species(cp1d; only_electrons_ions=:ions, only_thermal_fast=:fast)]
    ]

    # Check if evolve_densities contains all of dd thermal species
    @assert sort([specie for (specie, evolve) in evolve_densities]) == sort(dd_species) "Mismatch: dd species $(sort(dd_species)) VS evolve_densities species : $(sort([specie for (ispecie,j) in evolve_densities]))"

    # Check that either all species are fixed, or there is 1 quasi_neutrality specie when evolving densities
    if any(evolve == :zeff for (specie, evolve) in evolve_densities if specie != :electrons)
        txt = "When flux_matching densities, either none or all ion species must be :zeff"
        @assert all(evolve == :zeff for (specie, evolve) in evolve_densities if specie != :electrons)
    elseif all(evolve == :fixed for (specie, evolve) in evolve_densities if evolve != :quasi_neutrality)
        txt = "When flux_matching densities, no more than one species can be set to :quasi_neutrality"
        @assert length([specie for (specie, evolve) in evolve_densities if evolve == :quasi_neutrality]) <= 1 txt
    elseif any(evolve == :flux_match for (specie, evolve) in evolve_densities)
        txt = "When flux_matching densities, one an only one species must be set to :quasi_neutrality"
        @assert length([specie for (specie, evolve) in evolve_densities if evolve == :quasi_neutrality]) == 1 txt
    else
        error("Invalid evolve_densities = $(evolve_densities)")
    end
end

"""
    setup_density_evolution_electron_flux_match_rest_ne_scale(cp1d::IMAS.core_profiles__profiles_1d)

Sets up the evolve_density dict to evolve only ne and keep the rest matching the ne_scale lengths
"""
function setup_density_evolution_electron_flux_match_rest_ne_scale(cp1d::IMAS.core_profiles__profiles_1d)
    dd_thermal = Symbol[specie.name for specie in IMAS.species(cp1d; only_electrons_ions=:ions, only_thermal_fast=:thermal)]
    dd_fast = Symbol[specie.name for specie in IMAS.species(cp1d; only_electrons_ions=:all, only_thermal_fast=:fast)]
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
    dd_thermal = Symbol[specie.name for specie in IMAS.species(cp1d; only_electrons_ions=:ions, only_thermal_fast=:thermal)]
    dd_fast = Symbol[specie.name for specie in IMAS.species(cp1d; only_electrons_ions=:all, only_thermal_fast=:fast)]
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
    all_species = [specie.name for specie in IMAS.species(cp1d)]
    return evolve_densities_dict_creation(Symbol[]; fixed_species=all_species, quasi_neutrality_specie=false)
end

"""
    setup_density_evolution_electron_flux_match_fixed_zeff(cp1d::IMAS.core_profiles__profiles_1d)

Sets up the evolve_density dict to evolve only ne, and scale ion species in such a way to keep Zeff fixed
"""
function setup_density_evolution_electron_flux_match_fixed_zeff(cp1d::IMAS.core_profiles__profiles_1d)
    evolve_density = Dict{Symbol,Symbol}()
    evolve_density[:electrons] = :flux_match
    for ion_name in IMAS.species(cp1d; only_electrons_ions=:ions)
        evolve_density[ion_name.name] = :zeff
    end
    return evolve_density
end

"""
    evolve_densities_dict_creation(flux_match_species::Vector, fixed_species::Vector, match_ne_scale_species::Vector; quasi_neutrality_specie::Union{Symbol,Bool}=false)

Create the density_evolution dict based on input vectors: flux_match_species, fixed_species, match_ne_scale_species, quasi_neutrality_specie
"""
function evolve_densities_dict_creation(
    flux_match_species::Vector;
    fixed_species::Vector{Symbol}=Symbol[],
    match_ne_scale_species::Vector{Symbol}=Symbol[],
    quasi_neutrality_specie::Union{Symbol,Bool}=false)

    parse_list = vcat([[sym, :flux_match] for sym in flux_match_species], [[sym, :match_ne_scale] for sym in match_ne_scale_species], [[sym, :fixed] for sym in fixed_species])
    if typeof(quasi_neutrality_specie) <: Symbol
        parse_list = vcat(parse_list, [[quasi_neutrality_specie, :quasi_neutrality]])
    end
    return Dict(sym => evolve for (sym, evolve) in parse_list)
end

"""
    check_output_fluxes(output::Vector{Float64}, what::String)

Checks if there are any NaNs in the output
"""
function check_output_fluxes(output::Vector{Float64}, what::String)
    @assert !any(isnan, output) "The transport flux is NaN check your transport model fluxes in dd.core_transport ($(what)): $(output)"
end

function cp1d_copy_primary_quantities(cp1d::IMAS.core_profiles__profiles_1d{T}) where {T<:Real}
    to_cp1d = IMAS.core_profiles__profiles_1d{T}()
    to_cp1d.grid.rho_tor_norm = deepcopy(cp1d.grid.rho_tor_norm)
    to_cp1d.electrons.density_thermal = deepcopy(cp1d.electrons.density_thermal)
    to_cp1d.electrons.temperature = deepcopy(cp1d.electrons.temperature)
    resize!(to_cp1d.ion, length(cp1d.ion))
    for (initial_ion, ion) in zip(to_cp1d.ion, cp1d.ion)
        initial_ion.element = ion.element
        initial_ion.density_thermal = deepcopy(ion.density_thermal)
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
    for (initial_ion, ion) in zip(to_cp1d.ion, cp1d.ion)
        initial_ion.density_thermal .= ion.density_thermal
        initial_ion.temperature .= ion.temperature
    end
    to_cp1d.rotation_frequency_tor_sonic .= cp1d.rotation_frequency_tor_sonic
    return to_cp1d
end
