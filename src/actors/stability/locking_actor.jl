import Roots
using Plots
import ModeLocking
using ModeLocking: ODEparams, NNparams, LockingNNModel, LockingResults
#import FUSE: coordinates


#================== =#
#  ActorLockingProbability  #
#================== =#
Base.@kwdef mutable struct FUSEparameters__ActorLocking{T<:Real} <: ParametersActor{T}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    _time::Float64 = NaN
    m_pol::Entry{Int} = Entry{Int}("_", "poloidal mode number of the mode"; default=2)
    n_tor::Entry{Int} = Entry{Int}("_", "toroidal mode number of the mode"; default=1)
    q_surf::Entry{Float64} = Entry{Float64}("_", "rational surface of interest, usually 2.0"; default=2.)
    grid_size::Entry{Int} = Entry{Int}("-", "grid resolution for control space"; default=100)
    t_final::Entry{Float64} = Entry{Float64}("-", "Final integration time in units of tearing time (~ms)"; default=100.)
    time_steps::Entry{Int} = Entry{Int}("-", "number of time steps for the ODE integration"; default=200)
    overwrite_params::Entry{Bool} = Entry{Bool}("-", "Whether to overwrite ODE parameters to reproduce PoP2024 results"; default=false)
    control_type::Switch{Symbol} = Switch{Symbol}([:EF, :LinStab, :NLsaturation], # EF: error field
        "-",                                                            # LinStab: vary stability_index,
        "Use a user specified Control case to run the locking models"; default=:EF) # NLsaturation: vary NL saturation
    task::Switch{Symbol} = Switch{Symbol}(
        [:solve_system, :single_case, :calc_prob, :Monte_Carlo, :evaluate_probability, :transfer_learning],
        "-",
        "Choose whether to simulate system on full-grid, do a single case, retrain NN only (:calc_prob), or Monte Carlo (NOT implemented)",
        default = :solve_system
    )
    application::Switch{String} = Switch{String}(
        ["RP-RW", "RP-IW", "RP-RP-RW", "RP-RW-IW"],
        "-", 
        "Type of application: 'RP-RW' for resistive plasma with a single rational surface interacting with a resistive wall; 
                              'RP-IW' for resistive plasma with a single rational surface interacting with an ideal wall;
                              'RP-RP-RW' for resistive plasma with two rational surfaces interacting with a resistive wall;
                              'RP-RW-IW' for resistive plasma with a single rational surface interacting with both walls"; 
        default="RP-RW"
    )
    EF_phase::Entry{Float64} = Entry{Float64}("-", "Phase of the applied error field (degrees)"; default=0.)
    NL_saturation::Entry{Bool} = Entry{Bool}("-", "Nonlinear saturation parameter for the mode"; default=false)
    RPRW_stability_index::Entry{Float64} = Entry{Float64}(
        "-", 
        "Stability index of the system (set to Neg. value for now)", default=-0.5)
    b0::Entry{Float64} = Entry{Float64}(
        "Tesla", 
        "Scale for magnetic perturbations, usually ~10Gauss"; default=1.e-3)
    t0::Entry{Float64} = Entry{Float64}(
        "seconds", 
        "Characteristic time scale for normalization , usually TM/RW growth rate"; default=1.e-3)
    r0::Entry{Float64} = Entry{Float64}(
        "meter", 
        "Length scale for the integration , usually minor radius"; default=1.)
    source_torque::Entry{Float64} = Entry{Float64}("Newton.meter", "NBI torque"; default=5.)
    plot_orientation::Switch{Symbol} = Switch{Symbol}([:portrait, :landscape], "-",
        "Tile-plot layout: :portrait = 3×2 (paper), :landscape = 2×3 (slides)"; default=:portrait)
end


mutable struct ActorLocking{D,P} <: SingleAbstractActor{D,P}
    dd::IMAS.dd{D}
    par::OverrideParameters{P,FUSEparameters__ActorLocking{P}}
    ode_params::Union{Nothing, ODEparams}
    results::Union{Nothing, LockingResults}
    nn_params::NNparams

    function ActorLocking(
        dd::IMAS.dd{D},
        par::FUSEparameters__ActorLocking{P};
        ode_params = nothing,
        nn_params  = NNparams(),
        kw...
    ) where {D<:Real,P<:Real}

        logging_actor_init(ActorLocking)

        # Apply standard FUSE parameter overrides
        par = OverrideParameters(par; kw...)

        # Handle ODE params
        ode = if ode_params === nothing
            nothing
        elseif ode_params isa ODEparams
            ode_params
        elseif ode_params isa NamedTuple
            ODEparams(; ode_params...)
        else
            error("ode_params must be nothing, ODEparams, or NamedTuple")
        end

        return new{D,P}(dd, par, ode, nothing, nn_params)
    end
end


"""
    ActorLocking(dd::IMAS.dd, act::ParametersAllActors; kw...)

Provides a common interface to run different equilibrium actors
"""
function ActorLocking(dd::IMAS.dd, act::ParametersAllActors; kw...)
    #par = act.ActorLocking(kw...)
    actor = ActorLocking(dd, act.ActorLocking; kw...)
    step(actor)
    finalize(actor)
    return actor
end



function _step(actor::ActorLocking)
    dd = actor.dd
    par = actor.par
    task = par.task
    application = par.application

    # Populate the physical parameters needed to solve the ODEs
    if actor.ode_params === nothing
        actor.ode_params = ODEparams(;)
    end
    actor.ode_params = set_up_ode_params!(dd, par, actor.ode_params)

    # Main driver routine
    # Time evolve the ODEs, calculate locking probability, or load/evaluate model
    if task == :single_case
        @info "Solving one case for system"
        solve_one_case(par, actor.ode_params, application)
    
    elseif task == :solve_system
        @info "Solving the ODEs for $(application) system"
        # Solve the ODE system on the whole control grid, normalize, and
        # classify (prob=nothing; filled in-place by train_locking_nn below)
        actor.results = _solve_grid_and_classify(actor, application)

        # Train NN with hyperparameters stored on the actor
        #train_locking_nn(actor, actor.nn_params)

        # Checkpoint ODE results to disk so calc_prob can run in a future session
        save_ode_results(actor)

    elseif task == :calc_prob
        # Retrain NN only — ODEs are NOT re-solved.
        # Use in-memory results if available; otherwise load from disk.
        if actor.results === nothing
            @info "No in-memory ODE results — loading from disk"
            load_ode_results!(actor)
        end
        @info "Training NN classifier to calculate probability of locking as a function of Controls"
        train_locking_nn(actor, actor.nn_params)

    elseif task == :evaluate_probability
        # Load NN model from disk (cached after first load).
        # Also load ODE results if not already in memory (needed for plot_sols).
        if actor.results === nothing
            @info "No in-memory ODE results — loading from disk"
            load_ode_results!(actor)
        end
        prob_model = load_locking_nn(; control_type=par.control_type)
        actor.results.prob = prob_model
        @info "Locking NN model ready — call actor.results.prob(C1, C2) or plot_sols(actor)"

    elseif task == :transfer_learning
        # Solve the (typically focused/sparse, via ode_params.Control1_min/max
        # and Control2_min/max + a smaller grid_size) grid for the new dd,
        # normalize, and classify — same pipeline as :solve_system.
        actor.results = _solve_grid_and_classify(actor, application)

        # Fine-tune the saved base NN model on the new equilibrium's data,
        # freezing all but the last layer.
        @info "Loading base NN model for transfer learning"
        base_model = load_locking_nn(; control_type=par.control_type)

        X_new, y_new = ModeLocking.prepare_nn_data(actor.results.locking_labels,
                                        actor.ode_params.Control1, actor.ode_params.Control2)

        actor.results.prob = ModeLocking.transfer_learn_locking_nn(base_model, X_new, y_new; nn_params=actor.nn_params)


    elseif task == :Monte_Carlo
        error("Monte-Carlo not implemented yet")

    else
        error("Unknown task: $(task)")

    end
    
    return actor
end

function _finalize(actor::ActorLocking)
    dd = actor.dd
    par = actor.par

    @info "Finalizing ActorLocking at time"
    mhd = dd.mhd_linear
    mhd_ts = resize!(mhd.time_slice; wipe=false)
    mode = resize!(mhd_ts.toroidal_mode, "perturbation_type.name" => "Locked mode", "n_tor" => 1)
    mode.perturbation_type.description = "Locked mode limit from ActorLocking"
    mode.amplitude_multiplier = 1.0 / (par.b0*par.r0)
    #mode.plasma.b_field_perturbed.coordinate1.coefficients_real = AbstractArray(1.0, 0., 0.)
        
    return actor
end


# ─────────────────────────────────────────────────────────────────────────────
#  Setup: build ODE parameters from dd (called once per actor run)
# ─────────────────────────────────────────────────────────────────────────────

function set_up_ode_params!(dd::IMAS.dd, par, ode_params::ODEparams)
    """
    Initialize the ODE parameters for the locking actor.
    
    Args:
        dd: IMAS data structure
        par: Parameters for the simulation
    Returns:
        ode_params: Initialized ODE parameters
    """
    
    # find the normalized radius of the q=2 surface
    q_prof = dd.equilibrium.time_slice[].profiles_1d.q
    rho = dd.equilibrium.time_slice[].profiles_1d.rho_tor_norm
    ode_params.rat_surface = find_rat_surface(q_prof, rho, par.q_surf)
    
    # calculate the stability indices and mutual inductances
    ode_params = calculate_stability_index!(dd, par, ode_params)

    # Set physical parameters in dimensionless form
    ode_params = set_phys_params!(dd, par, ode_params)
    
    # Overwrite params to reproduce PoP2024
    if par.overwrite_params
        @info "Overwriting ODE parameters to reproduce PoP2024 results"
        ode_params.mu = 0.1
        ode_params.Inertia = 1
        ode_params.rat_surface = 0.67 # overwrite for now
    end

    
    # Prepare control parameters based on the control type
    ode_params = set_control_parameters!(dd, par, ode_params)

    return ode_params
end




function find_rat_surface(q_prof::Vector{Float64}, rho::Vector{Float64}, rat_surface::Float64)
    """
    Find the location of the q=rat_surface surface given the q profile
    """
    q_interp = IMAS.interp1d(rho, q_prof)
    x0 = findfirst(x -> abs(x) > rat_surface, q_prof)
    f = x -> abs(q_interp(x)) - rat_surface
    rho_rat = Roots.secant_method(f, (rho[x0-1], rho[x0]))
    @info "Found q=$rat_surface surface at: $rho_rat"
    return rho_rat
end



function calculate_stability_index!(dd::IMAS.dd, par, ode_params::ODEparams)
    rt = ode_params.rat_surface  
    rw = ode_params.res_wall
    rc = ode_params.control_surf
    m0 = par.m_pol
    
    rat21 = (rw / rt)^m0
    rat12 = rat21^(-1)
    rat32 = (rc / rw)^m0
    rat23 = rat32^(-1)

    ode_params.l12 = (2 * m0 / rw) / (rat21 - rat12)
    ode_params.l21 = (2 * m0 / rt) / (rat21 - rat12)
    ode_params.l32 = (2 * m0 / rw) / (rat32 - rat23)

    DeltaWl = (m0 / rw) * (rat21 + rat12) / (rat21 - rat12)
    DeltaWr = -(m0 / rw) * (rat32 + rat23) / (rat32 - rat23)

    ode_params.DeltaW = DeltaWr - DeltaWl
    ode_params.stability_index = par.RPRW_stability_index + ode_params.l21 * ode_params.l12 / ode_params.DeltaW

    ## In case control_type=:LinStab, adjust the upper Deltat to keep the
    ## RP-RW system weakly stable. Then, adjust the lower Deltat
    if par.control_type == :LinStab
        @info "Adjusting Deltat range for LinStab control to keep RW-RP system stable"
        if ode_params.Control2_min > 0.0
            ode_params.Control2_min = -3.5
            println("Setting Deltat minimum to $(ode_params.Control2_min) for RP-RW system to be stable")
        end

        ode_params.Control2_max = ode_params.l21 * ode_params.l12 / ode_params.DeltaW - 5.e-2
        if ode_params.Control2_max > 0.0
            error("Maximum Deltat for stability must be <0 for RP-RW system to be stable")
        end
        
        println(ode_params.Control2_min, ode_params.Control2_max)
    end
    

    return ode_params
end


function set_phys_params!(dd::IMAS.dd, par, ode_params::ODEparams)
    """
    Set the physical parameters in dimensionless form
    
    Args:
        dd: IMAS data structure
        par: Parameters for the simulation
        ode_params: ODE parameters to be set
    Returns:
        ode_params: Updated ODE parameters with physical constants set
    """

    # Define some constants
    #      Also NEED Zeff
    cp1d = dd.core_profiles.profiles_1d[1]
    mu0_val = IMAS.mks.μ_0
    mass_ion = cp1d.ion[1].element[1].a * IMAS.mks.m_p  # kg, mass of the main ion species

    rt = ode_params.rat_surface
    rho = dd.core_sources.source[1].profiles_1d[1].grid.rho_tor_norm
    mass_dens = mass_ion * cp1d.ion[1].density_thermal
    rho_interp = IMAS.interp1d(rho, mass_dens)
    mass_dens_atq2 = rho_interp(rt)  # Get the mass density at the q=2 surface
    #Zeff = dd.core_profiles.profiles_1d.zeff[40]
 
    # Set all the scales
    psi0 = par.b0 * par.r0
    U0 = psi0^2 * par.r0 / mu0_val
    R0 = dd.equilibrium.vacuum_toroidal_field.r0
    # approximate core rotation
    rot_core = dd.core_profiles.profiles_1d[1].rotation_frequency_tor_sonic[10]
    
    # NBI torque at the rational surface:
    # Prefer the FUSE NBI actor output (dd.core_sources) if it has been populated upstream;
    # fall back to the user-specified par.source_torque scalar if not.
    torque_at_rat_surf = try
        cp1d          = dd.core_profiles.profiles_1d[]
        src1d         = IMAS.total_sources(dd.core_sources, cp1d;
                            time0=dd.global_time, fields=[:torque_tor_inside])
        rho_src       = src1d.grid.rho_tor_norm
        torque_inside = src1d.torque_tor_inside
        if isempty(rho_src) || all(iszero, torque_inside)
            error("empty")
        end
        torque_val = IMAS.interp1d(rho_src, torque_inside)(rt)
        @info "NBI torque at rational surface from dd.core_sources: $(round(torque_val; sigdigits=4)) N·m"
        torque_val
    catch
        @warn "dd.core_sources torque not available — falling back to par.source_torque = $(par.source_torque) N·m"
        par.source_torque
    end

    # Calculate the drag coefficient in SI
    muSI = torque_at_rat_surf / rot_core
    # Calculate first the moment of inertia in the layer  
    inertia = (2*π)^2 * mass_dens_atq2 * R0 * rt^3 * ode_params.layer_width

    # set nonlinear saturation
    if par.NL_saturation == false
        ode_params.saturation_param = 0.
    end

    # Set the dimensionless quantities
    ode_params.mu = muSI / (U0 * par.t0)
    ode_params.Inertia = inertia / (U0 * par.t0^2)

    return ode_params
end


function set_control_parameters!(dd::IMAS.dd, par, ode_params::ODEparams) 
    """
    Prepare the control parameters based on the control type.
    
    Args:
        control_type: Symbol indicating the type of control
        ode_params: ODE parameters to be modified
        par: Parameters for the simulation
    Returns:
        ode_params: Updated ODE parameters with control settings
    """

    control_type = par.control_type
    N = par.grid_size
    M = par.grid_size
    b0 = par.b0

    l21 = ode_params.l21
    l12 = ode_params.l12
    DeltaW = ode_params.DeltaW
    rt = ode_params.rat_surface
    c2min = ode_params.Control2_min
    c2max = ode_params.Control2_max

    # figure out the rotation rate at the q=2 surface (informational only —
    # the actual Control1/Ω0 sweep range is set via ode_params.Control1_min/max,
    # which default to [1e-2, 10] but can be overridden, e.g. for a focused
    # transfer-learning sweep)
    rho = dd.core_sources.source[1].profiles_1d[1].grid.rho_tor_norm
    rot_core = dd.core_profiles.profiles_1d[1].rotation_frequency_tor_sonic
    rot_interp = IMAS.interp1d(rho, rot_core)
    Omega0_calc = rot_interp(rt) * 2 * π * par.t0 # in units of kHz
    @info("Calculated toroidal rotation frequency at rational surface: $Omega0_calc kHz")
    println("Using Control1 (Ω0) range: [$(ode_params.Control1_min), $(ode_params.Control1_max)]")

    # finally set Control1: currently HARD-coded to be momentum source
    Control1_vals = range(ode_params.Control1_min, ode_params.Control1_max, length=N) |> collect
    ode_params.Control1 = vec(repeat(Control1_vals, 1, M))

    # Initialize the other control parameter based on the control type
    if control_type == :EF
        # For :EF, Control2_min/max are interpreted directly as the EF amplitude
        # in Gauss. Convert to the flux-equivalent perturbation psi_eps (still
        # referred to as "Eps" in the code), which carries units of b0*r0 (same
        # as psi0):
        #   EF_Tesla = EF_Gauss * 1e-4
        #   psi_eps  = -i * r_c * EF_Tesla / m_pol   (magnitude = r_c*EF_Tesla/m_pol)
        # The -i indicates psi_eps is -90deg out of phase with the true EF — not
        # yet propagated as an actual phase offset in the RHS (TODO, left as-is).
        EF_Gauss_vals = range(c2min, c2max, length=M) |> collect
        EF_Tesla_vals = EF_Gauss_vals .* 1.e-4 / b0   # Gauss -> Tesla -> dimensionless
        rc    = ode_params.control_surf
        m_pol = par.m_pol
        Control2_vals = rc .* EF_Tesla_vals ./ m_pol  # psi_eps, units of b0*r0

        @info("Maximum error field is $(c2max) Gauss")
    else
        Control2_vals = range(c2min, c2max, length=M) |> collect

        EpsUp = ode_params.error_field * par.b0 / 1.e-4  # Convert to Gauss
        @info("Maximum error field is $(EpsUp) Gauss")
    end

    Control2 = vec(repeat(Control2_vals', N, 1))
    ode_params.Control2 = Control2

    return ode_params
end


# ─────────────────────────────────────────────────────────────────────────────
#  Thin wrappers around ModeLocking — keep backward-compatible call sites
# ─────────────────────────────────────────────────────────────────────────────

"""
    solve_one_case(par, ode_params::ODEparams, application::String)

Solve and plot a single trajectory (task = :single_case). Delegates the ODE
solve to `ModeLocking.simulate_one_case` and the plotting to
`ModeLocking.plot_time_traces`.
"""
function solve_one_case(par, ode_params::ODEparams, application::String)
    sol, norm_t = ModeLocking.simulate_one_case(ode_params, application, par.n_tor, par.control_type,
                                                  par.source_torque, par.t_final, par.time_steps)

    # Time-dependent figures: TM (ψ_tN), RWM (ψ_wN, RP-RW only), and Ω_tN vs time
    fig = ModeLocking.plot_time_traces(norm_t, sol.t)

    # Same traces, but in physical units (Gauss for magnetic amplitudes, kHz for rotation)
    fig_phys = ModeLocking.plot_time_traces(sol, [par.b0, par.t0, par.r0])

    # Stack both figures in a single combined plot so the second doesn't
    # overtake/replace the first on display
    fig_combined = plot(fig, fig_phys; layout=(2, 1), size=(800, 800))
    display(fig_combined)

    return fig
end


"""
    _solve_grid_and_classify(actor, application) → LockingResults

Solve the ODE system over the current `actor.ode_params.Control1/Control2`
grid, normalize the results, and classify them via k-means into
locked/unlocked. When NL saturation is off, also computes the analytic
bifurcation boundary.

Returns a `LockingResults` with `prob = nothing` — callers fill that in
(e.g. via `train_locking_nn` for `:solve_system`, or
`transfer_learn_locking_nn` for `:transfer_learning`).
"""
function _solve_grid_and_classify(actor::ActorLocking, application::String)
    par = actor.par
    return ModeLocking.solve_and_classify(actor.ode_params, application, par.n_tor, par.control_type,
                                           par.t_final, par.time_steps, par.NL_saturation, par.grid_size)
end


"""
    train_locking_nn(actor, nn_params=NNparams()) → LockingNNModel

Train a binary NN classifier (C1, C2) → P(locked) on the k-means labels in
`actor.results`.  Updates `actor.results.prob` in-place and returns the model.

The stored model is callable: `actor.results.prob(C1, C2)` ∈ [0, 1].
To search for better hyperparameters first, call `tune_locking_nn(actor)`.
"""
function train_locking_nn(actor::ActorLocking, nn_params::NNparams=NNparams())
    actor.results === nothing && error("No results — run the actor with task=:solve_system first")
    return ModeLocking.train_locking_nn(actor.results, actor.ode_params, nn_params)
end


"""
    tune_locking_nn(actor; n_trials=20, n_folds=3, rng=GLOBAL_RNG) → NNparams

Random hyperparameter search using k-fold CV on the current actor results.
After finding the best configuration, retrains the full model with those
hyperparameters and updates `actor.results.prob` in-place.

Returns the best `NNparams` (useful for Task 2 transfer learning).
"""
function tune_locking_nn(actor::ActorLocking; kwargs...)
    actor.results === nothing && error("No results — run the actor with task=:solve_system first")
    return ModeLocking.tune_locking_nn(actor.results, actor.ode_params; kwargs...)
end


"""
    save_ode_results(actor; filename="ode_results.bson", dir=ModeLocking.LOCKING_RESULTS_DIR) → path

Save the ODE grid results to disk so that `task=:calc_prob` can be run in a
future session without re-solving the ODEs.
"""
function save_ode_results(actor::ActorLocking; kwargs...)
    actor.results === nothing && error("No ODE results — run task=:solve_system first")
    return ModeLocking.save_ode_results(actor.results, actor.ode_params;
        control_type=actor.par.control_type, kwargs...)
end


"""
    load_ode_results!(actor; filename="ode_results.bson", dir=ModeLocking.LOCKING_RESULTS_DIR)

Load previously saved ODE results from disk into `actor.results` and
`actor.ode_params.Control1/Control2`.  Called by `task=:calc_prob` when
`actor.results` is nothing (fresh session).
"""
function load_ode_results!(actor::ActorLocking; kwargs...)
    results, Control1, Control2 = ModeLocking.load_ode_results(;
        control_type=actor.par.control_type, kwargs...)
    actor.ode_params === nothing && (actor.ode_params = ODEparams())
    actor.ode_params.Control1 = Control1
    actor.ode_params.Control2 = Control2
    actor.results = results
    return actor.results
end


"""
    save_locking_nn(actor; filename="nn_model.bson", dir=ModeLocking.LOCKING_RESULTS_DIR) → path

Save the trained LockingNNModel to disk.
"""
function save_locking_nn(actor::ActorLocking; kwargs...)
    actor.results === nothing      && error("No results — run task=:solve_system first")
    actor.results.prob === nothing && error("No trained model — run train_locking_nn first")
    return ModeLocking.save_locking_nn(actor.results.prob; control_type=actor.par.control_type, kwargs...)
end


"""
    load_locking_nn(; filename="nn_model.bson", dir=ModeLocking.LOCKING_RESULTS_DIR) → LockingNNModel

Load a saved LockingNNModel from disk. Cached in memory after first load.
"""
load_locking_nn(; kwargs...) = ModeLocking.load_locking_nn(; kwargs...)


# ─────────────────────────────────────────────────────────────────────────────
#  Plotting wrappers
# ─────────────────────────────────────────────────────────────────────────────

"""
    plot_sols(actor) → (fig1, fig2, fig3)

Calls plot_scatter, plot_phase_diagrams, and (when a trained NN model is
available) plot_probability.  Returns all three handles; fig3 is nothing
when no model has been trained yet.
"""
plot_sols(actor::ActorLocking) = ModeLocking.plot_sols(actor.results, actor.ode_params, actor.par.grid_size, actor.par.control_type;
    b0=actor.par.b0, t0=actor.par.t0, m_pol=Float64(actor.par.m_pol), orientation=actor.par.plot_orientation)

"plot_scatter(actor) → Figure 1 — multi-panel scatter of state variables"
plot_scatter(actor::ActorLocking) = ModeLocking.plot_scatter(actor.results; orientation=actor.par.plot_orientation)

"plot_phase_diagrams(actor) → Figure 2 — pcolor of Ω_n and ψ_tn over control space"
plot_phase_diagrams(actor::ActorLocking) = ModeLocking.plot_phase_diagrams(actor.results, actor.ode_params, actor.par.grid_size, actor.par.control_type;
    b0=actor.par.b0, t0=actor.par.t0, m_pol=Float64(actor.par.m_pol))

"plot_probability(actor) → Figure 3 — NN locking probability"
plot_probability(actor::ActorLocking) = ModeLocking.plot_probability(actor.results, actor.ode_params, actor.par.control_type;
    b0=actor.par.b0, t0=actor.par.t0, m_pol=Float64(actor.par.m_pol))
