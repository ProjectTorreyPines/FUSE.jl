using Distributed
using DifferentialEquations
import Roots
using Flux
using Random
import BSON
using Plots
using Clustering
#import FUSE: coordinates


#================== =#
#  ActorLockingProbability  #
#================== =#
Base.@kwdef mutable struct FUSEparameters__ActorLocking{T<:Real} <: ParametersActor{T}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    _time::Float64 = NaN
    n_mode::Entry{Int} = Entry{Int}("_", "toroidal mode number of the mode"; default=1)
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
    NL_saturation_ON::Entry{Bool} = Entry{Bool}("-", "Nonlinear saturation parameter for the mode"; default=false)
    RPRW_stability_index::Entry{Float64} = Entry{Float64}(
        "-", 
        "Stability index of the system (set to Neg. value for now)"; default=-0.5)
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
end


# Make sure workers know the ODEparams layout so they can receive an instance.
# (Same field order/types as your definition.)
#if !isdefined(Main, :ODEparams)
Base.@kwdef mutable struct ODEparams
    # --- user-set physical parameters (sensible defaults) ---
    saturation_param::Float64 = 0.2   # controls nonlinear island saturation
    error_field::Float64      = 0.5   # fixed error field when control_type is NOT :EF
    layer_width::Float64      = 1e-3  # resistive layer width from tearing theory (~1 mm)
    rat_surface::Float64      = 0.67  # q=2 surface location (dimensionless)
    res_wall::Float64         = 1.0   # resistive wall location (dimensionless)
    control_surf::Float64     = 1.25  # control surface location (dimensionless)
    mu::Float64               = 0.1   # anomalous perpendicular plasma viscosity
    Inertia::Float64          = 0.1   # moment of inertia of the layer
    Taut_Tauw::Float64        = 1.0   # ratio of tearing time to wall time
    hyper_cube_dims::Vector{Float64} = [1., 1., 1.] # initial condition hypercube dimensions
    Control1_min::Float64     = 1.0e-2
    Control1_max::Float64     = 10.0
    Control2_min::Float64     = 0.01
    Control2_max::Float64     = 1.0

    # --- computed at run time by calculate_stability_index! ---
    DeltaW::Float64          = NaN  # intrinsic RW stability — calculated at run time
    stability_index::Float64 = NaN  # tearing mode Delta' — calculated at run time
    l12::Float64             = NaN  # mutual inductance rational surface → RW — calculated at run time
    l21::Float64             = NaN  # mutual inductance RW → rational surface — calculated at run time
    l32::Float64             = NaN  # mutual inductance control surface → RW — calculated at run time

    # --- populated by set_control_parameters! ---
    Control1::Vector{Float64} = Float64[]  # rotation frequency grid — populated at run time
    Control2::Vector{Float64} = Float64[]  # swept control parameter grid — populated at run time
end

"Hyperparameters for the locking NN classifier"
Base.@kwdef struct NNparams
    hidden_sizes::Vector{Int}  = [100, 100, 100]  # neurons in each hidden layer
    activation::Symbol         = :relu             # :tanh, :relu, :sigmoid
    learning_rate::Float64     = 1e-3
    n_epochs::Int              = 1000
    batch_size::Int            = 200
    weight_decay::Float64      = 1e-8              # L2 regularisation coefficient
    val_fraction::Float64      = 0.0               # fraction held out for early stopping; 0 = disabled
    patience::Int              = 20                # early-stopping patience (epochs); ignored when val_fraction=0
end

"Trained NN model; callable as prob(C1, C2) → P(locked) ∈ [0,1]"
struct LockingNNModel
    model::Any          # Flux Chain — kept as Any for Task 2 transfer learning
    nn_params::NNparams
end

function (m::LockingNNModel)(C1::Real, C2::Real)
    return Float64(first(m.model(Float32[C1, C2])))
end

mutable struct LockingResults
    ode_sols::Matrix{Float64}                        # (N*M × n_states) raw final states
    prob::Any                                        # NN model once Task 1 is done; callable as prob(C1, C2)
    norm_sols::Matrix{Float64}                       # (N*M × n_states) normalized solutions
    locking_labels::Vector{Int}                      # k-means class assignments, one per grid point
    bifurcation_bounds::Union{Matrix{Float64}, Nothing}  # nothing when NL saturation is active
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
        @info "Solveing one case for system"
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
        prob_model = load_locking_nn()
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
        base_model = load_locking_nn()

        X_new, y_new = prepare_nn_data(actor.results.locking_labels,
                                        actor.ode_params.Control1, actor.ode_params.Control2)

        actor.results.prob = transfer_learn_locking_nn(base_model, X_new, y_new; nn_params=actor.nn_params)


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
    m0 = par.n_mode
    
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
    mu0_val = 4*pi*1.e-7  # Physical constant
    mass_ion = 1.67e-27  # kg, mass of ion (approx. for deuterium)

    rt = ode_params.rat_surface
    rho = dd.core_sources.source[1].profiles_1d[1].grid.rho_tor_norm
    mass_dens = mass_ion * dd.core_profiles.profiles_1d[1].ion[1].density_thermal
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
    if par.NL_saturation_ON == false
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

    Control1_vals = range(ode_params.Control1_min, ode_params.Control1_max, length=N) |> collect
    ode_params.Control1 = vec(repeat(Control1_vals, 1, M))

    # Initialize the other control parameter based on the control type
    Control2_vals = range(c2min, c2max, length=M) |> collect
    if control_type == :EF
        EpsUp = c2max * par.b0 * 1.e-4  # Convert to Gauss
        @info("Maximum error field is $(EpsUp) Gauss")

    elseif control_type == :LinStab
        #ode_params.error_field = 0.6  # Example value for error field
        ## Check to make sure the system is still weakly stable
        DeltatRW = Control2_vals .- l21*l12/DeltaW
        if any(DeltatRW .> 0)
            println("*** ALERT: You set up a case with an unstable RP-RW mode! ***")
            println("*** Maximum RP-RW stability set to ", maximum(DeltatRW))
            error("Deltat_RW > 0 for your range of TM stability values ***")
        end
        ode_params.stability_index = Control2_vals
    elseif control_type == :NLsaturation
        println("Do nothing for NL saturation control")
        
        ode_params.saturation_param = Control2
    end
    
    Control2 = vec(repeat(Control2_vals', N, 1))
    ode_params.Control2 = Control2

    return ode_params
end



# ─────────────────────────────────────────────────────────────────────────────
#  task = :single_case — solve and plot one trajectory
# ─────────────────────────────────────────────────────────────────────────────

function solve_one_case(par, ode_params::ODEparams, application::String)

    control1 = par.source_torque
   
    # harvest what's already under the hood to set these inputs 
    if par.control_type == :EF
        control2 = ode_params.error_field
    elseif par.control_type == :LinStab
        control2 = ode_params.stability_index - 0.5 # small adjustment to move away from marginality
    elseif par.control_type == :NLsaturation
        control2 = ode_params.saturation_param
    else
        @info "Control scenario NOT set, guessing the value of control2"
        control2 = 0.5
    end

    sol = solve_ODEs(par, ode_params, application, control1, control2; full_output=true)

    norm_t = reduce(vcat, (normalize_ode_results(u, ode_params, control2, control1, par.control_type)'
                            for u in sol.u))
    println("final normalized solution = ", norm_t[end, :])

    # Time-dependent figures: TM (ψ_tN), RWM (ψ_wN, RP-RW only), and Ω_tN vs time
    fig = plot_time_traces(norm_t, sol.t)
    display(fig)

    return fig
end



"""
    plot_time_traces(norm_t, t) → Plot

Plot the time-dependent, normalized tearing-mode amplitude (TM, ψ_tn) and
resistive-wall-mode amplitude (RWM, ψ_wn — only for the RP-RW system) on the
left y-axis, and the normalized rotation frequency (Ω_tN) vs time on a
secondary (right) y-axis, all on a single labeled plot.

`norm_t` is the (n_times × n_states) matrix of per-timestep normalized
states (as produced by `normalize_ode_results` applied to each `sol.u[i]`),
and `t` is the corresponding vector of times (`sol.t`).
"""
function plot_time_traces(norm_t::AbstractMatrix{<:Real}, t::AbstractVector{<:Real})

    psi_tN = norm_t[:, 1]
    Om_tN  = norm_t[:, 3]

    plt = plot(t, psi_tN; xlabel="time", ylabel="Normalized ψ", label="ψ_tN (TM)",
               lw=2, color=:steelblue, title="Time-dependent traces", legend=:topleft)

    if size(norm_t, 2) == 5 # RP-RW layout includes the resistive-wall-mode state
        psi_wN = norm_t[:, 4]
        plot!(plt, t, psi_wN; label="ψ_wN (RWM)", lw=2, color=:darkorange)
    end

    # Ω_tN on its own scale on the right y-axis
    plot!(twinx(plt), t, Om_tN; ylabel="Ω_tN", label="Ω_tN", lw=2,
          color=:firebrick, linestyle=:dash, legend=:topright)

    return plt
end



# ─────────────────────────────────────────────────────────────────────────────
#  task = :solve_system / :transfer_learning — grid solve + classification
# ─────────────────────────────────────────────────────────────────────────────

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
    dd  = actor.dd
    par = actor.par

    control1 = actor.ode_params.Control1
    control2 = actor.ode_params.Control2
    ode_sols = solve_system(actor, application)
    norm_sols = normalize_ode_results(ode_sols, actor.ode_params,
         control2, control1, par.control_type)

    ## classify normalized solutions
    R = hcat(norm_sols[:,1], norm_sols[:,3])
    kmc = kmeans(R', 2)
    locking_labels = kmc.assignments

    # Fix label polarity: the point with maximum OmN (col 3) is always unlocked.
    # k-means labels are {1,2}; we want that point to carry label 1 so that
    # after the {1,2}→{0,1} shift in prepare_nn_data it maps to 0 (unlocked).
    if locking_labels[argmax(norm_sols[:, 3])] != 1
        locking_labels = 3 .- locking_labels   # flip 1↔2
    end

    # Analytic bifurcation boundary — not defined when NL saturation is active
    bifurcation_bounds = par.NL_saturation_ON ? nothing :
        calculate_bifurcation_bounds(dd, par, actor.ode_params)

    return LockingResults(
        ode_sols,
        nothing,
        norm_sols,
        locking_labels,
        bifurcation_bounds
    )
end



function solve_system(actor::ActorLocking, task::String)

    println("Solving the FULL system, this may take a few seconds")

    par = actor.par
    ode_params = actor.ode_params

    n_mode = par.n_mode
    control_type = par.control_type
    t_final = par.t_final
    
    # Control1 = C1 (rotation, Y-axis), Control2 = C2 = EF/Δ′/α (X-axis)
    control1 = ode_params.Control1
    control2 = ode_params.Control2
    inputs = collect(zip(control1, control2))   # (C1, C2) pairs — natural order

    # Shrink what we ship to workers by clearing large control arrays
    ode_params_send = deepcopy(ode_params)
    ode_params_send.Control1 = Float64[]
    ode_params_send.Control2 = Float64[]

    # Parallel map over the grid, returning final states as a (N*M × n_states) matrix
    finals = pmap(inputs) do (C1, C2)
        solve_ODEs(par, ode_params_send, task, C1, C2)
    end

    return Matrix(reduce(hcat, finals)')  # Vector{Vector} → Matrix{Float64} (N*M × n_states)
end



"""
    calculate_bifurcation_bounds(dd, par, ode_params) → Matrix{Float64}

Compute the cubic-discriminant bifurcation boundary analytically over the full
control grid.  Returns a (grid_size × grid_size) matrix; negative entries mark
the parameter region where locking is possible.

Two wall configurations are supported via `par.application`:
  - "RP-RW"  (5th-order system): uses DeltatRW and the full RW geometry factors.
  - "RP-IW"  (3rd-order system): uses Deltat directly, no wall inductance factors.

Note: result is not meaningful when NL saturation is active.
"""
function calculate_bifurcation_bounds(dd::IMAS.dd, par, ode_params::ODEparams)
    m0     = par.n_mode
    mu     = ode_params.mu
    rt     = ode_params.rat_surface
    l21    = ode_params.l21
    l32    = ode_params.l32
    DeltaW = ode_params.DeltaW

    # resolve_control handles all control-type branching; wall branching is separate
    sc = resolve_control(ode_params, par)
    (; Y, Deltat, eps, DeltatRW) = sc

    # Cubic-discriminant coefficients — formula depends on wall configuration
    if par.application == "RP-RW"
        # 5th-order system: rational surface + resistive wall
        q = (DeltatRW ./ m0).^2 .+ rt .* (l32 .* l21 .* eps ./ DeltaW).^2 ./ (m0 * mu)
        r = -Y .* DeltatRW.^2 ./ m0^2
    elseif par.application == "RP-IW"
        # 3rd-order system: rational surface + ideal wall only
        q = (Deltat ./ m0).^2 .+ (l21 .* eps).^2 ./ mu
        r = -Y .* Deltat.^2 ./ m0^2
    else
        error("calculate_bifurcation_bounds not implemented for application: $(par.application)")
    end

    a = -(Y.^2) ./ 3.0 .+ q
    b = 2.0 .* (-Y).^3 ./ 27.0 .- q .* (-Y) ./ 3.0 .+ r

    bifurcation_bounds = reshape(b.^2 ./ 4 .+ a.^3 ./ 27, par.grid_size, par.grid_size)

    return bifurcation_bounds
end



# ─────────────────────────────────────────────────────────────────────────────
#  Shared low-level ODE solving and normalization
# ─────────────────────────────────────────────────────────────────────────────

function solve_ODEs(par, ode_params::ODEparams, task::String, C1::Float64, C2::Float64; full_output::Bool=false)
    n_mode       = par.n_mode
    control_type = par.control_type
    t_final      = par.t_final

    rhs!     = make_rhs_function(task)
    ode_rhs! = make_ode_func(rhs!)
    y0 = make_initial_condition(ode_params.hyper_cube_dims, task)

    # rhs! expects (C2, C1) as its first two positional args — swap once here so
    # all callers use the natural (C1, C2) convention
    p = (C2, C1, ode_params, n_mode, control_type)
    tspan = (0.0, t_final)

    prob = ODEProblem(ode_rhs!, y0, tspan, p)

    if full_output
        # Return the full time-resolved solution (used for time-dependent plots)
        tsave = range(0.0, t_final; length=par.time_steps)
        sol = solve(prob, Tsit5(); saveat=tsave, reltol=1e-8, abstol=1e-10)
        return sol
    else
        # Only the final state is needed (e.g. grid scans for classification/NN)
        sol = solve(prob, Tsit5(); saveat=t_final, reltol=1e-8, abstol=1e-10)
        return sol.u[end]
    end
end


"""
normalize_ode_results(results, ode_params, eps_vec, C1_vec, control_type)

Normalize the final solutions from ODE runs.

- `results` may be:
    * a single solution vector (Vector{Float64}), or
    * a Vector of solution vectors (Vector{Vector{Float64}}).
- `eps_vec` (Control2 values) and `C1_vec` (Control1 values) may be:
    * single Float64 values (if results is one vector), or
    * Vector{Float64} of the same length as results.

Normalization:
    psiN  = final_sol[1] * (Deltat * DeltaW - l12 * l21) / (l32 * l21 * eps)
    psiwN = final_sol[4] * (Deltat * DeltaW - l12 * l21) / (l32 * abs(Deltat) * eps)
    OmN   = final_sol[3] / C1
"""
function normalize_ode_results(results, ode_params::ODEparams, C2_vec, C1_vec, control_type)
    # Extract parameters
    l12    = ode_params.l12
    l21    = ode_params.l21
    l32    = ode_params.l32
    DeltaW = ode_params.DeltaW

    
    # Normalization for one solution — control branching handled by resolve_control
    function normalize_one(final_sol::AbstractVector{<:Real}, C2::Float64, C1::Float64)
        sc       = resolve_control(ode_params, control_type, C2)
        Deltat   = sc.Deltat
        eps      = sc.eps
        alpha    = sc.alpha
        DeltatRW = sc.DeltatRW

        if length(final_sol) == 5 # RP-RW layout
            psit    = final_sol[1]
            theta_t = mod(final_sol[2], 2π)
            OmN     = final_sol[3] / C1
            psiw    = final_sol[4]
            theta_w = mod(final_sol[5], 2π)
            rho = abs(DeltatRW / Deltat)

            if iszero(alpha)  # linear regime: saturation_param set to 0. when NL_saturation_ON=false
                num     = abs(Deltat * DeltaW) - l12 * l21
                psitMax = l32 * l21 * eps / num
                psiwMax = l32 * abs(Deltat) * eps / abs(DeltatRW * DeltaW)
                psiwMin = l32 * eps / abs(DeltaW)
            else              # NL saturation active
                psitMax = -(DeltatRW + sqrt(DeltatRW^2 + 4*alpha*l21*l32*eps*Deltat/DeltaW)) /
                           (2*alpha*Deltat)
                psiwMax = -(l32*eps + l12*psitMax) / DeltaW
                psiwMin = l32 * eps / abs(DeltaW)
            end

            psiN  = abs(psit / psitMax)
            psiwN = abs(psiw / psiwMax)
            psiwN = (psiwN - rho) / (1 - rho)
            #psiwN = abs((psiw - psiwMin) / (psiwMax - psiwMin))

            return [psiN, theta_t, OmN, psiwN, theta_w]

        elseif length(final_sol) == 3 # RP-IW layout
            psit    = final_sol[1]
            theta_t = mod(final_sol[2], 2π)
            OmN     = final_sol[3] / C1

            if iszero(alpha)  # linear regime: saturation_param set to 0. when NL_saturation_ON=false
                psitMax = l21 * eps / abs(Deltat)
            else              # NL saturation active
                psitMax = (-1.0 + sqrt(1.0 - 4*l21*alpha*eps/Deltat)) / (2*alpha)
            end

            psiN = abs(psit / psitMax)

            return [psiN, theta_t, OmN]

        else
            throw(ArgumentError("Unexpected final_sol length: $(length(final_sol))"))
        end
    end

    if isa(results, AbstractVector{<:Real}) && isa(C2_vec, Real) && isa(C1_vec, Real)
        # Single case
        return normalize_one(results, C2_vec, C1_vec)

    elseif isa(results, AbstractMatrix{<:Real}) &&
           isa(C2_vec, AbstractVector{<:Real}) &&
           isa(C1_vec, AbstractVector{<:Real})
        # Many cases — results is (N*M × n_states), iterate over rows
        size(results, 1) == length(C2_vec) == length(C1_vec) ||
            throw(ArgumentError("results, C2_vec, and C1_vec must all have the same length"))
        return reduce(vcat, (normalize_one(sol, e, c1)'
              for (sol, e, c1) in zip(eachrow(results), C2_vec, C1_vec)))

    else
        throw(ArgumentError("Input types do not match expected patterns"))
    end
end



# ─────────────────────────────────────────────────────────────────────────────
#  ODE right-hand-side machinery
# ─────────────────────────────────────────────────────────────────────────────

function make_rhs_function(phys_type::String)
    phys_type == "RP-RW" ? rhs_RW! :
    phys_type == "RP-IW"   ? rhs_basic! :
    error("Unknown application type: $phys_type")
   
end


function rhs_RW!(dydt, y, t, Control2::Float64, C1::Float64, ode_params::ODEparams, n_mode::Int, control_type::Symbol)
    m0 = n_mode
    DeltaW = ode_params.DeltaW
    rt = ode_params.rat_surface
    l21 = ode_params.l21
    l12 = ode_params.l12
    l32 = ode_params.l32
    Tt_Tw = ode_params.Taut_Tauw
    mu = ode_params.mu
    I = ode_params.Inertia

    alpha, errF, Deltat = control_adjustments(ode_params, Control2, control_type)

    psi, theta, Om, psiW, thW = y

    dydt[1] = Deltat * psi * (1.0 + alpha * abs(psi)) + l21 * psiW * cos(theta - thW)
    dydt[2] = -m0 * Om - l21 * psiW * sin(theta - thW) / psi
    dydt[3] = (rt * l21 * psiW * psi * sin(theta - thW) + mu * (C1 - Om)) / I
    dydt[4] = Tt_Tw * (DeltaW * psiW + l12 * psi * cos(theta - thW) + l32 * errF * cos(thW))
    dydt[5] = Tt_Tw * (l12 * psi * sin(theta - thW) - l32 * errF * sin(thW)) / psiW
end

"Right-hand side for basic system"

function rhs_basic!(dydt, y, t, Control2::Float64, C1::Float64, ode_params::ODEparams, n_mode::Int, control_type::Symbol)
    m0 = n_mode
    rt = ode_params.rat_surface
    l21 = ode_params.l21
    mu = ode_params.mu
    I = ode_params.Inertia

    alpha, errF, Deltat = control_adjustments(ode_params, Control2, control_type)

    psi, theta, Om = y

    dydt[1] = Deltat * psi * (1.0 + alpha * abs(psi)) + l21 * errF * cos(theta)
    dydt[2] = -m0 * Om - l21 * errF * sin(theta) / psi
    dydt[3] = (rt * l21 * errF * psi * sin(theta) + mu * (C1 - Om)) / I
end

"Return the correct RHS function for a given application"

function make_initial_condition(dims::Vector{Float64}, application::String)

    if application == "RP-RW"
        # 5D system
        return [
            rand() * (dims[1] - 0.001) + 0.001,   # y1
            rand() * 2π - π,                      # y2
            rand() * (dims[2] - 0.001) + 0.001,   # y3
            rand() * (dims[3] - 0.001) + 0.001,   # y4
            rand() * 2π - π                       # y5
        ]

    elseif application == "RP-IW"
        # 3D system
        return [
            rand() * (dims[1] - 0.001) + 0.001,   # y1
            rand() * 2π - π,                      # y2
            rand() * (dims[2] - 0.001) + 0.001    # y3
        ]

    else
       error("Unknown application: $application")
    end
end


function make_ode_func(rhs!)
    return function (du, u, p, t)
        # p layout: (C2, C1, ode_params, n_mode, control_type)
        # C2 = swept control param (EF / Δ′ / α); C1 = rotation frequency
        # This matches the argument order of rhs_RW! and rhs_basic!
        C2, C1, ode_params, n_mode, control_type = p
        rhs!(du, u, t, C2, C1, ode_params, n_mode, control_type)
    end
end



"""
    resolve_control(ode_params, control_type, C2) → NamedTuple

Scalar resolver: map a single raw Control2 value to concrete physical quantities.
This is the *only* place in the code where control-type branching lives.

Returns `(; Deltat, eps, alpha, DeltatRW)`.
- `Deltat`    : effective tearing stability index
- `eps`       : error-field amplitude
- `alpha`     : nonlinear saturation parameter
- `DeltatRW`  : effective resistive-wall stability  (Deltat − l21·l12/ΔW)
"""
function resolve_control(ode_params::ODEparams, control_type::Symbol, C2::Real)
    Deltat = control_type == :EF           ? ode_params.stability_index :
             control_type == :LinStab      ? Float64(C2) :
             control_type == :NLsaturation ? ode_params.stability_index :
             error("Unknown control_type: $control_type")

    eps    = control_type == :EF           ? Float64(C2) :
             control_type == :LinStab      ? ode_params.error_field :
             control_type == :NLsaturation ? ode_params.error_field :
             error("Unknown control_type: $control_type")

    alpha  = control_type == :NLsaturation ? Float64(C2) : ode_params.saturation_param

    DeltatRW = Deltat - ode_params.l21 * ode_params.l12 / ode_params.DeltaW

    return (; Deltat, eps, alpha, DeltatRW)
end


"""
    resolve_control(ode_params, par) → NamedTuple

Vector resolver: operates over the full Control1/Control2 grid stored in `ode_params`.
Returns `(; X, Y, Deltat, eps, alpha, DeltatRW)` — all `Vector{Float64}`.
"""
function resolve_control(ode_params::ODEparams, par)
    ctrl = par.control_type
    X    = ode_params.Control2   # swept control axis
    Y    = ode_params.Control1   # rotation-frequency axis
    n    = length(X)

    Deltat = ctrl == :EF           ? fill(ode_params.stability_index, n) :
             ctrl == :LinStab      ? copy(X) :
             ctrl == :NLsaturation ? fill(ode_params.stability_index, n) :
             error("Unknown control_type: $ctrl")

    eps    = ctrl == :EF           ? copy(X) :
             ctrl == :LinStab      ? fill(ode_params.error_field, n) :
             ctrl == :NLsaturation ? fill(ode_params.error_field, n) :
             error("Unknown control_type: $ctrl")

    alpha  = ctrl == :NLsaturation ? copy(X) : fill(ode_params.saturation_param, n)

    DeltatRW = @. Deltat - ode_params.l21 * ode_params.l12 / ode_params.DeltaW

    return (; X, Y, Deltat, eps, alpha, DeltatRW)
end

"Apply control_type adjustments — delegates to resolve_control"

function control_adjustments(ode_params::ODEparams, C2::Float64, control_type::Symbol)
    sc = resolve_control(ode_params, control_type, C2)
    return sc.alpha, sc.eps, sc.Deltat
end

"Right-hand side for RW system"

# ─────────────────────────────────────────────────────────────────────────────
#  Neural-network classifier (Task 1 / Task 2)
# ─────────────────────────────────────────────────────────────────────────────

"""
    train_locking_nn(actor, nn_params=NNparams()) → LockingNNModel

Train a binary NN classifier (C1, C2) → P(locked) on the k-means labels in
`actor.results`.  Updates `actor.results.prob` in-place and returns the model.

The stored model is callable: `actor.results.prob(C1, C2)` ∈ [0, 1].
To search for better hyperparameters first, call `tune_locking_nn(actor)`.
"""
function train_locking_nn(actor::ActorLocking, nn_params::NNparams=NNparams())
    r  = actor.results
    op = actor.ode_params
    r === nothing && error("No results — run the actor with task=:solve_system first")

    X, y = prepare_nn_data(r.locking_labels, op.Control1, op.Control2)

    @info "Training NN: architecture=$(nn_params.hidden_sizes)  epochs=$(nn_params.n_epochs)"
    model = _fit_nn_model(build_locking_nn(nn_params), X, y, nn_params)
    @info "NN training complete. Final loss=$(round(Flux.binarycrossentropy(model(X), y); digits=4))"

    prob_model = LockingNNModel(model, nn_params)
    actor.results.prob = prob_model
    return prob_model
end


"""
    prepare_nn_data(locking_labels, control1, control2) → (X, y)

Prepare (X, y) training data from ODE results.
  X : (2 × N) Float32 — raw [C1, C2] values
  y : (1 × N) Float32 — k-means labels shifted from {1,2} → {0,1}
"""
function prepare_nn_data(locking_labels::Vector{Int},
                          control1::Vector{Float64}, control2::Vector{Float64})
    X = Float32.(hcat(control1, control2)')   # (2 × N)
    y = reshape(Float32.(locking_labels .- 1), 1, :)   # (1 × N) Matrix{Float32}, {1,2} → {0,1}
    return X, y
end


function build_locking_nn(nn_params::NNparams)
    act   = get_activation(nn_params.activation)
    sizes = [2; nn_params.hidden_sizes; 1]
    layers = []
    for i in 1:length(sizes)-2
        push!(layers, Dense(sizes[i], sizes[i+1], act))
    end
    push!(layers, Dense(sizes[end-1], 1, sigmoid))
    return Chain(layers...)
end


# ─────────────────────────────────────────────────────────────────────────────
#  Neural-network classifier  (Task 1)
# ─────────────────────────────────────────────────────────────────────────────

"Map activation symbol to the corresponding Flux function"
function get_activation(act::Symbol)
    act == :tanh    ? tanh    :
    act == :relu    ? relu    :
    act == :sigmoid ? sigmoid :
    error("Unknown activation symbol: $act.  Choose :tanh, :relu, or :sigmoid")
end

"Build a Flux Chain from NNparams: 2 inputs → hidden layers → 1 sigmoid output"

"""
Train a Flux model on (X, y) with L2 regularisation.
When `nn_params.val_fraction > 0`, a validation split is held out and
early stopping is applied (patience = `nn_params.patience`).
When `val_fraction == 0` (default), all data is used and training runs
for the full `n_epochs` — appropriate for physics-driven classification
problems where the labels are deterministic.
Returns the best model found (or the final model when validation is off).

An optional pre-built `opt_state` can be supplied (e.g. with some layers
frozen via `Flux.freeze!`, as used by `transfer_learn_locking_nn`); if
omitted, a fresh `Adam` optimizer state is created for the whole model.
"""
function _fit_nn_model(model, X::Matrix{Float32}, y::Matrix{Float32}, nn_params::NNparams;
                        opt_state = Flux.setup(Adam(nn_params.learning_rate), model))
    use_val = nn_params.val_fraction > 0.0

    if use_val
        n_val   = max(1, round(Int, nn_params.val_fraction * size(X, 2)))
        n_train = size(X, 2) - n_val
        X_tr = X[:, 1:n_train];     y_tr = y[:, 1:n_train]
        X_va = X[:, n_train+1:end]; y_va = y[:, n_train+1:end]
    else
        X_tr = X;  y_tr = y
    end

    data       = Flux.DataLoader((X_tr, y_tr); batchsize=nn_params.batch_size, shuffle=true)
    wd         = Float32(nn_params.weight_decay)
    best_val   = Inf
    n_wait     = 0
    best_model = model

    for epoch in 1:nn_params.n_epochs
        for (xb, yb) in data
            _, grads = Flux.withgradient(model) do m
                l2 = sum(l -> sum(abs2, l.weight), m.layers)
                Flux.binarycrossentropy(m(xb), yb) + wd * l2
            end
            Flux.update!(opt_state, model, grads[1])
        end

        if use_val
            val = Flux.binarycrossentropy(model(X_va), y_va)
            if val < best_val
                best_val = val;  n_wait = 0;  best_model = deepcopy(model)
            else
                n_wait += 1
                if n_wait >= nn_params.patience
                    @info "  Early stopping at epoch $epoch  best_val=$(round(best_val; digits=4))"
                    break
                end
            end
            epoch % 200 == 0 && @info "  NN epoch $epoch/$(nn_params.n_epochs)  val=$(round(val; digits=4))"
        else
            epoch % 200 == 0 && @info "  NN epoch $epoch/$(nn_params.n_epochs)  loss=$(round(Flux.binarycrossentropy(model(X), y); digits=4))"
        end
    end
    return best_model
end

# Local directory for all locking actor results (ODE grid + NN model)

"""
    transfer_learn_locking_nn(base_model, X_new, y_new; nn_params=base_model.nn_params) → LockingNNModel

Fine-tune `base_model` (typically loaded via `load_locking_nn()`) on new
`(X_new, y_new)` data from a different equilibrium/`dd`, freezing every
layer except the last (output) `Dense` layer.

`nn_params` controls the fine-tuning run (learning rate, epochs, etc.) —
pass a separate, typically gentler, `NNparams` than the one used for the
original full training. Returns a new `LockingNNModel`; does not mutate
`base_model`.
"""
function transfer_learn_locking_nn(base_model::LockingNNModel,
                                    X_new::Matrix{Float32}, y_new::Matrix{Float32};
                                    nn_params::NNparams = base_model.nn_params)
    model     = deepcopy(base_model.model)
    opt_state = Flux.setup(Adam(nn_params.learning_rate), model)

    # Freeze every layer except the last (output) Dense layer
    for i in 1:length(model.layers)-1
        Flux.freeze!(opt_state.layers[i])
    end

    @info "Transfer learning: fine-tuning last layer only ($(length(model.layers)) total layers, $(nn_params.n_epochs) epochs)"
    fine_tuned = _fit_nn_model(model, X_new, y_new, nn_params; opt_state=opt_state)
    @info "Transfer learning complete. Final loss=$(round(Flux.binarycrossentropy(fine_tuned(X_new), y_new); digits=4))"

    return LockingNNModel(fine_tuned, nn_params)
end


"""
    tune_locking_nn(actor; n_trials=20, n_folds=3, rng=GLOBAL_RNG) → NNparams

Random hyperparameter search using k-fold CV on the current actor results.
After finding the best configuration, retrains the full model with those
hyperparameters and updates `actor.results.prob` in-place.

Returns the best `NNparams` (useful for Task 2 transfer learning).
"""
function tune_locking_nn(actor::ActorLocking; n_trials::Int=20, n_folds::Int=3,
                          rng::AbstractRNG=Random.GLOBAL_RNG)
    r  = actor.results
    op = actor.ode_params
    r === nothing && error("No results — run the actor with task=:solve_system first")

    X, y = prepare_nn_data(r.locking_labels, op.Control1, op.Control2)

    # Search space — mirrors Python's RandomizedSearchCV param_dist
    hidden_pool = [[10,10], [10,20,10], [100,100], [200,100,100], [200,100,100,200]]
    lr_pool     = exp10.(range(-4, -2; length=5))   # logspace(-4,-2,5)
    wd_pool     = exp10.(range(-10, -6; length=5))  # logspace(-10,-6,5)
    batch_pool  = [100, 200, 400, 800]
    epoch_pool  = [400, 1000]                        # CV uses half; final uses full

    N         = size(X, 2)
    fold_size = N ÷ n_folds
    best_params = NNparams()
    best_loss   = Inf

    @info "Hyperparameter search: $n_trials trials, $n_folds-fold CV"
    for trial in 1:n_trials
        params = NNparams(
            hidden_sizes  = rand(rng, hidden_pool),
            activation    = :relu,                       # Python fixes relu
            learning_rate = rand(rng, lr_pool),
            n_epochs      = rand(rng, epoch_pool) ÷ 2,  # half epochs for CV speed
            batch_size    = rand(rng, batch_pool),
            weight_decay  = rand(rng, wd_pool),
            patience      = 20,
        )

        # k-fold cross-validation loss
        cv_loss = 0.0
        for k in 1:n_folds
            val_idx   = ((k-1)*fold_size + 1):min(k*fold_size, N)
            train_idx = setdiff(1:N, val_idx)
            m = _fit_nn_model(build_locking_nn(params), X[:, train_idx], y[:, train_idx], params)
            cv_loss += Flux.binarycrossentropy(m(X[:, val_idx]), y[:, val_idx])
        end
        cv_loss /= n_folds

        @info "  Trial $trial/$n_trials  cv_loss=$(round(cv_loss; digits=4))  arch=$(params.hidden_sizes)  lr=$(round(params.learning_rate;sigdigits=2))  wd=$(round(params.weight_decay;sigdigits=2))"
        if cv_loss < best_loss
            best_loss   = cv_loss
            best_params = params
        end
    end

    @info "Best: arch=$(best_params.hidden_sizes)  act=$(best_params.activation)  lr=$(best_params.learning_rate)  CV_loss=$(round(best_loss; digits=4))"

    # Retrain on full data restoring full epoch count (CV used half)
    full_params = NNparams(
        hidden_sizes  = best_params.hidden_sizes,
        activation    = best_params.activation,
        learning_rate = best_params.learning_rate,
        n_epochs      = best_params.n_epochs * 2,
        batch_size    = best_params.batch_size,
        weight_decay  = best_params.weight_decay,
        patience      = best_params.patience,
    )
    @info "Retraining final model ($(full_params.n_epochs) epochs)..."
    train_locking_nn(actor, full_params)
    return full_params
end




# ─────────────────────────────────────────────────────────────────────────────
#  Persistence — ODE results & NN model checkpointing
# ─────────────────────────────────────────────────────────────────────────────

const LOCKING_RESULTS_DIR = joinpath(homedir(), ".julia", "locking_results")

# In-memory cache — avoids reloading from disk on every evaluate_probability call
const _locking_nn_cache = Dict{String, LockingNNModel}()


"""
    save_ode_results(actor; filename="ode_results.bson", dir=LOCKING_RESULTS_DIR) → path

Save the ODE grid results to disk so that `task=:calc_prob` can be run in a
future session without re-solving the ODEs.  Saves: ode_sols, norm_sols,
locking_labels, bifurcation_bounds, Control1, Control2.
"""
function save_ode_results(actor::ActorLocking;
                           filename::String = "ode_results.bson",
                           dir::String      = LOCKING_RESULTS_DIR)
    actor.results === nothing && error("No ODE results — run task=:solve_system first")
    mkpath(dir)
    path = joinpath(dir, filename)
    r    = actor.results
    op   = actor.ode_params
    ode_sols           = r.ode_sols
    norm_sols          = r.norm_sols
    locking_labels     = r.locking_labels
    bifurcation_bounds = r.bifurcation_bounds
    Control1           = op.Control1
    Control2           = op.Control2
    BSON.@save path ode_sols norm_sols locking_labels bifurcation_bounds Control1 Control2
    @info "Saved ODE results → $path"
    return path
end


"""
    load_ode_results(actor; filename="ode_results.bson", dir=LOCKING_RESULTS_DIR)

Load previously saved ODE results from disk into `actor.results` and
`actor.ode_params.Control1/Control2`.  Called by `task=:calc_prob` when
`actor.results` is nothing (fresh session).
"""
function load_ode_results!(actor::ActorLocking;
                            filename::String = "ode_results.bson",
                            dir::String      = LOCKING_RESULTS_DIR)
    path = joinpath(dir, filename)
    isfile(path) || error("No saved ODE results at $path — run task=:solve_system first")
    d = BSON.load(path, @__MODULE__)
    actor.ode_params === nothing && (actor.ode_params = ODEparams())
    actor.ode_params.Control1 = d[:Control1]
    actor.ode_params.Control2 = d[:Control2]
    actor.results = LockingResults(
        d[:ode_sols],
        nothing,
        d[:norm_sols],
        d[:locking_labels],
        d[:bifurcation_bounds],
    )
    @info "Loaded ODE results ← $path"
    return actor.results
end


"""
    save_locking_nn(actor; filename="nn_model.bson", dir=LOCKING_RESULTS_DIR) → path

Save the trained LockingNNModel to disk.
Default location: ~/.julia/locking_results/
"""
function save_locking_nn(actor::ActorLocking;
                          filename::String = "nn_model.bson",
                          dir::String      = LOCKING_RESULTS_DIR)
    actor.results === nothing        && error("No results — run task=:solve_system first")
    actor.results.prob === nothing   && error("No trained model — run train_locking_nn first")
    mkpath(dir)
    path      = joinpath(dir, filename)
    model     = actor.results.prob.model
    nn_params = actor.results.prob.nn_params
    BSON.@save path model nn_params
    @info "Saved locking NN model → $path"
    return path
end


"""
    load_locking_nn(; filename="nn_model.bson", dir=LOCKING_RESULTS_DIR) → LockingNNModel

Load a saved LockingNNModel from disk. Cached in memory after first load.
"""
function load_locking_nn(; filename::String = "nn_model.bson",
                           dir::String      = LOCKING_RESULTS_DIR)
    path = joinpath(dir, filename)
    haskey(_locking_nn_cache, path) && return _locking_nn_cache[path]
    isfile(path) || error("No saved model at $path — run save_locking_nn first")
    d = BSON.load(path, @__MODULE__)
    prob_model = LockingNNModel(d[:model], d[:nn_params])
    _locking_nn_cache[path] = prob_model
    @info "Loaded locking NN model ← $path"
    return prob_model
end


# ─────────────────────────────────────────────────────────────────────────────
#  Plotting
# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
#  Convenience wrapper
# ─────────────────────────────────────────────────────────────────────────────
"""
    plot_sols(actor) → (fig1, fig2, fig3)

Calls plot_scatter, plot_phase_diagrams, and (when a trained NN model is
available) plot_probability.  Returns all three handles; fig3 is nothing
when no model has been trained yet.
"""
function plot_sols(actor)
    fig1 = plot_scatter(actor);      display(fig1)
    fig2 = plot_phase_diagrams(actor); display(fig2)
    fig3 = (actor.results !== nothing && actor.results.prob !== nothing) ?
           (p = plot_probability(actor); display(p); p) : nothing
    return fig1, fig2, fig3
end




# ─────────────────────────────────────────────────────────────────────────────
#  Figure 1 — multi-panel scatter of state variables
# ─────────────────────────────────────────────────────────────────────────────
"""
    plot_scatter(actor) → Figure 1

Multi-panel scatter of raw and normalised ODE state variables.

RP-RW (5-state) — 2×2 layout:
  (a) raw ψ_t vs Ω_t
  (b) ψ_tn vs Ω_n, coloured by k-means class (unlocked = blue, locked = red)
  (c) ψ_tn vs ψ_wn
  (d) θ_t vs θ_w

RP-IW (3-state) — 2×1 layout:
  (a) raw ψ_t vs Ω_t
  (b) ψ_tn vs Ω_n, coloured by k-means class
"""
function plot_scatter(actor)
    r = actor.results
    r === nothing && error("No results — run the actor first")

    is_RPRW = size(r.norm_sols, 2) == 5

    psi_t   = r.ode_sols[:, 1]
    omega_t = r.ode_sols[:, 3]
    psi_tn  = r.norm_sols[:, 1]
    omega_n = r.norm_sols[:, 3]

    idx_U = findall(r.locking_labels .== 1)   # unlocked
    idx_L = findall(r.locking_labels .== 2)   # locked

    # ── (a) raw ψ_t vs Ω_t ──────────────────────────────────────────────────
    p_a = scatter(psi_t, omega_t;
        xlabel          = "ψ_t",
        ylabel          = "Ω_t",
        label           = false,
        alpha           = 0.3,
        markersize      = 3,
        markerstrokewidth = 0,
        grid            = true,
    )
    xra = extrema(psi_t); yra = extrema(omega_t)
    annotate!(p_a, xra[1] + 0.02*(xra[2]-xra[1]),
                   yra[1] + 0.04*(yra[2]-yra[1]),
                   Plots.text("(a)", 12, :left))

    # ── (b) ψ_tn vs Ω_n coloured by class ───────────────────────────────────
    p_b = scatter(psi_tn[idx_U], omega_n[idx_U];
        xlabel          = "ψ_tn",
        ylabel          = "Ω_n",
        label           = "Unlocked",
        color           = :steelblue,
        alpha           = 0.3,
        markershape     = :circle,
        markersize      = 4,
        markerstrokewidth = 0,
        grid            = true,
    )
    scatter!(p_b, psi_tn[idx_L], omega_n[idx_L];
        label       = "Locked",
        color       = :red,
        alpha       = 0.5,
        markershape = :xcross,
        markersize  = 5,
    )
    xlims!(p_b, -0.02, 1.06)
    ylims!(p_b, -0.02, 1.06)
    annotate!(p_b, 0.02, 0.02, Plots.text("(b)", 12, :left))

    if is_RPRW
        psi_w   = r.ode_sols[:, 4]
        psi_wn  = r.norm_sols[:, 4]
        theta_t = r.norm_sols[:, 2]
        theta_w = r.norm_sols[:, 5]

        # ── (c) raw ψ_w vs Ω_t ──────────────────────────────────────────────
        p_c = scatter(psi_w, omega_t;
            xlabel          = "ψ_w",
            ylabel          = "Ω_t",
            label           = false,
            color           = :steelblue,
            alpha           = 0.3,
            markersize      = 3,
            markerstrokewidth = 0,
            grid            = true,
        )
        xrc = extrema(psi_w); yrc = extrema(omega_t)
        annotate!(p_c, xrc[1] + 0.02*(xrc[2]-xrc[1]),
                       yrc[1] + 0.04*(yrc[2]-yrc[1]),
                       Plots.text("(c)", 12, :left))

        # ── (d) ψ_wn vs Ω_n coloured by class ──────────────────────────────
        p_d = scatter(psi_wn[idx_U], omega_n[idx_U];
            xlabel          = "ψ_wn",
            ylabel          = "Ω_n",
            label           = "Unlocked",
            color           = :steelblue,
            alpha           = 0.3,
            markershape     = :circle,
            markersize      = 4,
            markerstrokewidth = 0,
            grid            = true,
        )
        scatter!(p_d, psi_wn[idx_L], omega_n[idx_L];
            label       = "Locked",
            color       = :red,
            alpha       = 0.5,
            markershape = :xcross,
            markersize  = 5,
        )
        xlims!(p_d, -0.02, 1.06)
        ylims!(p_d, -0.02, 1.06)
        annotate!(p_d, 0.02, 0.02, Plots.text("(d)", 12, :left))

        # ── (e) ψ_tn vs ψ_wn ────────────────────────────────────────────────
        p_e = scatter(psi_tn, psi_wn;
            xlabel          = "ψ_tn",
            ylabel          = "ψ_wn",
            label           = false,
            color           = :steelblue,
            alpha           = 0.3,
            markersize      = 3,
            markerstrokewidth = 0,
            grid            = true,
        )
        xre = extrema(psi_tn); yre = extrema(psi_wn)
        annotate!(p_e, xre[1] + 0.02*(xre[2]-xre[1]),
                       yre[1] + 0.04*(yre[2]-yre[1]),
                       Plots.text("(e)", 12, :left))

        # ── (f) θ_t vs θ_w ──────────────────────────────────────────────────
        p_f = scatter(theta_t, theta_w;
            xlabel          = "θ_t (rad)",
            ylabel          = "θ_w (rad)",
            label           = false,
            color           = :steelblue,
            alpha           = 0.3,
            markersize      = 3,
            markerstrokewidth = 0,
            grid            = true,
        )
        xrf = extrema(theta_t); yrf = extrema(theta_w)
        annotate!(p_f, xrf[1] + 0.02*(xrf[2]-xrf[1]),
                       yrf[1] + 0.04*(yrf[2]-yrf[1]),
                       Plots.text("(f)", 12, :left))

        plt = plot(p_a, p_b, p_c, p_d, p_e, p_f; layout=(3, 2), size=(900, 1050))
    else
        plt = plot(p_a, p_b; layout=(2, 1), size=(600, 750))
    end

    return plt
end


# ─────────────────────────────────────────────────────────────────────────────
#  Figure 2 — phase diagrams (pcolor of Ω_n and ψ_tn over control space)
# ─────────────────────────────────────────────────────────────────────────────
"""
    plot_phase_diagrams(actor) → Figure 2

Two stacked pcolor panels of normalised solutions over the (C2, C1) control space.
  (a) Ω_n — normalised rotation      (RdBu colormap)
  (b) ψ_tn — normalised TM amplitude (RdBu_r colormap)
Analytic bifurcation boundary (D = 0) overlaid in black when NL saturation is off.
"""
function plot_phase_diagrams(actor)
    r   = actor.results
    par = actor.par
    r === nothing && error("No results — run the actor first")

    gs = par.grid_size
    x  = unique(actor.ode_params.Control2)   # Control2 values  (x-axis)
    y  = unique(actor.ode_params.Control1)   # Control1/Ω0 values (y-axis)

    OmN_grid   = reshape(r.norm_sols[:, 3], gs, gs)   # Ω_n   (rows=C1, cols=C2)
    PsiTN_grid = reshape(r.norm_sols[:, 1], gs, gs)   # ψ_tn  (rows=C1, cols=C2)

    xlabel_ctrl = _ctrl_xlabel(par)
    xr = extrema(x); yr = extrema(y)

    # ── (a) Ω_n ─────────────────────────────────────────────────────────────
    p_a = heatmap(x, y, OmN_grid;
        xlabel         = xlabel_ctrl,
        ylabel         = "Ω_0",
        title          = "Ω_n",
        color          = cgrad(:RdBu),
        colorbar_title = "Ω_n",
        clims          = (0.0, 1.0),
        left_margin    = 8Plots.mm,
    )
    _overlay_bifurcation!(p_a, x, y, r.bifurcation_bounds)
    annotate!(p_a, xr[1]+0.05*(xr[2]-xr[1]), yr[1]+0.85*(yr[2]-yr[1]),
              Plots.text("UNLOCKED", 14, :white, :left))
    annotate!(p_a, xr[1]+0.65*(xr[2]-xr[1]), yr[1]+0.05*(yr[2]-yr[1]),
              Plots.text("LOCKED",   14, :white, :left))
    annotate!(p_a, xr[1]+0.01*(xr[2]-xr[1]), yr[1]+0.05*(yr[2]-yr[1]),
              Plots.text("(a)", 12, :white, :left))

    # ── (b) ψ_tn ────────────────────────────────────────────────────────────
    p_b = heatmap(x, y, PsiTN_grid;
        xlabel         = xlabel_ctrl,
        ylabel         = "Ω_0",
        title          = "ψ_tn",
        color          = cgrad(:RdBu, rev=true),
        colorbar_title = "ψ_tn",
        clims          = (0.0, 1.0),
        left_margin    = 8Plots.mm,
    )
    _overlay_bifurcation!(p_b, x, y, r.bifurcation_bounds)
    annotate!(p_b, xr[1]+0.05*(xr[2]-xr[1]), yr[1]+0.85*(yr[2]-yr[1]),
              Plots.text("UNLOCKED", 14, :white, :left))
    annotate!(p_b, xr[1]+0.65*(xr[2]-xr[1]), yr[1]+0.05*(yr[2]-yr[1]),
              Plots.text("LOCKED",   14, :white, :left))
    annotate!(p_b, xr[1]+0.01*(xr[2]-xr[1]), yr[1]+0.05*(yr[2]-yr[1]),
              Plots.text("(b)", 12, :white, :left))

    plt = plot(p_a, p_b; layout=(2, 1), size=(650, 1100))
    return plt
end


# ─────────────────────────────────────────────────────────────────────────────
#  Figure 3 — NN locking probability
# ─────────────────────────────────────────────────────────────────────────────
"""
    plot_probability(actor) → Figure 3

Contourf of NN locking probability P(locked) over the (C2, C1) control space.
  • dashed black  : P = 0.5 decision boundary
  • dashed yellow : analytic bifurcation boundary (D = 0), when available
"""
function plot_probability(actor)
    r   = actor.results
    par = actor.par
    r === nothing      && error("No results — run the actor first")
    r.prob === nothing && error("No trained NN model — run train_locking_nn first")

    x = unique(actor.ode_params.Control2)   # Control2 (x-axis)
    y = unique(actor.ode_params.Control1)   # Control1/Ω0 (y-axis)

    prob_grid   = [r.prob(c1, c2) for c1 in y, c2 in x]
    xlabel_ctrl = _ctrl_xlabel(par)
    xr = extrema(x); yr = extrema(y)

    plt = contourf(x, y, prob_grid;
        xlabel         = xlabel_ctrl,
        ylabel         = "Ω_0",
        title          = "Locking probability P(locked) — NN",
        colorbar_title = "P(locked)",
        clims          = (0.0, 1.0),
        levels         = 20,
        color          = cgrad(:RdBu, rev=true),
    )
    contour!(plt, x, y, prob_grid;
        levels    = [0.5],
        linecolor = :black,
        linestyle = :dash,
        linewidth = 2,
        colorbar  = false,
        label     = "P = 0.5",
    )
    _overlay_bifurcation!(plt, x, y, r.bifurcation_bounds; color=:yellow, style=:dash)

    annotate!(plt, xr[1]+0.05*(xr[2]-xr[1]), yr[1]+0.85*(yr[2]-yr[1]),
              Plots.text("UNLOCKED", 14, :white, :left))
    annotate!(plt, xr[1]+0.70*(xr[2]-xr[1]), yr[1]+0.05*(yr[2]-yr[1]),
              Plots.text("LOCKED",   14, :white, :left))

    return plt
end


# ─────────────────────────────────────────────────────────────────────────────
#  Plotting helpers
# ─────────────────────────────────────────────────────────────────────────────

"Return the x-axis label for the swept control parameter C2"
function _ctrl_xlabel(par)
    par.control_type == :EF           ? "Error Field"         :
    par.control_type == :LinStab      ? "Linear Stability Δ′" :
    par.control_type == :NLsaturation ? "NL Saturation α"     : "Control 2"
end

"""
    _zero_isoline(x, y, z) → (xs, ys)

Compute the D=0 isoline of a 2-D scalar field `z` defined on grid `(x, y)`
by linear interpolation along cell edges (simplified marching squares).
Returns flat vectors suitable for `plot!`; NaN separates disjoint segments.
`z` must be (length(y) × length(x)) — the same convention as Plots.jl heatmap.
"""
function _zero_isoline(x::AbstractVector, y::AbstractVector, z::AbstractMatrix)
    xs = Float64[]
    ys = Float64[]
    nx, ny = length(x), length(y)   # x → columns, y → rows

    for i in 1:ny-1, j in 1:nx-1
        pts = NTuple{2,Float64}[]

        # bottom edge: row i,   col j → j+1
        v1, v2 = z[i,j], z[i,j+1]
        if v1 * v2 < 0
            t = v1 / (v1 - v2)
            push!(pts, (x[j] + t*(x[j+1]-x[j]), y[i]))
        end
        # top edge:    row i+1, col j → j+1
        v1, v2 = z[i+1,j], z[i+1,j+1]
        if v1 * v2 < 0
            t = v1 / (v1 - v2)
            push!(pts, (x[j] + t*(x[j+1]-x[j]), y[i+1]))
        end
        # left edge:   col j,   row i → i+1
        v1, v2 = z[i,j], z[i+1,j]
        if v1 * v2 < 0
            t = v1 / (v1 - v2)
            push!(pts, (x[j], y[i] + t*(y[i+1]-y[i])))
        end
        # right edge:  col j+1, row i → i+1
        v1, v2 = z[i,j+1], z[i+1,j+1]
        if v1 * v2 < 0
            t = v1 / (v1 - v2)
            push!(pts, (x[j+1], y[i] + t*(y[i+1]-y[i])))
        end

        if length(pts) >= 2
            push!(xs, pts[1][1], pts[2][1], NaN)
            push!(ys, pts[1][2], pts[2][2], NaN)
        end
    end
    return xs, ys
end

"""
Overlay analytic bifurcation boundary (D=0 isoline) — no-op when bb is nothing.

Drawn as a plain `plot!` line series from manually-computed crossing points,
rather than a `contour!` series: heatmap + contour! on one subplot share a
single color/z-scale in GR, and `bb`'s native range (e.g. up to ~700) versus
the heatmap's `[0,1]` range makes that shared scale unworkable either way
(heatmap goes flat, or the D=0 level gets clipped away). A line series has no
z/colormap at all, so it cannot disturb the heatmap's color scale.
"""
function _overlay_bifurcation!(p, x, y, bb; color=:black, style=:solid)
    bb === nothing && return
    xs, ys = _zero_isoline(x, y, bb)
    isempty(xs) && return
    plot!(p, xs, ys; linecolor=color, linestyle=style, linewidth=2.5, label=false)
end


# ─────────────────────────────────────────────────────────────────────────────
#  Legacy plotting helpers
# ─────────────────────────────────────────────────────────────────────────────

function make_contour(X::AbstractArray, Y::AbstractArray, Z::AbstractMatrix)
    # Determine target shape
    m, n = size(Z)

    # If C1 and C2 are vectors, try to reshape them
    if ndims(X) == 1
        X = unique(X)
    end
    if ndims(Y) == 1
        Y = unique(Y)
    end
    
    plt = plot()
    plt = heatmap!(plt, X,Y, Z)#; linewidth=2)
    display(plt)   # explicitly display
    return plt
end

function make_contour(X::AbstractArray, Y::AbstractArray, Z::AbstractMatrix, levels::Vector{Float64}, control_type)
    lblsz = 16
    
    # Ensure unique 1D grids
    if ndims(X) == 1 
        X = unique(X)
    end
    if ndims(Y) == 1 
        Y = unique(Y)
    end

    # Axis label
    xlabel = control_type == :EF          ? "Error Field" :
             control_type == :LinStab     ? "Linear Stability" :
             control_type == :NLsaturation ? "NL saturation" : "Control1"

    # Contour plot
    plt = contour(X, Y, Z; levels=levels, linewidth=2, 
                  xlabel=xlabel, ylabel="Normalized Torque", 
                  clabel=false)

    # Compute axis ranges for relative placement
    xmin, xmax = extrema(X)
    ymin, ymax = extrema(Y)

    annotate!(xmin + 0.5*(xmax-xmin), ymin + 0.85*(ymax-ymin), text("UNLOCKED", lblsz))
    annotate!(xmin + 0.80*(xmax-xmin), ymin + 0.02*(ymax-ymin), text("LOCKED", lblsz))

    if any(Z .< 0)
        ind = argmin(Z)
        i, j = Tuple(ind)   # row, col indices
        xloc = 0.9*X[j]
        yloc = 0.9*Y[i]
        annotate!(xloc, yloc, text("Locking\n(Possible)", lblsz, :red))
    end

    display(plt)
    return plt
end

