using Distributed
using DifferentialEquations
import Roots
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
        [:solve_system, :single_case, :Monte_Carlo, :evaluate_probability, :transfer_learning],
        "-",  
        "Choose whether to simulate system on full-grid, do a single case, or Monte Carlo (NOT implmented) which does 1e4 simulations on the full grid",  # docstring
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
    MaxErrorField::Entry{Float64} = Entry{Float64}(
        "Tesla*meter", 
        "Maximum error field perturbation for control, usually ~10 Gauss.meter"; default=1e-3)
    b0::Entry{Float64} = Entry{Float64}(
        "Tesla", 
        "Scale for magnetic perturbations, usually ~10Gauss"; default=1.e-3)
    t0::Entry{Float64} = Entry{Float64}(
        "seconds", 
        "Characteristic time scale for normalization , usually TM/RW growth rate"; default=1.e-3)
    r0::Entry{Float64} = Entry{Float64}(
        "meter", 
        "Length scale for the integration , usually minor radius"; default=1.)
    NBItorque::Entry{Float64} = Entry{Float64}("Newton.meter", "NBI torque"; default=5.)    
end


# Make sure workers know the ODEparams layout so they can receive an instance.
# (Same field order/types as your definition.)
#if !isdefined(Main, :ODEparams)
Base.@kwdef mutable struct ODEparams
    #sim_time::Vector{Float64} = Float64[]# simulation time vector
    saturation_param::Float64 = 0.2 # controls nonlinear island saturation
    error_field::Float64 = 0.5 # Set this when control_type is NOT :EF
    DeltaW::Float64 = 0.0 # intrinsic stability of RW mode
    stability_index::Float64 = 0.0 # intrinsic stability of the mode: Delta'
    layer_width::Float64 = 1.e-3 # linear layer width from tearing theory, ususally ~1 mm
    rat_surface::Float64 = 0.67 # q=2 surface location in dimensionless units"; default=0.67)
    res_wall::Float64 = 1.0  # Resistive wall location in dimensionless units"; default=1.0)
    control_surf::Float64 = 1.25    # Control surface location in dimensionless units"; default=1.25)
    mu::Float64 = 0.1              # anomalous perp. plasma viscosity
    Inertia::Float64 = 0.1          # moment of inertia of the layer
    Control1::Vector{Float64} = Float64[] # control parameter 1, e.g. error field
    Control2::Vector{Float64} = Float64[] # control parameter 2,
    l12::Float64 = 1.0        # mutual inductance between rational surface and RW
    l21::Float64 = 1.0        # mutual inductance between RW and rational surface
    l32::Float64 = 1.0        # mutual inductance between the control surface and RW
    Taut_Tauw::Float64 = 1.0  # ratio of tearing time to wall time
    hyper_cube_dims::Vector{Float64} = [1., 1., 1.] # dimensions of the initial condtion hypercube
    Control2_min::Float64 = 0.01
    Control2_max::Float64 = 1.0
    # keep adding default and magic things
end

Base.@kwdef struct ContourData
    x::Vector{Float64}
    y::Vector{Float64}
    z::Matrix{Float64}
end

mutable struct LockingResults
    ode_sols::Matrix{Float64}            # (N*M × n_states) raw final states
    prob::Any                            # NN model once Task 1 is done; callable as prob(C1, C2)
    norm_sols::Matrix{Float64}           # (N*M × n_states) normalized solutions
    locking_labels::Vector{Int}          # k-means class assignments, one per grid point
    contour_data::ContourData
end

mutable struct ActorLocking{D,P} <: SingleAbstractActor{D,P}
    dd::IMAS.dd{D}
    par::OverrideParameters{P,FUSEparameters__ActorLocking{P}}
    ode_params::Union{Nothing, ODEparams}
    results::Union{Nothing, LockingResults}

    function ActorLocking(
        dd::IMAS.dd{D},
        par::FUSEparameters__ActorLocking{P};
        ode_params = nothing,
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

        return new{D,P}(dd, par, ode, nothing)
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
    actor.ode_params = ODEparams(;)
    #pressure = dd.equilibrium.time_slice[].pressure

    actor.ode_params = set_up_ode_params!(dd, par, actor.ode_params)
    
    #bifurcation_bounds = calculate_bifurcation_bounds(dd, par, actor.ode_params)

    ## Time evolve the ODEs
    if task == :single_case
        control1 = 1. # normalized rotation frequency.
        solve_one_case(par, actor.ode_params, application, control1)
        return actor
    elseif task == :solve_system
        # Solve the ODE system on the whole control grid):
        control1 = actor.ode_params.Control1
        control2 = actor.ode_params.Control2
        ode_sols = solve_system(actor, application)
        norm_sols = normalize_ode_results(ode_sols, actor.ode_params,
             control2, control1, par.control_type)

        # ## plot normalize solution scatters
        #plot_sols_scatter(norm_sols; xcol=1,ycol=3)
        #inputs = [(c1, c2) for (c1, c2) in zip(control1, control2)]
        #inputs = vec(inputs)  # flatten

        # If you want them back on a 2D grid (N x M), reshape here:
        N, M = size(norm_sols)
        norm_sols_2D= reshape(norm_sols, (par.grid_size, par.grid_size, M))
        psi_tN = norm_sols_2D[:,:,1]
        

        ## classify normalized solutions
        R = hcat(norm_sols[:,1], norm_sols[:,3])
        kmc = kmeans(R', 2)
        locking_labels = kmc.assignments 
        #plot_sols_scatter(norm_sols; labels=locking_labels)

        contour = ContourData(
            x = control2,
            y = control1,
            z = psi_tN
)

        # store back in the actor
        prob = nothing # placeholder until Task 1 (NN classifier) is implemented
        actor.results = LockingResults(
            ode_sols,
            prob,
            norm_sols,
            locking_labels,
            contour
        )

    elseif task == :evaluate_probability
        @info   "Evaluating saved Locking probability (not implemented yet)"
        # C1 = .5; C2 = .5 # place holders for now, need to implement a way to specify which case(s) to evaluate probability for
        #probability = actor.results.prob(C1, C2) # placeholder for now, need to implement a way to evaluate probability from the ODE solutions
    elseif task == :Monte_Carlo
        error("Monte-Carlo not implemented yet")
    elseif task == :transfer_learning
        @info   "Transfer Learning training (not implemented yet)"    

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


function plot_sols(actor)
    r   = actor.results
    par = actor.par

    xlabel_ctrl = par.control_type == :EF          ? "Error Field (a.u.)"  :
                  par.control_type == :LinStab      ? "Linear Stability Δ′" :
                  par.control_type == :NLsaturation ? "NL Saturation α"     : "Control 2"

    # Unique 1-D axes (contour_data stores the full N*M vectors)
    x = unique(r.contour_data.x)   # Control2 axis  (length M)
    y = unique(r.contour_data.y)   # Control1 axis  (length N)
    z = r.contour_data.z           # N × M matrix of ψ at t_final

    # Locking boundary: reshape 1-D label vector back to (N, M) grid
    labels_2D = reshape(r.locking_labels, par.grid_size, par.grid_size)

    p1 = heatmap(x, y, z;
        xlabel         = xlabel_ctrl,
        ylabel         = "Rotation Frequency (a.u.)",
        #title          = "ψ at t_final",
        colorbar_title = "ψ_N",
    )
    # Overlay the locking boundary as a white contour at the class transition
    contour!(p1, x, y, float.(labels_2D);
        levels    = [1.5],
        linecolor = :white,
        linewidth = 2,
        colorbar  = false,
        label     = "locking boundary",
    )

    # Scatter: ψ_N vs Ω_N coloured by locking class
    xs      = r.norm_sols[:, 1]
    ys      = r.norm_sols[:, 3]
    classes = sort(unique(r.locking_labels))
    shapes  = [:circle, :x]

    p2 = plot(;
        xlabel = "NormalizedTM amplitude ψ_N",
        ylabel = "Normalized Rotation Ω_N",
        title  = "Locking classification",
    )
    for (cl, sh) in zip(classes, shapes)
        idx = findall(r.locking_labels .== cl)
        scatter!(p2, xs[idx], ys[idx]; markershape=sh, label="Class $cl", markersize=5)
    end

    plt = plot(p1, p2; layout=(1, 2), size=(1100, 450))
    display(plt)
    return plt
end



function solve_one_case(par, ode_params::ODEparams, task::String, control1::Float64)    
    
    if par.control_type == :EF
        control2 = 0.5
    elseif par.control_type == :LinStab
        control2 = -6.
    else
        control2 = 0.2 
    end

    final_sol = solve_ODEs(par, ode_params, task, control1, control2)
    println("final raw solution = ", final_sol)

    sols_norm = normalize_ode_results(final_sol, ode_params, control2, control1, par.control_type)
    println("final normalized solution = ", sols_norm)

    return
end


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

# function set_sim_time(time0::Float64, tfinal::Float64, time_steps::Int64) 
#     """
#     initialize time for the ODE integration
#     """
#     sim_time = collect(range(time0, tfinal, length=time_steps))
    
#     return sim_time
# end


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
    
    # get the actual torque from dd
    #cp1d = dd.core_profiles.profiles_1d[]
    #total_source1d = IMAS.total_sources(dd.core_sources, cp1d; time0=dd.global_time, fields=[:torque_tor_inside])
    #dd.core_sources.source[15].global_quantities[1].torque_tor

    # Calculate the drag coefficient in SI
    muSI = par.NBItorque / rot_core
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

    # figure out the rotation rate at the q=2 surface
    rho = dd.core_sources.source[1].profiles_1d[1].grid.rho_tor_norm
    rot_core = dd.core_profiles.profiles_1d[1].rotation_frequency_tor_sonic
    rot_interp = IMAS.interp1d(rho, rot_core)
    Omega0 = rot_interp(rt) * 2 * π * par.t0 # in units of kHz
    @info("Calculated toroidal rotation frequency at rational surface: $Omega0 kHz")
    Omega0 = 10  # overwrite for now to DEBUG
    println("Toroidal rotation requency at rational surface: ", Omega0, " kHz")

    Om0Vals = range(1.0e-2, Omega0, length=N) |> collect
    Om0s = repeat(Om0Vals, 1, M)
    Control1 = vec(Om0s)
    ode_params.Control1 = Control1

    # Initialize the other control parameter based on the control type
    Control2_vals = range(c2min, c2max, length=M) |> collect
    if control_type == :EF
        println("Overwriting default Control2 range with Max Error Field")
        EpsUp = par.MaxErrorField / (par.b0 * par.r0)  # Convert to dimensionless units
        Control2_vals = range(1e-2, EpsUp, length=M) |> collect

    elseif control_type == :LinStab
        #ode_params.error_field = 0.6  # Example value for error field
        #Control2_vals = range(ode_params.Delta_lower, ode_params.Delta_upper, length=M) |> collect
        ## Check to make sure the system is still weakly stable
        DeltatRW = Control2_vals .- l21*l12/DeltaW
        if any(DeltatRW .> 0)
            println("*** ALERT: You set up a case with an unstable RP-RW mode! ***")
            println("*** Maximum RP-RW stability set to ", maximum(DeltatRW))
            error("Deltat_RW > 0 for your range of TM stability values ***")
        end

    elseif control_type == :NLsaturation
        println("Do nothing for NL saturation control")
        #Control2_vals = range(ode_params.alpha_lower, ode_params.alpha_upper, length=M) |> collect
        
        #ode_params.saturation_param = Control2
    end
    
    Control2 = vec(repeat(Control2_vals', N, 1))
    ode_params.Control2 = Control2

    return ode_params
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
function make_rhs_function(phys_type::String)
    phys_type == "RP-RW" ? rhs_RW! :
    phys_type == "RP-IW"   ? rhs_basic! :
    error("Unknown application type: $phys_type")
   
end

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
        q = @. (DeltatRW / m0)^2 + rt * (l32 * l21 * eps / DeltaW)^2 / (m0 * mu)
        r = @. -Y * DeltatRW^2 / m0^2
    elseif par.application == "RP-IW"
        # 3rd-order system: rational surface + ideal wall only
        q = @. (Deltat / m0)^2 + (l21 * eps)^2 / mu
        r = @. -Y * Deltat^2 / m0^2
    else
        error("calculate_bifurcation_bounds not implemented for application: $(par.application)")
    end

    a = @. -(Y^2) / 3.0 + q
    b = @. 2.0 * (-Y)^3 / 27.0 - q * (-Y) / 3.0 + r

    bifurcation_bounds = reshape(@. b^2 / 4 + a^3 / 27, par.grid_size, par.grid_size)

    make_contour(ode_params.Control2, ode_params.Control1, bifurcation_bounds, [0.0, 1e-10], par.control_type)

    return bifurcation_bounds
end


function solve_ODEs(par, ode_params::ODEparams, task::String, C1::Float64, C2::Float64)
    n_mode       = par.n_mode
    control_type = par.control_type
    t_final      = par.t_final

    rhs!     = make_rhs_function(task)
    ode_rhs! = make_ode_func(rhs!)
    ode_params.hyper_cube_dims = [3., 3., 3.]
    y0 = make_initial_condition(ode_params.hyper_cube_dims, task)

    # rhs! expects (C2, C1) as its first two positional args — swap once here so
    # all callers use the natural (C1, C2) convention
    p = (C2, C1, ode_params, n_mode, control_type)
    tspan = (0.0, t_final)

    prob = ODEProblem(ode_rhs!, y0, tspan, p)
    sol  = solve(prob, Tsit5(); saveat=t_final, reltol=1e-8, abstol=1e-10)

    return sol.u[end]
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

            if iszero(alpha)  # linear regime: saturation_param set to 0. when NL_saturation_ON=false
                num     = abs(Deltat * DeltaW) - l12 * l21
                psitMax = l32 * l21 * eps / num
                psiwMax = l32 * abs(DeltaW) * eps / num
            else              # NL saturation active
                psitMax = -(DeltatRW + sqrt(DeltatRW^2 + 4*alpha*l21*l32*eps*Deltat/DeltaW)) /
                           (2*alpha*Deltat)
                psiwMax = -(l32*eps + l12*psitMax) / DeltaW
            end

            psiN  = abs(psit / psitMax)
            psiwN = abs(psiw / psiwMax)

            return [psiN, theta_t, OmN, psiwN, theta_w]

        elseif length(final_sol) == 3 # RP-IW layout
            psiN    = final_sol[1] / (l21 * eps)
            theta_t = mod(final_sol[2], 2π)
            OmN     = final_sol[3] / C1
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
plot_sols_scatter(norm_sol; xcol=1, ycol=3)

Make a scatter plot from normalized solutions.

Arguments:
- `norm_sol`: 
    * Vector{Vector{Float64}} or 
    * Vector{NTuple{N,Float64}} (N=3 or 5 typically).
  Each entry is one normalized solution.
- `xcol`, `ycol`: column indices (1-based) to plot.

Example:
    plot_normalized_scatter(norm_sol; xcol=1, ycol=3)
"""
function plot_sols_scatter(
        norm_sol; 
        xcol::Int=1, 
        ycol::Int=3,
        xlabel::AbstractString="TM amplitude",
        ylabel::AbstractString="Rotation at the rat. surf.",
        title::AbstractString="Normalized solution scatter"
    )
    
    if eltype(norm_sol) <: AbstractVector{<:Real}
        # Convert vector-of-vectors to a matrix
        data = reduce(vcat, (x' for x in norm_sol))
    elseif eltype(norm_sol) <: NTuple
        # Convert vector-of-tuples to a matrix
        data = reduce(vcat, (collect(x)' for x in norm_sol))
    else
        throw(ArgumentError("norm_sol must be Vector{Vector{Float64}} or Vector{NTuple{N,Float64}}"))
    end

    plt = plot()
    plt = scatter!(plt, data[:, xcol], data[:, ycol],
                  xlabel=xlabel,
                  ylabel=ylabel,
                  title=title)

    display(plt)
    return plt
end


function plot_sols_scatter(
        norm_sols::AbstractMatrix{<:Real};
        xcol::Int = 1,
        ycol::Int = 3,
        labels::Union{Nothing,AbstractVector{<:Integer}} = nothing,
        xlabel::AbstractString = "TM amplitude",
        ylabel::AbstractString = "Rotation at the rat. surf.",
        title::AbstractString = "Normalized solution scatter"
    )

    N, M = size(norm_sols)

    @assert 1 ≤ xcol ≤ M "xcol out of bounds"
    @assert 1 ≤ ycol ≤ M "ycol out of bounds"

    x = norm_sols[:, xcol]
    y = norm_sols[:, ycol]

    plt = plot()

    if labels === nothing
        scatter!(
            plt,
            x, y;
            xlabel = xlabel,
            ylabel = ylabel,
            title  = title,
            legend = false,
            markersize = 4
        )
    else
        @assert length(labels) == N "labels must match number of rows"
        # Ensure two unique classes
        classes = unique(labels)
        @assert length(classes) == 2 "This version expects exactly 2 classes"

        class1, class2 = classes

        idx1 = findall(labels .== class1)
        idx2 = findall(labels .== class2)

        scatter!(
            plt,
            x[idx1], y[idx1];
            markershape = :circle,
            label = "Class $class1",
            markersize = 5
        )

        scatter!(
            plt,
            x[idx2], y[idx2];
            markershape = :x,
            label = "Class $class2",
            markersize = 6
        )
        # scatter!(
        #     plt,
        #     x, y;
        #     group = labels,
        #     xlabel = xlabel,
        #     ylabel = ylabel,
        #     title  = title,
        #     legend = :topright,
        #     markersize = 4
        # )
    end

    display(plt)
    return plt
end



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

