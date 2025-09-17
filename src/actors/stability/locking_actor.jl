#import DifferentialEquations as DiffEqs
using Distributed
using DifferentialEquations
import Roots

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
    control_type::Switch{Symbol} = Switch{Symbol}([:EF, :LinStab, :NLsaturation], # EF: error field
        "-",                                                            # LinStab: vary stability_index,
        "Use a user specified Control case to run the locking models"; default=:EF) # NLsaturation: vary NL saturation
    NL_saturation_ON::Entry{Bool} = Entry{Bool}("-", "Nonlinear saturation parameter for the mode"; default=true)
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
        "Time scale for the integration , usually TM/RW growth rate"; default=1.e-3)
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
    rat_index::Int64 = 50        # index of the q=2 surface in the rho_tor_norm profile
    mu::Float64 = 0.1              # shear modulus
    Inertia::Float64 = 0.1          # moment of inertia of the layer
    Control1::Vector{Float64} = Float64[] # control parameter 1, e.g. error field
    Control2::Vector{Float64} = Float64[] # control parameter 2,
    l12::Float64 = 1.0        # mutual inductance between rational surface and RW
    l21::Float64 = 1.0        # mutual inductance between RW and rational surface
    l32::Float64 = 1.0        # mutual inductance between the control surface and RW
    Taut_Tauw::Float64 = 1.0  # ratio of tearing time to wall time
    hyper_cube_dims::Vector{Float64} = [1., 1., 1.] # dimensions of the initial condtion hypercube
    #alpha_upper::Float64 = 0.6
    #alpha_lower::Float64 = 0.1
    # keep adding default and magic things
end

mutable struct ActorLocking{D,P} <: SingleAbstractActor{D,P}
    dd::IMAS.dd{D}
    par::OverrideParameters{P,FUSEparameters__ActorLocking{P}}
    ode_params::Union{Nothing, ODEparams}
    # Add the ODEstruct thing here so you can access it when you debug it
    function ActorLocking(dd::IMAS.dd{D}, par::FUSEparameters__ActorLocking{P}; kw...) where {D<:Real,P<:Real}
        logging_actor_init(ActorLocking)
        par = OverrideParameters(par; kw...)
        return new{D,P}(dd, par, nothing)
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

    # actor.ode_params :: nothing
    actor.ode_params = ODEparams(;)
    #dd.global_time,par.t_final, par.time_steps)
    #pressure = dd.equilibrium.time_slice[].pressure

    actor.ode_params = set_up_ode_params!(dd, par, actor.ode_params)
    
    bifurcation_bounds = calculate_bifurcation_bounds(dd, par, actor.ode_params)

    # Solve a single case to debug
    # pick Control1 and Control 2
    C1 = .1; C2 = 5.; task = "solveRW"
    solve_one_case(par, actor.ode_params, task, C1, C2 )

    #all_sols = solve_system(actor, task)
    #control1 = actor.ode_params.Control1
    #control2 = actor.ode_params.Control2
    #inputs = [(c1, c2) for (c1, c2) in zip(control1, control2)]
    #inputs = vec(inputs)  # flatten

    # NEW (distributed over the whole control grid):
    #finals = solve_system(actor, task)

    return actor
end

function solve_one_case(par, ode_params::ODEparams, task::String, eps::Float64, Om0::Float64)    
    final_sol = solve_ODEs(par, ode_params, task, eps, Om0)
    println("final raw solution = ", final_sol)
    
    sols_norm = normalize_ode_results(final_sol, ode_params, eps, Om0, par.control_type)
    println("final normalized solution", sols_norm)

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
    
    # SIM time is obsolete since the integrator uses adaptive stepping
    #ode_params.sim_time = set_sim_time(dd.global_time, par.t_final, par.time_steps)
    
    # find the normalized radius of the q=2 surface
    q_prof = dd.equilibrium.time_slice[].profiles_1d.q
    rho = dd.equilibrium.time_slice[].profiles_1d.rho_tor_norm
    ode_params.rat_surface = find_rat_surface(q_prof, rho, par.q_surf)
    
    # calculate the stability indices and mutual inductances
    ode_params = calculate_stability_index!(dd, par, ode_params)

    # Set physical parameters in dimensionless form
    ode_params = set_phys_params!(dd, par, ode_params)

    # Prepare control parameters based on the control type
    ode_params = set_control_parameters!(dd, par, ode_params)

    return ode_params #ODEparams(;sim_time, rat_surface=rat_surf,rat_index,stability_index=Deltat, DeltaW=Deltaw)
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
    @time rho_rat = Roots.secant_method(f, (rho[x0-1], rho[x0]))
    println("Found q=$rat_surface surface at: ", rho_rat)
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

    return ode_params
end

function set_phys_params!(dd::IMAS.DD, par, ode_params::ODEparams)
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

    # figure out the rotation rate at the q=2 surface
    rt = ode_params.rat_surface
    rho = dd.core_sources.source[1].profiles_1d[1].grid.rho_tor_norm
    rot_core = dd.core_profiles.profiles_1d[1].rotation_frequency_tor_sonic
    rot_interp = IMAS.interp1d(rho, rot_core)
    Omega0 = rot_interp(rt) * 2 * π * par.t0 # Convert to
    println("Dimensionless Omega0 at q=2 surface: ", Omega0)

    Om0Vals = range(1.0e-2, Omega0, length=N) |> collect
    Om0s = repeat(Om0Vals, 1, M)
    Control1 = vec(Om0s)
    ode_params.Control1 = Control1

    # Initialize the other control parameter based on the control type
    if control_type == :EF
        EpsUp = par.MaxErrorField / (par.b0 * par.r0)  # Convert to dimensionless units
        Control2_vals = range(1e-2, EpsUp, length=M) |> collect

    elseif control_type == :LinStab
        ode_params.error_field = 0.5  # Example value for error field
        Control2_vals = range(DeltaLow, DeltaUp, length=M) |> collect
        
    elseif control_type == :NLsaturation
        Control2_vals = range(ode_params.alpha_lower, ode_params.alpha_upper, length=M) |> collect
        
        #ode_params.saturation_param = Control2
    end
    
    Control2 = vec(repeat(Control2_vals', N, 1))
    ode_params.Control2 = Control2

    return ode_params
end


"Apply control_type adjustments"
function control_adjustments(ode_params::ODEparams, Control1::Float64, control_type::Symbol)
    alpha = ode_params.saturation_param
    errF = ode_params.error_field
    Deltat = ode_params.stability_index

    if control_type == :NLsaturation
        alpha = Control1
    elseif control_type == :EF
        errF = Control1
    elseif control_type == :StabIndex
        Deltat = Control1
    end

    return alpha, errF, Deltat
end

"Right-hand side for RW system"
function rhs_RW!(dydt, y, t, eps::Float64, Om0::Float64, ode_params::ODEparams, n_mode::Int, control_type::Symbol)
    m0 = n_mode
    DeltaW = ode_params.DeltaW
    rt = ode_params.rat_surface
    l21 = ode_params.l21
    l12 = ode_params.l12
    l32 = ode_params.l32
    Tt_Tw = ode_params.Taut_Tauw
    mu = ode_params.mu
    I = ode_params.Inertia

    alpha, errF, Deltat = control_adjustments(ode_params, eps, control_type)

    psi, theta, Om, psiW, thW = y

    dydt[1] = Deltat * psi * (1.0 + alpha * abs(psi)) + l21 * psiW * cos(theta - thW)
    dydt[2] = -m0 * Om - l21 * psiW * sin(theta - thW) / psi
    dydt[3] = (rt * l21 * psiW * psi * sin(theta - thW) + mu * (Om0 - Om)) / I
    dydt[4] = Tt_Tw * (DeltaW * psiW + l12 * psi * cos(theta - thW) + l32 * errF * cos(thW))
    dydt[5] = Tt_Tw * (l12 * psi * sin(theta - thW) - l32 * errF * sin(thW)) / psiW
end

"Right-hand side for basic system"
function rhs_basic!(dydt, y, t, eps, Om0, ode_params::ODEparams, n_mode::Int, control_type::Symbol)
    m0 = n_mode
    rt = ode_params.rat_surface
    l21 = ode_params.l21
    mu = ode_params.mu
    I = ode_params.Inertia

    alpha, errF, Deltat = control_adjustments(ode_params, eps, control_type)

    psi, theta, Om = y

    dydt[1] = Deltat * psi * (1.0 + alpha * abs(psi)) + l21 * errF * cos(theta)
    dydt[2] = -m0 * Om - l21 * errF * sin(theta) / psi
    dydt[3] = (rt * l21 * errF * psi * sin(theta) + mu * (Om0 - Om)) / I
end

"Return the correct RHS function for a given task"
function make_rhs_function(task::String)
    #task == "solveRW" ? rhs_RW! :
    #task == "solve"   ? rhs_basic! :
    #error("Unknown task: $task")
    if task == "solveRW"
        return rhs_RW!
    elseif task == "solve"
        return rhs_basic!
    else
        error("Unknown task: $task")
    end
end

function make_initial_condition(dims::Vector{Float64}, task::String)
    
    if task == "solveRW"
        # 5D system
        return [
            rand() * (dims[1] - 0.001) + 0.001,   # y1
            rand() * 2π - π,                      # y2
            rand() * (dims[2] - 0.001) + 0.001,   # y3
            rand() * (dims[3] - 0.001) + 0.001,   # y4
            rand() * 2π - π                       # y5
        ]

    elseif task == "solve"
        # 3D system
        return [
            rand() * (dims[1] - 0.001) + 0.001,   # y1
            rand() * 2π - π,                      # y2
            rand() * (dims[2] - 0.001) + 0.001    # y3
        ]

    else
       error("Unknown task: $task")
    end
end

function make_ode_func(rhs!)
    return function (du, u, p, t)
        # p must match this tuple layout:
        # (ode_params, eps, Om0, n, control_type)
        X1, X2, ode_params, n_mode, control_type = p
        rhs!(du, u, t, X1, X2, ode_params, n_mode, control_type)
    end
end


function calculate_bifurcation_bounds(dd::IMAS.DD, par, ode_params::ODEparams)
    """
    Calculate the bifurcation boundaries analytically
    Note: Output will NOT be correct if NLsatON is true
    """

    N = par.grid_size
    M = par.grid_size

    Y = ode_params.Control1
    m0 = par.n_mode
    mu = ode_params.mu
    DeltatRW = par.RPRW_stability_index
    DeltaW = ode_params.DeltaW
    rt = ode_params.rat_surface
    l21 = ode_params.l21
    l12 = ode_params.l12
    l32 = ode_params.l32
    
    RWon = false
    if par.RPRW_stability_index !== nothing
        RWon = true
    end
    
    if par.control_type == :EF
        X = ode_params.Control2
        Deltat = ode_params.stability_index
        if RWon
            q = (DeltatRW / m0)^2 .+ rt * (l32 * l21 * X / DeltaW).^2 / (m0 * mu)
            r = -Y * DeltatRW^2 / m0^2
        else
            q = (Deltat / m0)^2 .+ (l21 * X).^2 / mu
            r = -Y * Deltat^2 / m0^2
        end
    elseif par.control_type == :Stab
        Eps = ode_params.error_field
        X = ode_params.Control2
        DeltatRW = X - l21*l21/DeltaW # overwrite the stability index since Deltat is control
        if RWon
            q = (DeltatRW / m0)^2 .+ rt * (l32 * l21 * Eps / DeltaW).^2 / (m0 * mu)
            r = -Y * DeltatRW^2 / m0^2
        else
            q = (Deltat / m0)^2 + (l21 * Eps).^2 / mu
            r = -Om0s * Deltat^2 / m0^2
        end
    end
    
    # remaining coefficients
    a = -(Y.^2) / 3.0 .+ q
    b = 2.0 * (-Y).^3 / 27.0 - q .* (-Y) / 3.0 .+ r
    bifurcation_bounds = b.^2 / 4 + a.^3 / 27 
    bifurcation_bounds = reshape(bifurcation_bounds, M, N)
    
    make_contour(ode_params.Control2, ode_params.Control1, bifurcation_bounds, [0.0, 1e-10], par.control_type)
    #make_contour(ode_params.Control2, ode_params.Control1, bifurcation_bounds)


    return bifurcation_bounds
end


function solve_ODEs(par, ode_params::ODEparams, task::String,
                        control1::Float64, control2::Float64)
        
    n_mode = par.n_mode
    control_type = par.control_type
    t_final = par.t_final

    rhs! = make_rhs_function(task)
    ode_rhs! = make_ode_func(rhs!)
    y0 = make_initial_condition(ode_params.hyper_cube_dims, task)
    p = (control1, control2, ode_params, n_mode, control_type)
    tspan = (0.0, t_final)

    prob = ODEProblem(ode_rhs!, y0, tspan, p)
    sol = solve(prob, Tsit5(); saveat=t_final, reltol=1e-8, abstol=1e-10)

    final = sol.u[end]

    return final
end

"""
normalize_results(results, ode_params, eps, Om0)

Normalize the final solutions from ODE runs.

- `results` may be:
    * a single solution vector (Vector{Float64}), or
    * a Vector of solution vectors (Vector{Vector{Float64}}).
- `eps` and `Om0` may be:
    * single Float64 values (if results is one vector), or
    * Vector{Float64} of the same length as results.

Normalization:
    psiN  = final_sol[1] * (Deltat * DeltaW - l12 * l21) / (l32 * l21 * eps)
    psiwN = final_sol[4] * (Deltat * DeltaW - l12 * l21) / (l32 * abs(Deltat) * eps)
    OmN   = final_sol[3] / Om0
"""
function normalize_ode_results(results, ode_params::ODEparams, Control1, Control2, control_type)
    # Extract parameters
    l12    = ode_params.l12
    l21    = ode_params.l21
    l32    = ode_params.l32
    DeltaW = ode_params.DeltaW

    # check controls
    if control_type != :Stab 
        Deltat = ode_params.stability_index
    else
        Deltat = Control1
    end

    
    # Normalization for one solution
    function normalize_one(final_sol::AbstractVector{<:Real}, eps_val::Float64, Om0_val::Float64)
        if length(final_sol) == 5 # assumes solveRW layout
            num = Deltat * DeltaW - l12 * l21

            psiN   = final_sol[1] * num / (l32 * l21 * eps_val)
            theta_t = mod(final_sol[2], 2π)  # Normalize phase
            OmN    = final_sol[3] / Om0_val
            psiwN = final_sol[4] * num / (l32 * abs(Deltat) * eps_val)  
            theta_w = mod(final_sol[5], 2π)  # Normalize phase
            return (psiN, theta_t, OmN, psiwN, theta_w)

        elseif length(final_sol) == 3
            psiN   = final_sol[1] / (l21 * eps_val)
            theta_t = mod(final_sol[2], 2π)  # Normalize phase
            OmN    = final_sol[3] / Om0_val
            return (psiN, theta_t, OmN)

        else
            throw(ArgumentError("Unexpected final_sol length: $(length(final_sol))"))
        end
    end  

    if isa(results, AbstractVector{<:Real}) && isa(Control1, Real) && isa(Control2, Real)
        # Single case
        return normalize_one(results, Control1, Control2)

    elseif isa(results, AbstractVector{<:AbstractVector{<:Real}}) &&
           isa(Control1, AbstractVector{<:Real}) &&
           isa(Control2, AbstractVector{<:Real})
        # Many cases
        length(results) == length(Control1) == length(Control2) ||
            throw(ArgumentError("results, eps, and Om0 must all have the same length"))
        return [normalize_one(sol, e, o) for (sol, e, o) in zip(results, Control1, Control2)]

    else
        throw(ArgumentError("Input types do not match expected patterns"))
    end
end


function solve_system(actor::ActorLocking, task::String)
    par = actor.par
    ode_params = actor.ode_params

    n_mode = par.n_mode
    control_type = par.control_type
    t_final = par.t_final
    
    # Get the control params
    control1 = ode_params.Control1
    control2 = ode_params.Control2

    # NOTE: rhs expects (eps, Om0), and in your EF case Control2 = eps, Control1 = Om0  
    inputs = collect(zip(control1, control2))

    # Optional: shrink what we ship by clearing large control arrays
    ode_params_send = deepcopy(ode_params)
    ode_params_send.Control1 = Float64[]
    ode_params_send.Control2 = Float64[]

    # Parallel map over the grid, returning final states
    finals = pmap(inputs) do (C1, C2)
        solve_ODEs(par, ode_params_send, task, C2, C1)
    end

    # finals[i] is the final state vector for inputs[i].
    # If you want them back on a 2D grid (N x M), reshape here:
    # finals_mat = reshape(finals, (par.grid_size, par.grid_size))

    return finals
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
function plot_sols_scatter(norm_sol; xcol::Int=1, ycol::Int=3)
    if eltype(norm_sol) <: AbstractVector{<:Real}
        # Convert vector-of-vectors to a matrix
        data = reduce(vcat, (x' for x in norm_sol))
    elseif eltype(norm_sol) <: NTuple
        # Convert vector-of-tuples to a matrix
        data = reduce(vcat, (collect(x)' for x in norm_sol))
    else
        throw(ArgumentError("norm_sol must be Vector{Vector{Float64}} or Vector{NTuple{N,Float64}}"))
    end

    scatter(data[:, xcol], data[:, ycol],
            xlabel="Component $xcol",
            ylabel="Component $ycol",
            title="Normalized solution scatter")
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
    
    
    plt = contour(X,Y, Z; linewidth=2)
    display(plt)   # explicitly display
    return plt
end


function make_contour(X::AbstractArray, Y::AbstractArray, Z::AbstractMatrix, levels::Vector{Float64}, control_type)
    # Determine target shape
    m, n = size(Z)

    # If C1 and C2 are vectors, try to reshape them
    if ndims(X) == 1 
        X = unique(X)
    end
    if ndims(Y) == 1 
        Y = unique(Y)
    end

    if control_type==:EF
        xlabel = "Error Field"
    elseif control_type==:LinStab
        xlabel = "Linear Stability"
    elseif control_type==:NLsaturation
        xlabel = "NL saturation"
    end

    ## Create the contour plot
    #pythonplot()
    plt = contour(X, Y, Z; levels=levels, linewidth=2, xlabel=xlabel, ylabel="Normalized Torque")
    #plt.xlabel("Control1")
    display(plt)  # show plot

    return plt
end



