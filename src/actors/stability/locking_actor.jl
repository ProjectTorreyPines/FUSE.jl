using DifferentialEquations
#import DifferentialEquations as DiffEqs
using Distributed

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


#= ================= =#
#  ActorLockingProbability  #
#= ================= =#
Base.@kwdef mutable struct FUSEparameters__ActorLocking{T<:Real} <: ParametersActor{T}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    _time::Float64 = NaN
    n_mode::Entry{Int} = Entry{Int}("_", "toroidal mode number of the mode"; default=1)
    q_surf::Entry{Float64} = Entry{Float64}("_", "rational surface of interest, usually 2.0"; default=2.)
    grid_size::Entry{Int} = Entry{Int}("-", "grid resolution for control space"; default=100)
    t_final::Entry{Float64} = Entry{Float64}("-", "Final integration time in units of tearing time (~ms)")
    time_steps::Entry{Int} = Entry{Int}("-", "number of time steps for the ODE integration"; default=200)
    control_type::Switch{Symbol} = Switch{Symbol}([:EF, :LinStab, :NLsaturation], # EF: error field
        "-",                                                           # StabIndex: vary stability_index,
        "Use a user specified Control case to run the locking models") # SatParam: vary NL saturation
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

    task = "solveRW"  #MAKE this an ODEparam or param
    sols = solve_ODEs(par, actor.ode_params, task, 5., 0.5)

    #all_sols = solve_system(actor, task)
    #control1 = actor.ode_params.Control1
    #control2 = actor.ode_params.Control2
    #inputs = [(c1, c2) for (c1, c2) in zip(control1, control2)]
    #inputs = vec(inputs)  # flatten

    ## run distributed ODE solves
    #results = pmap(inputs) do (eps, Om0)
    #    solve_ODEs(par, ode_params, "LockingTask", eps, Om0)
    #end
    
    
    return actor
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
    ode_params.rat_surface = find_q2_surface(q_prof, rho, par.q_surf)
    
    # calculate the stability indices and mutual inductances
    ode_params = calculate_stability_index!(dd, par, ode_params)

    # Set physical parameters in dimensionless form
    ode_params = set_phys_params!(dd, par, ode_params)

    # Prepare control parameters based on the control type
    ode_params = set_control_parameters!(dd, par, ode_params)

    return ode_params #ODEparams(;sim_time, rat_surface=rat_surf,rat_index,stability_index=Deltat, DeltaW=Deltaw)
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


# function set_sim_time(time0::Float64, tfinal::Float64, time_steps::Int64) 
#     """
#     initialize time for the ODE integration
#     """
#     sim_time = collect(range(time0, tfinal, length=time_steps))
    
#     return sim_time
# end


function find_q2_surface(q_prof::Vector{Float64}, rho::Vector{Float64}, rat_surface::Float64)
    """
    Find the location of the q=2 surface given the q profile
    """

    # Interpolate rho(q) to find the q=2 surface
    # Use absolute values to ensure we find the correct surface
    rho_interp = IMAS.interp1d(abs.(q_prof), rho)
    rho_rat = rho_interp.(2.0)  
    rat_indx = findfirst(x-> x > rat_surface, abs.(q_prof))
    println("Found q=2 surface at: ", rho_rat)
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
    # Calculate the drag coefficient in SI
    #cp1d = dd.core_profiles.profiles_1d[]
    #total_source1d = IMAS.total_sources(dd.core_sources, cp1d; time0=dd.global_time, fields=[:torque_tor_inside])
    #dd.core_sources.source[15].global_quantities[1].torque_tor
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

    alpha, errF, Deltat = control_adjustments(ode_params, Control1, control_type)

    psi, theta, Om = y

    dydt[1] = Deltat * psi * (1.0 + alpha * abs(psi)) + l21 * errF * cos(theta)
    dydt[2] = -m0 * Om - l21 * errF * sin(theta) / psi
    dydt[3] = (rt * l21 * errF * psi * sin(theta) + mu * (Om0 - Om)) / I
end

"Return the correct RHS function for a given task"
function make_rhs_function(task::String)
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

function make_problem_factory(prob_template, ode_params, task, dt)
    return function(input)
        eps, Om0 = input
        y0 = make_initial_condition(ode_params.hyper_cube_dims, task)
        p = (ode_params, eps, Om0)
        prob = remake(prob_template; u0=y0, p=p)
        solve(prob, Tsit5(); saveat=dt, reltol=1e-9)
    end
end


## ---- Define worker solve (runs on each worker)
#@everywhere function SolveODE_RW(inputs, prob_template, task, ode_params)
#    control1, control2 = inputs
#    y0 = make_initial_condition(ode_params.hyper_cube_dims, task)
#    println(y0)
#    p = (ode_params, control1, control2)
#    prob = DiffEqs.remake(prob_template; u0=y0)
#    return DiffEqs.solve(prob, DiffEqs.Tsit5(), reltol=1e-9)
#end


function make_ode_func(rhs!)
    return function (du, u, p, t)
        # p must match this tuple layout:
        # (ode_params, eps, Om0, n, control_type)
        X1, X2, ode_params, n_mode, control_type = p
        rhs!(du, u, t, X1, X2, ode_params, n_mode, control_type)
    end
end


function solve_ODEs(par, ode_params::ODEparams, task::String, eps::Float64=1.0, Om0::Float64=1.0)
    
    n = par.n_mode
    control_type = par.control_type
    t_final = par.t_final
    DeltaW = ode_params.DeltaW
    Deltat = ode_params.stability_index
    l21 = ode_params.l21
    l12 = ode_params.l12
    l32 = ode_params.l32 

    # Create ODE problem
    tspan = (0.0, t_final)
    dt = t_final/100.
    #p = (ode_params, eps, Om0)  # Parameters
    p = (eps, Om0, ode_params, n, control_type)
    
    # --- pick correct RHS once
    rhs! = make_rhs_function(task)
    ode_rhs! = make_ode_func(rhs!)                  # wraps it to f!(du,u,p,t)
    
    # choose the inifial condition
    y0 = make_initial_condition(ode_params.hyper_cube_dims, task)
    
    # Define ODE function for DifferentialEquations.jl
    #function ode_func!(dydt, y, p, t)
    #    ode_params, eps, Om0 = p
    #    rhs!(dydt, y, t, eps, Om0, ode_params, n, control_type)
    #end

    prob = ODEProblem(ode_rhs!, y0, tspan, p)
    sol = solve(prob, Tsit5(), saveat=dt, reltol=1e-9)

    if task == "solveRW"
        # Initial conditions: [psiMag, theta, Om, psiW, thetaW]
        
        # Extract final solution
        final_sol = sol.u[end]
        final_sol[2] = mod(final_sol[2], 2π)  # Normalize phase
        final_sol[5] = mod(final_sol[5], 2π)  # Normalize phase

        println("Raw sol = ", final_sol)

        
        # Normalized solutions
        psiN = final_sol[1] * (Deltat * DeltaW - l12 * l21) / (l32 * l21 * eps)
        psiwN = final_sol[4] * (Deltat * DeltaW - l12 * l21) / (l32 * abs(Deltat) * eps)
        
        println("Normalized psit, psiW, Omt = ", psiN, ", ", psiwN, ", ", final_sol[3]/Om0)

        #plt = plot(sol.u[:, 1]); display(plt)

    elseif task == "solve"
        # Define ODE function for DifferentialEquations.jl
        # function ode_func!(dydt, y, p, t)
        #     ode_params, eps, Om0 = p
        #     rhs!(dydt, y, t, eps, Om0, ode_params, control_type)
        # end

        
        # Implement 3rd order system solution
        println("3rd order system solution not yet implemented")
        return nothing
    else
        throw(ArgumentError("Unknown task: $task"))
    end

    
    return sol.u
end

function solve_system(actor::ActorLocking, task::String)
    par = actor.par
    ode_params = actor.ode_params

    n = par.n_mode
    control_type = par.control_type
    t_final = par.t_final
    
    # Get the control params
    control1 = ode_params.Control1
    control2 = ode_params.Control2


    # Create ODE problem
    tspan = (0.0, t_final)
    dt = t_final/100.
    
    # --- pick correct RHS once
    rhs! = make_rhs_function(task)
    
    # ---- Build a template problem once (dummy y0)
    y0_dummy = make_initial_condition(ode_params.hyper_cube_dims, task)
    prob_template = ODEProblem(ode_func!, y0_dummy, tspan, (ode_params, 0.1, 0.1))

    # Run distributed solve
    results = pmap(SolveODE_RW, inputs)

    return results
    #sols = pmap(inp -> SolveODE_RW(inp, prob_template, task, ode_params), inp)
    # return prob_template
end
