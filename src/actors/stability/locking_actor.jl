Base.@kwdef mutable struct ODEparams
    sim_time::Vector{Float64} = Float64[]# simulation time vector
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
    control_type::Switch{Symbol} = Switch{Symbol}([:EF, :StabIndex, :SatParam], # EF: error field
        "-",                                                           # StabIndex: vary stability_index,
        "Use a user specified Control case to run the locking models") # SatParam: vary NL saturation
    NL_saturation_ON::Entry{Bool} = Entry{Bool}("-", "Nonlinear saturation parameter for the mode"; default=true)
    RPRW_stability_index::Entry{Float64} = Entry{Float64}(
        "-", 
        "Stability index of the system (set to Neg. value for now)"; default=-0.5)
    max_tor_rot::Entry{Float64} = Entry{Float64}(
        "2pi*Hertz", 
        "Maximum angular rotation rate observed in the shot")
    b0::Entry{Float64} = Entry{Float64}(
        "Tesla", 
        "Scale for magnetic perturbations, usually ~10Gauss"; default=1.e-3)
    t0::Entry{Float64} = Entry{Float64}(
        "seconds", 
        "Time scale for the integration , usually TM/RW growth rate"; default=1.e3)
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

    actor.ode_params = set_up_ode_params(dd, par, actor.ode_params)
    #actor.ode_params = calculate_inductances(dd, par, actor.ode_params)
    
    #actor.ode_params = set_phys_params(dd, par, actor.ode_params)

    
    return actor
end

function _finalize(actor::ActorLocking)
    print("WHAT?")
    return actor
end

function set_up_ode_params(dd, par, ode_params)
    """
    Initialize the ODE parameters for the locking actor.
    
    Args:
        dd: IMAS data structure
        par: Parameters for the simulation
    Returns:
        ode_params: Initialized ODE parameters
    """
    
    #ode_params = ODEparams
    ode_params.sim_time = set_sim_time(dd.global_time, par.t_final, par.time_steps)
    
    # find the normalized radius of the q=2 surface
    q_prof = dd.equilibrium.time_slice[].profiles_1d.q
    rho = dd.equilibrium.time_slice[].profiles_1d.rho_tor_norm
    ode_params.rat_surface, ode_params.rat_index = find_q2_surface(q_prof, rho, par.q_surf)
    println(ode_params.rat_surface, ode_params.rat_index)

    # Set physical parameters in dimensionless form
    ode_params.stability_index, ode_params.DeltaW = calculate_inductances(dd, par, ode_params)

    #set_phys_params(dd, par, ODEparams(;sim_time, rat_surface=rat_surf, rat_index=rat_index))

    # Prepare control parameters based on the control type
    #ode_params = prepare_control_parameters(par.control_type, ode_params, par)

    return ode_params #ODEparams(;sim_time, rat_surface=rat_surf,rat_index,stability_index=Deltat, DeltaW=Deltaw)
end

function set_sim_time(time0::Float64, tfinal::Float64, time_steps::Int64) 
    """
    initialize time for the ODE integration
    """
    sim_time = collect(range(time0, tfinal, length=time_steps))
    
    return sim_time
end


function find_q2_surface(q_prof, rho, rat_surface)
    """
    Find the location of the q=2 surface given the q profile
    """

    # Interpolate rho(q) to find the q=2 surface
    # Use absolute values to ensure we find the correct surface
    rho_interp = IMAS.interp1d(abs.(q_prof), rho)
    rho_rat = rho_interp.(2.0)  
    rat_indx = findfirst(x-> x > rat_surface, abs.(q_prof))
    println("Found q=2 surface at: ", rho_rat)
    return rho_rat, rat_indx
end

function calculate_inductances(dd, par, ode_params)
    rt = ode_params.rat_surface  
    rw = ode_params.res_wall
    rc = ode_params.control_surf
    m0 = par.n_mode
    
    rat21 = (rw / rt)^m0
    rat12 = rat21^(-1)
    rat32 = (rc / rw)^m0
    rat23 = rat32^(-1)

    l12 = (2 * m0 / rw) / (rat21 - rat12)
    l21 = (2 * m0 / rt) / (rat21 - rat12)
    l32 = (2 * m0 / rw) / (rat32 - rat23)

    DeltaWl = (m0 / rw) * (rat21 + rat12) / (rat21 - rat12)
    DeltaWr = -(m0 / rw) * (rat32 + rat23) / (rat32 - rat23)
    DeltaW = DeltaWr - DeltaWl
    Deltat = par.RPRW_stability_index + l21 * l12 / DeltaW

    #ode_params.DeltaW = DeltaW
    #ode_params.stability_index = Deltat
    
    return Deltat, DeltaW
end

function set_phys_params(dd, par, ode_params)
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
    mass_dens = mass_ion * dd.core_profiles.profiles_1d[1].ion[1].density_thermal[83]
    #Zeff = dd.core_profiles.profiles_1d.zeff[40]
 
    # Set all the scales
    psi0 = par.b0 * par.r0
    U0 = psi0^2 * par.r0 / mu0_val
    R0 = dd.equilibrium.vacuum_toroidal_field.r0
    # approximate core rotation
    rot_core = dd.core_profiles.profiles_1d[1].rotation_frequency_tor_sonic[10]

    # Calculate the drag coefficient in SI
    muSI = par.NBItorque / rot_core
    
    # Calculate first the moment of inertia in the layer  
    inertia = (2*π)^2 * mass_dens * R0 * rt^3 * ode_params.layer_width

    # Set the dimensionless quantities
    ode_params.mu = muSI / (U0 * par.t0)
    ode_params.Inertia = inertia / (U0 * par.t0^2)

    return ode_params
end


function prepare_control_parameters(control_type::Symbol, ode_params::ODEparams, par::FUSEparameters__ActorLocking)
    """
    Prepare the control parameters based on the control type.
    
    Args:
        control_type: Symbol indicating the type of control
        ode_params: ODE parameters to be modified
        par: Parameters for the simulation
    Returns:
        ode_params: Updated ODE parameters with control settings
    """

    Om0Vals = range(1.0e-2, par.Omega0_max, length=par.N) |> collect

    if control_type == :EF
        Eps = par.eps
        EpsUp = par.MaxErrorField
        psiWvals = range(1e-2, EpsUp, length=par.M) |> collect
        Om0s = repeat(Om0Vals, 1, par.M)
        Wflux = repeat(psiWvals', par.N, 1)
        par.Control2 = vec(Wflux')

        ode_params.error_field = par.b0 * 1.e-3  # Example value for error field
    elseif control_type == :StabIndex
        ode_params.stability_index = par.RPRW_stability_index
    elseif control_type == :SatParam
        ode_params.saturation_param = par.NL_saturation_ON ? 0.2 : 0.0
    end
    
    return ode_params
end


