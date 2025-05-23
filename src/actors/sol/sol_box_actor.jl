#= ===================== =#
#       ActorSolBox       #
#= ===================== =#

# Part of act structure (act.ActorSolBox here)
Base.@kwdef mutable struct FUSEparameters__ActorSolBox{T<:Real} <: ParametersActor{T}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    _time::Float64 = NaN
    # Define parameters in act.ActorSolBox here: - ONLY INPUTS
    # Syntax is: name::Entry{Type} = Entry{Type}("unit", "Description of parameter"; default = default value)
    Te_t::Entry{T}   = Entry{T}("eV","Input electron temperature at the target"; default = 10.0)
    Ti_t::Entry{T}   = Entry{T}("eV","Input ion temperature at the target"; default = 10.0)
    frac_cond::Entry{T}   = Entry{T}("-","Fraction of power carried by electron conduction"; default = 0.7)
    frac_mom::Entry{T}   = Entry{T}("-","Fraction of momentum lost due to collisions with neutrals, atomic processes and viscous forces"; default = 0.5)
    κ0_e::Entry{T}   = Entry{T}("-","Coefficient of electron conductivity"; default = 2000.0)
    κ0_i::Entry{T}   = Entry{T}("-","Coefficient of ion conductivity"; default = 60.0)
    qpar_i::Entry{T}   = Entry{T}("W m^-2","Upstream parallel ion heat flux"; default = 10.0e+09)
    qpar_e::Entry{T}   = Entry{T}("W m^-2","Upstream parallel electron heat flux"; default = 10.0e+09)
    Γ_i::Entry{T}   = Entry{T}("pcles m^-2 s^-1","Upstream ion particle flux"; default = 1.0e+22)
    Γ_e::Entry{T}   = Entry{T}("pcles m^-2 s^-1","Upstream electron particle flux"; default =1.0e+22)
    mass_ion::Entry{T}   = Entry{T}("kg","Ion mass in multiples of amu"; default = 2.5)
    recycling_coeff_i::Entry{T}   = Entry{T}("-","Ion particle recycling coefficient"; default = 0.98)
    recycling_coeff_e::Entry{T}   = Entry{T}("-","Electron particle recycling coefficient"; default = 0.98)
    # Plot switch
    do_plot::Entry{Bool} = act_common_parameters(; do_plot=false)
end

# This struct is can be passed around by FUSE - do not want to pass Act around all the time
# create struct for the actor - INPUTS and anything that should be saved
mutable struct ActorSolBox{D,P} <: SingleAbstractActor{D,P}
    dd::IMAS.dd{D}
    par::OverrideParameters{P,FUSEparameters__ActorSolBox{P}}
    # insert here the structs or variables that must be passed in the actor while running
    # i.e. the outputs of the step function must be passed to the finalize function to populate the dd or saved in general
    Te_u::Union{Nothing, Float64} # Electron temperature at the upstream (separatrix) [eV]
    Ti_u::Union{Nothing, Float64} # Ion temperature at the upstream (separatrix) [eV]
    ne_t::Union{Nothing, Float64} # Electron density at the target [pcles m^-2 s^-1]
    ni_t::Union{Nothing, Float64} # Ion density at the target [pcles m^-2 s^-1]
    ne_u::Union{Nothing, Float64} # Electron density at the upstream (separatrix) [pcles m^-2 s^-1]
    ni_u::Union{Nothing, Float64} # Ion density at the upstream (separatrix) [pcles m^-2 s^-1]
    cst::Union{Nothing, Float64} # Sound speed at the sheath entrance [m s^-1]
    sol_connection_length::Union{Nothing, Float64} # Parallel connection length [m]
    sol_total_Fx::Union{Nothing, Float64} # Total flux expansion [-]
end

# define actor and its dispatced functions (leave as is)
function ActorSolBox(dd::IMAS.dd{D}, par::FUSEparameters__ActorSolBox{P}; kw...) where {D<:Real,P<:Real}
    logging_actor_init(ActorSolBox)
    par = OverrideParameters(par; kw...)
    return ActorSolBox(dd, par, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing) # should be the same number of inputs as struct above
end

"""
    ActorSolBox(dd::IMAS.dd, act::ParametersAllActors; kw...)

    Box model for the Scrape-Off layer. It returns plasma parameters at the divertor targets with the heat load deposited at the strike points 
    (change this description!)
"""
# FUSE.ActorSolBox(dd,act) will trigger: Step -> Fianlize 
function ActorSolBox(dd::IMAS.dd, act::ParametersAllActors; kw...)
    actor = ActorSolBox(dd, act.ActorSolBox; kw...)
    step(actor)
    finalize(actor)
    return actor
end

# Step function
function _step(actor::ActorSolBox{D,P}) where {D<:Real, P<:Real}
    # retrieve variables defined in the struct ActorSolBox
    dd = actor.dd
    par = actor.par

    # Define physical constants
    amu = 1.660538921 * 1.0e-27 # Atomic mass unit [kg]
    e_charge = 1.60e-19 # Magnitude of the electron charge [C]

    # Write here all your models
    # To call the act parameters use: par.parameter
    # Example: par.Te_t,
    
    # Step 1 - get connection length and total flux expansion for the separatrix flux tube between midplane and outer target

    # Trace field lines in the SOL
    sol = IMAS.sol(dd; levels=1)

    # The first (and only) traced flux surface in the LFS region is the primary separatrix
    separatrix = sol[:lfs][1]

    # Parallel connection length from midplane to target (last point)
    actor.sol_connection_length = separatrix.s[end]

    # Total flux expansion from midplane to target (last point)
    actor.sol_total_Fx = separatrix.total_flux_expansion[end]

    # Step 2 - calculate the sound speed
    actor.cst = sqrt(e_charge*(par.Te_t + par.Ti_t)/(par.mass_ion*amu))

    # Step 3 - calculate the upstream temperature(s)

    # Intermediate variable. Calculate once for reduced computation.
    alpha = par.frac_cond*1.75*actor.sol_connection_length*(1.0-0.5*(actor.sol_total_Fx - 1.0))

    # Upstream ion temperature
    actor.Ti_u = ( (par.Ti_t)^(7.0/2.0) + alpha*(par.qpar_i/par.κ0_i) )^(2.0/7.0)
    
    # Upstream electron temperature
    actor.Te_u = ( (par.Te_t)^(7.0/2.0) + alpha*(par.qpar_e/par.κ0_e) )^(2.0/7.0)

    # Step 4 - calculate the target densities

    # Calculate the particle transmission coefficients

    # Ion recycling coefficient
    γi = 1.0 - par.recycling_coeff_i

    # Electron recycling coefficient
    γe = 1.0 - par.recycling_coeff_e

    # Intermediate variable. Calculate once for reduced computation.
    beta = (actor.sol_total_Fx + 1)/(2.0*actor.sol_total_Fx)

    # Target ion density
    actor.ni_t = beta * (par.Γ_i/(γi*cst))

    # Target electron density
    actor.ne_t = beta * (par.Γ_e/(γe*cst))

    # Step 5 - calculate the upstream densities

    # Upstream ion density
    actor.ni_u = ((2.0*ni_t)/(par.frac_mom)) * ((par.Ti_t + par.Te_t)/(actor.Ti_u + actor.Te_u))
    
    # Upstream electron density
    actor.ne_u = ((2.0*ne_t)/(par.frac_mom)) * ((par.Ti_t + par.Te_t)/(actor.Ti_u + actor.Te_u))

    #plot
    if par.do_plot
        #probably you will plot something
    end

    # # save output in struct actor to be passed to finalize
    # # Example:
    # actor.wall_heat_flux = Float64[] # using the example I used before
    return actor
end

# Finalize function
function _finalize(actor::ActorSolBox)
    # Finalization after computation is done in step
    # (Probably populate dd; for details of dd: https://fuse.help/dev/dd.html)
    return actor
end