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
    Te_u::Entry{T}   = Entry{T}("eV","Input electron temperature at the upstream"; default = 100.0)
    Ti_u::Entry{T}   = Entry{T}("eV","Input ion temperature at the upstream"; default = 100.0)
    frac_cond::Entry{T}   = Entry{T}("-","Fraction of power carried by electron conduction"; default = 0.7)
    κ0_e::Entry{T}   = Entry{T}("-","Coefficient of electron conductivity"; default = 2000.0)
    κ0_i::Entry{T}   = Entry{T}("-","Coefficient of ion conductivity"; default = 60.0)
    qpar_i::Entry{T}   = Entry{T}("W m^-2","Upstream parallel ion heat flux"; default = 1.0e+09)
    qpar_e::Entry{T}   = Entry{T}("W m^-2","Upstream parallel electron heat flux"; default = 1.0e+09)
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
    Te_t::Union{Nothing, Float64} # Electron temperature at the target [eV]
    Ti_t::Union{Nothing, Float64} # Ion temperature at the target [eV]
    sol_connection_length::Union{Nothing, Float64} # Parallel connection length [m]
    sol_total_Fx::Union{Nothing, Float64} # Total flux expansion [-]
end

# define actor and its dispatced functions (leave as is)
function ActorSolBox(dd::IMAS.dd{D}, par::FUSEparameters__ActorSolBox{P}; kw...) where {D<:Real,P<:Real}
    logging_actor_init(ActorSolBox)
    par = OverrideParameters(par; kw...)
    return ActorSolBox(dd, par, nothing, nothing, nothing, nothing)
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

    # Write here all your models
    # To call the act parameters use: par.parameter
    # Example: par.Te_t,
    
    # Step 1 - get connection length and total flux expansion for the separatrix flux tube between midplane and outer target

    # Trace field lines in the SOL
    sol = IMAS.sol(dd; levels=1)

    # The first (and only) traced flux surface in the LFS region is the primary separatrix
    separatrix = sol[:lfs][1]

    # Parallel connection length from midplane to target (last point)
    actor.sol_connection_length = separatrix[s][end]

    # Total flux expansion from midplane to target (last point)
    actor.sol_total_Fx = separatrix[total_flux_expansion][end]

    # Step 2 - calculate the target temperature(s)

    # Intermediate variable. Calculate once for reduced computation
    alpha = par.frac_cond*1.75*actor.sol_connection_length*(1.0-0.5*(actor.sol_total_Fx - 1.0))

    # Target ion temperature
    actor.Ti_t = ( par.Ti_u^3.5 - (par.qpar_i/par.κ0_i)*alpha)^(1.0/3.5)

    # Target electron temperature
    actor.Te_t = ( par.Te_u^3.5 - (par.qpar_e/par.κ0_e)*alpha)^(1.0/3.5)

    #plot
    if par.do_plot
        #probably you will plot something
    end

    # save output in struct actor to be passed to finalize
    # Example:
    actor.wall_heat_flux = Float64[] # using the example I used before
    return actor
end

# Finalize function
function _finalize(actor::ActorSolBox)
    # Finalization after computation is done in step
    # (Probably populate dd; for details of dd: https://fuse.help/dev/dd.html)
    return actor
end