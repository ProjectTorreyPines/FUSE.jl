#= ===================== =#
#  ActorAnalyticPedestal  #
#= ===================== =#
Base.@kwdef mutable struct FUSEparameters__ActorAnalyticPedestal{T<:Real} <: ParametersActor{T}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    _time::Float64 = NaN
    #== common pedestal parameters==#
    rho_nml::Entry{T} = Entry{T}("-", "Defines rho at which the no man's land region starts")
    rho_ped::Entry{T} = Entry{T}("-", "Defines rho at which the pedestal region starts") # rho_nml < rho_ped
    T_ratio_pedestal::Entry{T} =
        Entry{T}("-", "Ratio of ion to electron temperatures (or rho at which to sample for that ratio, if negative; or rho_nml-(rho_ped-rho_nml) if 0.0)"; default=1.0)
    Te_sep::Entry{T} = Entry{T}("-", "Separatrix electron temperature"; default=80.0, check=x -> @assert x > 0 "Te_sep must be > 0")
    ip_from::Switch{Symbol} = switch_get_from(:ip)
    βn_from::Switch{Symbol} = switch_get_from(:βn)
    ne_from::Switch{Symbol} = switch_get_from(:ne_ped)
    zeff_from::Switch{Symbol} = switch_get_from(:zeff_ped)
    #== actor parameters ==#
	model::Switch{Symbol} = Switch{Symbol}([:MAST, :NSTX], "-", "Pedestal width model, w_ped~beta_p,ped^gamma, where gamma=1.0 for :NSTX and 0.5 for :MAST"; default=:MAST)
    width_coefficient::Entry{T} = Entry{T}("-", "Pedestal width coefficient, C2, w_ped=C2*beta_p,ped^gamma ")
    height_coefficient::Entry{T} = Entry{T}("-", "Pedestal height coefficient, C1, beta_p,ped=C1*w_ped/IN^1/3 ")
    #ped_to_core_fraction::Entry{T} = Entry{T}("-", "Ratio of edge (@rho=0.9) to core stored energy [0.05 for L-mode, 0.3 for neg-T plasmas, missing keeps original ratio]")
    #== display and debugging parameters ==#
    do_plot::Entry{Bool} = act_common_parameters(; do_plot=false)
end

mutable struct ActorAnalyticPedestal{D,P} <: SingleAbstractActor{D,P}
    dd::IMAS.dd{D}
    #par::FUSEparameters__ActorAnalyticPedestal{P}
	par::OverrideParameters{P,FUSEparameters__ActorAnalyticPedestal{P}}
	act::ParametersAllActors{P}

    function ActorAnalyticPedestal(dd::IMAS.dd, par::FUSEparameters__ActorAnalyticPedestal; kw...)
        logging_actor_init(ActorAnalyticPedestal)
        par = par(kw...)
        return new{D,P}(dd, par)
    end
    
end

"""
    ActorAnalyticPedestal(dd::IMAS.dd, act::ParametersAllActors; kw...)

Evaluates the pedestal boundary condition (height and width)
"""
function ActorAnalyticPedestal(dd::IMAS.dd, act::ParametersAllActors; kw...)
    actor = ActorAnalyticPedestal(dd, act.ActorAnalyticPedestal; kw...)
    step(actor)
    finalize(actor)
    return actor
end

# Step function
function _step(actor::ActorAnalyticPedestal{D,P}) where {D<:Real, P<:Real}

    # Retrieve variables defined in the struct ActorAnalyticPedestal
    dd = actor.dd
    par = actor.par

    # Define physical constants
    amu = 1.660538921 * 1.0e-27 # Atomic mass unit [kg]
    e_charge = 1.60e-19 # Magnitude of the electron charge [C]
   
	#println("\nHello")

    if par.do_debug
		# Debug output goes here - WIP
    end

    # Plotting
    if par.do_plot
        # Plotting goes here - WIP
    end

    return actor
end

# Finalize function
function _finalize(actor::ActorAnalyticPedestal)
    # Finalization after computation is done in step
    # (Probably populate dd; for details of dd: https://fuse.help/dev/dd.html)
    return actor
end
