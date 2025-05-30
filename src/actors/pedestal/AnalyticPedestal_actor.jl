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
    width_coefficient::Entry{T} = Entry{T}("-", "Pedestal width coefficient, C2, w_ped=C2*beta_p,ped^gamma"; default=0.0)
    height_coefficient::Entry{T} = Entry{T}("-", "Pedestal height coefficient, C1, beta_p,ped=C1*w_ped/IN^1/3"; default=0.0)
    #== display and debugging parameters ==#
    do_plot::Entry{Bool} = act_common_parameters(; do_plot=false)
end

mutable struct ActorAnalyticPedestal{D,P} <: SingleAbstractActor{D,P}
    dd::IMAS.dd{D}
	par::OverrideParameters{P,FUSEparameters__ActorAnalyticPedestal{P}}
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

function ActorAnalyticPedestal(dd::IMAS.dd, par::FUSEparameters__ActorAnalyticPedestal; kw...)
    logging_actor_init(ActorAnalyticPedestal)
	par = OverrideParameters(par; kw...)
    return ActorAnalyticPedestal(dd, par)
end

# Step function
function _step(actor::ActorAnalyticPedestal{D,P}) where {D<:Real, P<:Real}
    dd = actor.dd
    par = actor.par

	eqt = dd.equilibrium.time_slice[]

	IN = (eqt.global_quantities.ip/1e6)/(eqt.boundary.minor_radius*eqt.global_quantities.magnetic_axis.b_field_tor) # Normalised plasma current IN=IP[MA]/(a[m]*BT[T])

	println(par.height_coefficient)

	# NSTX like width - w_ped~beta_p,ped^0.5
	if par.model == :MAST
        if par.height_coefficient==0.0
            par.height_coefficient=4.0
        end
        if par.width_coefficient==0.0
            par.width_coefficient=0.11
        end
        w_ped = par.height_coefficient^0.8*par.width_coefficient^1.6/IN^0.25 
        p_ped = 32*par.height_coefficient^1.6*par.width_coefficient^1.2*IN^1.5*eqt.global_quantities.magnetic_axis.b_field_tor^2/(1+eqt.boundary.elongation^2) 
    # NSTX like width - w_ped~beta_p,ped^1.0
    elseif par.model == :NSTX
        if par.height_coefficient==0.0
            par.height_coefficient=2.0
        end
        if par.width_coefficient==0.0
            par.width_coefficient=0.4
        end
        w_ped = par.height_coefficient^4.0*par.width_coefficient^4.0/IN^(4/3) 
        p_ped = 32*par.height_coefficient^4.0*par.width_coefficient^3.0*IN^(2/3)*eqt.global_quantities.magnetic_axis.b_field_tor^2/(1+eqt.boundary.elongation^2)
    else
        error("Undefined model! Model should be one of :MAST, :NSTX")
    end

	print("\nmodel=",par.model)
	print("\nw_ped=",w_ped)
	print("\np_ped=",p_ped,"\n")

    #if par.do_debug
		# Debug output goes here - WIP
    #end

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
