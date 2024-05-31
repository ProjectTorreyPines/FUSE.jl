#= ============= =#
#  ActorWPED  #
#= ============= =#
Base.@kwdef mutable struct FUSEparameters__ActorWPED{T<:Real} <: ParametersActorPlasma{T}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    _time::Float64 = NaN
    #== actor parameters ==#
    rho_nml::Entry{T} = Entry{T}("-", "Defines rho at which the no man's land region starts")
    rho_ped::Entry{T} = Entry{T}("-", "Defines rho at which the pedestal region starts") # rho_nml < rho_ped
    T_ratio_pedestal::Entry{T} =
        Entry{T}("-", "Ratio of ion to electron temperatures (or rho at which to sample for that ratio, if negative; or rho_nml-(rho_ped-rho_nml) if 0.0)"; default=1.0)
    stored_energy_fraction::Entry{T} = Entry{T}("-", "ratio of pedestal stored energy to core stored energy"; default = 0.08)
end

mutable struct ActorWPED{D,P} <: SingleAbstractActor{D,P}
    dd::IMAS.dd{D}
    par::FUSEparameters__ActorWPED{P}
    result::Nothing
end

"""
    ActorWPED(dd::IMAS.dd, act::ParametersAllActors; kw...)

Evaluates the pedestal boundary condition (height and width)
"""
function ActorWPED(dd::IMAS.dd, act::ParametersAllActors; kw...)
    actor = ActorWPED(dd, act.ActorWPED; kw...)
    step(actor)
    finalize(actor)
    return actor
end


function ActorWPED(dd::IMAS.dd, par::FUSEparameters__ActorWPED; kw...)
    logging_actor_init(ActorWPED)
    par = par(kw...)
    return ActorWPED(dd, par, nothing)
end

"""
    _step(actor::ActorWPED)

Runs pedestal actor to evaluate pedestal width and height
"""
function _step(actor::ActorWPED{D,P}) where {D<:Real,P<:Real}
    return actor
end

"""
    _finalize(actor::ActorWPED)

Writes results to dd.summary.local.pedestal and possibly updates core_profiles
"""
function _finalize(actor::ActorWPED)
    dd = actor.dd
    par = actor.par

    return actor
end