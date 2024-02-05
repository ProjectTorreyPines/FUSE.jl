#= ================ =#
#  ActorCostingGASC  #
#= ================ =#

# FUSEparameters__ActorCostingGASC must match name and types of FUSEparameters__ActorCostingSheffield
Base.@kwdef mutable struct FUSEparameters__ActorCostingGASC{T<:Real} <: ParametersActor{T}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    _time::Float64 = NaN
    # NOTE: parameters below should reflect the parameters in FUSEparameters__ActorCostingSheffield except for capitalize_blanket and capitalize_divertor
    construction_lead_time::Entry{T} = Entry{T}("year", "Duration of construction"; default=8.0)
    fixed_charge_rate::Entry{T} = Entry{T}("-", "Constant dollar fixed charge rate"; default=0.078)
    capitalize_blanket::Entry{Bool} = Entry{Bool}("-", "If true, include cost of 1st blanket in direct captial cost"; default=true)
    capitalize_divertor::Entry{Bool} = Entry{Bool}("-", "If true, include cost of 1st divertor in direct captial cost"; default=true)
    divertor_fluence_lifetime::Entry{T} = Entry{T}("MW*yr/m^2", "Divertor fluence over its lifetime"; default=10.0)
    blanket_fluence_lifetime::Entry{T} = Entry{T}("MW*yr/m^2", "Blanket fluence over its lifetime"; default=15.0)
end

mutable struct ActorCostingGASC{D,P} <: SingleAbstractActor{D,P}
    dd::IMAS.dd{D}
    par::FUSEparameters__ActorCostingGASC{P}
    sheffield_actor::ActorCostingSheffield
    function ActorCostingGASC(dd::IMAS.dd{D}, par::FUSEparameters__ActorCostingGASC{P}; kw...) where {D<:Real,P<:Real}
        logging_actor_init(ActorCostingGASC)
        par = par(kw...)
        # convert GASC costing actor parameters to Sheffield costing actor parameters
        sheffield_par = FUSEparameters__ActorCostingSheffield{P}()
        for field in fieldnames(typeof(par))
            setfield!(sheffield_par, field, getfield(par, field))
        end
        return new{D,P}(dd, par, ActorCostingSheffield(dd, sheffield_par))
    end
end

"""
    ActorCostingGASC(dd::IMAS.dd, act::ParametersAllActors; kw...)

Like ActorCostingSheffield but by default blanket and divertor are capitalized
"""
function ActorCostingGASC(dd::IMAS.dd, act::ParametersAllActors; kw...)
    par = act.ActorCostingGASC(kw...)
    actor = ActorCostingGASC(dd, par)
    step(actor)
    finalize(actor)
    return actor
end

function _step(actor::ActorCostingGASC)
    step(actor.sheffield_actor)
    return actor
end

function _finalize(actor::ActorCostingGASC)
    finalize(actor.sheffield_actor)
    return actor
end
