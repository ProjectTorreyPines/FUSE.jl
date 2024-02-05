#= ================ =#
#  ActorGASCCosting  #
#= ================ =#

# FUSEparameters__ActorGASCCosting must match name and types of FUSEparameters__ActorSheffieldCosting
Base.@kwdef mutable struct FUSEparameters__ActorGASCCosting{T<:Real} <: ParametersActor{T}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    _time::Float64 = NaN
    # NOTE: parameters below should reflect the parameters in FUSEparameters__ActorSheffieldCosting except for capitalize_blanket and capitalize_divertor
    construction_lead_time::Entry{T} = Entry{T}("year", "Duration of construction"; default=8.0)
    fixed_charge_rate::Entry{T} = Entry{T}("-", "Constant dollar fixed charge rate"; default=0.078)
    capitalize_blanket::Entry{Bool} = Entry{Bool}("-", "If true, include cost of 1st blanket in direct captial cost"; default=true)
    capitalize_divertor::Entry{Bool} = Entry{Bool}("-", "If true, include cost of 1st divertor in direct captial cost"; default=true)
    divertor_fluence_lifetime::Entry{T} = Entry{T}("MW*yr/m^2", "Divertor fluence over its lifetime"; default=10.0)
    blanket_fluence_lifetime::Entry{T} = Entry{T}("MW*yr/m^2", "Blanket fluence over its lifetime"; default=15.0)
end

mutable struct ActorGASCCosting{D,P} <: FacilityAbstractActor{D,P}
    dd::IMAS.dd{D}
    par::FUSEparameters__ActorGASCCosting{P}
    sheffield_actor::ActorSheffieldCosting
    function ActorGASCCosting(dd::IMAS.dd{D}, par::FUSEparameters__ActorGASCCosting{P}; kw...) where {D<:Real,P<:Real}
        logging_actor_init(ActorGASCCosting)
        par = par(kw...)
        # convert GASC costing actor parameters to Sheffield costing actor parameters
        sheffield_par = FUSEparameters__ActorSheffieldCosting{P}(par)
        for field in fieldnames(typeof(par))
            setfield!(sheffield_par, field, getfield(par, field))
        end
        return new{D,P}(dd, par, ActorSheffieldCosting(dd, sheffield_par))
    end
end

"""
    ActorGASCCosting(dd::IMAS.dd, act::ParametersAllActors; kw...)

Like ActorSheffieldCosting but by default blanket and divertor are capitalized
"""
function ActorGASCCosting(dd::IMAS.dd, act::ParametersAllActors; kw...)
    par = act.ActorGASCCosting(kw...)
    actor = ActorGASCCosting(dd, par)
    step(actor)
    finalize(actor)
    return actor
end

function _step(actor::ActorGASCCosting)
    step(actor.sheffield_actor)
    return actor
end

function _finalize(actor::ActorGASCCosting)
    finalize(actor.sheffield_actor)
    return actor
end
