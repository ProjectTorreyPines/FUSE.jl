#= ================ =#
#  ActorGASCCosting  #
#= ================ =#

# FUSEparameters__ActorGASCCosting must match name and types of FUSEparameters__ActorSheffieldCosting
mutable struct FUSEparameters__ActorGASCCosting{T} <: ParametersActor where {T<:Real}
    _parent::WeakRef
    _name::Symbol
    construction_lead_time::Entry{T}
    fixed_charge_rate::Entry{T}
    capitalize_blanket::Entry{Bool}
    capitalize_divertor::Entry{Bool}
    divertor_fluence_lifetime::Entry{T}
    blanket_fluence_lifetime::Entry{T}
end

function FUSEparameters__ActorGASCCosting{T}() where {T<:Real}
    # copy default value parameters from Sheffield costing actor
    # but set capitalize_blanket=true and capitalize_divertor=true by default
    sheffield_par = FUSEparameters__ActorSheffieldCosting{T}()
    values = []
    for field in fieldnames(typeof(sheffield_par))
        push!(values, getfield(sheffield_par, field))
        if field in [:capitalize_blanket, :capitalize_divertor]
            values[end].value = values[end].base = values[end].default = true
        end
    end
    return FUSEparameters__ActorGASCCosting{T}(values...)
end

function FUSEparameters__ActorSheffieldCosting{T}(par::FUSEparameters__ActorGASCCosting) where {T<:Real}
    sheffield_par = FUSEparameters__ActorSheffieldCosting{T}()
    for field in fieldnames(typeof(par))
        setfield!(sheffield_par, field, getfield(par, field))
    end
    return sheffield_par
end

mutable struct ActorGASCCosting{D,P} <: FacilityAbstractActor
    dd::IMAS.dd{D}
    par::FUSEparameters__ActorGASCCosting{P}
    sheffield_actor::ActorSheffieldCosting
    function ActorGASCCosting(dd::IMAS.dd{D}, par::FUSEparameters__ActorGASCCosting{P}; kw...) where {D<:Real,P<:Real}
        logging_actor_init(ActorGASCCosting)
        par = par(kw...)
        # convert GASC costing actor parameters to Sheffield costing actor parameters
        sheffield_par = FUSEparameters__ActorSheffieldCosting{P}(par)
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
