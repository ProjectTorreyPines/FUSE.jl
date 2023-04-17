# #= ==================== =#
# #  ActorStabilityLimits  #
# #= ==================== =#

# Base.@kwdef mutable struct FUSEparameters__ActorStabilityLimits{T} <: ParametersActor where {T<:Real}
#     _parent::WeakRef = WeakRef(nothing)
#     _name::Symbol = :not_set
#    model_ids::Entry{Vector{Symbol}} = Entry(Vector{Symbol}, "-", "Models for the limit calculation"; default=[:force_fail])
# end

# mutable struct ActorStabilityLimits <: PlasmaAbstractActor
#     dd::IMAS.dd
#     par::FUSEparameters__ActorStabilityLimits
#     models::Vector{IMAS.stability__model}
# end

# """
# ActorStabilityLimits(dd::IMAS.dd, act::ParametersAllActors; kw...)

# Runs all the limit actors. 
# """
# function ActorStabilityLimits(dd::IMAS.dd, act::ParametersAllActors; kw...)
#     par = act.ActorStabilityLimits(kw...)
#     actor = ActorStabilityLimits(dd, par)
#     step(actor)
#     finalize(actor)
#     return actor
# end

# function ActorStabilityLimits(dd::IMAS.dd, par::FUSEparameters__ActorStabilityLimits; kw...)
#     logging_actor_init(ActorStabilityLimits)
#     par = par(kw...)

#     @ddtime(dd.stability.time = dd.global_time)

#     number_of_models = length(par.model_ids)
#     models = Vector{IMAS.stability__model}(undef, number_of_models)
#     for (idx, model_id) in enumerate(par.model_ids)
#         model_index = IMAS.name_2_index(dd.stability.model)[model_id]
#         model = resize!(dd.stability.model, "identifier.index" => model_index)
#         models[idx] = model        
#     end    
#     return ActorStabilityLimits(dd, par, models)

# end

# # function ActorStabilityLimits(dd::IMAS.dd, par::FUSEparameters__ActorStabilityLimits; kw...)
# #     logging_actor_init(ActorStabilityLimits)
# #     par = par(kw...)

# #     @ddtime(dd.stability.time = dd.global_time)

# #     model_index = IMAS.name_2_index(dd.stability.model)[par.model_id]
# #     model = resize!(dd.stability.model, "identifier.index" => model_index)

# #     return ActorStabilityLimits(dd, par, model, model_index)
# # end

# """
#     step(actor::ActorStabilityLimits)

# Runs through the selected stability actor's step
# """
# function _step(actor::ActorStabilityLimits)
#     dd = actor.dd
#     par = actor.par

#     println("Step: ", actor.models)

#     for model in actor.models
#         limit_models[model.identifier.index](dd, par, model)
#     end

#     return actor
# end

# """
#     finalize(actor::ActorStabilityLimits)

#     Finalizes the selected stability actor
# """
# function _finalize(actor::ActorStabilityLimits)
#     # put a sort here by index
#     return actor
# end





# #######################################################


