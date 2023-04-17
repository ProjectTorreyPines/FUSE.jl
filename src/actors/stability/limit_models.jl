##### QUESTIONS #####
# 1) Putting the models in a constant Dictionary of functions
#    > should match the constant dataset
#    > Is there a better way to link them? 
#    > something more automatic?
# 2) Add actual links to references in the docstring?
#    > How to do this? 
#    > Can we clone the docstring to the id.description?
# 3) Best way to implement collections?
#    > Just run all the model functions?
#    >>> This wouldn't do the `dd` correctly
#    > Run new actors for each? 
# 4) Best way to handle inputs to keep compatability with function Dict?
#    > new mutable struct? function params? 

##### MODEL COLLECTIONS #####
function force_fail(dd::IMAS.dd, par::FUSEparameters__ActorStabilityLimits, model::IMAS.stability__model) 
    error(raw"¯\_(ツ)_/¯")
end

function default(dd::IMAS.dd, par::FUSEparameters__ActorStabilityLimits, model::IMAS.stability__model)
    #error("Not implemented")
    logging(Logging.Error, :actors, "ActorStabilityLimits: default model may not be sufficent check on stability.")
end

function collection_cleared()
    all_cleared = Vector{Int64}([])
    for (time_index, time) in enumerate(stability.time)
        cleared = []
        for model in stability.model
            append!(cleared, model.cleared[time_index])
        end
        append!(all_cleared,prod(cleared))
    end
    stability.all_cleared = all_cleared
end


##### BETA LIMIT MODELS #####

function beta_limits(dd::IMAS.dd, par::FUSEparameters__ActorStabilityLimits, model::IMAS.stability__model)
    model.identifier.name = "Beta Limits"
    model.identifier.description = "Checks all beta limits models"
    
    println(model)

    par.model_ids = [:troyon_1984, :troyon_1985]
    step(ActorStabilityLimits(dd, par))

    


end

"""
    beta_troyon_1984(dd::IMAS.dd, par::FUSEparameters__ActorStabilityLimits, model::IMAS.stability__model)

Limit in normalized beta using Troyon scaling
Model Formulation: Beta_{N} < 3.5
Citation:  F Troyon et al 1984 Plasma Phys. Control. Fusion 26 209
"""
#funciton model_101(dd::IMAS.dd, par::FUSEparameters__ActorStabilityLimits, model::IMAS.stability__model)
function beta_troyon_1984(dd::IMAS.dd, par::FUSEparameters__ActorStabilityLimits, model::IMAS.stability__model)
    model.identifier.name = "Troyon 1984"
    model.identifier.description = "Beta_{N} < 3.5"
 
    beta_normal = dd.equilibrium.time_slice[].global_quantities.beta_normal
    model_value = beta_normal 
    target_value = 3.5

    @ddtime(model.fraction = model_value / target_value)
    # RHS might need to be a vector?
    #     @ddtime(lim.fraction = [model_value / target_value])
end

"""
    beta_troyon1985(dd::IMAS.dd, par::FUSEparameters__ActorStabilityLimits, model::IMAS.stability__model)

Limit in normalized beta using classical scaling using combined kink and ballooning stability
Model Formulation: Beta_{N} < 2.8
Citation: 
"""
function beta_troyon_1985(dd::IMAS.dd, par::FUSEparameters__ActorStabilityLimits, model::IMAS.stability__model)
    model.identifier.name = "Troyon 1985"
    model.identifier.description = "Beta_{N} < 2.8"

    beta_normal = dd.equilibrium.time_slice[].global_quantities.beta_normal
    model_value = beta_normal 
    target_value = 2.8

    @ddtime(model.fraction = model_value / target_value)
end


const limit_models = Dict(
    0 => force_fail,
    1 => default,
    100 => beta_limits,
    101 => beta_troyon_1984,
    102 => beta_troyon_1985
)


### Everything below has par::LimFUSE instead of just FUSE?



# """
#     beta_tuda1985(dd::IMAS.dd, par::FUSEparameters__ActorStabilityLimits, lim::IMAS.stability__model)

# Limit in beta_normal using classical scaling using only kink stability
# Model Formulation: Beta_{N} < 3.2
# Citation: 
# """
# function beta_tuda1985(dd::IMAS.dd, par::FUSEparameters__ActorStabilityLimits, lim::IMAS.stability__model)
#     lim.model.name = "Standard::KinkOnly"
#     lim.model.formula = "Beta_{N} < 3.2"

#     beta_normal = dd.equilibrium.time_slice[].global_quantities.beta_normal
#     model_value = beta_normal 
#     target_value = 3.2

#     lim.model.fraction = model_value / target_value
# end

# """
# beta_bernard1983(dd::IMAS.dd, par::FUSEparameters__ActorStabilityLimits, lim::IMAS.stability__model)

# Limit in normalized beta using classical scaling using only ballooning stability
# Model Formulation: Beta_{N} < 2.8
# Citation: 
# """
# function beta_bernard1983(dd::IMAS.dd, par::FUSEparameters__ActorStabilityLimits, lim::IMAS.stability__model)
#     lim.model.name = "Standard::BallooningOnly"
#     lim.model.formula = "Beta_{N} < 4.4"

#     beta_normal = dd.equilibrium.time_slice[].global_quantities.beta_normal
#     model_value = beta_normal 
#     target_value = 4.4

#     lim.model.fraction = model_value / target_value
# end



# """
#     beta_betali_a(dd::IMAS.dd, par::FUSEparameters__ActorStabilityLimits, lim::IMAS.stability__model)

# Modern limit in normlaized beta normalized by plasma inductance
# Model Formulation: beta_{N} / li < C_{beta}
# Citation: 
# """
# function beta_betali(dd::IMAS.dd, par::FUSEparameters__ActorStabilityLimits, lim::IMAS.stability__model)
#     lim.model.name = "BetaLi::Troyon"
#     lim.model.formula = "Beta_{N} / Li < 4.0"

#     beta_normal = dd.equilibrium.time_slice[].global_quantities.beta_normal
#     plasma_inductance =  dd.equilibrium.time_slice[].global_quantities.li_3
    
#     model_value = beta_normal / plasma_inductance
#     target_value = 4.4

#     lim.model.fraction = model_value / target_value
# end


# ##### CURRENT LIMIT MODELS #####

# """
#     current_standard(dd::IMAS.dd)

# Standard limit in edge current via the safety factor
# Model Formulation: q95 < C
# Citation: 
# """
# function current_standard(dd::IMAS.dd)
    
#     q95 = dd.equilibrium.time_slice[].global_quantities.q_95 
    
#     model_value = 1/abs(q95)

#     return model_value
# end

# """
#     current_standard_a(dd::IMAS.dd, par::FUSEparameters__ActorStabilityLimits, lim::IMAS.stability__model)

# Standard limit in edge current via the safety factor
# Model Formulation: q95 < 2
# Citation: 
# """
# function current_standard_a(dd::IMAS.dd, par::FUSEparameters__ActorStabilityLimits, lim::IMAS.stability__model)
#     lim.model.name = "Standard::q95"
#     lim.model.formula = "q_95 > 2.0"

#     model_value = current_standard(dd)
#     target_value = 0.5

#     lim.model.fraction = model_value / target_value
# end


# ##### DENSITY LIMIT MODELS #####


# """
#     density_standard(dd::IMAS.dd)

# Standard limit in density
# Model Formulation: f_GW < 1
# Citation: 
# """
# function density_standard(dd::IMAS.dd)
 
#     eqt = dd.equilibrium.time_slice[]
#     cp1d = dd.core_profiles.profiles_1d[]

#     model_value = IMAS.greenwald_fraction(eqt, cp1d)

#     return model_value
# end

# """
#     density_standard_a(dd::IMAS.dd, par::FUSEparameters__ActorStabilityLimits, lim::IMAS.stability__model)

# Standard limit in density using IMAS greenwald fraction
# Model Formulation: f_{GW,IMAS} < 1.0
# Citation: 
# """
# function density_standard_a(dd::IMAS.dd, par::FUSEparameters__ActorStabilityLimits, lim::IMAS.stability__model)
#     lim.model.name = "Standard::IMAS_Greenwald"
#     lim.model.formula = "IMAS.greenwald_fraction < 1.0"

#     model_value = density_standard(dd)
#     target_value = 1.0

#     lim.model.fraction = model_value / target_value
# end


# function model_list(limits::Vector{Symbol})
#     model_list = []
#     :BetaLimit in limits && append!(model_list,beta_model_list())
#     :CurrentLimit in limits && append!(model_list,current_model_list())
#     :DensityLimit in limits && append!(model_list,density_model_list())
#     return model_list
# end    
# function beta_model_list()
#     model_list = [
#         :None,
#         :Troyon1983,
#         :Troyon1985,
#         :Tuda1985,
#         :Bernard1983,
#         :BetaLi
#     ]
#     return model_list
# end
# function current_model_list()
#     model_list = [
#         :None,
#         :Standard__q95
#     ]
#     return model_list
# end
# function density_model_list()
#     model_list = [
#         :None,
#         :IMAS_Greenwald
#     ]
#     return model_list
# end
