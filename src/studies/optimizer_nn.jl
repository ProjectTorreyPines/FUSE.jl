import Flux

#= ================= =#
#   StudyOptimizerNN  #
#= ================= =#

"""
    study_parameters(::Type{Val{:OptimizerNN}})::Tuple{FUSEparameters__ParametersStudyOptimizerNN,ParametersAllActors}

Uses NN surrogate of an optimization database for uncertainty quantification
"""
function study_parameters(::Type{Val{:OptimizerNN}})::Tuple{FUSEparameters__ParametersStudyOptimizerNN,ParametersAllActors}
    sty = FUSEparameters__ParametersStudyOptimizerNN{Real}()
    act = ParametersActors()

    set_new_base!(sty)
    set_new_base!(act)

    return sty, act
end

Base.@kwdef mutable struct FUSEparameters__ParametersStudyOptimizerNN{T<:Real} <: ParametersStudy{T}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :StudyMultiObjectiveOptimizer
    preloaded_interpolator_model::Entry{Union{Missing,<:Flux.Chain}} = Entry{Union{Missing,<:Flux.Chain}}("-", "Pre-trained interpolator model to load"; default = missing)
    preloaded_classifier_model::Entry{Union{Missing,<:Flux.Chain}} = Entry{Union{Missing,<:Flux.Chain}}("-", "Pre-trained classifier model to load"; default = missing)
    actuator_list::Entry{Union{Vector{Symbol},Missing}} = Entry{Union{Vector{Symbol},Missing}}("-", "List of all actuators in original optimization run"; default = missing)
    actuator_variation::Entry{Float64} = Entry{Float64}("-", "Fraction by which to vary each of the actuators"; default = 0.1)
    n_test_points::Entry{Int} = Entry{Int}("-", "Number of test points for sensitivty analysis"; default = 10)
end

mutable struct StudyOptimizerNN{T<:Real} <: AbstractStudy
    sty::OverrideParameters{T,FUSEparameters__ParametersStudyOptimizerNN{T}}
    ini::ParametersAllInits
    act::ParametersAllActors
    constraint_functions::Vector{IMAS.ConstraintFunction}
    objective_functions::Vector{IMAS.ObjectiveFunction}
    interpolator_nn_inputs::Union{Matrix{Float64},Missing}
    interpolator_nn_outputs::Union{Matrix{Float64},Missing}
    classifier_nn_inputs::Union{Matrix{Float64},Missing}
    classifier_nn_outputs::Union{Matrix{Float64},Missing}
    dataframe::Union{DataFrame,Missing}
    model_interpolator::Union{<:Flux.Chain,Missing}
    model_classifier::Union{<:Flux.Chain,Missing}
end

function StudyOptimizerNN(
    sty::ParametersStudy,
    ini::ParametersAllInits,
    act::ParametersAllActors,
    constraint_functions::Vector{IMAS.ConstraintFunction},
    objective_functions::Vector{IMAS.ObjectiveFunction},
    interpolator_nn_inputs::Union{Matrix{Float64},Missing},
    interpolator_nn_outputs::Union{Matrix{Float64},Missing},
    classifier_nn_inputs::Union{Matrix{Float64},Missing},
    classifier_nn_outputs::Union{Matrix{Float64},Missing},
    dataframe::Union{DataFrame,Missing},
    model_interpolator::Union{<:Flux.Chain,Missing},
    model_classifier::Union{<:Flux.Chain,Missing};
    kw...
)
    sty = OverrideParameters(sty; kw...)
    study = StudyOptimizerNN(sty, ini, act, constraint_functions, objective_functions, missing, missing, missing, missing, dataframe, missing, missing)
    return setup(study)
end

function _setup(study::StudyOptimizerNN)
    sty = study.sty
    return study
end

function _run(study::StudyOptimizerNN)
    sty = study.sty
    df = study.dataframe

    if !ismissing(sty.preloaded_classifier_model) && !ismissing(sty.preloaded_interpolator_model)
        study.model_interpolator = sty.preloaded_interpolator_model
        study.model_classifier = sty.preloaded_classifier_model
    else 
        # do training of the models and maybe hyperparamemter tuning 
    end 

    # this originally was the filtered df, need to filter for the constraints 
    inputs = hcat(df.R0, df.B0, df.ip, df.δ, df.fGW, df.Pec, df.Pic, df."<zeff>", df.κ, df.a)
    outputs = hcat(df.q95, df.capital_cost, df.Pelectric_net)

    valid_rows = map(is_valid_row, eachrow(inputs)) .& map(is_valid_row, eachrow(outputs)) 
    
    study.interpolator_nn_inputs = permutedims(inputs[valid_rows, :])
    study.interpolator_nn_outputs = permutedims(outputs[valid_rows, :])

    ########################################################################

    cl_inputs = hcat(df.R0, df.B0, df.ip, df.δ, df.fGW, df.Pec, df.Pic, df."<zeff>", df.κ, df.a, df.Pelectric_net)

    valid_rows = map(is_valid_row, eachrow(cl_inputs))

    study.classifier_nn_inputs = permutedims(cl_inputs[valid_rows, :])

    # for actuator in sty.actuator_list
    #     # do the constrained sensitivity for each actuator 
    idx_base = 10002
    test_vals, interpolator_outputs, classifier_outputs = constrained_sensitivity(study, "R0", idx_base)
    @show interpolator_outputs 
    # end

    # check sensitivity of outputs to inputs 
    # show a summary of how e.g. a default 10% variation in inputs affects each objective output
    return study
end 


# instead of idx_base, user should be able to select the design with params closest to
# some chosen value 
# then vary around that 
function constrained_sensitivity(study::StudyOptimizerNN, input_label::String, idx_base::Int)
    sty = study.sty 

    interpolator_nn_inputs = study.interpolator_nn_inputs
    interpolator_nn_outputs = study.interpolator_nn_outputs
    classifier_nn_inputs = study.classifier_nn_inputs

    percent_diff = sty.actuator_variation
    n_test = sty.n_test_points

    labelsx = ["R0", "B0", "ip", "δ", "fGW", "Pec", "Pic", "zeff", "kappa", "minor radius"]
    labelsy = ["q95", "capital cost", "Pelectric net"]

    cl_labelsx = ["R0", "B0", "ip", "δ", "fGW", "Pec", "Pic", "zeff", "kappa", "minor radius", "net electric power"]

    x, x_min, x_max = minmax_normalize(interpolator_nn_inputs)
    y, y_min, y_max = minmax_normalize(interpolator_nn_outputs)

    cl_x, cl_x_min, cl_x_max = minmax_normalize(classifier_nn_inputs)

    input = findfirst(isequal(input_label), labelsx)

    mean_val = sum(interpolator_nn_inputs[input,:]) / length(interpolator_nn_inputs[input,:])
    test_vals = collect(range(mean_val - (mean_val*percent_diff), mean_val + (mean_val*percent_diff), n_test))

    base_input = interpolator_nn_inputs[:,idx_base]

    outputs = []
    classifier_inputs = []
    inp = deepcopy(base_input)
    for i in 1:length(test_vals)
        
        inp[1] = test_vals[i]
    
        input_vec = vec(Float32.(inp))
        input_norm = minmax_normalize(input_vec, x_min, x_max)
        output_norm = study.model_interpolator(input_norm)
        output = vec(minmax_unnormalize(output_norm, y_min, y_max))
        push!(outputs, output)

        classifier_input = vcat(input_vec, output[3]) # tack net electric power onto inputs for classifier
        push!(classifier_inputs, classifier_input)  
    end

    classifier_inputs = reduce(hcat, classifier_inputs)'
    
    classifier_outputs = []
    for i in 1:length(test_vals)
        norm_classifier_input = minmax_normalize(classifier_inputs[i,:], cl_x_min, cl_x_max)
        classifier_output = only(round.(study.model_classifier(norm_classifier_input)))
        push!(classifier_outputs, classifier_output)
    end

    outputs = hcat(outputs...)'
    return (test_vals, outputs, classifier_outputs)
end

function minmax_normalize(x)
    min_x = minimum(x, dims=2)
    max_x = maximum(x, dims=2)
    x_norm = (x .- min_x) ./ (max_x .- min_x .+ eps())
    return x_norm, min_x, max_x
end

function minmax_normalize(x, min_x, max_x)
    x_norm = (x .- min_x) ./ (max_x .- min_x .+ eps())
    return x_norm
end

function minmax_unnormalize(x_norm, min_x, max_x)
    return x_norm .* (max_x .- min_x .- eps()) .+ min_x
end

function is_valid_row(row)
    all(x -> isfinite(x) && !ismissing(x), row)
end


