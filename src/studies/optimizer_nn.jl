import Flux
import PrettyTables

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

    idx_base = 10002
    results = complete_sensitivity_analysis(study, idx_base)
    # store the results in study 
    return study
end 
# instead of idx_base, user should be able to select the design with params closest to
# some chosen value then vary around that 

function constrained_sensitivity2(study::StudyOptimizerNN, input_label::String, idx_base::Int)
    sty = study.sty

    interpolator_nn_inputs  = study.interpolator_nn_inputs
    interpolator_nn_outputs = study.interpolator_nn_outputs
    classifier_nn_inputs    = study.classifier_nn_inputs

    percent_diff = sty.actuator_variation
    n_test       = sty.n_test_points

    labelsx = ["R0","B0","ip","δ","fGW","Pec","Pic","zeff","kappa","minor radius"] # this eventually needs to be a variable 
    
    j = findfirst(==(input_label), labelsx)
    j === nothing && error("input_label '$input_label' not found in labelsx")

    x, x_min, x_max         = minmax_normalize(interpolator_nn_inputs)
    y, y_min, y_max         = minmax_normalize(interpolator_nn_outputs)
    cl_x, cl_x_min, cl_x_max = minmax_normalize(classifier_nn_inputs)

    row = view(interpolator_nn_inputs, j, :)
    mean_val = sum(row) / length(row)

    test_vals = collect(range(mean_val * (1 - percent_diff), mean_val * (1 + percent_diff), n_test))

    base_input = interpolator_nn_inputs[:, idx_base]

    outputs = Vector{Vector{Float32}}()
    classifier_inputs = Vector{Vector{Float32}}()

    for v in test_vals
        inp = copy(base_input)
        inp[j] = v

        input_vec   = vec(Float32.(inp))
        input_norm  = minmax_normalize(input_vec, x_min, x_max)
        output_norm = study.model_interpolator(input_norm)
        output      = vec(minmax_unnormalize(output_norm, y_min, y_max))
        push!(outputs, output)

        # append net electric power to inputs for classifier
        push!(classifier_inputs, vcat(input_vec, output[3]))
    end

    classifier_inputs_mat = reduce(hcat, classifier_inputs)'

    classifier_outputs = Vector{Int}()
    for i in 1:length(test_vals)
        norm_cl_in = minmax_normalize(classifier_inputs_mat[i, :], cl_x_min, cl_x_max)
        cl_out = only(round.(study.model_classifier(norm_cl_in)))
        push!(classifier_outputs, cl_out)
    end

    outputs_mat = hcat(outputs...)'
    return (test_vals, outputs_mat, classifier_outputs)
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

function calculate_output_variations(interpolator_outputs)
    # Get the range (max - min) for each output column
    q95_range = maximum(interpolator_outputs[:, 1]) - minimum(interpolator_outputs[:, 1])
    cost_range = maximum(interpolator_outputs[:, 2]) - minimum(interpolator_outputs[:, 2])
    pelectric_range = maximum(interpolator_outputs[:, 3]) - minimum(interpolator_outputs[:, 3])
    
    # Calculate baseline as the first value (or you could use middle value)
    q95_baseline = interpolator_outputs[1, 1]
    cost_baseline = interpolator_outputs[1, 2]  
    pelectric_baseline = interpolator_outputs[1, 3]
    
    q95_percent_var = (q95_range / abs(q95_baseline)) * 100
    cost_percent_var = (cost_range / abs(cost_baseline)) * 100
    pelectric_percent_var = (pelectric_range / abs(pelectric_baseline)) * 100
    
    return q95_percent_var, cost_percent_var, pelectric_percent_var
end

function run_sensitivity_analysis(study, input_params, idx_base)
    results = []
    
    for param in input_params
        println("Running sensitivity analysis for: $param")
        
        # Run your existing sensitivity function
        test_vals, interpolator_outputs, classifier_outputs = constrained_sensitivity2(study, param, idx_base)
        
        # Calculate variations
        q95_var, cost_var, pelectric_var = calculate_output_variations(interpolator_outputs)
        
        # Calculate total sensitivity (sum of absolute variations)
        total_sensitivity = abs(q95_var) + abs(cost_var) + abs(pelectric_var)
        
        # Store results
        push!(results, (
            parameter = param,
            q95_variation = q95_var,
            cost_variation = cost_var,
            pelectric_variation = pelectric_var,
            total_sensitivity = total_sensitivity,
            interpolator_outputs = interpolator_outputs  # Store for reference if needed
        ))
    end
    
    # Sort by total sensitivity (descending)
    sort!(results, by = x -> x.total_sensitivity, rev = true)
    
    return results
end

function display_sensitivity_results(results)
    n_params = length(results)
    data = Matrix{Any}(undef, n_params, 5)
    
    for (i, result) in enumerate(results)
        data[i, 1] = result.parameter
        data[i, 2] = round(result.q95_variation, digits=2)
        data[i, 3] = round(result.cost_variation, digits=2)
        data[i, 4] = round(result.pelectric_variation, digits=2)
        data[i, 5] = round(result.total_sensitivity, digits=2)
    end
    
    header = ["Parameter", "q95 Variation (%)", "Cost Variation (%)", "P_electric Variation (%)", "Total Sensitivity (%)"]
    
    pretty_table(data; 
                header = header,
                header_crayon = crayon"bold blue",
                alignment = [:l, :r, :r, :r, :r],
                formatters = ft_printf("%.2f", [2,3,4,5]))
    
    return data
end

# Main function to run everything
function complete_sensitivity_analysis(study, idx_base)
    # Define your input parameters
    input_params = ["R0", "B0", "ip", "δ", "fGW", "Pec", "Pic", "zeff", "kappa", "minor radius"] # should be sty.input_params instead 
    
    println("="^80)
    println("RUNNING SENSITIVITY ANALYSIS")
    println("="^80)
    
    # Run analysis for all parameters
    results = run_sensitivity_analysis(study, input_params, idx_base)
    
    println("\n" * "="^80)
    println("SENSITIVITY ANALYSIS RESULTS")
    println("Parameters ranked by total output sensitivity")
    println("="^80)
    
    # Display results
    display_sensitivity_results(results)
    
    return results
end