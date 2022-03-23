import DataFrames
import CSV
import Statistics

function load_dataframe(filename::String)
    return CSV.read(filename, DataFrames.DataFrame)
end

function R_squared(x, y)
    return 1 - sum((x .- y) .^ 2) / sum((x .- Statistics.mean(x)) .^ 2)
end

function mean_relative_error(x, y)
    return 100 - 100 * sum(abs.((x .- y)) ./ abs.(x)) / length(x)
end


"""
    plot_x_y_regression(dataframe::DataFrame
     x_name :: Union{String, Symbol},
     y_name :: Union{String, Symbol})
Plots x vs y regression plot in log-form on x_name and y_name in the dataframe
"""
function plot_x_y_regression(dataframe::DataFrames.DataFrame, x_name::Union{String,Symbol}, y_name::Union{String,Symbol})
    if x_name == "TAUTH"
        x_ylim = [5e-3, 1e1]
    else
        x_ylim = [minimum(abs.(dataframe[:, x_name])) / 1e1, maximum(dataframe[:, x_name]) * 1e1]
    end
    DataFrames.dropmissing!(dataframe)
    R² = round(R_squared(dataframe[:, x_name], dataframe[:, y_name]), digits = 2)
    MRE = round(mean_relative_error(dataframe[:, x_name], dataframe[:, y_name]), digits = 2)
    plot(dataframe[:, x_name], dataframe[:, y_name], seriestype = :scatter, xaxis = :log, yaxis = :log, ylim = x_ylim, xlim = x_ylim, xlabel = x_name, ylabel = y_name, label = "mean_relative_error = $MRE %")
    plot!([0.5 * x_ylim[1], 0.5 * x_ylim[2]], [2 * x_ylim[1], 2 * x_ylim[2]], linestyle = :dash, label = "+50%")
    plot!([2 * x_ylim[1], 2 * x_ylim[2]], [0.5 * x_ylim[1], 0.5 * x_ylim[2]], linestyle = :dash, label = "-50%", legend = :topleft)
    display(plot!([x_ylim[1], x_ylim[2]], [x_ylim[1], x_ylim[2]], label = nothing))
    println("R² = $(R²), mean_relative_error = $MRE)")
end