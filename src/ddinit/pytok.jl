import JSON

"""
    PyTok(filename::String; save_aggregated_json::Bool=false)

Parses PyTok output from directory with option to save aggregated json
"""
function PyTok(filename::String; save_aggregated_json::Bool=false)
    pytok = JSON.parsefile(filename)
    pytok["EXTRAS"] = Dict{String,Any}()

    dir = @show(joinpath(splitpath(filename)[1:end-1]))
    @show device = splitext(basename(filename))[1]

    for extra_filename in readdir(dir)
        extra_name = splitext(extra_filename)[1]
        if extra_name != device && match(device * r".*\.json", extra_filename) !== nothing
            data = JSON.parsefile(joinpath(dir, extra_filename))
            pytok["EXTRAS"][extra_name[length(device)+2:end]] = data
        end
    end

    if save_aggregated_json
        open(joinpath(dir, "FUSE_aggregated_$device.json"), "w") do io
            return JSON.print(io, pytok)
        end
    end

    return pytok
end
