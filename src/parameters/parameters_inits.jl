
include("parameters_inits_plasma.jl")

include("parameters_inits_build.jl")

# DEFAULT init structure to use
function ParametersInits(args...; kw...)
    return ParametersInitsBuild{Float64}(args...; kw...)
end

"""
    ini2json(ini::ParametersAllInits, filename::AbstractString; kw...)

Save the FUSE parameters to a JSON file with give `filename`

`kw` arguments are passed to the JSON.print function
"""
function ini2json(ini::ParametersAllInits, filename::AbstractString; kw...)
    return SimulationParameters.par2json(ini, filename; kw...)
end

"""
    json2ini(filename::AbstractString, ini::ParametersAllInits=ParametersInits())

Load the FUSE parameters from a JSON file with given `filename`
"""
function json2ini(filename::AbstractString, ini::ParametersAllInits=ParametersInits())
    return SimulationParameters.json2par(filename, ini)
end

"""
    ini2yaml(ini::ParametersAllInits, filename::AbstractString; kw...)

Save the FUSE parameters to a YAML file with given `filename`

`kw` arguments are passed to the YAML.print function
"""
function ini2yaml(ini::ParametersAllInits, filename::AbstractString; kw...)
    return SimulationParameters.par2yaml(ini, filename; kw...)
end

"""
    yaml2ini(filename::AbstractString, ini::ParametersAllInits=ParametersInits())

Load the FUSE parameters from a YAML file with given `filename`
"""
function yaml2ini(filename::AbstractString, ini::ParametersAllInits=ParametersInits())
    return SimulationParameters.yaml2par(filename, ini)
end

"""
    plot_ini(ini::ParametersAllInits; time0=global_time(ini))

Plots ini time dependent time traces including plasma boundary
"""
@recipe function plot_ini(ini::ParametersAllInits; time0=global_time(ini))
    @assert typeof(time0) <: Float64

    # count number of time-dependent parameters
    N = 0
    for par in SimulationParameters.leaves(ini)
        if typeof(par.value) <: Function
            N += 1
        end
    end
    layout := @layout [N + 1]

    time_bkp = ini.time.simulation_start
    try
        ini.time.simulation_start = time0

        # plot equilibrium including x-points
        nx = n_xpoints(ini.equilibrium.xpoints)
        mxh = IMAS.MXH(ini)
        mxhb = fitMXHboundary(mxh, nx)
        wr = wall_radii(ini.build.layers, mxh.R0, mxh.ϵ * mxh.R0, ini.build.plasma_gap)
        @series begin
            label := ""
            seriestype := :vline
            subplot := 1
            colour := :black
            [wr.r_hfs, wr.r_lfs]
        end
        @series begin
            label := ""
            subplot := 1
            aspectratio := :equal
            xlim := (0, mxh.R0 * 2)
            mxhb
        end

        # plot time dependent parameters
        k = 1
        for par in SimulationParameters.leaves(ini)
            if typeof(par.value) <: Function
                k += 1
                @series begin
                    label := ""
                    subplot := k
                    time0 := time0
                    par
                end
            end
        end

    finally
        ini.time.simulation_start = time_bkp
    end
end
