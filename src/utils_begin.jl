import ForwardDiff
import Distributed
import ClusterManagers
import TimerOutputs

# ====== #
# Timing #
# ====== #
const timer = TimerOutputs.TimerOutput()

function TimerOutputs.reset_timer!(to::TimerOutputs.TimerOutput, section::String)
    pop!(to.inner_timers, section, nothing)
    to.prev_timer_label = ""
    to.prev_timer = nothing
end

function TimerOutputs.reset_timer!(section::String)
    pop!(timer.inner_timers, section, nothing)
    timer.prev_timer_label = ""
    timer.prev_timer = nothing
end

# ==== #
# Math #
# ==== #
function unwrap(v, inplace=false)
    unwrapped = inplace ? v : copy(v)
    for i = 2:length(v)
        while (unwrapped[i] - unwrapped[i-1] >= pi)
            unwrapped[i] -= 2pi
        end
        while (unwrapped[i] - unwrapped[i-1] <= -pi)
            unwrapped[i] += 2pi
        end
    end
    return unwrapped
end

function IMAS.force_float(x::ForwardDiff.Dual)
    ## we purposly do not do it recursively since generally
    ## ForwardDiff.Dual of ForwardDiff.Dual is an indication of someghing going wrong
    # return force_float(x.value)
    return x.value
end

"""
    same_length_vectors(args...)

Returns scalars and vectors as vectors of the same lengths
For example:

    same_length_vectors(1, [2], [3,3,6], [4,4,4,4,4,4])
    
    4-element Vector{Vector{Int64}}:
    [1, 1, 1, 1, 1, 1]
    [2, 2, 2, 2, 2, 2]
    [3, 3, 6, 3, 3, 6]
    [4, 4, 4, 4, 4, 4]
"""
function same_length_vectors(args...)
    function length2(x)
        if ismissing(x)
            return 1
        else
            return length(x)
        end
    end
    n = maximum(map(length2, args))
    args = collect(map(x -> isa(x, Vector) ? x : [x], args))
    args = map(x -> vcat([x for k = 1:n]...)[1:n], args)
end

"""
    mirror_bound(x, l, u)

Return tuple with value of x bounded between l and u
The bounding is done by mirroring the value at the bound limits.
"""
function mirror_bound(x, l, u)
    d = (u - l) / 2.0
    c = (u + l) / 2.0
    x0 = (x .- c) / d
    while abs(x0) > 1.0
        if x0 < 1.0
            x0 = -2.0 - x0
        else
            x0 = 2.0 - x0
        end
    end
    return x0 * d + c
end

"""
    mirror_bound(x, l, u)

Return tuple with value of x bounded between l and u and error (cost) for going out of bounds
The bounding is done by mirroring the value at the bound limits.
"""
function mirror_bound_w_cost(x, l, u)
    y = mirror_bound.(x, l, u)
    return y, abs.((x .- l) .* (x .< l) .+ (x .- u) .* (x .> u))
end

# =========== #
# Convex Hull #
# =========== #
import VacuumFields: convex_hull

# ======== #
# TraceCAD #
# ======== #
import FileIO, JpegTurbo

"""
    dd_build_layers_to_ini(dd::IMAS.dd)

Utility function to convert layers in dd.build to layers in `ini.build.layers = layers = OrderedCollections.OrderedDict()`
"""
function dd_build_layers_to_ini(dd::IMAS.dd)
    for layer in dd.build.layer
        name = replace(layer.name, " " => "_")
        println("layers[:$name] = $(layer.thickness)")
    end
end

struct TraceCAD
    name::Symbol
    x_length::Real
    x_offset::Real
    y_offset::Real
end

function TraceCAD(device::Symbol)
    return TraceCAD(Val{device})
end

@recipe function plot_trace_cad(cad::TraceCAD)
    img = FileIO.load(joinpath(dirname(@__DIR__), "cases", "$(cad.name).jpg"))
    img .= img[end:-1:1, 1:end]
    x = LinRange(0, cad.x_length, size(img)[1]) .+ cad.x_offset
    y = LinRange(-0.5, 0.5, size(img)[2]) .* (size(img)[1] / size(img)[2]) * (x[end] - x[1]) .- cad.y_offset
    @series begin
        flip --> false
        x, y, img
    end
end

# ==== #
# fuse #
# ==== #
function fuse()
    return """
███████╗██╗   ██╗███████╗███████╗
██╔════╝██║   ██║██╔════╝██╔════╝
█████╗  ██║   ██║███████╗█████╗  
██╔══╝  ██║   ██║╚════██║██╔══╝  
██║     ╚██████╔╝███████║███████╗
╚═╝      ╚═════╝ ╚══════╝╚══════╝
"""
end

# ======== #
# parallel #
# ======== #
"""
    parallel_environment(cluster::String="localhost", nprocs_max::Integer=0, kw...) 

Start multiprocessing environment

kw arguments are passed to the Distributed.addprocs
"""
function parallel_environment(cluster::String="localhost", nprocs_max::Integer=0, kw...)
    if cluster == "saga"
        if gethostname() == "saga.cluster"
            nodes = 4
            np = 30 * nodes
            if nprocs_max > 0
                np = min(np, nprocs_max)
            end
            ENV["JULIA_WORKER_TIMEOUT"] = "180"
            if nprocs() < np
                Distributed.addprocs(ClusterManagers.SlurmManager(np - nprocs()), exclusive="", topology=:master_worker, kw...)
            end
            println("Working with $(nprocs()) distributed processes on $(gethostname())")
        else
            error("Not running on saga cluster")
        end

    elseif cluster == "localhost"
        np = length(Sys.cpu_info())
        if nprocs_max > 0
            np = min(np, nprocs_max)
        end
        if nprocs() < np + 1
            Distributed.addprocs(np - nprocs() + 1, topology=:master_worker)
        end
        println("Working with $(nprocs()-1) processes on $(gethostname())")

    else
        error("Cluster $server is unknown. Add it to the FUSE.parallel_environment")
    end
end