import ForwardDiff
import Distributed
import ClusterManagers
import TimerOutputs
import Dates
using SimulationParameters

# =========== #
# global_time #
# =========== #
# this is used to disambiguate between global_time of SimulationParameters and IMAS
function global_time(par::Union{<:SimulationParameters.AbstractParameters,<:SimulationParameters.AbstractParameter})::Float64
    return SimulationParameters.global_time(par)
end

function global_time(par::Union{<:SimulationParameters.AbstractParameters,<:SimulationParameters.AbstractParameter}, time0::Float64)::Float64
    return SimulationParameters.global_time(par, time0)
end

function global_time(@nospecialize(ids::Union{IMAS.IDS,IMAS.IDSvector}))::Float64
    return IMAS.global_time(ids)
end

function global_time(@nospecialize(ids::Union{IMAS.IDS,IMAS.IDSvector}), time0::Float64)
    return IMAS.global_time(ids, time0)
end

# ====== #
# Timing #
# ====== #
const timer = TimerOutputs.TimerOutput()

function TimerOutputs.reset_timer!(to::TimerOutputs.TimerOutput, section::String)
    pop!(to.inner_timers, section, nothing)
    to.prev_timer_label = ""
    return to.prev_timer = nothing
end

function TimerOutputs.reset_timer!(section::String)
    pop!(timer.inner_timers, section, nothing)
    timer.prev_timer_label = ""
    return timer.prev_timer = nothing
end

# ==== #
# Math #
# ==== #
function unwrap(v, inplace=false)
    unwrapped = inplace ? v : copy(v)
    for i in 2:length(v)
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
    return args = map(x -> vcat([x for k in 1:n]...)[1:n], args)
end

"""
    mirror_bound(x::T, l::T, u::T) where {T<:Real}

Return tuple with value of x bounded between l and u
The bounding is done by mirroring the value at the bound limits.
"""
function mirror_bound(x::T, l::T, u::T) where {T<:Real}
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

# =========== #
# Convex Hull #
# =========== #
import VacuumFields: convex_hull!, convex_hull

# ======== #
# TraceCAD #
# ======== #
import FileIO, JpegTurbo

"""
    dd_build_layers_to_ini(dd::IMAS.dd)

Utility function to convert layers in dd.build to layers in `ini.build.layers = layers = OrderedCollections.OrderedDict{Symbol,Float64}()`
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

"""
TIP! use it with Interact.@manipulate, like this:

    import Interact
    Interact.@manipulate for x_length in range(0.5,5,100), x_offset in range(-1,1,100), y_offset in range(-1,1,100)
        plot(FUSE.TraceCAD(:D3D, x_length, x_offset, y_offset))
        plot!(dd.pf_active)
    end
"""
@recipe function plot_trace_cad(cad::TraceCAD)
    img = FileIO.load(joinpath(@__DIR__, "cases", "$(cad.name).jpg"))
    img .= img[end:-1:1, 1:end]
    x = range(0, cad.x_length, size(img)[1]) .+ cad.x_offset
    y = range(-0.5, 0.5, size(img)[2]) .* (size(img)[1] / size(img)[2]) * (x[end] - x[1]) .- cad.y_offset
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
    function parallel_environment(cluster::String="localhost", nworkers::Integer=0, cpus_per_task::Int=1,memory_usage_fraction::Float64=0.5, kw...)

Start multiprocessing environment

  - kw arguments are passed to the Distributed.addprocs

  - nworkers == 0 uses as many workers as the number of available CPUs
  - cpus_per_task can be used to control memory usage
  - memory_usage_fraction is the fraction of peak memory that can be used
"""
function parallel_environment(cluster::String="localhost", nworkers::Integer=0, cpus_per_task::Int=1; memory_usage_fraction::Float64=0.5, kw...)
    if cluster == "omega"
        if occursin("omega", gethostname())
            gigamem_per_node = 512
            cpus_per_node = 128
            if nworkers > 0
                nodes = 4 # omega has 12 ga-ird nodes
                nprocs_max = cpus_per_node * nodes
                nworkers = min(nworkers, nprocs_max)
            end
            np = nworkers + 1
            gigamem_per_cpu = Int(round(memory_usage_fraction * gigamem_per_node / cpus_per_node * cpus_per_task))
            ENV["JULIA_WORKER_TIMEOUT"] = "360"
            if Distributed.nprocs() < np
                Distributed.addprocs(
                    ClusterManagers.SlurmManager(np - Distributed.nprocs());
                    partition="ga-ird",
                    exclusive="",
                    topology=:master_worker,
                    time="99:99:99",
                    cpus_per_task,
                    exeflags=["--threads=$(cpus_per_task)", "--heap-size-hint=$(gigamem_per_cpu)G"],
                    kw...
                )
            end
        else
            error("Not running on omega cluster")
        end

    elseif cluster == "saga"
        if occursin("saga", gethostname())
            gigamem_per_node = 192
            cpus_per_node = 48
            if nworkers > 0
                nodes = 4  # saga has 6 nodes
                nprocs_max = cpus_per_node * nodes
                nworkers = min(nworkers, nprocs_max)
            end
            np = nworkers + 1
            gigamem_per_cpu = Int(round(memory_usage_fraction * gigamem_per_node / cpus_per_node * cpus_per_task))
            ENV["JULIA_WORKER_TIMEOUT"] = "180"
            if Distributed.nprocs() < np
                Distributed.addprocs(
                    ClusterManagers.SlurmManager(np - Distributed.nprocs());
                    exclusive="",
                    topology=:master_worker,
                    cpus_per_task,
                    exeflags=["--threads=$(cpus_per_task)", "--heap-size-hint=$(gigamem_per_cpu)G"],
                    kw...
                )
            end
        else
            error("Not running on saga cluster")
        end

    elseif cluster == "localhost"
        mem_size = Int(ceil(localhost_memory() * memory_usage_fraction))

        if nworkers > 0
            nprocs_max = length(Sys.cpu_info())
            nworkers = min(nworkers, nprocs_max)
        end
        np = nworkers + 1
        if Distributed.nprocs() < np
            Distributed.addprocs(np - Distributed.nprocs(); topology=:master_worker, exeflags=["--heap-size-hint=$(mem_size)G"])
        end

    else
        error("Cluster `$cluster` is unknown. Use `localhost` or add `$cluster` to the FUSE.parallel_environment")
    end

    return println("Working with $(Distributed.nprocs()-1) workers on $(gethostname())")
end

"""
    localhost_memory()

Determines what the maximum memory is based on the device type (apple, windows, unix,linux)
"""
function localhost_memory()
    if Sys.isapple()
        cmd = `sysctl hw.memsize` # for OSX
        mem_size = parse(Int, match(r"\d+", readchomp(cmd)).match) / 1024^3
    elseif Sys.isunix()
        # General Unix command (including macOS and Linux)
        cmd = `free -b` # get memory in bytes
        mem_size = parse(Int, match(r"\d+", readchomp(cmd)).match) / 1024^3
    elseif Sys.iswindows()
        # Windows command
        cmd = `wmic ComputerSystem get TotalPhysicalMemory`
        mem_size = parse(Int, match(r"\d+", readchomp(cmd)).match) / 1024^3
    elseif Sys.islinux()
        # Linux-specific command
        cmd = `grep MemTotal /proc/meminfo`
        mem_size = parse(Int, match(r"\d+", readchomp(cmd)).match) / 1024^2 # Linux reports in KB
    else
        error("couldn't determine the mem_size")
    end
    return mem_size
end