import ForwardDiff
import Distributed
import ClusterManagers
import TimerOutputs
import Dates
using InteractiveUtils: summarysize, format_bytes, Markdown

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

# ====== #
# Memory #
# ====== #
Base.@kwdef struct MemTrace
    data::Vector{Tuple{Dates.DateTime,String,Int}} = Tuple{Dates.DateTime,String,Int}[]
end

const memtrace = MemTrace()

function memory_time_tag(txt::String)
    return push!(memtrace.data, (Dates.now(), txt, get_julia_process_memory_usage()))
end

"""
    plot_memtrace(memtrace::MemTrace, n_big_jumps::Int=5, ignore_first_seconds::Int=0)

Plot the memory usage over time from a `MemTrace` object.

# Arguments

  - `n_big_jumps`: number of significant memory jumps to highlight in the plot.
  - `ignore_first_seconds`: number of initial seconds to ignore in the plot. Memory usage will be plotted relative to the memory after this cutoff. Default is `0` (no seconds are ignored).

# Returns

A plot with the following characteristics:

  - Time is displayed on the x-axis as delta seconds since the first recorded time (after ignoring the specified initial seconds).
  - Memory usage is displayed on the y-axis in MB.
  - Memory usage is shown as a scatter plot.
  - The `n_big_jumps` largest jumps in memory are highlighted in red with annotations indicating the action causing each jump.
"""
@recipe function plot_memtrace(memtrace::MemTrace, n_big_jumps::Int=5, ignore_first_seconds::Int=0)
    if isempty(memtrace.data)
        cutoff = 0.0
    else
        cutoff = memtrace.data[1][1] + Dates.Second(ignore_first_seconds)
    end
    filtered_data = filter(point -> point[1] > cutoff, memtrace.data)

    dates = [point[1] for point in filtered_data]
    if !isempty(memtrace.data)
        dates = Dates.value.(dates .- dates[1]) ./ 1000
    end
    action = [point[2] for point in filtered_data]
    mem = [point[3] / 1024 / 1024 for point in filtered_data]
    if !isempty(memtrace.data) && ignore_first_seconds > 0
        mem = mem .- mem[1]
    end

    @series begin
        seriestype := scatter
        label --> ""
        dates, mem
    end

    index = sortperm(diff(mem))[end-min(n_big_jumps - 1, length(mem) - 2):end]
    for i in index
        @series begin
            primary := false
            series_annotations := Plots.text.([action[i], action[i+1]], 6, :red)
            color := :red
            [dates[i], dates[i+1]], [mem[i], mem[i+1]]
        end
    end
    @series begin
        primary := false
        xlabel := "ΔTime [s]"
        ylabel := (ignore_first_seconds > 0 ? "Δ" : "") * "Memory [MB]"
        [], []
    end
end

function get_julia_process_memory_usage()
    pid = getpid()
    mem_info = read(`ps -p $pid -o rss=`, String)
    mem_usage_kb = parse(Int, strip(mem_info))
    return mem_usage_kb * 1024
end

"""
    save(memtrace::MemTrace, filename::String="memtrace.txt")

Save a FUSE memory trace to file
"""
function save(memtrace::MemTrace, filename::String="memtrace.txt")
    open(filename, "w") do file
        for (date, txt, kb) in memtrace.data
            println(file, "$date $kb \"$txt\"")
        end
    end
end

"""
    load(memtrace::MemTrace, filename::String="memtrace.txt")

Load a FUSE memory trace from file
"""
function load(memtrace::MemTrace, filename::String="memtrace.txt")
    open(filename, "r") do file
        for line in eachline(file)
            m = match(r"(\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}\.\d{3})\s+(\d+)\s+\"(.*?)\"", line)
            if m !== nothing
                date = Dates.DateTime(m[1], "yyyy-mm-ddTHH:MM:SS.sss")
                kb = parse(Int, m[2])
                txt = m[3]
                push!(memtrace.data, (date, txt, kb))
            end
        end
    end
    return memtrace
end

"""
    varinfo(m::Module=Main, pattern::Regex=r""; all=false, imported=false, recursive=false, sortby::Symbol=:name, minsize::Int=0)

Return a markdown table giving information about exported global variables in a module, optionally restricted
to those matching `pattern`.

The memory consumption estimate is an approximate lower bound on the size of the internal structure of the object.

  - `all` : also list non-exported objects defined in the module, deprecated objects, and compiler-generated objects.
  - `imported` : also list objects explicitly imported from other modules.
  - `recursive` : recursively include objects in sub-modules, observing the same settings in each.
  - `sortby` : the column to sort results by. Options are `:name` (default), `:size`, and `:summary`.
  - `minsize` : only includes objects with size at least `minsize` bytes. Defaults to `0`.

The output of `varinfo` is intended for display purposes only.  See also [`names`](@ref) to get an array of symbols defined in
a module, which is suitable for more general manipulations.
"""
function varinfo(m::Module=Main, pattern::Regex=r""; all::Bool=false, imported::Bool=false, recursive::Bool=false, sortby::Symbol=:name, minsize::Int=0)
    sortby in (:name, :size, :summary) || throw(ArgumentError("Unrecognized `sortby` value `:$sortby`. Possible options are `:name`, `:size`, and `:summary`"))
    rows = Vector{Any}[]
    workqueue = [(m, "")]
    parents = Module[m]
    while !isempty(workqueue)
        m2, prep = popfirst!(workqueue)
        for v in names(m2; all, imported)
            if !isdefined(m2, v) || !occursin(pattern, string(v))
                continue
            end
            value = getfield(m2, v)
            isbuiltin = value === Base || value === Main || value === Core
            if recursive && !isbuiltin && isa(value, Module) && value !== m2 && nameof(value) === v && value ∉ parents
                push!(parents, value)
                push!(workqueue, (value, "$v."))
            end
            ssize_str, ssize = if isbuiltin
                ("", typemax(Int))
            else
                ss = summarysize(value)
                (format_bytes(ss), ss)
            end
            if ssize >= minsize
                push!(rows, Any[string(prep, v), ssize_str, summary(value), ssize])
            end
        end
    end
    let (col, rev) = if sortby === :name
            1, false
        elseif sortby === :size
            4, true
        elseif sortby === :summary
            3, false
        else
            @assert "unreachable"
        end
        sort!(rows; by=r -> r[col], rev)
    end
    pushfirst!(rows, Any["name", "size", "summary"])

    return Markdown.MD(Any[Markdown.Table(map(r -> r[1:3], rows), Symbol[:l, :r, :l])])
end

"""
    malloc_trim_if_glibc()

Check if the underlying system uses the GNU C Library (GLIBC) and, if so, calls `malloc_trim(0)`
to attempt to release memory back to the system. This function is primarily useful on Linux systems
where GLIBC is the standard C library.

On non-Linux systems, or Linux systems not using GLIBC, the function does nothing

## Notes

  - This method uses a heuristic by checking for the presence of "glibc" in `/proc/version`.
    While this is indicative of GLIBC's presence on many typical systems, it's not a guaranteed check.
  - The actual `malloc_trim` function is specific to the glibc implementation of the C standard library.
"""
function malloc_trim_if_glibc()
    # Check if on Linux
    if Sys.islinux()
        # Check for glibc by reading /proc/version (this is specific to Linux)
        try
            version_info = read("/proc/version", String)
            if occursin("glibc", version_info)
                ccall(:malloc_trim, Cvoid, (Cint,), 0)
                #println("malloc_trim called.")
            else
                #println("Not using glibc.")
            end
        catch e
            #println("Error reading /proc/version: ", e)
        end
    else
        #println("Not on a Linux system.")
    end
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
    parallel_environment(cluster::String="localhost", nworkers::Integer=0, cpus_per_task::Int=1, kw...)

Start multiprocessing environment

  - kw arguments are passed to the Distributed.addprocs

  - nworkers == 0 uses as many workers as the number of available CPUs

  - cpus_per_task can be used to control memory usage
"""
function parallel_environment(cluster::String="localhost", nworkers::Integer=0, cpus_per_task::Int=1, kw...)
    if cluster == "omega"
        if gethostname() ∈ ("omega-a.gat.com", "omega-b.gat.com")
            gigamem_per_node = 512
            cpus_per_node = 128
            if nworkers > 0
                nodes = 4 # omega has 12 ga-ird nodes
                nprocs_max = cpus_per_node * nodes
                nworkers = min(nworkers, nprocs_max)
            end
            np = nworkers + 1
            gigamem_per_cpu = Int(round(gigamem_per_node / cpus_per_node * cpus_per_task))
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
        if gethostname() == "saga.cluster"
            gigamem_per_node = 192
            cpus_per_node = 48
            if nworkers > 0
                nodes = 4  # saga has 6 nodes
                nprocs_max = cpus_per_node * nodes
                nworkers = min(nworkers, nprocs_max)
            end
            np = nworkers + 1
            gigamem_per_cpu = Int(round(gigamem_per_node / cpus_per_node * cpus_per_task))
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
        if nworkers > 0
            nprocs_max = length(Sys.cpu_info())
            nworkers = min(nworkers, nprocs_max)
        end
        np = nworkers + 1
        if Distributed.nprocs() < np
            Distributed.addprocs(np - Distributed.nprocs(); topology=:master_worker)
        end

    else
        error("Cluster `$cluster` is unknown. Use `localhost` or add `$cluster` to the FUSE.parallel_environment")
    end

    return println("Working with $(Distributed.nprocs()-1) workers on $(gethostname())")
end
