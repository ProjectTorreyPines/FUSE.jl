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
    parallel_environment(cluster::String="localhost", nworkers::Integer=0, cpus_per_task::Int=1,memory_usage_fraction::Float64=0.5, kw...)

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
                nodes = 4  # don't use more than 4 nodes (omega has 12 ird nodes)
                nprocs_max = cpus_per_node * nodes
                nworkers = min(nworkers, nprocs_max)
            end
            np = nworkers + 1
            gigamem_per_cpu = Int(ceil(memory_usage_fraction * gigamem_per_node / cpus_per_node * cpus_per_task))
            ENV["JULIA_WORKER_TIMEOUT"] = "360"
            if Distributed.nprocs() < np
                Distributed.addprocs(
                    ClusterManagers.SlurmManager(np - Distributed.nprocs());
                    partition="ga-ird",
                    exclusive="",
                    topology=:master_worker,
                    time="99:99:99",
                    cpus_per_task,
                    exeflags=["--threads=$(cpus_per_task)", "--heap-size-hint=$(gigamem_per_cpu)G", "--project=$(Base.active_project())"],
                    kw...
                )
            end
        else
            error("Not running on omega cluster")
        end

    elseif cluster == "stellar"
        if occursin("stellar", gethostname())
            gigamem_per_node = 192
            cpus_per_node = 96
            if nworkers > 0
                nodes = 4 # don't use more than 4 nodes
                nprocs_max = cpus_per_node * nodes
                nworkers = min(nworkers, nprocs_max)
            end
            np = nworkers + 1
            gigamem_per_cpu = Int(ceil(memory_usage_fraction * gigamem_per_node / cpus_per_node * cpus_per_task))
            ENV["JULIA_WORKER_TIMEOUT"] = "360"
            if Distributed.nprocs() < np
                Distributed.addprocs(
                    ClusterManagers.SlurmManager(np - Distributed.nprocs());
                    partition="pppl-medium",
                    exclusive="",
                    topology=:master_worker,
                    time="00:48:00",
                    cpus_per_task,
                    exeflags=["--threads=$(cpus_per_task)", "--heap-size-hint=$(gigamem_per_cpu)G", "--project=$(Base.active_project())"],
                    kw...
                )
            end
        else
            error("Not running on stellar cluster")
        end

    elseif cluster == "saga"
        if occursin("saga", gethostname())
            gigamem_per_node = 192
            cpus_per_node = 48
            if nworkers > 0
                nodes = 4 # don't use more than 4 nodes (saga has 6 nodes)
                nprocs_max = cpus_per_node * nodes
                nworkers = min(nworkers, nprocs_max)
            end
            np = nworkers + 1
            gigamem_per_cpu = Int(ceil(memory_usage_fraction * gigamem_per_node / cpus_per_node * cpus_per_task))
            ENV["JULIA_WORKER_TIMEOUT"] = "180"
            if Distributed.nprocs() < np
                Distributed.addprocs(
                    ClusterManagers.SlurmManager(np - Distributed.nprocs());
                    exclusive="",
                    topology=:master_worker,
                    cpus_per_task,
                    exeflags=["--threads=$(cpus_per_task)", "--heap-size-hint=$(gigamem_per_cpu)G", "--project=$(Base.active_project())"],
                    kw...
                )
            end
        else
            error("Not running on saga cluster")
        end

    elseif cluster == "feynman"
        if occursin("feynman", gethostname())
            gigamem_per_node = 800
            cpus_per_node = 30
            if nworkers > 0
                nodes = 1 # don't use more than 1 node
                nprocs_max = cpus_per_node * nodes
                nworkers = min(nworkers, nprocs_max)
            end
            np = nworkers + 1
            gigamem_per_cpu = Int(ceil(memory_usage_fraction * gigamem_per_node / cpus_per_node * cpus_per_task))
            ENV["JULIA_WORKER_TIMEOUT"] = "360"
            if Distributed.nprocs() < np
                Distributed.addprocs(
                    ClusterManagers.SlurmManager(np - Distributed.nprocs());
                    partition="LocalQ",
                    topology=:master_worker,
                    time="99:99:99",
                    cpus_per_task,
                    exeflags=["--threads=$(cpus_per_task)", "--heap-size-hint=$(gigamem_per_cpu)G", "--project=$(Base.active_project())"],
                    kw...
                )
            end
        else
            error("Not running on feynman cluster")
        end

    elseif cluster == "engaging"
        if occursin("eofe", gethostname())
            gigamem_per_node = 512
            cpus_per_node = 64
            if nworkers > 0
                nodes = 4 # don't use more than 4 nodes
                nprocs_max = cpus_per_node * nodes
                nworkers = min(nworkers, nprocs_max)
            end
            np = nworkers + 1
            gigamem_per_cpu = Int(ceil(memory_usage_fraction * gigamem_per_node / cpus_per_node * cpus_per_task))
            ENV["JULIA_WORKER_TIMEOUT"] = "360"
            if Distributed.nprocs() < np
                Distributed.addprocs(
                    ClusterManagers.SlurmManager(np - Distributed.nprocs());
                    partition="sched_mit_psfc_r8",
                    exclusive="",
                    topology=:master_worker,
                    time="7:59:59",
                    cpus_per_task,
                    exeflags=["--threads=$(cpus_per_task)", "--heap-size-hint=$(gigamem_per_cpu)G", "--project=$(Base.active_project())"],
                    kw...
                )
            end
        else
            error("Not running on engaging cluster")
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
        mem_size = parse(Int, match(r"\d+", readchomp(cmd)).match) / 1024^3 # GiB
    elseif Sys.iswindows()
        # Windows command
        cmd = `powershell -Command "Get-CimInstance Win32_ComputerSystem | Select-Object -ExpandProperty TotalPhysicalMemory"`
        mem_bytes = parse(Int, readchomp(cmd))
        mem_size = mem_bytes / 1024^3 # GiB
    elseif Sys.islinux()
        # Linux-specific command
        cmd = `grep MemTotal /proc/meminfo`
        mem_size = parse(Int, match(r"\d+", readchomp(cmd)).match) / 1024^2 # GiB
    else
        error("couldn't determine the mem_size")
    end
    return mem_size
end

# ====== #
# errors #
# ====== #
struct MissingExtensionError <: Exception
    actor_name::String
    package_name::String
end

Base.showerror(io::IO, e::MissingExtensionError) =
    print(io, "The FUSE actor $(e.actor_name) cannot be run because the Julia package $(e.package_name).jl is not loaded. Please load it to enable this feature.")


"""
    take_last!(ch::Distributed.RemoteChannel)

Drains all available messages from `ch` and returns the last one.
If no messages are present, returns `nothing`.
"""
function take_last!(ch::Distributed.RemoteChannel)
    last_msg = nothing
    while isready(ch)
        last_msg = take!(ch)
    end
    return last_msg
end

"""
    prepare_file_lock_channels()

Creates file lock channels for distributed file I/O control.
Returns a named tuple `(w2m, m2w)` where:

  - `w2m` is a Dict mapping worker IDs to RemoteChannels (created on master)
    for worker-to-master signals.
  - `m2w` is a Dict mapping worker IDs to RemoteChannels (created on each worker)
    for master-to-worker signals.
"""
function prepare_file_lock_channels()
    # worker-to-master
    w2m = Dict{Int,Distributed.RemoteChannel}()
    for w in Distributed.workers()
        w2m[w] = Distributed.RemoteChannel(() -> Distributed.Channel{Symbol}(100), 1)
        put!(w2m[w], :unlock)  # Initialize with unlock signal
    end

    # master-to-worker
    m2w = Dict{Int,Distributed.RemoteChannel}()
    for w in Distributed.workers()
        m2w[w] = Distributed.RemoteChannel(() -> Distributed.Channel{Symbol}(100), w)
        Distributed.remotecall_wait(() -> put!(m2w[w], :unlock), w)  # Initialize with unlock signal
    end

    channels = (w2m=w2m, m2w=m2w)
    return channels
end

"""
    database_IO_manager(save_folder, file_lock_channels, study_status, interval; verbose=false)

Periodically manages file I/O for a distributed study.
It waits for `interval` seconds, sends lock signals to workers via `m2w`,
then polls `w2m` until all workers report `:unlock`,
merges temporary files, and finally sends unlock signals.
Continues while `study_status[] == :running`.
"""
function database_IO_manager(
    save_folder::AbstractString,
    file_lock_channels::NamedTuple{(:w2m, :m2w),Tuple{Dict{Int,Distributed.RemoteChannel},Dict{Int,Distributed.RemoteChannel}}},
    study_status::Base.RefValue{Symbol},
    interval::Integer;
    verbose::Bool=false)

    @assert interval > 0 "interval must be a positive integer"

    start_time = time()
    function elapsed()
        seconds = time() - start_time
        return @sprintf("[%d s] ", seconds)
    end

    printstyled("\n" * elapsed() * "Master: database_IO_manager starts\n"; bold=true, color=:blue)

    while study_status[] == :running
        NN = interval*5
        for _ in 1:NN
            sleep(0.2)
            if study_status[] == :finished
                printstyled("\n" * elapsed() * "Master: database_IO_manager finished\n"; bold=true, color=:blue)
                return
            end
        end

        for (w, ch) in file_lock_channels.m2w
            @async put!(ch, :lock)
        end
        verbose && println(elapsed() * "Master: lock signals are sent")

        nworkers = Distributed.nworkers()
        worker_status = Dict{Int,Symbol}()
        while true
            for wid in Distributed.workers()
                if isready(file_lock_channels.w2m[wid])
                    status = take_last!(file_lock_channels.w2m[wid])
                    worker_status[wid] = status

                    if status == :unlock
                        put!(file_lock_channels.w2m[wid], :unlock)
                    elseif status == :lock
                        verbose && println(elapsed() * "Master: received signal from worker $wid: $status")
                    end
                end
            end
            # Only break when we have received a response from all workers and they are all :unlock.
            if length(worker_status) == nworkers && all(==(:unlock), values(worker_status))
                break
            end
            sleep(0.1)  # Avoid busy waiting
        end
        verbose && println(elapsed() * "Master: all workers unlocked files")

        printstyled("\n"*elapsed() * "Master: merges temporary study files\n"; bold=true, color=:blue)
        _merge_tmp_study_files(save_folder; cleanup=true)
        verbose && println(elapsed() * "Master: merge finished")

        for (w, ch) in file_lock_channels.m2w
            @async put!(ch, :unlock)
        end
        verbose && println(elapsed() * "Master: unlock signals are sent")
    end

    printstyled("\n" * elapsed() * "Master: database_IO_manager finished\n"; bold=true, color=:blue)

    return
end


"""
    wait_for_unlock!(ch::Distributed.RemoteChannel; timeout=30, poll_interval=0.1, stable_period=0.5, verbose=false)

Continuously polls the channel `ch` until an `:unlock` signal is observed for at least `stable_period` seconds,
or a timeout is reached.
If the unlock state persists for the stable period, reinserts `:unlock` into `ch` and returns true.
If a lock signal is received during that period or the timeout is reached, it logs a warning and returns false.
"""
function wait_for_unlock!(ch::Distributed.RemoteChannel; timeout::Real=300.0, poll_interval::Real=0.1, stable_period::Real=0.5, verbose::Bool=false)
    myid = Distributed.myid()
    start_time = time()

    stable_start = nothing  # time when we first see a stable unlock

    # Continuously poll the chanel for a message.
    while true
        if isready(ch)
            # Drain the channel and get the most recent message.
            latest = take_last!(ch)
            if latest == :unlock
                put!(ch, :unlock)  # Reinsert unlock so state persists.
                # If we haven't started the stable timer, start it.
                if isnothing(stable_start)
                    stable_start = time()
                    verbose && println("Worker $(myid) observed unlock; starting stable period.")
                end
                # If unlock has persisted for the required stable period, return true.
                if time() - stable_start >= stable_period
                    verbose && println("Worker $(myid) maintained unlock for stable period; resuming.")
                    return true
                end
            elseif latest == :lock
                verbose && println("Worker $(myid) sees lock signal, still waiting")
                stable_start = nothing  # Reset stable timer if a lock appears.
            else
                verbose && println("Worker $(myid) received unknown signal: ", latest)
            end
        end

        sleep(poll_interval)
        if time() - start_time > timeout
            @warn "Worker $(myid): Timeout reached waiting for unlock signal, proceeding anyway"
            return false
        end
    end
end
