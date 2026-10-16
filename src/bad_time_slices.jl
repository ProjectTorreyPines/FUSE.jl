# ================ #
# bad time slices  #
# ================ #
"""
    BadTimeSlices()

Record of the experimental time slices a case could not use: one entry per slice with its `time`, the
`stage` that failed on it, the `reason` (the error) and the `action` taken. Passed as `bad_time_slices`
to the steps that work one time slice at a time (`case_parameters(:D3D)`, `set_ini_act_from_ods!`,
`ActorFitProfiles`, `StudyPostdictive`), it turns a slice that fails from an error of the whole case into
a dropped slice, following [`droppable_time_slices`](@ref). Without it (`nothing`) those steps throw as before.
"""
const BadTimeSliceRecord = @NamedTuple{time::Float64, stage::String, reason::String, action::String}

struct BadTimeSlices
    records::Vector{BadTimeSliceRecord}
end

BadTimeSlices() = BadTimeSlices(BadTimeSliceRecord[])

function record!(bad::BadTimeSlices, time::Real, stage::AbstractString, reason, action::AbstractString)
    reason_str = reason isa Exception ? first(sprint(showerror, reason), 300) : string(reason)
    push!(bad.records, BadTimeSliceRecord((time, stage, reason_str, action)))
    return bad
end

"""
    droppable_time_slices(bad::AbstractVector{<:Integer}, n::Integer) -> (drop::Vector{Int}, cut::Union{Nothing,Int})

Of `n` time slices of which `bad` (sorted) cannot be used: runs of them at the start or at the end, and
single ones in between, are dropped. The first run of two or more in between ends the shot there: it and
every later slice are dropped, and `cut` is its first index. Errors if fewer than two slices are left.
"""
function droppable_time_slices(bad::AbstractVector{<:Integer}, n::Integer)
    runs = UnitRange{Int}[]
    for k in bad
        !isempty(runs) && last(runs[end]) == k - 1 ? (runs[end] = first(runs[end]):k) : push!(runs, k:k)
    end
    k_cut = findfirst(r -> first(r) > 1 && last(r) < n && length(r) > 1, runs)
    cut = k_cut === nothing ? nothing : first(runs[k_cut])
    drop = cut === nothing ? collect(Int, bad) : vcat(filter(<(cut), bad), cut:n)
    n - length(drop) >= 2 || error("only $(n - length(drop)) of $n time slices usable")
    return drop, cut
end

"""
    each_time_slice!(f, bad_time_slices::Union{Nothing,BadTimeSlices}, times::AbstractVector{<:Real}, stage::AbstractString) -> drop::Vector{Int}

`f(k)` for every slice `k` of `times`, in order. Without `bad_time_slices` a failure throws. With it,
the failing slices are recorded under `stage` and the indices to drop are returned (see
[`droppable_time_slices`](@ref)); slices dropped only because the shot ends before them are recorded
too. Removing them is the caller's job, e.g. [`drop_time_slices!`](@ref).
"""
function each_time_slice!(f::Function, bad_time_slices::Union{Nothing,BadTimeSlices}, times::AbstractVector{<:Real}, stage::AbstractString)
    if bad_time_slices === nothing
        foreach(f, eachindex(times))
        return Int[]
    end
    errs = Dict{Int,Any}()
    for k in eachindex(times)
        try
            f(k)
        catch e
            isa(e, InterruptException) && rethrow()
            errs[k] = e
        end
    end
    isempty(errs) && return Int[]

    bad = sort!(collect(keys(errs)))
    drop, cut = try
        droppable_time_slices(bad, length(times))
    catch e
        throw(ErrorException("$stage: $(sprint(showerror, e)); first failure at $(times[bad[1]]) s: $(first(sprint(showerror, errs[bad[1]]), 300))"))
    end
    for k in drop
        if haskey(errs, k)
            record!(bad_time_slices, times[k], stage, errs[k], cut !== nothing && k >= cut ? "dropped, shot ends here" : "dropped")
        else
            record!(bad_time_slices, times[k], stage, "after the shot ended at $(times[cut]) s", "dropped")
        end
    end
    @warn "$stage: dropped $(length(drop)) of $(length(times)) time slices ($(length(bad)) failed" *
          (cut === nothing ? ")" : "; the shot ends at $(times[cut]) s)") * "; first failure at $(times[bad[1]]) s: $(first(sprint(showerror, errs[bad[1]]), 200))"
    return drop
end

"""
    drop_time_slices!(dd::IMAS.DD, times::AbstractVector{<:Real}, drop::AbstractVector{<:Integer})

Remove the `drop` slices of the time basis `times` from the time-sliced experimental data a case is built
on: the equilibrium time slices (with `equilibrium.time` and `vacuum_toroidal_field.b0`) and the
core_profiles profiles. When the first or last slice of `times` is dropped, everything before the first
or after the last slice kept goes too, so a basis other than these IDSs' own still ends them there.
"""
function drop_time_slices!(dd::IMAS.DD, times::AbstractVector{<:Real}, drop::AbstractVector{<:Integer})
    isempty(drop) && return dd
    keep = setdiff(eachindex(times), drop)
    t_first = first(eachindex(times)) in drop ? times[first(keep)] : -Inf
    t_last = last(eachindex(times)) in drop ? times[last(keep)] : Inf
    dropped(t) = t < t_first || t > t_last || t in view(times, drop)

    eq = dd.equilibrium
    if !isempty(eq.time_slice)
        gone = findall(eqt -> dropped(eqt.time), eq.time_slice)
        kept = setdiff(eachindex(eq.time_slice), gone)
        b0 = ismissing(eq.vacuum_toroidal_field, :b0) ? Float64[] : eq.vacuum_toroidal_field.b0
        length(b0) == length(eq.time_slice) && (eq.vacuum_toroidal_field.b0 = b0[kept])
        eq.time = [eq.time_slice[k].time for k in kept]
        foreach(k -> deleteat!(eq.time_slice, k), reverse(gone))
    end
    cp = dd.core_profiles
    if !isempty(cp.profiles_1d)
        foreach(k -> deleteat!(cp.profiles_1d, k), reverse(findall(cp1d -> dropped(cp1d.time), cp.profiles_1d)))
        cp.time = [cp1d.time for cp1d in cp.profiles_1d]
    end
    return dd
end

"""
    flux_surfaces!(dd::IMAS.DD, bad_time_slices::Union{Nothing,BadTimeSlices})

`IMAS.flux_surfaces` on every equilibrium time slice; with `bad_time_slices`, the slices it fails on are
dropped from `dd` (see [`each_time_slice!`](@ref)).
"""
function flux_surfaces!(dd::IMAS.DD, bad_time_slices::Union{Nothing,BadTimeSlices})
    eq = dd.equilibrium
    wall_r, wall_z = IMAS.first_wall(dd.wall)
    if bad_time_slices === nothing
        IMAS.flux_surfaces(eq, wall_r, wall_z)
        return dd
    end
    # in order: each slice starts its magnetic-axis search from the previous slice's axis, so a failed
    # slice's axis (wherever its search diverged to) is cleared before the next slice reads it
    drop = each_time_slice!(bad_time_slices, eq.time, "equilibrium flux surfaces") do k
        eqt = eq.time_slice[k]
        findfirst(:rectangular, eqt.profiles_2d) === nothing && return
        try
            IMAS.flux_surfaces(eqt, wall_r, wall_z)
        catch
            empty!(eqt.global_quantities.magnetic_axis, :r)
            empty!(eqt.global_quantities.magnetic_axis, :z)
            rethrow()
        end
    end
    return drop_time_slices!(dd, eq.time, drop)
end

const MAX_START_DROPS = 5

"""
    init_dropping_bad_start!(init::Function, dd::IMAS.DD, ini::ParametersAllInits, bad_time_slices::Union{Nothing,BadTimeSlices};
        prepare::Function=() -> nothing)

`prepare()` then `init()`: what follows from `ini.time.simulation_start` (e.g. the time window), then the
initialization of `dd` from it. With `bad_time_slices`, an `init()` failure drops the experimental slice at
the start from `ini.general.dd`, moves the start to the next equilibrium slice, empties `dd` and does both
again, recording each slice dropped; after `MAX_START_DROPS` (5) of them, or with no later slice, the
failure is thrown. A `prepare()` failure is thrown as is.
"""
function init_dropping_bad_start!(init::Function, dd::IMAS.DD, ini::ParametersAllInits, bad_time_slices::Union{Nothing,BadTimeSlices};
    prepare::Function=() -> nothing)
    if bad_time_slices === nothing
        prepare()
        return init()
    end
    for n_dropped in 0:MAX_START_DROPS
        prepare()
        try
            return init()
        catch e
            isa(e, InterruptException) && rethrow()
            t0 = ini.time.simulation_start
            eq_time = ini.general.dd.equilibrium.time
            k0 = something(findlast(<=(t0), eq_time), 1)
            (n_dropped == MAX_START_DROPS || k0 + 2 > length(eq_time)) && rethrow()
            record!(bad_time_slices, eq_time[k0], "init at the simulation start", e, "dropped, start moved to $(eq_time[k0+1]) s")
            ini.time.simulation_start = eq_time[k0+1]
            drop_time_slices!(ini.general.dd, eq_time, [k0])
            empty!(dd)
        end
    end
end
