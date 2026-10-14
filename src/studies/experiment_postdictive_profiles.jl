#= ================================================== =#
#  Externally supplied kinetic profiles in postdiction  #
#= ================================================== =#
# Postdictive runs replay `ini.general.dd.core_profiles` (`act.ActorReplay.replay_dd`): profiles from another
# reconstruction can be written there before `init!`. The replay reads `profiles_1d[t]`, the slice at or before t,
# so the core_profiles time basis decides which profile each time step sees. These are the pieces to build that
# basis, assign source slices to it, and report which source slice each step replayed.

"""
    postdictive_instants(t0::Real, Δt::Real, Nt::Int)

ActorDynamicPlasma's half-step instants `t0 + kk * δt / 2`, `kk in 0:2Nt`, with its own arithmetic (`δt = Δt / Nt`).
Every other entry (`[1:2:end]`) is a full step.
"""
function postdictive_instants(t0::Real, Δt::Real, Nt::Int)
    δt = Δt / Nt
    return [t0 + kk * (δt / 2) for kk in 0:(2*Nt)]
end

"""
    core_profiles_time_base(eq_times::AbstractVector{<:Real}, source_times::AbstractVector{<:Real}, sim_times::AbstractVector{<:Real}; offset::Float64=1e-9)

A core_profiles time basis for replaying profiles from another source: the union of the equilibrium times, the
source times and a tick `offset` below every simulation instant, clipped to the equilibrium span (core_profiles
grids are not defined outside it), strictly increasing.

On the equilibrium basis alone the zero-order hold of the replay can serve a profile one whole equilibrium
interval old; the tick below each instant makes it land on a slice made for that instant. An equilibrium or
source time that falls between a tick and its instant is dropped, since the hold would land on it instead.
"""
function core_profiles_time_base(
    eq_times::AbstractVector{<:Real},
    source_times::AbstractVector{<:Real},
    sim_times::AbstractVector{<:Real};
    offset::Float64=1e-9
)
    @assert !isempty(eq_times) "core_profiles_time_base: empty equilibrium time base"
    t_lo, t_hi = Float64(first(eq_times)), Float64(last(eq_times))
    in_span(t) = t_lo <= t <= t_hi
    sim_ticks = sort!(unique!(Float64[t - offset for t in sim_times if in_span(t - offset)]))
    function shadows_sim_tick(x::Float64)
        isempty(sim_ticks) && return false
        j = searchsortedlast(sim_ticks, x)
        return j >= firstindex(sim_ticks) && x <= sim_ticks[j] + offset
    end
    extra = Float64[]
    for t in Iterators.flatten((Float64.(eq_times), Float64.(source_times)))
        if in_span(t) && !shadows_sim_tick(t)
            push!(extra, t)
        end
    end
    out = sort!(unique!(vcat(sim_ticks, extra)))
    @assert issorted(out) && allunique(out) "core_profiles_time_base: base is not strictly increasing"
    return out
end

"""
    reseed_core_profiles_time!(dd::IMAS.DD, act::ParametersAllActors, time_base::AbstractVector{<:Real}; time_averaging::Float64=0.05, rho_averaging::Float64=0.25)

Refit `dd.core_profiles` from the diagnostics on `time_base` (`ActorFitProfiles` with `time_basis_ids=:core_profiles`).
Thomson is not rescaled to the interferometer again: the `case_parameters` fit already scaled it in place, and rescales compound.
Returns the number of slices, or 0 (nothing done) if fewer than 5 CER channels can be fitted.
"""
function reseed_core_profiles_time!(dd::IMAS.DD, act::ParametersAllActors, time_base::AbstractVector{<:Real};
    time_averaging::Float64=0.05, rho_averaging::Float64=0.25)
    # as case_parameters(:D3D) counts them (a channel with one time sample cannot be fitted)
    n_cer = count(ch -> !isempty(ch.ion) && IMAS.hasdata(ch.ion[1].t_i, :data) && length(ch.ion[1].t_i.time) >= 2, dd.charge_exchange.channel)
    if n_cer < 5
        @warn "reseed_core_profiles_time!: only $n_cer CER channel(s), keeping the core_profiles time basis"
        return 0
    end
    dd.core_profiles.time = Float64.(time_base)
    ActorFitProfiles(dd, act; time_averaging, rho_averaging, time_basis_ids=:core_profiles, use_interferometer=false)
    return length(dd.core_profiles.profiles_1d)
end

"""
    assign_profile_slices(cp_times::AbstractVector{<:Real}, source_times::AbstractVector{<:Real})

Nearest source slice for each core_profiles time (`k[i]` indexes `source_times`), decided for all times at once.

Distances are compared rounded to 1 μs, so floating-point noise cannot decide a tie between two equidistant
slices; a tie goes to the candidate not yet used. Returns `(; k, n_deadheat, n_recovered)`.
"""
function assign_profile_slices(cp_times::AbstractVector{<:Real}, source_times::AbstractVector{<:Real})
    @assert issorted(cp_times) "assign_profile_slices: cp_times must be sorted"
    n = length(source_times)
    cand = Vector{Tuple{Int,Int}}(undef, length(cp_times))
    for (i, t) in enumerate(cp_times)
        j = searchsortedfirst(source_times, t)
        if j <= 1
            cand[i] = (1, 1)
        elseif j > n
            cand[i] = (n, n)
        else
            dl = round(t - source_times[j-1]; digits=6)
            dr = round(source_times[j] - t; digits=6)
            cand[i] = dl < dr ? (j - 1, j - 1) : dr < dl ? (j, j) : (j - 1, j)
        end
    end
    k = Vector{Int}(undef, length(cp_times))
    used = falses(n)
    n_deadheat = 0
    n_recovered = 0
    for i in eachindex(cp_times)
        lo, hi = cand[i]
        lo == hi || (n_deadheat += 1)
        if lo != hi && used[lo]
            k[i] = hi
            used[hi] || (n_recovered += 1)
        else
            k[i] = lo
        end
        used[k[i]] = true
    end
    return (; k, n_deadheat, n_recovered)
end

"""
    replay_pairing(cp_times::AbstractVector{<:Real}, sim_times::AbstractVector{<:Real}, k_of::AbstractDict, source_times::AbstractVector{<:Real}, replaced::AbstractVector{<:Real}; t_lo::Real, t_hi::Real)

Which source slice the replay serves at each simulation time: the core_profiles slice at or before t (`tick`), and
the source slice `k_of[tick]` if that slice was replaced. `outcome` is `:source`, `:fit_skipped` (inside
`[t_lo, t_hi]` but not replaced), `:fit_out_of_band` or `:no_slice`. One `(; t_sim, tick, source_time, outcome)` per time.
"""
function replay_pairing(cp_times::AbstractVector{<:Real}, sim_times::AbstractVector{<:Real}, k_of::AbstractDict,
    source_times::AbstractVector{<:Real}, replaced::AbstractVector{<:Real}; t_lo::Real, t_hi::Real)
    @assert issorted(cp_times) "replay_pairing: cp_times must be sorted"
    is_replaced = Set{Float64}(Float64.(replaced))
    rows = NamedTuple{(:t_sim, :tick, :source_time, :outcome),Tuple{Float64,Float64,Union{Nothing,Float64},Symbol}}[]
    for t in sim_times
        j = searchsortedlast(cp_times, t)
        if j < firstindex(cp_times)
            push!(rows, (; t_sim=Float64(t), tick=NaN, source_time=nothing, outcome=:no_slice))
            continue
        end
        tick = Float64(cp_times[j])
        if tick in is_replaced
            push!(rows, (; t_sim=Float64(t), tick, source_time=Float64(source_times[k_of[tick]]), outcome=:source))
        else
            push!(rows, (; t_sim=Float64(t), tick, source_time=nothing, outcome=(t_lo <= tick <= t_hi ? :fit_skipped : :fit_out_of_band)))
        end
    end
    return rows
end
