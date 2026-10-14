# Periodic sawtooth crash for ActorQED (`crash_period > 0`): a reset of the q implied by the current in the
# latest equilibrium geometry, applied where q dips below `crash_q`, as in OpenFUSIONToolkit's `saw_reset_1d` (rule 2)

"""
    sawtooth_q_times_I(ρ::AbstractVector, eqt::IMAS.equilibrium__time_slice)

q·|I| = dΦ/dρ dV/dρ gm2 / (4π² μ0) of `eqt` on `ρ` (rho_tor_norm).

The psi-uniform equilibrium grid has no surface between the axis and ρ ≈ 0.05, so H/ρ² (regular at the axis) is
least-squares fitted in s = ρ² over ρ ≤ 0.7 with a cubic plus exp(-s/0.08²), blended into the plain interpolation over ρ 0.6-0.7.
"""
function sawtooth_q_times_I(ρ::AbstractVector, eqt::IMAS.equilibrium__time_slice)
    e1 = eqt.profiles_1d
    x = e1.rho_tor_norm
    H = abs.(2π * eqt.global_quantities.vacuum_toroidal_field.b0 .* e1.rho_tor .* e1.dvolume_drho_tor .* e1.gm2) ./ (4π^2 * 4e-7π)
    k = findall(>(0), x)
    s, h = x[k] .^ 2, H[k] ./ x[k] .^ 2
    basis(s) = [s^0 s s^2 s^3 exp(-s / 0.08^2)]
    m = x[k] .<= 0.7
    c = reduce(vcat, basis.(s[m])) \ h[m]
    hf = [(basis(r^2)*c)[1] for r in ρ]
    hi = IMAS.interp1d(x[k], h, :cubic).(clamp.(ρ, x[k[1]], x[end]))
    b = [(u = clamp((r - 0.6) / 0.1, 0.0, 1.0); u^2 * (3 - 2u)) for r in ρ]
    return ρ .^ 2 .* ((1 .- b) .* hf .+ b .* hi)
end

"""
    sawtooth_true_q(dd::IMAS.DD, q_s::Real, dq::Real)

q implied by `cp1d.j_total` in the latest equilibrium geometry, and `im`: the first index beyond the last q < `q_s`
where q ≥ `q_s + dq` (the mixing radius), or `nothing`.
"""
function sawtooth_true_q(dd::IMAS.DD, q_s::Real, dq::Real)
    cp1d = dd.core_profiles.profiles_1d[]
    eqt = dd.equilibrium.time_slice[]
    ρ, area = cp1d.grid.rho_tor_norm, cp1d.grid.area
    I = IMAS.cumtrapz(area, IMAS.Jpar_2_Jtor(ρ, cp1d.j_total, true, eqt))
    q = sawtooth_q_times_I(ρ, eqt) ./ abs.(I)
    q[1] = q[2]
    is = findlast(<(q_s), q)
    im = findnext(>=(q_s + dq), q, is === nothing ? 1 : is + 1)
    return (; ρ, area, I, q, im)
end

sawtooth_mixing_radius(r) = r.im === nothing || r.im == 1 ? 0.0 : r.ρ[r.im]

# reset shape p(x): p(0) = p'(0) = 0, p(1) = 1, p'(1) = m; the cubic is monotone only for m ≤ 3, so x^m beyond
sawtooth_shape(x, m) = m <= 3 ? (3 - 2x) * x^2 + m * (x^3 - x^2) : x^m

# ρ of the minimum q over nodes `r`, refined by the parabola through it and its neighbours, clamped to [lo, hi]
function sawtooth_dip_centre(ρ, q, r, lo, hi)
    k = r[argmin(q[r])]
    k == 1 && return 0.0
    d0, d2 = ρ[k-1] - ρ[k], ρ[k+1] - ρ[k]
    s0, s2 = (q[k-1] - q[k]) / d0, (q[k+1] - q[k]) / d2
    c = (s2 - s0) / (d2 - d0)
    return min(max(c > 0 ? ρ[k] - 0.5 * (s0 - c * d0) / c : ρ[k], lo), hi)
end

# (ρ, q) of the top of the hump at node `i`: vertex of the parabola through it and its neighbours
function sawtooth_hump_top(ρ, q, i)
    d0, d2 = ρ[i-1] - ρ[i], ρ[i+1] - ρ[i]
    s0, s2 = (q[i-1] - q[i]) / d0, (q[i+1] - q[i]) / d2
    c = (s2 - s0) / (d2 - d0)
    c >= 0 && return ρ[i], q[i]
    b = s0 - c * d0
    t = min(max(-0.5 * b / c, d0), d2)
    return ρ[i] + t, max(q[i] + b * t + c * t^2, q[i])
end

# two-sided reset: q_s at ρc, C1 onto (ρo, q_s + Δo, slope so) outside it and onto (ρi, q_s + Δi, slope si) inside it
function sawtooth_target(x, ρi, Δi, si, ρc, ρo, Δo, so, q_s)
    x >= ρc && return Δo > 0 && ρo > ρc ? q_s + Δo * sawtooth_shape((x - ρc) / (ρo - ρc), max((ρo - ρc) * so / Δo, 0.0)) : q_s
    return Δi > 0 ? q_s + Δi * sawtooth_shape((ρc - x) / (ρc - ρi), max(-(ρc - ρi) * si / Δi, 0.0)) : q_s
end

"""
    sawtooth_reset_q(ρ::AbstractVector, q::AbstractVector, q_s::Real, dq::Real; ramp::Real)

Target q of a sawtooth crash on nodes `ρ` (axis first), as OpenFUSIONToolkit's `saw_reset_1d` rule 2.

Each maximal run of nodes with q < q_s + dq that holds a q < q_s (and does not reach the last node) is a region between
its q_s + dq crossings ρ_in and ρ_out (ρ_in = 0 on the axis). Inside a region q is reset to q_s at its minimum and rises C1
onto q at both crossings; q outside the regions, and on the axis when it is outside them, is unchanged, so the current moves
only inside each region (Ip conserved). Each reset is weighted by `clamp((q_s - min q) / ramp, 0, 1)` (1 for `ramp <= 0`), so
a dip that appears is reset continuously. Dips separated by a hump blend into separate resets as the hump nears q_s + dq.

Returns `(; q_new, regions)`, one `(; ρ_in, ρ_c, ρ_out, q_min, w, s)` per region.
"""
function sawtooth_reset_q(ρ::AbstractVector, q::AbstractVector, q_s::Real, dq::Real; ramp::Real)
    n = length(q)
    q_new = copy(q)
    regions = NamedTuple{(:ρ_in, :ρ_c, :ρ_out, :q_min, :w, :s),NTuple{6,Float64}}[]
    qm = q_s + dq
    g = IMAS.gradient(ρ, q; method=:central)

    # reset weight of the nodes `r`
    weight(r) = ramp > 0 ? min(max((q_s - minimum(q[r])) / ramp, 0.0), 1.0) : 1.0

    # q_s + dq crossing between nodes k and k+1, and q' there
    function crossing(k)
        x = ρ[k] + (qm - q[k]) / (q[k+1] - q[k]) * (ρ[k+1] - ρ[k])
        return x, g[k] + (x - ρ[k]) / (ρ[k+1] - ρ[k]) * (g[k+1] - g[k])
    end

    # local maximum of q at node i; node 2 is none when q[1] == q[2] (axis value copied from the first surface)
    ismax(i) = i >= 2 && !(i == 2 && q[1] == q[2]) && q[i] >= q[i-1] && q[i] >= q[i+1]

    # edge hump: the highest local max in `maxima`, measured from its base max(q_s, min q over `span(ih)`)
    function edge_hump(maxima, span)
        ih = maxima[argmax(q[maxima])]
        ρh, h = sawtooth_hump_top(ρ, q, ih)
        base = max(q_s, minimum(q[span(ih)]))
        return ρh, h, qm > base ? min(max((h - base) / (qm - base), 0.0), 1.0) : 0.0
    end

    a = 1
    while a <= n
        if q[a] >= qm
            a += 1
            continue
        end
        # region: maximal run of nodes a0:b below q_s + dq
        a0 = a
        b = something(findnext(>=(qm), q, a0), n + 1) - 1
        r = a0:b
        a = b + 2
        (b == n || all(>=(q_s), q[r])) && continue

        ρo, so = crossing(b)
        if a0 == 1
            ρi, Δi, si = 0.0, max(q[1] - q_s, 0.0), 0.0
        else
            ρi, si = crossing(a0 - 1)
            Δi = dq
        end
        ρc = sawtooth_dip_centre(ρ, q, r, ρi, ρo)
        wm = weight(r)

        # runs of nodes below q_s inside the region
        d = filter(j -> q[j] < q_s, r)
        cut = findall(diff(d) .> 1)
        runs = [d[i1]:d[i2] for (i1, i2) in zip([1; cut .+ 1], [cut; length(d)])]

        # segments (0: run-free edge, j: run j) and the hump tops (ρ_h, h, s) between them
        segs, humps = Int[], NTuple{3,Float64}[]
        maxima = filter(ismax, a0:first(runs[1])-1)
        if !isempty(maxima)
            push!(segs, 0)
            push!(humps, edge_hump(maxima, ih -> a0:ih))
        end
        for j in eachindex(runs)
            push!(segs, j)
            j < length(runs) || continue
            ρh, h = sawtooth_hump_top(ρ, q, last(runs[j]) + argmax(q[last(runs[j])+1:first(runs[j+1])-1]))
            push!(humps, (ρh, h, min(max((h - q_s) / dq, 0.0), 1.0)))
        end
        maxima = filter(ismax, last(runs[end])+1:b)
        if !isempty(maxima)
            push!(segs, 0)
            push!(humps, edge_hump(maxima, ih -> ih:b))
        end

        ρh, sh, ns = first.(humps), last.(humps), length(segs)
        ends = [ρi; ρh; ρo]
        heights = [0.0; [max(h[2] - q_s, 0.0) for h in humps]; 0.0]
        # split target of each run: sawtooth_target's (ρi, Δi, si, ρc, ρo, Δo, so) and the run's weight
        run_reset = map(1:ns) do m
            segs[m] == 0 && return nothing
            nodes = runs[segs[m]]
            target = (ends[m], m == 1 ? Δi : heights[m], m == 1 ? si : 0.0, sawtooth_dip_centre(ρ, q, nodes, ends[m], ends[m+1]),
                ends[m+1], m == ns ? dq : heights[m+1], m == ns ? so : 0.0)
            return target, weight(nodes)
        end

        # each node: the whole-region reset, blended toward its run's own reset by the hump fraction s
        for k in r
            x = ρ[k]
            qt = sawtooth_target(x, ρi, Δi, si, ρc, ρo, dq, so, q_s)
            if isempty(humps)
                q_new[k] = q[k] + wm * (qt - q[k])
                continue
            end
            i = searchsortedlast(ρh, x)
            s = i == 0 ? sh[1] : i == length(ρh) ? sh[end] : sh[i] + (x - ρh[i]) * (sh[i+1] - sh[i]) / (ρh[i+1] - ρh[i])
            qsp, wk = run_reset[i+1] === nothing ? (q[k], 0.0) : (sawtooth_target(x, run_reset[i+1][1]..., q_s), run_reset[i+1][2])
            q_new[k] = q[k] + (1 - s) * wm * (qt - q[k]) + s * wk * (qsp - q[k])
        end
        push!(regions, (; ρ_in=ρi, ρ_c=ρc, ρ_out=ρo, q_min=minimum(q[r]), w=wm, s=isempty(sh) ? 0.0 : maximum(sh)))
    end
    return (; q_new, regions)
end

"""
    sawtooth_source_jpar!(dd::IMAS.DD, ρ0::Real; flat_factor::Real, exclude_indexes::AbstractVector{Int})

Rebuild the j_parallel of the sawteeth source (701) at this time as `IMAS.sawteeth_source!` builds it for inversion radius `ρ0`,
and move `cp1d.j_non_inductive` by the same change, so that j_ohmic does not carry the flattening.
"""
function sawtooth_source_jpar!(dd::IMAS.DD, ρ0::Real; flat_factor::Real, exclude_indexes::AbstractVector{Int})
    cp1d = dd.core_profiles.profiles_1d[]
    ρ, area = cp1d.grid.rho_tor_norm, cp1d.grid.area
    src = findfirst(:sawteeth, dd.core_sources.source)
    @assert src !== nothing "sawtooth_source_jpar!: no sawteeth source"
    s1d = src.profiles_1d[]
    @assert s1d.time == dd.global_time && length(s1d.grid.rho_tor_norm) == length(ρ) "sawtooth_source_jpar!: sawteeth source at t=$(s1d.time) is not at t=$(dd.global_time) on the same grid"
    j = IMAS.total_sources(dd.core_sources, cp1d; time0=dd.global_time, exclude_indexes, fields=[:j_parallel]).j_parallel
    jsaw = ρ0 > 0 ? (IMAS.flatten_profile!(copy(j), ρ, area, ρ0, min(ρ0 / 4, 0.05)) .- j) .* flat_factor : zero(j)
    old = ismissing(s1d, :j_parallel) ? zero(j) : s1d.j_parallel
    s1d.j_parallel = jsaw
    setproperty!(s1d, IMAS._core_sources_value_keys[:j_parallel], IMAS.cumtrapz(s1d.grid.area, jsaw))
    # set, not refrozen: the j_non_inductive expression reads 701 at cp1d.time, which may not be this slice
    cp1d.j_non_inductive = cp1d.j_non_inductive .+ (jsaw .- old)
    return jsaw
end

"""
    qed_sawtooth_crash!(actor::ActorQED)

After the current diffusion of a step ending at `dd.global_time`: on the first step past each multiple of `crash_period`,
resets q per [`sawtooth_reset_q`](@ref), adding the change of current to `cp1d.j_total` (so j_ohmic carries it). Then
rebuilds the sawteeth source current at the mixing radius and stores that radius in `actor.saw_radius`.
"""
function qed_sawtooth_crash!(actor::ActorQED)
    dd = actor.dd
    par = actor.par
    cp1d = dd.core_profiles.profiles_1d[]
    t = dd.global_time
    q_s = par.crash_q > 0 ? par.crash_q : par.qmin_desired
    r = sawtooth_true_q(dd, q_s, par.crash_dq)
    ρm = sawtooth_mixing_radius(r)
    (; q_new, regions) = sawtooth_reset_q(r.ρ, r.q, q_s, par.crash_dq; ramp=par.crash_ramp)
    if !isempty(regions) && floor((t + 1e-9) / par.crash_period) > floor((t - par.Δt + 1e-9) / par.crash_period)
        Δj = IMAS.Jtor_2_Jpar(r.ρ, IMAS.gradient(r.area, r.I .* (r.q ./ q_new .- 1)), false, dd.equilibrium.time_slice[])
        cp1d.j_total = cp1d.j_total .+ Δj
        cp1d.q = q_new
    end
    saw = actor.act.ActorSawteethSource
    sawtooth_source_jpar!(dd, ρm; saw.flat_factor, exclude_indexes=sawteeth_exclude_indexes(saw))
    actor.saw_radius = ρm
    return actor
end
