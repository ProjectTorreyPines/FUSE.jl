#= ============================================================================ =#
#  fuse_efit inputs: what the host must compute itself                            #
#                                                                                #
#  The model repo's docs/INPUT_CONTRACT.md (branch `fuse-efit`) is the authority #
#  for every definition below; section numbers refer to it.                      #
#                                                                                #
#  The model was trained on EFIT01 (the between-shot, magnetics-only DIII-D      #
#  reconstruction). FUSE feeds the same quantities evaluated on whichever        #
#  equilibrium `nn_state_source` selects:                                        #
#    :simulation  FUSE's own evolving equilibrium (FRESCO, or the flux map       #
#                 GSLite sends over ZMQ) — self-consistent and predictive, but   #
#                 not EFIT01, which is where an accuracy loss against fuse29     #
#                 would come from (the authors measured +31 % loss when feeding  #
#                 real-time EFIT instead of EFIT01, mostly through the Shafranov #
#                 shift, the squareness and the strike points; §8)               #
#    :replay      the experimental equilibrium kept in `act.ActorReplay.replay_dd`#
#                 — isolates the model's own error from FUSE's equilibrium error #
#= ============================================================================ =#

const FUSE_EFIT_DRSEP_CLIP = 0.1         # m, §5: drsep is clipped to ±0.1 m per slice
const FUSE_EFIT_DRSEP_SATURATION = 0.4   # m, §5: EFIT saturates drsep at ±0.4 m when no second separatrix is in range
const FUSE_EFIT_MIN_BOUNDARY_POINTS = 10 # §4.3: a boundary slice counts only with at least 10 points
# An X-point whose normalised flux is within this of 1 lies on the boundary: it is the boundary's lowest
# (highest) point, which the traced outline rounds off.
const FUSE_EFIT_ACTIVE_XPOINT_DPSIN = 2e-3
# COCOS 11: ψ is in Wb and rises from the axis to the boundary for a positive (counter-clockwise) plasma
# current; the ohmic transformer that sustains a positive Ip raises the boundary flux, so the loop voltage
# is +dψ_boundary/dt. Checked on DIII-D 199055 (Ip > 0: ψ_boundary rises by 0.1 Wb/s in flat-top, and
# IMAS.vloop is positive) against DIII-D's VLOOP, which has the sign of Ip (training median +0.28 V at
# Ip median +0.98 MA; reference shot 193190: +0.4 V at +1.2 MA).
const FUSE_EFIT_VLOOP_SIGN_DPSI_DT = 1.0

# ── Geometry helpers ────────────────────────────────────────────────────────

# Luce squareness of one quadrant (§4.3; T.C. Luce, PPCF 55 (2013) 095009), same construction as
# `IMAS.luce_squareness` and the model repo's reference `external/efit01_shape.py`: O is the "ellipse
# centre" of the quadrant, E its rectangle corner, D where the ray O→E crosses the (closed) boundary
# polygon, C where it crosses the inscribed ellipse; zeta = (|OD| − |OC|) / (|OE| − |OC|).
# NaN when the quadrant is undefined (the contract's "node hole"), where IMAS returns 0.0.
function _fuse_efit_luce(pr::AbstractVector{Float64}, pz::AbstractVector{Float64}, PO::NTuple{2,Float64}, PE::NTuple{2,Float64})
    OE = hypot(PE[1] - PO[1], PE[2] - PO[2])
    OE > 0.0 || return NaN
    crossings = IMAS.intersection([PO[1], PE[1]], [PO[2], PE[2]], pr, pz).crossings
    isempty(crossings) && return NaN
    OD = minimum(hypot(c[1] - PO[1], c[2] - PO[2]) for c in crossings)
    OC = OE / sqrt(2.0)
    return (OD - OC) / (OE - OC)
end

# First intersection of the segment P→Q with a closed polygon (wall); `nothing` if there is none.
function _fuse_efit_segment_hit(p::NTuple{2,Float64}, q::NTuple{2,Float64}, wr::AbstractVector{Float64}, wz::AbstractVector{Float64})
    best_t = Inf
    dx, dy = q[1] - p[1], q[2] - p[2]
    n = length(wr)
    @inbounds for k in 1:n
        k2 = k == n ? 1 : k + 1
        ex, ey = wr[k2] - wr[k], wz[k2] - wz[k]
        den = dx * ey - dy * ex
        abs(den) < 1e-14 && continue
        t = ((wr[k] - p[1]) * ey - (wz[k] - p[2]) * ex) / den
        s = ((wr[k] - p[1]) * dy - (wz[k] - p[2]) * dx) / den
        if 0.0 <= t <= 1.0 && 0.0 <= s <= 1.0 && t < best_t
            best_t = t
        end
    end
    isfinite(best_t) || return nothing
    return (p[1] + best_t * dx, p[2] + best_t * dy)
end

# Saddle points of ψ inside `in_domain`, located on IMAS's cubic interpolant of ψ(R,Z): grid minima of |∇ψ|²
# refined by Newton iterations on ∇ψ = 0, kept when the Hessian is indefinite (a saddle, not the axis).
function _fuse_efit_saddles(r::AbstractRange, z::AbstractRange, PSI, in_domain::Function)
    nr, nz = length(r), length(z)
    grad(x, y) = IMAS.FI.gradient(PSI, (x, y))
    h = 1e-3 * min(Base.step(r), Base.step(z))   # the Hessian by central differences of the analytic gradient
    g2 = Matrix{Float64}(undef, nr, nz)
    @inbounds for j in 1:nz, i in 1:nr
        g = grad(Float64(r[i]), Float64(z[j]))
        g2[i, j] = g[1]^2 + g[2]^2
    end
    saddles = NTuple{2,Float64}[]
    @inbounds for j in 2:nz-1, i in 2:nr-1
        v = g2[i, j]
        is_min = true
        for dj in -1:1, di in -1:1
            (di == 0 && dj == 0) && continue
            if g2[i+di, j+dj] < v
                is_min = false
                break
            end
        end
        is_min || continue
        x, y = Float64(r[i]), Float64(z[j])
        ok = false
        for _ in 1:30
            g = grad(x, y)
            gxp, gxm, gyp, gym = grad(x + h, y), grad(x - h, y), grad(x, y + h), grad(x, y - h)
            Hxx = (gxp[1] - gxm[1]) / 2h
            Hyy = (gyp[2] - gym[2]) / 2h
            Hxy = ((gyp[1] - gym[1]) + (gxp[2] - gxm[2])) / 4h
            det = Hxx * Hyy - Hxy^2
            abs(det) < 1e-30 && break
            dx = -(Hyy * g[1] - Hxy * g[2]) / det
            dy = -(-Hxy * g[1] + Hxx * g[2]) / det
            x += dx
            y += dy
            (r[1] <= x <= r[end] && z[1] <= y <= z[end]) || break
            if hypot(dx, dy) < 1e-8
                ok = det < 0.0
                break
            end
        end
        ok || continue
        # the refinement must not have wandered off to another critical point
        hypot(x - r[i], y - z[j]) <= 2.0 * hypot(Base.step(r), Base.step(z)) || continue
        in_domain(x, y) || continue
        any(hypot(x - s[1], y - s[2]) < 1e-3 for s in saddles) && continue
        push!(saddles, (x, y))
    end
    return saddles
end

# Major radius where normalised flux `psin_target` is first reached going outward along the midplane
# z = ZA from the magnetic axis; `nothing` if it is not reached inside the grid.
function _fuse_efit_midplane_radius(r::AbstractRange, PSI, RA::Float64, ZA::Float64, psi_axis::Float64, psi_bnd::Float64, psin_target::Float64)
    n = 2000
    R1 = Float64(r[end])
    psin(R) = (PSI(R, ZA) - psi_axis) / (psi_bnd - psi_axis)
    R_prev, f_prev = RA, psin(RA)
    for k in 1:n
        R = RA + (R1 - RA) * k / n
        f = psin(R)
        if f >= psin_target
            f == f_prev && return R
            return R_prev + (psin_target - f_prev) / (f - f_prev) * (R - R_prev)
        end
        R_prev, f_prev = R, f
    end
    return nothing
end

# Wall hits of the flux surface `psi_level` (just outside the separatrix of X-point `X`) on the X-point's
# side of the midplane: the strike points of that X-point's divertor legs. Each contour that passes close
# to X is walked from its vertex nearest to X in both directions until it first crosses the wall.
function _fuse_efit_strike_hits(eqt::IMAS.equilibrium__time_slice, psi_level::Float64, X::NTuple{2,Float64}, ZA::Float64,
                                wr::AbstractVector{Float64}, wz::AbstractVector{Float64}; near::Float64)
    hits = NTuple{2,Float64}[]
    surfaces = IMAS.flux_surface(eqt, psi_level, :any, Float64[], Float64[])
    for (cr, cz) in surfaces
        length(cr) < 2 && continue
        d2, k0 = findmin((cr .- X[1]) .^ 2 .+ (cz .- X[2]) .^ 2)
        sqrt(d2) <= near || continue
        for dir in (1, -1)
            k = k0
            while 1 <= k + dir <= length(cr)
                hit = _fuse_efit_segment_hit((Float64(cr[k]), Float64(cz[k])), (Float64(cr[k+dir]), Float64(cz[k+dir])), wr, wz)
                if hit !== nothing
                    if sign(hit[2] - ZA) == sign(X[2] - ZA) && !any(hypot(hit[1] - h[1], hit[2] - h[2]) < 1e-3 for h in hits)
                        push!(hits, hit)
                    end
                    break
                end
                k += dir
            end
        end
    end
    return hits
end

"""
    fuse_efit_shape(eqt::IMAS.equilibrium__time_slice, wall_r, wall_z) -> Union{Nothing,NamedTuple}

The 14 boundary-shape inputs of fuse_efit and the 3 strike-point presence flags
for ONE equilibrium slice, with the definitions of the model's
`docs/INPUT_CONTRACT.md` (§2.3, §4, §5). `nothing` when the slice cannot be
used (no boundary outline of at least 10 points, no ψ(R,Z) map, non-positive
minor or major radius: §3.2 rule 1, "failed reconstructions are excluded").

| field | definition |
|---|---|
| `aminor`, `rsurf`, `zsurf` | half the radial extent, and the geometric centre, of the boundary [m] |
| `kappa` | boundary height / width |
| `tritop`, `tribot` | (rsurf − R at the highest / lowest boundary point) / aminor |
| | the highest / lowest boundary point is the X-point when one lies on the boundary (ψ_N within 2e-3 of 1, no further than 0.2 `aminor` from the outline), as in EFIT; this enters `zsurf`, `kappa`, `tritop`, `tribot` |
| `shafranov` | magnetic-axis R − `rsurf` [m] (§4.2) |
| `zeta_ou`, `zeta_ol` | Luce squareness of the outer-upper / outer-lower quadrant from the discrete boundary vertices, no refinement (§4.3); NaN when undefined |
| `drsep` | outer-midplane distance between the separatrices of the lower and the upper X-point [m]: positive for upper single null (and double null), negative for lower single null; ±0.4 m when one of the two X-points does not exist; **clipped to ±0.1 m** (§5). NaN when there is no X-point at all |
| `gapin` | radial gap between the innermost boundary point and the wall [m]; NaN without a wall |
| `rvsout` | major radius of the lower outer strike point [m] |
| `rvsou`, `rvsiu` | major radius of the upper outer / upper inner strike point [m] |
| `rvsout_present`, `rvsou_present`, `rvsiu_present` | whether that strike point exists; when `false` its radius is NaN — never a zero or a placeholder radius (§5, §6.3) |

Also returned for diagnostics: `drsep_raw` (before the clip) and the X-points
`x_lower`, `x_upper` (`(r, z, psi_norm)` or `nothing`).

X-points are the saddle points of ψ inside the wall (inside a box around the
plasma when there is no wall). The strike points of an X-point are where the
flux surface just outside its separatrix first meets the wall on that X-point's
side of the midplane; the X-point that is not the active one contributes strike
points only when its separatrix is within EFIT's ±0.4 m range at the midplane,
which is when EFIT01 reports them.

The squareness carries a resolution bias of order 0.02 that depends on the
number of boundary vertices (by design, §4.3); EFIT01 boundaries have a median
of 92 points.
"""
function fuse_efit_shape(eqt::IMAS.equilibrium__time_slice, wall_r::AbstractVector{<:Real}=Float64[], wall_z::AbstractVector{<:Real}=Float64[])
    (ismissing(eqt.boundary.outline, :r) || ismissing(eqt.boundary.outline, :z)) && return nothing
    pr = Float64.(eqt.boundary.outline.r)
    pz = Float64.(eqt.boundary.outline.z)
    (length(pr) >= FUSE_EFIT_MIN_BOUNDARY_POINTS && length(pr) == length(pz)) || return nothing
    if pr[1] != pr[end] || pz[1] != pz[end]   # the polygon must be closed (§4.3)
        push!(pr, pr[1])
        push!(pz, pz[1])
    end

    (_, _, _, _, r_at_max_z, max_z, r_at_min_z, min_z, z_at_max_r, max_r, z_at_min_r, min_r) = IMAS.fluxsurface_extrema(pr, pz)
    aminor = 0.5 * (max_r - min_r)
    rsurf = 0.5 * (max_r + min_r)
    zsurf = 0.5 * (max_z + min_z)
    (aminor > 0.0 && rsurf > 0.0 && isfinite(aminor) && isfinite(rsurf)) || return nothing
    kappa = 0.5 * (max_z - min_z) / aminor
    tritop = (rsurf - r_at_max_z) / aminor
    tribot = (rsurf - r_at_min_z) / aminor
    # the squareness is defined on the outline's own vertices (§4.3), so it is taken before the X-point
    # correction below: restoring the X-point moves it by 0.002 on DIII-D 193190, a tenth of its resolution bias
    zeta_ou = _fuse_efit_luce(pr, pz, (r_at_max_z, z_at_max_r), (max_r, max_z))
    zeta_ol = _fuse_efit_luce(pr, pz, (r_at_min_z, z_at_max_r), (max_r, min_z))

    RA = Float64(eqt.global_quantities.magnetic_axis.r)
    ZA = Float64(eqt.global_quantities.magnetic_axis.z)
    shafranov = RA - rsurf

    wr = Float64.(wall_r)
    wz = Float64.(wall_z)
    has_wall = length(wr) >= 3 && length(wr) == length(wz)

    # inner gap: from the innermost boundary point radially inward to the wall
    gapin = NaN
    if has_wall
        hit = _fuse_efit_segment_hit((min_r, z_at_min_r), (0.0, z_at_min_r), wr, wz)
        hit === nothing || (gapin = max(min_r - hit[1], 0.0))
    end

    drsep_raw = NaN
    x_lower = nothing
    x_upper = nothing
    rvsout = rvsou = rvsiu = NaN
    rvsout_present = rvsou_present = rvsiu_present = false

    eqt2d = findfirst(:rectangular, eqt.profiles_2d)
    if eqt2d !== nothing && !ismissing(eqt2d, :psi)
        r, z, PSI = IMAS.ψ_interpolant(eqt2d)
        psi_axis = Float64(PSI(RA, ZA))
        psi_bnd = Float64(eqt.profiles_1d.psi[end])
        if psi_bnd != psi_axis
            psin_of(x, y) = (PSI(x, y) - psi_axis) / (psi_bnd - psi_axis)
            if has_wall
                wall_polygon = collect(zip(wr, wz))
                in_domain = (x, y) -> PolygonOps.inpolygon((x, y), wall_polygon) == 1
            else
                z_mid, half_height = zsurf, kappa * aminor   # of the outline as traced
                in_domain = (x, y) -> abs(x - rsurf) <= 1.5 * aminor && abs(y - z_mid) <= 2.0 * half_height
            end
            saddles = [(s[1], s[2], psin_of(s[1], s[2])) for s in _fuse_efit_saddles(r, z, PSI, in_domain)]
            filter!(s -> s[3] > 0.9, saddles)
            lower = filter(s -> s[2] < ZA, saddles)
            upper = filter(s -> s[2] > ZA, saddles)
            isempty(lower) || (x_lower = lower[argmin([s[3] for s in lower])])
            isempty(upper) || (x_upper = upper[argmin([s[3] for s in upper])])

            # EFIT's ZSURF, KAPPA, TRITOP and TRIBOT take the X-point as the boundary's vertical extremum. The
            # outline IMAS traces stops just inside the separatrix and rounds that corner off: on DIII-D 193190
            # its lowest vertex sits 7 mm above and 28 mm outboard of the X-point, which lowers `tribot` by
            # 0.03-0.05. With the X-point restored the three scalars agree with EFIT01's to 3e-5 (median).
            # (an X-point further than 0.2 a from the outline's own extremum does not belong to this outline)
            on_boundary(X, r_ext, z_ext) = abs(X[3] - 1.0) <= FUSE_EFIT_ACTIVE_XPOINT_DPSIN && hypot(X[1] - r_ext, X[2] - z_ext) <= 0.2 * aminor
            if x_lower !== nothing && x_lower[2] < min_z && on_boundary(x_lower, r_at_min_z, min_z)
                r_at_min_z, min_z = x_lower[1], x_lower[2]
            end
            if x_upper !== nothing && x_upper[2] > max_z && on_boundary(x_upper, r_at_max_z, max_z)
                r_at_max_z, max_z = x_upper[1], x_upper[2]
            end
            zsurf = 0.5 * (max_z + min_z)
            kappa = 0.5 * (max_z - min_z) / aminor
            tritop = (rsurf - r_at_max_z) / aminor
            tribot = (rsurf - r_at_min_z) / aminor

            R_lower = x_lower === nothing ? nothing : _fuse_efit_midplane_radius(r, PSI, RA, ZA, psi_axis, psi_bnd, x_lower[3])
            R_upper = x_upper === nothing ? nothing : _fuse_efit_midplane_radius(r, PSI, RA, ZA, psi_axis, psi_bnd, x_upper[3])
            if R_lower !== nothing && R_upper !== nothing
                drsep_raw = clamp(R_lower - R_upper, -FUSE_EFIT_DRSEP_SATURATION, FUSE_EFIT_DRSEP_SATURATION)
            elseif R_lower !== nothing
                drsep_raw = -FUSE_EFIT_DRSEP_SATURATION   # lower single null, no second separatrix in range
            elseif R_upper !== nothing
                drsep_raw = FUSE_EFIT_DRSEP_SATURATION    # upper single null, no second separatrix in range
            end

            if has_wall && isfinite(drsep_raw)
                dpsin = 2e-3   # just outside the separatrix, in normalised flux (~0.5 mm at the midplane)
                # an X-point's legs exist for EFIT when its separatrix is the boundary or within range of it
                lower_in_range = R_lower !== nothing && (drsep_raw <= 0.0 || abs(drsep_raw) < FUSE_EFIT_DRSEP_SATURATION)
                upper_in_range = R_upper !== nothing && (drsep_raw >= 0.0 || abs(drsep_raw) < FUSE_EFIT_DRSEP_SATURATION)
                near = 0.5 * aminor
                if lower_in_range
                    level = psi_axis + (x_lower[3] + dpsin) * (psi_bnd - psi_axis)
                    hits = _fuse_efit_strike_hits(eqt, level, (x_lower[1], x_lower[2]), ZA, wr, wz; near)
                    outer = filter(h -> h[1] >= x_lower[1], hits)
                    if length(hits) >= 2
                        rvsout = maximum(h[1] for h in hits)
                        rvsout_present = true
                    elseif !isempty(outer)
                        rvsout = outer[1][1]
                        rvsout_present = true
                    end
                end
                if upper_in_range
                    level = psi_axis + (x_upper[3] + dpsin) * (psi_bnd - psi_axis)
                    hits = _fuse_efit_strike_hits(eqt, level, (x_upper[1], x_upper[2]), ZA, wr, wz; near)
                    if length(hits) >= 2
                        rvsou = maximum(h[1] for h in hits)
                        rvsiu = minimum(h[1] for h in hits)
                        rvsou_present = rvsiu_present = true
                    elseif length(hits) == 1
                        if hits[1][1] >= x_upper[1]
                            rvsou = hits[1][1]
                            rvsou_present = true
                        else
                            rvsiu = hits[1][1]
                            rvsiu_present = true
                        end
                    end
                end
            end
        end
    end

    drsep = isfinite(drsep_raw) ? clamp(drsep_raw, -FUSE_EFIT_DRSEP_CLIP, FUSE_EFIT_DRSEP_CLIP) : NaN

    return (aminor=aminor, rsurf=rsurf, zsurf=zsurf, kappa=kappa, tritop=tritop, tribot=tribot, drsep=drsep,
            rvsout=rvsout, gapin=gapin, rvsou=rvsou, rvsiu=rvsiu, shafranov=shafranov, zeta_ou=zeta_ou, zeta_ol=zeta_ol,
            rvsout_present=rvsout_present, rvsou_present=rvsou_present, rvsiu_present=rvsiu_present,
            drsep_raw=drsep_raw, x_lower=x_lower, x_upper=x_upper)
end

# The limiter outline exactly as stored: IMAS.first_wall's default simplification to an inscribed polygon
# moves vertices (and collapses a coarse outline), which shifts strike points and the inner gap.
function _fuse_efit_wall(dd::IMAS.DD)
    try
        fw = IMAS.first_wall(dd.wall; simplify_to_inscribed_fractional_area=1.0)
        return Float64.(fw.r), Float64.(fw.z)
    catch
        return Float64[], Float64[]
    end
end

# fuse_efit_shape with the result of immutable (replay) slices cached in dd._aux
function _fuse_efit_shape_cached(aux, eqt::IMAS.equilibrium__time_slice, wr, wz; cache::Bool)
    cache || return _fuse_efit_shape_or_nothing(eqt, wr, wz)
    store = get!(aux, :fuse_efit_shape_cache) do
        IdDict{Any,Any}()
    end
    return get!(store, eqt) do
        _fuse_efit_shape_or_nothing(eqt, wr, wz)
    end
end

function _fuse_efit_shape_or_nothing(eqt::IMAS.equilibrium__time_slice, wr, wz)
    try
        return fuse_efit_shape(eqt, wr, wz)
    catch e
        e isa InterruptException && rethrow()
        @warn "fuse_efit: could not evaluate the boundary shape of the equilibrium at t=$(eqt.time) s" exception = e maxlog = 3
        return nothing
    end
end

# Mean over slices of each shape channel, finite values only; a strike point counts as present when at
# least one slice has it (§3.2, §5).
function _fuse_efit_mean_shape(shapes::AbstractVector)
    isempty(shapes) && return nothing
    function mean_finite(field::Symbol)
        acc = 0.0
        n = 0
        for s in shapes
            v = getproperty(s, field)
            if isfinite(v)
                acc += v
                n += 1
            end
        end
        return n == 0 ? NaN : acc / n
    end
    vals = NamedTuple{Symbol.(FUSE_EFIT_SHAPE_NAMES)}(Tuple(mean_finite(Symbol(name)) for name in FUSE_EFIT_SHAPE_NAMES))
    return (; vals...,
            rvsout_present=any(s.rvsout_present for s in shapes),
            rvsou_present=any(s.rvsou_present for s in shapes),
            rvsiu_present=any(s.rvsiu_present for s in shapes))
end

# ── Plasma-state inputs from an equilibrium source ──────────────────────────

# Shape, plasma current [A] and loop voltage [V] of the experimental (replay) equilibrium for the tick
# (time − period, time].
function _fuse_efit_state_replay(aux, replay_dd::IMAS.DD, wr, wz, time::Float64, period::Float64)
    eq = replay_dd.equilibrium
    times = [Float64(eqt.time) for eqt in eq.time_slice]
    # §3.2: mean of the slices in (t − 50 ms, t], left edge open, right edge closed, 1 µs tolerance; no hold-last
    in_window = [k for k in eachindex(times) if time - period + 1e-6 < times[k] <= time + 1e-6]
    shapes = Any[]
    ips = Float64[]
    for k in in_window
        eqt = eq.time_slice[k]
        s = _fuse_efit_shape_cached(aux, eqt, wr, wz; cache=true)
        s === nothing || push!(shapes, s)
        ismissing(eqt.global_quantities, :ip) || push!(ips, Float64(eqt.global_quantities.ip))
    end
    shape = _fuse_efit_mean_shape(shapes)

    # plasma current: window mean, else linear interpolation between the neighbouring slices
    ip = nothing
    if !isempty(ips)
        ip = sum(ips) / length(ips)
    else
        ks = [k for k in eachindex(times) if !ismissing(eq.time_slice[k].global_quantities, :ip)]
        if length(ks) >= 2 && times[ks[1]] <= time <= times[ks[end]]
            ip = IMAS.interp1d(times[ks], [Float64(eq.time_slice[k].global_quantities.ip) for k in ks], :linear)(time)
        end
    end

    # loop voltage: there is no measured VLOOP in the dd, so it is the rate of change of the boundary flux,
    # as a causal least-squares slope over the slices of the last three ticks (a two-point difference of a
    # magnetics-only reconstruction is noisy at 50 ms)
    vloop = nothing
    ks = [k for k in eachindex(times) if time - 3 * period + 1e-6 < times[k] <= time + 1e-6 &&
                                          !ismissing(eq.time_slice[k].global_quantities, :psi_boundary)]
    if length(ks) >= 2
        tt = times[ks]
        pp = [Float64(eq.time_slice[k].global_quantities.psi_boundary) for k in ks]
        tm = sum(tt) / length(tt)
        pm = sum(pp) / length(pp)
        den = sum((tt .- tm) .^ 2)
        den > 0.0 && (vloop = FUSE_EFIT_VLOOP_SIGN_DPSI_DT * sum((tt .- tm) .* (pp .- pm)) / den)
    end
    return (shape=shape, ip=ip, vloop=vloop)
end

# Shape, plasma current [A] and loop voltage [V] of FUSE's own equilibrium and current diffusion at `time`.
function _fuse_efit_state_simulation(aux, dd::IMAS.DD, wr, wz, time::Float64, period::Float64, max_age::Float64)
    eq = dd.equilibrium
    isempty(eq.time_slice) && return (shape=nothing, ip=nothing, vloop=nothing)
    times = [Float64(eqt.time) for eqt in eq.time_slice]
    k = findlast(τ -> τ <= time + 1e-6, times)
    (k === nothing || time - times[k] > max_age + 1e-6) && return (shape=nothing, ip=nothing, vloop=nothing)
    eqt = eq.time_slice[k]
    shape1 = _fuse_efit_shape_cached(aux, eqt, wr, wz; cache=false)
    shape = shape1 === nothing ? nothing : _fuse_efit_mean_shape([shape1])
    ip = ismissing(eqt.global_quantities, :ip) ? nothing : Float64(eqt.global_quantities.ip)

    # Plasma current and loop voltage from the current-diffusion solution, which advances every step. The
    # equilibrium does not: when its solver fails to converge the previous slice is kept under the new time
    # (ActorEquilibrium), and its `ip` then lags the ramp (0.25 MA held for 0.2 s on DIII-D 199055).
    # The loop voltage is the local Ohm's law integrated over the cross-section, the estimate ActorZMQ sends
    # to GSLite; it needs only the present state.
    vloop = nothing
    cp = dd.core_profiles.profiles_1d
    if !isempty(cp)
        kc = findlast(c -> c.time <= time + 1e-6, cp)
        if kc !== nothing && time - cp[kc].time <= max_age + 1e-6
            ip_cp = try
                Float64(IMAS.Ip(cp[kc], eqt))
            catch e
                e isa InterruptException && rethrow()
                nothing
            end
            (ip_cp === nothing || !isfinite(ip_cp)) || (ip = ip_cp)
            vloop = try
                Float64(IMAS.vloop(cp[kc], eqt))
            catch e
                e isa InterruptException && rethrow()
                nothing
            end
        end
    end
    # fallback: rate of change of the boundary flux between the two latest slices
    if (vloop === nothing || !isfinite(vloop)) && k >= 2 && times[k] > times[k-1] &&
       !ismissing(eqt.global_quantities, :psi_boundary) && !ismissing(eq.time_slice[k-1].global_quantities, :psi_boundary)
        vloop = FUSE_EFIT_VLOOP_SIGN_DPSI_DT * (eqt.global_quantities.psi_boundary - eq.time_slice[k-1].global_quantities.psi_boundary) / (times[k] - times[k-1])
    end
    return (shape=shape, ip=ip, vloop=vloop)
end

"""
    build_fuse_efit_inputs(nn::FuseEfitNN, dd::IMAS.DD; source=:zmq, state_source=:simulation,
                           time=dd.global_time, replay_dd=nothing, max_age=2*nn.period_s)
        -> (x::Vector{Float32}, indicators::Vector{Float32}, live::BitVector)

Assemble one fuse_efit tick — the 25 inputs in **physical units** and
`FUSE_EFIT_INPUT_NAMES` order, and the 5 per-tick indicators — for the 50 ms
interval ending at `time`. `live[j]` says whether channel `j` got a value; a
channel without one is left at NaN, which the graph ignores when the channel's
indicator or availability flag is 0 and otherwise reads as the training mean.

Actuators (`source`, as for fuse29; see [`_pedestal_nn_common_actuators`](@ref)):

| channel | units | `:zmq` | `:dd` |
|---|---|---|---|
| `pinj` | MW | `dd._aux[:zmq_Pnbi]` / 1e6 | Σ `dd.pulse_schedule.nbi.unit[].power.reference` / 1e6; 0 when the dd has no beams |
| `tinj` | N·m | Σ `dd._aux[:zmq_tinj_per_beam]` | NBI `core_sources` torque, else the beam-geometry estimate; 0 when the dd has no beams |
| `ech_mw` | MW | `dd._aux[:zmq_Pech]` / 1e6 | Σ `dd.ec_launchers.beam[].power_launched` / 1e6; 0 when the dd has no launchers |
| `gasa_cal..gase_cal` | Torr·L/s | `dd._aux[:zmq_gas[a-e]_cal]` | `dd.gas_injection.valve[]` `GASA..GASE` |
| `bt` | T (signed) | `dd.equilibrium.vacuum_toroidal_field.b0` | same |

Plasma state (`state_source`):

| channel | `:simulation` | `:replay` |
|---|---|---|
| `ip` [MA] | `dd._aux[:zmq_Ip_avg]` when `source = :zmq` (the measured current GSLite streams), else `IMAS.Ip(cp1d, eqt)` from the current-diffusion solution (`global_quantities.ip` of FUSE's equilibrium slice without current profiles) | `global_quantities.ip` of the experimental slices in the tick window |
| 14 shape channels | [`fuse_efit_shape`](@ref) of FUSE's equilibrium slice at `time` (at most `max_age` old) | mean of [`fuse_efit_shape`](@ref) over the experimental slices in (time − 50 ms, time] |
| `vloop` [V] | `IMAS.vloop(cp1d, eqt)` (Ohm's law on the current-diffusion solution) | dψ_boundary/dt of the experimental equilibrium |

Indicators: `efit_valid` is 0 when the state source has no usable slice for the
tick (§3.2: no hold-last), `vloop_valid` is 0 when no loop voltage could be
computed, each `*_present` is 0 when the equilibrium has no such strike point.

`replay_dd` is the experimental dd kept by `init!` (`act.ActorReplay.replay_dd`);
`state_source = :replay` without one yields `efit_valid = vloop_valid = 0`.
"""
function build_fuse_efit_inputs(nn::FuseEfitNN, dd::IMAS.DD; source::Symbol=:zmq, state_source::Symbol=:simulation,
                                time::Float64=dd.global_time, replay_dd::Union{Nothing,IMAS.DD}=nothing,
                                max_age::Float64=2 * nn.period_s)
    state_source in (:simulation, :replay) || error("build_fuse_efit_inputs: state_source must be :simulation or :replay, got $(repr(state_source))")
    x = fill(NaN32, FUSE_EFIT_N_INPUTS)
    indicators = zeros(Float32, FUSE_EFIT_N_INDICATORS)
    live = falses(FUSE_EFIT_N_INPUTS)
    aux = getfield(dd, :_aux)
    period = nn.period_s

    function _set!(name::AbstractString, value)
        value === nothing && return
        isfinite(value) || return
        idx = findfirst(==(name), nn.input_names)
        idx === nothing && error("build_fuse_efit_inputs: channel \"$name\" is not in the bundle contract")
        x[idx] = Float32(value)
        live[idx] = true
        return
    end
    _flag!(name::AbstractString, on::Bool) = (indicators[findfirst(==(name), nn.indicator_names)] = on ? 1.0f0 : 0.0f0)

    # -- heating, torque, fuelling, field --
    common = _pedestal_nn_common_actuators(dd; source, time, period)
    for (name, value) in common
        _set!(name == "ech" ? "ech_mw" : name, value)
    end
    if source == :dd
        # A dd that carries no beams / no gyrotrons describes a shot without them: 0 MW is the measurement
        # (training has "ECH shown off at 0 MW" as a record, §6.1), not a missing channel.
        no_beams = isempty(dd.pulse_schedule.nbi.unit) && isempty(dd.nbi.unit)
        haskey(common, "pinj") || !no_beams || _set!("pinj", 0.0)
        haskey(common, "tinj") || !no_beams || _set!("tinj", 0.0)
        haskey(common, "ech") || !isempty(dd.ec_launchers.beam) || _set!("ech_mw", 0.0)
    end

    # -- plasma state: shape, current, loop voltage --
    wr, wz = _fuse_efit_wall(dd)
    if state_source == :replay
        state = replay_dd === nothing ? (shape=nothing, ip=nothing, vloop=nothing) :
                _fuse_efit_state_replay(aux, replay_dd, wr, wz, time, period)
    else
        state = _fuse_efit_state_simulation(aux, dd, wr, wz, time, period, max_age)
    end

    ip = state.ip
    if source == :zmq && state_source == :simulation
        v = _nn_aux_value_at(aux, :zmq_Ip_avg, time)   # raw PCS measurement [A], machine sign convention
        v === nothing || (ip = v)
    end
    ip === nothing || _set!("ip", ip / 1e6)

    shape = state.shape
    if shape !== nothing
        _flag!("efit_valid", true)
        for name in FUSE_EFIT_SHAPE_NAMES
            _set!(name, getproperty(shape, Symbol(name)))
        end
        _flag!("rvsout_present", shape.rvsout_present)
        _flag!("rvsou_present", shape.rvsou_present)
        _flag!("rvsiu_present", shape.rvsiu_present)
    end

    if state.vloop !== nothing && isfinite(state.vloop)
        _set!("vloop", state.vloop)
        _flag!("vloop_valid", true)
    end

    return x, indicators, live
end

# ── Per-tick record for validation ──────────────────────────────────────────

"""
    PedestalNNHistory

Every model tick of a run, kept in `dd._aux[:pedestal_nn_history]` for
validation and for comparing fuse29 with fuse_efit across two runs of a shot:
the inputs the model was given (`inputs`, plus `indicators` for fuse_efit), its
raw nine heads (`predictions`, ungated) and the hysteretic H-mode decision.
One entry per host evaluation time; a repeated evaluation at the same time
replaces the entry.
"""
mutable struct PedestalNNHistory
    model::Symbol
    input_names::Vector{String}
    indicator_names::Vector{String}
    head_names::Vector{String}
    time::Vector{Float64}
    inputs::Vector{Vector{Float32}}
    indicators::Vector{Vector{Float32}}
    predictions::Vector{Vector{Float32}}
    is_h_mode::Vector{Bool}
    n_ticks::Vector{Int}
end

function pedestal_nn_record!(aux, model::Symbol, input_names, indicator_names, head_names, time::Float64,
                             inputs::AbstractVector, indicators::AbstractVector, pred::AbstractVector, is_h_mode::Bool, n_ticks::Int)
    h = get(aux, :pedestal_nn_history, nothing)
    if !(h isa PedestalNNHistory) || h.model != model || (!isempty(h.time) && time < h.time[end] - 1e-9)
        h = PedestalNNHistory(model, collect(String, input_names), collect(String, indicator_names), collect(String, head_names),
                              Float64[], Vector{Float32}[], Vector{Float32}[], Vector{Float32}[], Bool[], Int[])
        aux[:pedestal_nn_history] = h
    end
    if !isempty(h.time) && abs(time - h.time[end]) <= 1e-9
        # repeated evaluation of the same host step; a zero-order-hold call (n_ticks == 0) keeps the inputs of the tick it holds
        n_ticks == 0 || (h.inputs[end] = Float32.(inputs); h.indicators[end] = Float32.(indicators))
        h.predictions[end] = Float32.(pred)
        h.is_h_mode[end] = is_h_mode
        h.n_ticks[end] = n_ticks
    else
        held = n_ticks == 0 && !isempty(h.inputs)   # host stepped faster than 50 ms: the model saw nothing new
        push!(h.time, time)
        push!(h.inputs, held ? copy(h.inputs[end]) : Float32.(inputs))
        push!(h.indicators, held ? copy(h.indicators[end]) : Float32.(indicators))
        push!(h.predictions, Float32.(pred))
        push!(h.is_h_mode, is_h_mode)
        push!(h.n_ticks, n_ticks)
    end
    return h
end

"""
    pedestal_nn_history(dd::IMAS.DD) -> Union{Nothing,Dict{String,Any}}

The run's [`PedestalNNHistory`](@ref) as a JSON-ready dictionary (`nothing` when
no pedestal predictor ran): `model`, `time`, and one time trace per input,
indicator and head name under `inputs`, `indicators`, `predictions`, plus
`is_h_mode` and `n_ticks`. NaN (an input the model was told to ignore, or had no
source for) is written as `nothing`.
"""
function pedestal_nn_history(dd::IMAS.DD)
    h = get(getfield(dd, :_aux), :pedestal_nn_history, nothing)
    h isa PedestalNNHistory || return nothing
    clean(v) = isfinite(v) ? Float64(v) : nothing
    traces(names, rows) = Dict{String,Any}(name => Any[clean(row[j]) for row in rows] for (j, name) in enumerate(names))
    return Dict{String,Any}(
        "model" => String(h.model),
        "time" => copy(h.time),
        "inputs" => traces(h.input_names, h.inputs),
        "indicators" => traces(h.indicator_names, h.indicators),
        "predictions" => traces(h.head_names, h.predictions),
        "is_h_mode" => copy(h.is_h_mode),
        "n_ticks" => copy(h.n_ticks))
end

# ── Driving the model from ActorPedestal ────────────────────────────────────

"""
    fuse_efit_predict!(actor::ActorPedestal)

Advance the fuse_efit pedestal predictor to `dd.global_time` and store the
result in `actor.nn_prediction` (same fields as fuse29's, see
[`fuse_efit_prediction`](@ref)). Same time-step bookkeeping as
[`fuse29_predict!`](@ref): the per-shot recurrent state lives in
`dd._aux[:fuse_efit]` and is advanced exactly once per 50 ms of simulated time.

- Actuators come from `par.fpe_source`, the plasma state (Ip, boundary shape,
  loop voltage) from `par.nn_state_source`. Warm-up ticks before the first call
  (`par.nn_warmup_time`) read the stored experimental traces and the
  experimental equilibrium, whatever the live sources are.
- The per-shot availability flags are fixed at the first call: a channel with
  no live source (typically a gas valve the dd does not carry) is flagged
  unavailable for the whole shot — the model was trained with such flags — where
  fuse29 has to fall back on the training median.
- The gate is the same Schmitt trigger; without `nn_hmode_enter` /
  `nn_hmode_exit` it uses the form's own validation thresholds (enter at the
  balanced one, leave below the conservative one: 0.54 / 0.50 for `:single`).
"""
function fuse_efit_predict!(actor::ActorPedestal)
    dd = actor.dd
    par = actor.par

    if !(actor.nn_predictor isa FuseEfitNN) || actor.nn_predictor.form != par.nn_efit_form
        actor.nn_predictor = load_pedestal_efit_nn(; form=par.nn_efit_form)
        @info "ActorPedestal: loaded fuse_efit pedestal predictor $(actor.nn_predictor)"
    end
    nn = actor.nn_predictor::FuseEfitNN
    aux = getfield(dd, :_aux)
    replay_dd = ismissing(actor.act.ActorReplay, :replay_dd) ? nothing : actor.act.ActorReplay.replay_dd
    t_now = dd.global_time
    eps = 1e-3 * nn.period_s

    tr = get(aux, :fuse_efit, nothing)
    if !(tr isa PedestalNNTracker) || t_now < tr.committed_time - eps
        if tr isa PedestalNNTracker
            @warn "ActorPedestal: fuse_efit time rewound from $(tr.committed_time) s to $t_now s; restarting the recurrent state from zero"
        end
        t_start = t_now - nn.period_s
        if !ismissing(par, :nn_warmup_time)
            t_start = min(par.nn_warmup_time, t_start)
        end
        tr = PedestalNNTracker(nn, t_start)
        aux[:fuse_efit] = tr
        # Hand-off time: ticks before this are replay, ticks from here on are live.
        aux[:fuse_efit_handoff] = t_now
        # Per-shot availability, from what the live sources can supply now. The plasma-state channels are
        # always available to FUSE (it has an equilibrium); their per-tick indicators say when a value is missing.
        _, _, live = build_fuse_efit_inputs(nn, dd; source=par.fpe_source, state_source=par.nn_state_source, time=t_now, replay_dd)
        available = Float32.(live)
        for (j, name) in enumerate(nn.input_names)
            (name == "ip" || name == "vloop" || name in FUSE_EFIT_SHAPE_NAMES) && (available[j] = 1.0f0)
        end
        aux[:fuse_efit_available] = available
        missing_ch = [name for (name, a) in zip(nn.input_names, available) if a < 0.5]
        if !isempty(missing_ch)
            @warn "ActorPedestal: fuse_efit channels with no `$(par.fpe_source)` source, flagged unavailable for this shot: $(join(missing_ch, ", "))"
        end
    end
    handoff = get(aux, :fuse_efit_handoff, t_now)
    available = aux[:fuse_efit_available]::Vector{Float32}
    # the host's own step: FUSE's equilibrium slice may be that old when a catch-up tick asks for it
    max_age = max(2 * nn.period_s, t_now - tr.committed_time)

    last_x = fill(NaN32, FUSE_EFIT_N_INPUTS)
    last_ind = zeros(Float32, FUSE_EFIT_N_INDICATORS)

    function step_at(t::Float64, s::PedestalNNState)
        # Warm-up (replay) ticks read FUSE's stored experimental traces and the experimental equilibrium:
        # GSLite pushes current values only and the simulation has no equilibrium before it starts.
        warmup = t < handoff - eps
        x, ind, _ = build_fuse_efit_inputs(nn, dd;
            source=warmup ? :dd : par.fpe_source,
            state_source=warmup ? :replay : par.nn_state_source,
            time=t, replay_dd, max_age)
        if !tr.units_checked && !warmup
            tr.units_checked = true
            for w in fuse_efit_check_units(nn, x, ind, available)
                @warn "ActorPedestal: fuse_efit unit check — $w"
            end
            ind[1] > 0.5 || @warn "ActorPedestal: fuse_efit has no boundary shape from `nn_state_source = $(repr(par.nn_state_source))` at t=$t s (efit_valid = 0); the prediction is not validated without geometry"
            ind[2] > 0.5 || @warn "ActorPedestal: fuse_efit has no loop voltage at t=$t s (vloop_valid = 0)"
        end
        exc = fuse_efit_in_distribution(nn, x, ind, available)
        isempty(exc) || @debug "ActorPedestal: fuse_efit inputs outside training p1..p99 at t=$t s" excursions = exc
        last_x .= x
        last_ind .= ind
        return fuse_efit_step(nn, x, ind, available, s)
    end

    enter = ismissing(par, :nn_hmode_enter) ? nn.hmode_balanced : Float64(par.nn_hmode_enter)
    exit = ismissing(par, :nn_hmode_exit) ? min(nn.hmode_conservative, enter) : Float64(par.nn_hmode_exit)
    pred, is_h, n = pedestal_nn_advance!(tr, nn.period_s, t_now, step_at; hmode_enter=enter, hmode_exit=exit)
    actor.nn_prediction = fuse_efit_prediction(nn, pred, is_h; time=t_now, n_ticks=n)
    pedestal_nn_record!(aux, :fuse_efit, nn.input_names, nn.indicator_names, nn.head_names, t_now, last_x, last_ind, pred, is_h, n)
    return actor
end

"""
    pedestal_nn_predict!(actor::ActorPedestal)

Run the pedestal predictor selected by `act.ActorPedestal.nn_model`
(`:fuse29` or `:fuse_efit`) up to `dd.global_time`; the result is in
`actor.nn_prediction` with the same fields for both models.
"""
function pedestal_nn_predict!(actor::ActorPedestal)
    model = actor.par.nn_model
    if model == :fuse29
        return fuse29_predict!(actor)
    elseif model == :fuse_efit
        return fuse_efit_predict!(actor)
    end
    return error("act.ActorPedestal.nn_model must be :fuse29 or :fuse_efit, got $(repr(model))")
end
