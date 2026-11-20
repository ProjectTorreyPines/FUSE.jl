#= ========================================== =#
#  DIII-D IDA-lite kinetic profiles for postdiction  #
#= ========================================== =#
# IDA-lite (one netCDF4 file per shot, IDA_<shot>_.cdf) holds GP-fitted kinetic profiles on a psi_n grid at ~20 ms
# cadence, either as one profile per slice (direct layout) or as a posterior ensemble. `apply_ida_lite_profiles!`
# writes them into `ini.general.dd.core_profiles` before `init!`, so postdiction replays them instead of the fit.
# The reading and the Zeff/ni route ladder follow bouquet's `io/ida.py` `read_ida`.

const IDA_LITE_CORE_DATASETS = ("time", "psi_n", "n_e", "T_e", "T_12C6", "Zeff")
const IDA_LITE_PROFILE_DATASETS = ("n_e", "T_e", "T_12C6", "Zeff", "n_12C6", "omega_tor_12C6", "rho_tor", "q")
const IDA_LITE_CARBON_FILL = 1e30

"""
    ida_lite_file(dirs::AbstractString, shot::Int)

`IDA_<shot>_.cdf` in the first of the ':'-separated directories `dirs` that has it, else `nothing`
"""
function ida_lite_file(dirs::AbstractString, shot::Int)
    for dir in split(dirs, ':'; keepempty=false)
        path = joinpath(dir, "IDA_$(shot)_.cdf")
        isfile(path) && return path
    end
    return nothing
end

"""
    ida_lite_bad_times(path::AbstractString)

IDA-lite slice times [s] an `ida_lite_filtering_provenance.json` (IDA cleaning) removes, or excluded before flagging
"""
function ida_lite_bad_times(path::AbstractString)
    prov = JSON.parsefile(path)
    return sort!(unique!(Float64[prov["removed_times"]; prov["excluded_before_flagging"]["times"]]))
end

"""
    ida_lite_time(path::AbstractString)

Time basis of an IDA-lite file [s]
"""
ida_lite_time(path::AbstractString) = HDF5.h5open(f -> Float64.(vec(read(f["time"]))) ./ 1e3, path, "r")

"""
    ida_lite_layout(f, path::AbstractString)

Layout of an open IDA-lite file from the shape of `n_e`: `:direct` (n_radial, n_time) or `:ensemble`
(n_radial, n_samples, n_time) as HDF5.jl reports the netCDF dimensions. Anything else is an error.
"""
function ida_lite_layout(f, path::AbstractString)
    datasets = sort!(String[string(k) for k in keys(f)])
    for key in IDA_LITE_CORE_DATASETS
        @assert haskey(f, key) "$path has no '$key' dataset. Datasets present: $(join(datasets, ", "))"
    end
    shp = size(f["n_e"])
    if length(shp) == 2
        layout, n_samples = :direct, 0
        n_radial, n_time = shp
    elseif length(shp) == 3
        layout = :ensemble
        n_radial, n_samples, n_time = shp
    else
        error("$path: 'n_e' has size $shp; supported are (n_radial, n_time) and (n_radial, n_samples, n_time)")
    end
    @assert n_radial >= 2 && n_time >= 1 && (layout === :direct || n_samples >= 1) "$path: 'n_e' has size $shp"
    present = Set{String}()
    for key in IDA_LITE_PROFILE_DATASETS
        haskey(f, key) || continue
        push!(present, key)
        @assert size(f[key]) == shp "$path: '$key' is $(size(f[key])) but 'n_e' is $shp"
    end
    return (; layout, n_radial, n_time, n_samples, datasets, present)
end

"""
    ida_lite_radial_grid(raw::AbstractArray, lay, path::AbstractString; grid_mode::Symbol=:shared)

The one radial grid of `psi_n` (stored 1-D, per slice or per sample and slice) and the spread of the stored copies
around it. `grid_mode=:shared` refuses copies that differ; `:index` takes the first copy anyway (as bouquet does).
"""
function ida_lite_radial_grid(raw::AbstractArray, lay, path::AbstractString; grid_mode::Symbol=:shared)
    @assert grid_mode in (:shared, :index) "grid_mode must be :shared or :index"
    a = Float64.(raw)
    nd = ndims(a)
    if nd == 1
        @assert length(a) == lay.n_radial "$path: 'psi_n' has $(length(a)) points, expected $(lay.n_radial)"
        grid = copy(a)
    elseif nd == 2 && lay.layout === :direct
        @assert size(a) == (lay.n_radial, lay.n_time) "$path: 'psi_n' is $(size(a))"
        grid = a[:, 1]
    elseif nd == 3 && lay.layout === :ensemble
        @assert size(a) == (lay.n_radial, lay.n_samples, lay.n_time) "$path: 'psi_n' is $(size(a))"
        grid = a[:, 1, 1]
    else
        error("$path: 'psi_n' has size $(size(a)), which does not match the $(lay.layout) layout")
    end
    worst, at, lo, hi = 0.0, 0, count(<=(1.0), grid), count(<=(1.0), grid)
    if nd > 1
        flat = reshape(a, lay.n_radial, :)
        for j in axes(flat, 2)
            col = @view flat[:, j]
            d = maximum(abs, col .- grid)
            d > worst && ((worst, at) = (d, j))
            c = count(<=(1.0), col)
            lo, hi = min(lo, c), max(hi, c)
        end
    end
    dispersion = (; worst_abs=worst, worst_index=at, n_in_lcfs_min=lo, n_in_lcfs_max=hi, n_copies=(nd == 1 ? 1 : div(length(a), lay.n_radial)))
    if grid_mode === :shared && worst > 1e-12 * max(1.0, maximum(abs, grid))
        error("$path: 'psi_n' is not the same grid at every time/sample (worst |Δ| $worst at copy $at; " *
              "count(psi_n <= 1) spans $lo..$hi). Use grid_mode=:index to take the first copy anyway.")
    end
    return grid, dispersion
end

"""
    ida_lite_rho_grid(path::AbstractString; grid_mode::Symbol=:shared)

Number of IDA-lite psi_n points inside the separatrix (a natural `ActorFitProfiles.rho_grid`); reads psi_n only
"""
function ida_lite_rho_grid(path::AbstractString; grid_mode::Symbol=:shared)
    return HDF5.h5open(path, "r") do f
        lay = ida_lite_layout(f, path)
        psi_n, _ = ida_lite_radial_grid(read(f["psi_n"]), lay, path; grid_mode)
        return count(<=(1.0), psi_n)
    end
end

# numpy.percentile(v, q) (linear) on a sorted vector
function _ida_lite_percentile(sorted_v::AbstractVector{Float64}, q::Float64)
    n = length(sorted_v)
    n == 1 && return sorted_v[1]
    idx = (q / 100) * (n - 1)
    lo = floor(Int, idx)
    hi = min(lo + 1, n - 1)
    return sorted_v[lo+1] + (idx - lo) * (sorted_v[hi+1] - sorted_v[lo+1])
end

"""
    ida_lite_reduce_ensemble(a::AbstractArray{<:Real,3}; estimator::Symbol, band::Symbol)

(n_radial, n_samples, n_time) to a central value and a symmetric band, both (n_radial, n_time), as bouquet:
`estimator` `:mean` or `:median`, `band` `:percentile` (0.5 (p84 - p16)) or `:std`.
"""
function ida_lite_reduce_ensemble(a::AbstractArray{<:Real,3}; estimator::Symbol, band::Symbol)
    @assert estimator in (:mean, :median) "estimator must be :mean or :median"
    @assert band in (:percentile, :std) "band must be :percentile or :std"
    nr, ns, nt = size(a)
    value = Matrix{Float64}(undef, nr, nt)
    sigma = Matrix{Float64}(undef, nr, nt)
    buf = Vector{Float64}(undef, ns)
    for t in 1:nt, r in 1:nr
        @inbounds for j in 1:ns
            buf[j] = Float64(a[r, j, t])
        end
        m = sum(buf) / ns
        if estimator === :mean
            value[r, t] = m
        else
            s = sort(buf)
            value[r, t] = isodd(ns) ? s[(ns+1)÷2] : (s[ns÷2] + s[ns÷2+1]) / 2
        end
        if band === :std
            sigma[r, t] = ns > 1 ? sqrt(sum(abs2, buf .- m) / ns) : 0.0
        else
            sort!(buf)
            sigma[r, t] = 0.5 * (_ida_lite_percentile(buf, 84.0) - _ida_lite_percentile(buf, 16.0))
        end
    end
    return value, sigma
end

"""
    ida_lite_route_flags(f, lay, path::AbstractString)

Per slice, whether each ni route has a usable envelope (bouquet `read_ida`): VB (Zeff) always on the ensemble
layout, else if `Zeff_err` is present; CER if `n_12C6` (and on the direct layout `n_12C6_err`) is usable at every radius.
"""
function ida_lite_route_flags(f, lay, path::AbstractString)
    nt = lay.n_time
    vb = lay.layout === :ensemble ? trues(nt) : (haskey(f, "Zeff_err") ? trues(nt) : falses(nt))
    cer = falses(nt)
    haskey(f, "n_12C6") || return vb, cer
    lay.layout === :direct && !haskey(f, "n_12C6_err") && return vb, cer
    ok(x) = isfinite(x) && x > 0 && abs(x) < IDA_LITE_CARBON_FILL
    nc, ne = read(f["n_12C6"]), read(f["n_e"])
    if lay.layout === :ensemble
        for k in 1:nt
            cer[k] = all(ok, selectdim(nc, 3, k)) && all(ok, selectdim(ne, 3, k))
        end
        return vb, cer
    end
    snc = read(f["n_12C6_err"])
    sne = haskey(f, "n_e_err") ? read(f["n_e_err"]) : nothing
    for (key, a) in (("n_12C6_err", snc), ("n_e_err", sne))
        @assert a === nothing || size(a) == size(ne) "$path: '$key' is $(size(a)), but 'n_e' is $(size(ne))"
    end
    ok_s(x) = isfinite(x) && x >= 0 && abs(x) < IDA_LITE_CARBON_FILL
    for k in 1:nt
        cer[k] = all(ok, @view(nc[:, k])) && all(ok, @view(ne[:, k])) && all(ok_s, @view(snc[:, k])) &&
                 (sne === nothing || all(isfinite, @view(sne[:, k])))
    end
    return vb, cer
end

"""
    read_ida_lite(path::AbstractString; ensemble_estimator::Symbol=:mean, ensemble_band::Symbol=:percentile, need_carbon::Bool=false, need_rotation::Bool=true, need_q::Bool=true, grid_mode::Symbol=:shared)

Kinetic profiles of one IDA-lite file, either layout: `psi_n`, `time` [s], and `n_e`, `T_e`, `T_12C6` (Ti), `Zeff`,
`n_12C6`, `omega_tor_12C6`, `q` as (n_radial, n_time) matrices in SI units (`nothing` if absent and not needed).
Also the layout, the dataset list, the ensemble band `sigma`, and the per-slice ni route flags `vb_live`/`cer_live`.
"""
function read_ida_lite(path::AbstractString;
    ensemble_estimator::Symbol=:mean,
    ensemble_band::Symbol=:percentile,
    need_carbon::Bool=false,
    need_rotation::Bool=true,
    need_q::Bool=true,
    grid_mode::Symbol=:shared)
    return HDF5.h5open(path, "r") do f
        lay = ida_lite_layout(f, path)
        for (key, needed) in (("n_12C6", need_carbon), ("omega_tor_12C6", need_rotation), ("q", need_q))
            @assert !needed || key in lay.present "$path has no '$key' dataset. Datasets present: $(join(lay.datasets, ", "))"
        end
        time = Float64.(vec(read(f["time"]))) ./ 1e3
        @assert length(time) == lay.n_time "$path: 'time' has $(length(time)) entries for n_time = $(lay.n_time)"
        @assert length(time) >= 2 && issorted(time) "$path: need at least 2 sorted time slices"
        psi_n, grid_dispersion = ida_lite_radial_grid(read(f["psi_n"]), lay, path; grid_mode)
        @assert issorted(psi_n) "$path: psi_n is not sorted"

        sigma = lay.layout === :ensemble ? Dict{Symbol,Matrix{Float64}}() : nothing
        function channel(key::AbstractString, name::Symbol)
            raw = read(f[key])
            if lay.layout === :direct
                a = Float64.(raw)
                @assert size(a) == (lay.n_radial, lay.n_time) "$path: '$key' is $(size(a))"
                return a
            end
            v, s = ida_lite_reduce_ensemble(raw; estimator=ensemble_estimator, band=ensemble_band)
            sigma[name] = s
            return v
        end
        optional(key, name) = key in lay.present ? channel(key, name) : nothing

        n_e = channel("n_e", :n_e)
        T_e = channel("T_e", :T_e)
        T_12C6 = channel("T_12C6", :T_12C6)
        Zeff = channel("Zeff", :Zeff)
        n_12C6 = optional("n_12C6", :n_12C6)
        omega_tor_12C6 = optional("omega_tor_12C6", :omega_tor_12C6)
        q = optional("q", :q)
        vb_live, cer_live = ida_lite_route_flags(f, lay, path)

        return (; path=String(path), psi_n, time, q, n_e, T_e, T_12C6, Zeff, n_12C6, omega_tor_12C6, vb_live, cer_live,
            layout=lay.layout, n_samples=lay.n_samples, datasets=lay.datasets,
            ensemble_estimator=(lay.layout === :ensemble ? ensemble_estimator : :none),
            ensemble_band=(lay.layout === :ensemble ? ensemble_band : :none),
            grid_mode, grid_dispersion, sigma)
    end
end

"""
    ida_lite_keep_times(ida, keep::AbstractVector{Bool})

`ida` restricted to the slices where `keep` is true
"""
function ida_lite_keep_times(ida, keep::AbstractVector{Bool})
    sub(a) = a === nothing ? nothing : a[:, keep]
    sigma = ida.sigma === nothing ? nothing : Dict{Symbol,Matrix{Float64}}(k => v[:, keep] for (k, v) in ida.sigma)
    return merge(ida, (; time=ida.time[keep], q=sub(ida.q), n_e=sub(ida.n_e), T_e=sub(ida.T_e), T_12C6=sub(ida.T_12C6),
        Zeff=sub(ida.Zeff), n_12C6=sub(ida.n_12C6), omega_tor_12C6=sub(ida.omega_tor_12C6),
        vb_live=ida.vb_live[keep], cer_live=ida.cer_live[keep], sigma))
end

"""
    ida_lite_routes(ni_source::Symbol, vb_live::Bool, cer_live::Bool, impurity_Z::Float64)

bouquet `read_ida`'s ni route ladder, VB+CER > CER > VB: returns `(use_vb, use_cer, tier)`.
`:Zeff` and `:CER` are strict; `:all` takes the live routes, else VB. CER needs `impurity_Z == 6`.
"""
function ida_lite_routes(ni_source::Symbol, vb_live::Bool, cer_live::Bool, impurity_Z::Float64)
    @assert ni_source in (:Zeff, :CER, :all) "ni_source must be :Zeff, :CER or :all"
    @assert impurity_Z > 1.0 "impurity_Z must exceed 1"
    @assert ni_source !== :CER || impurity_Z == 6.0 "ni_source=:CER subtracts the carbon density and requires impurity_Z=6.0"
    cer_live &= impurity_Z == 6.0
    if ni_source === :Zeff
        @assert vb_live "ni_source=:Zeff needs the file's 'Zeff_err' dataset (direct layout); use :all"
        return true, false, "VB"
    elseif ni_source === :CER
        return false, cer_live, cer_live ? "CER" : "none"
    end
    (vb_live || cer_live) || return true, false, "VB (ne-fraction sigma)"
    return vb_live, cer_live, (vb_live && cer_live) ? "VB+CER" : cer_live ? "CER" : "VB"
end

"""
    ida_lite_ni(ida; ni_source::Symbol, impurity_Z::Float64)

Main-ion density on the IDA grid per slice, ni = ne (Z - Zeff)/(Z - 1), with Zeff from the routes of
[`ida_lite_routes`](@ref) (each clamped to [1, Z], weighted equally); and each slice's route tier.
"""
function ida_lite_ni(ida; ni_source::Symbol, impurity_Z::Float64)
    Z = impurity_Z
    ni = fill(NaN, size(ida.n_e))
    tier = Vector{String}(undef, length(ida.time))
    n_c0 = zeros(size(ida.n_e, 1))
    for k in eachindex(ida.time)
        use_vb, use_cer, tier[k] = ida_lite_routes(ni_source, ida.vb_live[k], ida.cer_live[k], impurity_Z)
        (use_vb || use_cer) || continue
        ne = @view(ida.n_e[:, k])
        n_c = ida.n_12C6 === nothing ? n_c0 : @view(ida.n_12C6[:, k])
        w_v = (use_vb && use_cer) ? 0.5 : use_vb ? 1.0 : 0.0
        w_c = (use_vb && use_cer) ? 0.5 : use_cer ? 1.0 : 0.0
        zres = w_v .* clamp.(@view(ida.Zeff[:, k]), 1.0, Z)
        if w_c > 0
            zcer = 1.0 .+ Z * (Z - 1.0) .* n_c ./ max.(ne, 1e10)
            zres = w_v > 0 ? zres .+ w_c .* clamp.(zcer, 1.0, Z) : w_c .* clamp.(zcer, 1.0, Z)
        end
        zeff = clamp.(zres, 1.0, Z)
        ni[:, k] = ne .* (Z .- zeff) ./ (Z - 1.0)
    end
    return ni, tier
end

"""
    ida_lite_phi_n(psi_n::AbstractVector{<:Real}, q::AbstractVector{<:Real})

IDA's normalized toroidal flux at its `psi_n` points from its own q, ∫q dψ / ∫₀¹ q dψ (trapezoid, q linear across
psi_n = 1), as bouquet `coords.phi_n_from_q`. `nothing` unless psi_n starts at 0 and rises and q is finite and non-zero.
"""
function ida_lite_phi_n(psi_n::AbstractVector{<:Real}, q::AbstractVector{<:Real})
    x = Float64.(psi_n)
    y = abs.(Float64.(q))
    (length(x) >= 2 && abs(first(x)) <= 1e-9 && all(>(0), diff(x)) && all(isfinite, y) && all(>(0), y)) || return nothing
    phi = zeros(length(x))
    for i in 2:length(x)
        phi[i] = phi[i-1] + 0.5 * (y[i] + y[i-1]) * (x[i] - x[i-1])
    end
    j = findfirst(>=(1.0), x)
    j === 1 && return nothing
    i0 = j === nothing ? length(x) : j - 1
    q1 = j === nothing ? y[end] : y[i0] + (y[j] - y[i0]) * (1.0 - x[i0]) / (x[j] - x[i0])
    return phi ./ (phi[i0] + 0.5 * (y[i0] + q1) * (1.0 - x[i0]))
end

# slices whose q gives a Phi_N map, with finite positive ne/Te/Ti over the mapped points (and a live CER route for :CER)
function ida_lite_usable_slices(ida, n_psi::Int; ni_source::Symbol, impurity_Z::Float64)
    ok(a, k) = all(x -> isfinite(x) && x > 0.0, @view(a[1:n_psi, k]))
    return BitVector(map(eachindex(ida.time)) do k
        ni_source === :CER && !(ida.cer_live[k] && impurity_Z == 6.0) && return false
        ida_lite_phi_n(ida.psi_n[1:n_psi], @view(ida.q[1:n_psi, k])) === nothing && return false
        return ok(ida.n_e, k) && ok(ida.T_e, k) && ok(ida.T_12C6, k)
    end)
end

# IDA slice k onto cp1d's rho_tor_norm (Phi_N = rho_tor_norm²), written as ActorFitProfiles writes its fit
function ida_lite_map_slice!(cp1d::IMAS.core_profiles__profiles_1d, ida, k::Int, n_psi::Int; radial_interp::Symbol, rotation::Bool)
    skipped(why) = (; outcome=:skipped, notes=String[why])
    phi_k = ida_lite_phi_n(ida.psi_n[1:n_psi], @view(ida.q[1:n_psi, k]))
    phi_k === nothing && return skipped(@sprintf("SKIPPED: IDA slice %.3f s has no usable q, kept the fit", ida.time[k]))
    x_nodes = clamp.(cp1d.grid.rho_tor_norm .^ 2, first(phi_k), last(phi_k))
    on_rho(a) = IMAS.interp1d(phi_k, a[1:n_psi, k], radial_interp).(x_nodes)

    ne = on_rho(ida.n_e)
    te = on_rho(ida.T_e)
    ti = on_rho(ida.T_12C6)
    if any(x -> !isfinite(x) || x <= 0.0, ne) || any(x -> !isfinite(x) || x <= 0.0, te) || any(x -> !isfinite(x) || x <= 0.0, ti)
        return skipped("SKIPPED: non-finite or non-positive ne/Te/Ti, kept the fit")
    end
    ni = on_rho(ida.ni)

    bulk_ion = nothing
    imp_ion = nothing
    for ion in cp1d.ion
        if IMAS.is_hydrogenic(ion)
            bulk_ion === nothing && (bulk_ion = ion)
        else
            imp_ion === nothing && (imp_ion = ion)
        end
    end
    if bulk_ion === nothing || imp_ion === nothing
        return (; outcome=:no_pair, notes=String["SKIPPED: no hydrogenic/impurity ion pair on this slice, kept the fit"])
    end

    # as ActorFitProfiles: the impurity carries the dilution and the main ion follows from quasi-neutrality
    cp1d.electrons.density_thermal = ne
    IMAS.unfreeze!(cp1d.electrons, :density)
    cp1d.electrons.temperature = te
    for ion in cp1d.ion
        ion.temperature = ti
    end
    Zavg = IMAS.avgZ(imp_ion)
    bulk_ion.density_thermal = zero(ne)
    imp_ion.density_thermal = clamp.((ne .- ni) ./ Zavg, 0.0, ne ./ Zavg .* 0.99)
    IMAS.unfreeze!(bulk_ion, :density)
    IMAS.unfreeze!(imp_ion, :density)
    IMAS.enforce_quasi_neutrality!(cp1d, :D)

    notes = String[]
    if rotation && ida.omega_tor_12C6 !== nothing
        ωtor = on_rho(ida.omega_tor_12C6)
        if all(isfinite, ωtor)
            for ion in cp1d.ion
                ion.rotation_frequency_tor = ωtor
            end
            IMAS.ωtor2sonic!(cp1d)
        else
            push!(notes, "non-finite rotation, kept the fitted rotation")
        end
    end
    return (; outcome=:replaced, notes)
end

"""
    apply_ida_lite_profiles!(dd::IMAS.DD, path::AbstractString;
        ni_source::Symbol=:all, impurity_Z::Float64=6.0, rotation::Bool=true, time_pad::Float64=0.02, radial_interp::Symbol=:pchip,
        ensemble_estimator::Symbol=:mean, ensemble_band::Symbol=:percentile, grid_mode::Symbol=:shared,
        bad_times::AbstractVector{<:Real}=Float64[], sim_times::AbstractVector{<:Real}=Float64[], strict_pairing::Bool=false)

Overwrite `dd.core_profiles` (in a postdiction, `ini.general.dd` before `init!`) with the IDA-lite profiles of `path`.

Every core_profiles slice within `time_pad` of the IDA-lite coverage takes the nearest IDA slice
([`assign_profile_slices`](@ref)) verbatim; outside it the fit is kept. IDA slices at `bad_times` and slices
without a usable q or ne/Te/Ti are dropped first. The radial map is IDA's Phi_N from its q onto rho_tor_norm²,
as bouquet places IDA profiles. ni follows bouquet's Zeff route ladder ([`ida_lite_ni`](@ref)), the impurity carries
the dilution and the main ion follows from quasi-neutrality. Rotation is replaced if `rotation`.

[`replay_pairing`](@ref) over `sim_times` must find IDA at every step inside coverage (`strict_pairing`: else error, otherwise warn).
Returns what was done: `(; ida, assign, m, pairing, bad_dropped, unusable_times, t_lo, t_hi, n_psi, cp_times, in_band)`.
"""
function apply_ida_lite_profiles!(dd::IMAS.DD, path::AbstractString;
    ni_source::Symbol=:all,
    impurity_Z::Float64=6.0,
    rotation::Bool=true,
    time_pad::Float64=0.02,
    radial_interp::Symbol=:pchip,
    ensemble_estimator::Symbol=:mean,
    ensemble_band::Symbol=:percentile,
    grid_mode::Symbol=:shared,
    bad_times::AbstractVector{<:Real}=Float64[],
    sim_times::AbstractVector{<:Real}=Float64[],
    strict_pairing::Bool=false)

    ida = read_ida_lite(path; ensemble_estimator, ensemble_band, need_carbon=(ni_source === :CER), need_rotation=rotation, need_q=true, grid_mode)
    bad_dropped = Float64[]
    if !isempty(bad_times)
        bad = [round(Float64(t); digits=6) for t in bad_times]
        keep = [!(round(t; digits=6) in bad) for t in ida.time]
        @assert any(keep) "$path: no IDA-lite slice left after dropping bad_times"
        if !all(keep)
            bad_dropped = ida.time[.!keep]
            ida = ida_lite_keep_times(ida, keep)
        end
    end

    # IDA points up to and including the first psi_n >= 1
    n_psi = something(findfirst(>=(1.0), ida.psi_n), length(ida.psi_n))
    usable = ida_lite_usable_slices(ida, n_psi; ni_source, impurity_Z)
    @assert any(usable) "$path: no usable IDA-lite slice (q, ne/Te/Ti or CER route)"
    unusable_times = ida.time[.!usable]
    if !all(usable)
        @warn "apply_ida_lite_profiles!: dropped $(length(unusable_times)) unusable IDA-lite slice(s) before matching: $(unusable_times)"
        ida = ida_lite_keep_times(ida, usable)
    end
    ni, ni_tier = ida_lite_ni(ida; ni_source, impurity_Z)
    ida = merge(ida, (; ni, ni_tier))
    t_lo, t_hi = first(ida.time) - time_pad, last(ida.time) + time_pad

    in_band = [cp1d.time for cp1d in dd.core_profiles.profiles_1d if t_lo <= cp1d.time <= t_hi]
    first_in = findfirst(cp1d -> t_lo <= cp1d.time <= t_hi, dd.core_profiles.profiles_1d)
    if first_in !== nothing
        ions = dd.core_profiles.profiles_1d[first_in].ion
        @assert any(IMAS.is_hydrogenic, ions) && !all(IMAS.is_hydrogenic, ions) "IDA-lite profiles need a hydrogenic and an impurity ion in core_profiles"
    end
    assign = assign_profile_slices(in_band, ida.time)
    k_of = Dict(zip(in_band, assign.k))

    replaced, dt_used, notes, n_no_pair = Float64[], Float64[], Tuple{Float64,String}[], 0
    for cp1d in dd.core_profiles.profiles_1d
        time0 = cp1d.time
        t_lo <= time0 <= t_hi || continue
        k = k_of[time0]
        r = ida_lite_map_slice!(cp1d, ida, k, n_psi; radial_interp, rotation)
        for why in r.notes
            push!(notes, (time0, why))
        end
        r.outcome === :no_pair && (n_no_pair += 1)
        r.outcome === :replaced || continue
        push!(replaced, time0)
        push!(dt_used, time0 - ida.time[k])
    end
    m = (; replaced, dt_used, notes, n_no_pair)

    cp_times = [cp1d.time for cp1d in dd.core_profiles.profiles_1d]
    @assert !isempty(replaced) "no core_profiles slice inside the IDA-lite coverage [$t_lo, $t_hi] s was replaced (core_profiles times: $(extrema(cp_times; init=(NaN, NaN))))"
    pairing = replay_pairing(cp_times, sim_times, k_of, ida.time, replaced; t_lo, t_hi)
    bad_pairs = [r for r in pairing if t_lo + 1e-6 <= r.t_sim <= t_hi && r.outcome !== :source]
    if !isempty(bad_pairs)
        msg = "apply_ida_lite_profiles!: $(length(bad_pairs)) step(s) inside IDA-lite coverage would replay a non-IDA profile: " *
              join(("t=$(r.t_sim) -> tick $(r.tick) ($(r.outcome))" for r in bad_pairs[1:min(end, 5)]), "; ")
        strict_pairing ? error(msg) : @warn(msg)
    end
    dd.core_profiles.code.name = "IDA-lite"
    return (; ida, assign, m, pairing, bad_dropped, unusable_times, t_lo, t_hi, n_psi, cp_times, in_band)
end

"""
    ida_lite_provenance(r)

What [`apply_ida_lite_profiles!`](@ref) returned (`r`), as written to `ida_provenance.json`: the file, its coverage, the slices dropped,
and which IDA-lite slice the replay serves at each simulation step (`replay_pairing`)
"""
function ida_lite_provenance(r)
    pairing = [Dict("t_sim" => p.t_sim, "tick" => isnan(p.tick) ? nothing : p.tick, "ida_time" => p.source_time, "outcome" => p.outcome === :source ? "ida" : String(p.outcome))
               for p in r.pairing]
    return Dict("ida_file" => r.ida.path, "t_lo" => r.t_lo, "t_hi" => r.t_hi, "bad_dropped" => r.bad_dropped, "unusable_times" => r.unusable_times,
        "n_replaced" => length(r.m.replaced), "sim_times" => [p.t_sim for p in r.pairing], "replay_pairing" => pairing)
end
