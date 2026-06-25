import TJLFEP
import ALPHA
import TurbulentTransport

#= =========== =#
#  ActorTJLFEP  #
#= =========== =#
@actor_parameters_struct ActorTJLFEP{T} begin
    rho_scan::Entry{AbstractVector{T}} =
        Entry{AbstractVector{T}}("-", "rho_tor_norm radii at which TGLF-EP finds the critical EP-density scale SFmin"; default=[0.01, 0.21, 0.41, 0.61, 0.81, 0.95])
    is_ep::Entry{Int} = Entry{Int}("-", "ion index of the energetic-particle (EP) driver species in dd.core_profiles.ion"; default=3)
    process_in::Entry{Int} = Entry{Int}("-", "TGLF-EP PROCESS_IN (4 or 5: critical-EP-density-gradient scan)"; default=5)
    solver::Switch{Symbol} =
        Switch{Symbol}([:grid, :ad, :robust_ad, :truth], "-", "critical-factor engine: :ad (default; fast autodiff AE-onset Newton + IFT (kyhat,width) descent, width-extended via extend_mode=:locate — tracks :robust_ad's accuracy closely at a fraction of the cost: the speed/accuracy sweet spot), :grid (Fortran-equivalent kwscale_scan sweep; SFmin matches Fortran — use for Fortran equivalence or QL-diffusivity coupling), :robust_ad (robust autodiff: global-min of the all-filter onset over the (kyhat,width) grid + refine_rounds of window narrowing; most accurate AD path), or :truth (extends width below WIDTH_MIN for the genuine narrow-width EP-driven AEs; NOT Fortran-faithful, returns the most-unstable physical threshold)"; default=:ad)
    refine_rounds::Entry{Int} =
        Entry{Int}("-", "accuracy/speed knob for solver=:robust_ad: rounds of (kyhat,width) window narrowing around the running best (0 = coarse-grid min; higher = better resolution of off-node binding points such as the plasma edge, at proportionally higher cost). Ignored by :grid and :ad."; default=1)
    extend_mode::Switch{Symbol} =
        Switch{Symbol}([:locate, :wide], "-", "solver=:ad width-extension strategy: :locate (dense log-grid + multistart descents + grid-floor guard; tracks :robust_ad closely) or :wide (fast single-pass log-seeded multistart, ~2x cheaper than :locate and conservative — always >= robust_ad, within ~1-2x; recommended for bulk NN-database generation). Ignored by :grid/:robust_ad/:truth."; default=:locate)
    wide_kdesc::Entry{Int} =
        Entry{Int}("-", "solver=:ad extend_mode=:wide multistart breadth (number of well-separated descents). Higher closes the residual over-prediction gap to :robust_ad at proportionally higher cost; 2 is the accuracy/cost sweet spot. Ignored unless solver=:ad and extend_mode=:wide."; default=2)
    faithful_confirm::Entry{Bool} =
        Entry{Bool}("-", "solver=:ad: confirm the located onset with the expensive all-filter (IFLUX=true) keep-onset (true, production: conservative and matches :robust_ad's filtering). false = cheap 'pure AD' AE-band onset only (faster but can dangerously under-predict; not recommended). Ignored by :grid/:robust_ad/:truth."; default=true)
    threshold_flag::Entry{Int} = Entry{Int}("-", "TGLF-EP THRESHOLD_FLAG"; default=0)
    scan_method::Entry{Int} = Entry{Int}("-", "TGLF-EP SCAN_METHOD"; default=1)
    n_basis::Entry{Int} = Entry{Int}("-", "TGLF-EP N_BASIS (number of Hermite basis functions)"; default=2)
    ky_model::Entry{Int} = Entry{Int}("-", "TGLF-EP KY_MODEL"; default=2)
    nmodes::Entry{Int} = Entry{Int}("-", "number of eigenmodes"; default=4)
    nn::Entry{Int} = Entry{Int}("-", "TGLF-EP nn (number of scalefactor iterations)"; default=5)
    width_min::Entry{T} = Entry{T}("-", "minimum Gaussian width"; default=1.0)
    width_max::Entry{T} = Entry{T}("-", "maximum Gaussian width"; default=2.0)
    factor_in::Entry{T} = Entry{T}("-", "initial EP-density scalefactor"; default=1.0)
    reject_i_pinch::Entry{Bool} = Entry{Bool}("-", "reject ion-pinch modes"; default=false)
    reject_e_pinch::Entry{Bool} = Entry{Bool}("-", "reject electron-pinch modes"; default=false)
    reject_th_pinch::Entry{Bool} = Entry{Bool}("-", "reject thermal-pinch modes"; default=true)
    reject_ep_pinch::Entry{Bool} = Entry{Bool}("-", "reject EP-pinch modes"; default=false)
    reject_tearing::Entry{Bool} = Entry{Bool}("-", "reject tearing-parity modes"; default=true)
    rotational_suppression::Entry{Bool} = Entry{Bool}("-", "apply rotational suppression"; default=true)
    alpha_method::Switch{Symbol} =
        Switch{Symbol}([:density, :pressure], "-", "ALPHA critical-gradient variable (density or pressure threshold)"; default=:density)
    alpha_solver::Switch{Symbol} =
        Switch{Symbol}([:stiff, :marginal], "-", "ALPHA solver (:stiff = Fortran stiff-CGM; :marginal = analytic min(classical,marginal))"; default=:stiff)
    alpha_use_ql::Entry{Bool} =
        Entry{Bool}("-", "use TGLF-EP-derived gamma_star/diff_star in stiff-CGM QL diffusivity"; default=false)
    E_alpha::Entry{T} = Entry{T}("MeV", "EP birth energy (3.5 MeV for fusion alphas)"; default=3.5)
    use_gpu::Entry{Bool} = Entry{Bool}("-", "run TJLF on GPU when available"; default=false)
    inner::Switch{Symbol} =
        Switch{Symbol}([:threads, :mps_team], "-", "within-radius parallelism (SPMD per-radius layout only; the in-process actor run always uses :threads): :threads (serial-per-GPU baseline; fastest for solver=:ad) or :mps_team (MPS clients sharing each GPU; fastest for solver=:grid)"; default=:threads)
    mps_team::Entry{Int} = Entry{Int}("-", "MPS-team size per GPU when inner=:mps_team (0 disables); recommended only for solver=:grid (the :ad path is faster with :threads)"; default=8)
end

mutable struct ActorTJLFEP{D,P} <: SingleAbstractActor{D,P}
    dd::IMAS.dd{D}
    par::OverrideParameters{P,FUSEparameters__ActorTJLFEP{P}}
    rho_grid::Vector{D}            # full rho_tor_norm grid the critical gradients live on
    SFmin::Vector{D}               # critical EP-density scalefactor per scan radius (stability metric)
    width::Vector{D}               # Gaussian width per scan radius
    kymark::Vector{D}              # marginal ky per scan radius
    dndr_crit::Vector{D}           # critical EP density gradient on rho_grid [10^19 m^-3 / m]
    dpdr_crit::Vector{D}           # critical EP pressure gradient on rho_grid [10 kPa / m]
    alpha::Union{Nothing,ALPHA.AlphaResult{D}}
end

"""
    ActorTJLFEP(dd::IMAS.dd, act::ParametersAllActors; kw...)

Energetic-particle (EP) transport actor combining TGLF-EP and ALPHA.

TGLF-EP (`TJLFEP.runTHD`) finds the critical EP-density scale `SFmin` per radius —
the EP density at which the linear Alfvén-eigenmode (AE) drive turns on — and from
it the critical EP density- and pressure-gradient profiles. ALPHA (`ALPHA.run_alpha`)
then integrates those critical gradients (a steady-state critical-gradient model) to
obtain the EP radial profiles: density `n_EP`, pressure `p_EP`, temperature
`T_EP = p_EP / n_EP`, and the associated EP particle/energy flux.

`_finalize` writes the integrated EP profiles to `dd.core_profiles` as a fast-ion
population on the EP species and the EP flux to `dd.core_transport`. AE-stability
metrics (`SFmin`, `width`, `kymark`) are kept on the actor struct.
"""
function ActorTJLFEP(dd::IMAS.dd, act::ParametersAllActors; kw...)
    actor = ActorTJLFEP(dd, act.ActorTJLFEP; kw...)
    step(actor)
    finalize(actor)
    return actor
end

function ActorTJLFEP(dd::IMAS.dd{D}, par::FUSEparameters__ActorTJLFEP; kw...) where {D<:Real}
    logging_actor_init(ActorTJLFEP)
    par = OverrideParameters(par; kw...)
    return ActorTJLFEP(dd, par, D[], D[], D[], D[], D[], D[], nothing)
end

"""
    _optionsdict(par) -> Dict{String,Any}

Build the TGLF-EP `OptionsDict` consumed by `TJLFEP.runTHD` from the actor parameters.
"""
function _optionsdict(par::OverrideParameters{P,FUSEparameters__ActorTJLFEP{P}}) where {P<:Real}
    scan_n = length(par.rho_scan)
    return Dict{String,Any}(
        "nn" => par.nn, "jtscale_max" => 1, "nmodes" => par.nmodes,
        "PROCESS_IN" => par.process_in, "THRESHOLD_FLAG" => par.threshold_flag, "N_BASIS" => par.n_basis,
        "SCAN_METHOD" => par.scan_method,
        "REJECT_I_PINCH_FLAG" => Int(par.reject_i_pinch), "REJECT_E_PINCH_FLAG" => Int(par.reject_e_pinch),
        "REJECT_TH_PINCH_FLAG" => Int(par.reject_th_pinch), "REJECT_EP_PINCH_FLAG" => Int(par.reject_ep_pinch),
        "REJECT_TEARING_FLAG" => Int(par.reject_tearing), "ROTATIONAL_SUPPRESSION_FLAG" => Int(par.rotational_suppression),
        "QL_RATIO_THRESH" => 0.001, "THETA_SQ_THRESH" => 100.0, "Q_SCALE" => 1.0,
        "WRITE_WAVEFUNCTION" => 0, "KY_MODEL" => par.ky_model, "SCAN_N" => scan_n, "IRS" => 2,
        "FACTOR_IN_PROFILE" => false, "FACTOR_IN" => par.factor_in,
        "WIDTH_IN_FLAG" => false, "WIDTH_MIN" => par.width_min, "WIDTH_MAX" => par.width_max,
        "INPUT_PROFILE_METHOD" => 2, "N_ION" => 3, "IS_EP" => par.is_ep, "REAL_FREQ" => 1)
end

"""
    _step(actor::ActorTJLFEP)

Runs TGLF-EP (`TJLFEP.runTHD`) to obtain the critical EP density/pressure gradients,
then integrates them with ALPHA (`ALPHA.run_alpha`) into the EP profiles and flux.
"""
function _step(actor::ActorTJLFEP{D,P}) where {D<:Real,P<:Real}
    dd = actor.dd
    par = actor.par

    rho_scan = collect(Float64, par.rho_scan)
    OptionsDict = _optionsdict(par)

    # The AD solvers (:ad, :faithful_grid) return the critical factor/marking only (no
    # per-mode QL buffers), so QL-diffusivity coupling is unavailable on those paths.
    if par.alpha_use_ql && par.alpha_solver == :stiff && par.solver != :grid
        @warn "ActorTJLFEP: alpha_use_ql is not supported with solver=$(par.solver) (no QL buffers from the AD path); falling back to non-QL stiff-CGM."
    end
    use_ql = par.alpha_use_ql && par.alpha_solver == :stiff && par.solver == :grid
    width, kymark_out, SFmin, dpdr_crit_out, dndr_crit_out, marginal_ql =
        TJLFEP.runTHD(dd, rho_scan, OptionsDict; use_gpu=par.use_gpu, ql_flux_scan=use_ql,
            inner=par.inner, mps_team=par.mps_team, solver=par.solver, refine_rounds=par.refine_rounds,
            extend_mode=par.extend_mode, wide_kdesc=par.wide_kdesc, faithful_confirm=par.faithful_confirm)

    actor.rho_grid = D.(dd.core_profiles.profiles_1d[].grid.rho_tor_norm)
    actor.SFmin = D.(SFmin)
    actor.width = D.(width)
    actor.kymark = D.(kymark_out)
    actor.dndr_crit = D.(dndr_crit_out)
    actor.dpdr_crit = D.(dpdr_crit_out)

    # Where TGLF-EP found no AE limit the critical gradient is NaN / >=9000 sentinel;
    # map those to a large gradient so the integrated marginal profile stays well above
    # the classical slowing-down profile and `min(classical, marginal)` keeps the
    # (un-flattened) classical EP density there.
    dndr = _sanitize_crit_grad(dndr_crit_out)
    dpdr = _sanitize_crit_grad(dpdr_crit_out) ./ 0.16022   # 10 kPa/m -> 10^19 m^-3·keV /m (see runTHD)

    crit_grad = (; dndr_crit=dndr, dpdr_crit=dpdr)

    ql_modes = nothing
    transport_params = nothing
    if use_ql
        km = par.nmodes
        raw = TJLFEP.build_alpha_ql_modes(marginal_ql, rho_scan, actor.rho_grid, dndr; km_max=km)
        ql_modes = [
            ALPHA.QLModeInput{D}(;
                gamma_star=m.gamma_star,
                diff_star=m.diff_star,
                rg_n_crit=m.rg_n_crit,
                crit_index_shift=m.crit_index_shift,
                crit_scale=m.crit_scale,
            ) for m in raw
        ]
        transport_params = ALPHA.AlphaTransportParams{D}(; use_ql_diffusivity=true)
    end

    actor.alpha = ALPHA.run_alpha(dd, actor.rho_grid, crit_grad;
        solver=par.alpha_solver, method=par.alpha_method, E_alpha=Float64(par.E_alpha),
        transport_params=transport_params, ql_modes=ql_modes)

    return actor
end

"""Replace NaN / `>=9000` "no AE limit" sentinels with a large finite gradient (10^4)."""
function _sanitize_crit_grad(g::AbstractVector)
    return [(!isfinite(x) || abs(x) >= 9000.0) ? 1.0e4 : x for x in g]
end

"""
    _finalize(actor::ActorTJLFEP)

Writes the integrated EP profiles to `dd.core_profiles` (fast-ion population on the
EP species: `density_fast`, `pressure_fast_parallel`, `pressure_fast_perpendicular`)
and the EP particle/energy flux to `dd.core_transport`.
"""
function _finalize(actor::ActorTJLFEP{D,P}) where {D<:Real,P<:Real}
    dd = actor.dd
    par = actor.par
    res = actor.alpha
    res === nothing && return actor

    cp1d = dd.core_profiles.profiles_1d[]
    rho_cp = cp1d.grid.rho_tor_norm

    # keV (10^19 m^-3) -> Pa : n[m^-3]·T[keV]·1e3·e
    keV19_to_Pa = 1.0e19 * 1.0e3 * IMAS.mks.e

    n_fast = IMAS.interp1d(res.rho, res.n_EP).(rho_cp) .* 1e19           # m^-3
    p_fast = IMAS.interp1d(res.rho, res.p_EP).(rho_cp) .* keV19_to_Pa   # Pa

    ep_ion = cp1d.ion[par.is_ep]
    ep_ion.density_fast = max.(n_fast, 0.0)
    # isotropic fast pressure: p_par + 2·p_perp = p_fast  =>  p_par = p_perp = p_fast/3
    ep_ion.pressure_fast_parallel = max.(p_fast ./ 3.0, 0.0)
    ep_ion.pressure_fast_perpendicular = max.(p_fast ./ 3.0, 0.0)

    # EP flux to core_transport. These are local flux *densities* (Γ [m^-2 s^-1],
    # q [W/m^2]); unlike the GACODE->IMAS anomalous fluxes they carry no V' Miller
    # weighting, since ALPHA is a standalone steady-state EP solver (diagnostic
    # output, not iterated by the FUSE flux-matcher).
    model = resize!(dd.core_transport.model, :anomalous; wipe=false)
    model.identifier.name = "TJLFEP-ALPHA EP ($(par.alpha_method))"
    m1d = resize!(model.profiles_1d)
    m1d.grid_flux.rho_tor_norm = collect(res.rho)
    ion = resize!(m1d.ion, 1)[1]
    ion.label = ismissing(ep_ion, :label) ? "EP" : ep_ion.label
    ion.particles.flux = res.flux_particle .* 1e19              # m^-2 s^-1
    ion.energy.flux = res.flux_energy .* keV19_to_Pa           # W/m^2 (keV·10^19 m^-2 s^-1 -> W/m^2)

    return actor
end

#= ================= =#
#  TGLF-EP run recipe  #
#= ================= =#
"""
    plot(state::TurbulentTransport.TGLFEPRunState; isep=nothing, E_alpha=3.5)

Before/after energetic-particle profile plot for a completed `run_tjlfep` scan
(analogous to `plot(dd.core_profiles)`). Reads the post-run `dd_out.json` from the run
directory and overlays the fast-alpha content that `ActorTJLFEP`+ALPHA modify, in three
panels: EP fast-ion density, EP fast-ion pressure, and total core pressure.

  - `after`  = the AE-limited profile written to `dd_out.json`.
  - `before` = the classical (unlimited) slowing-down alpha profile reconstructed from the
    `dd_out` background plasma, so the plot works even when `dd_in.json` was not saved.

The EP/alpha species is auto-detected by element (helium, `z_n ≈ 2`); override with `isep`.
"""
@recipe function plot_TGLFEPRunState(state::TurbulentTransport.TGLFEPRunState; isep=nothing, E_alpha=3.5)
    dd_out = joinpath(state.basedir, "dd_out.json")
    isfile(dd_out) || error("plot(TGLFEPRunState): missing $(dd_out) — is the run complete?")
    ddo = IMAS.json2imas(dd_out)
    cp1d = ddo.core_profiles.profiles_1d[]
    rho = cp1d.grid.rho_tor_norm
    nz = zeros(length(rho))
    vget(o, f) = try (v = getproperty(o, f); (v === missing || isempty(v)) ? nz : v) catch; nz end
    dfast(ion) = vget(ion, :density_fast)
    pfast(ion) = vget(ion, :pressure_fast_parallel) .+ 2 .* vget(ion, :pressure_fast_perpendicular)

    # EP/alpha species: prefer helium (z_n ≈ 2), else the ion with the largest fast density
    znum(e) = try e.z_n catch; missing end
    is_he(ion) = any(e -> (z = znum(e); z !== missing && abs(z - 2) < 0.5), ion.element)
    if isep === nothing
        isep = something(findfirst(k -> is_he(cp1d.ion[k]), eachindex(cp1d.ion)),
            argmax([maximum(dfast(ion)) for ion in cp1d.ion]))
    end
    lab = cp1d.ion[isep].label

    # classical (unlimited) slowing-down alpha profile from the dd_out background plasma
    inp = ALPHA.AlphaInput(ddo, rho; E_alpha=E_alpha)
    n_cl, T_eq, _, _ = ALPHA.slowing_down(inp.ne, inp.Te, inp.Ti, inp.ni;
        E_alpha=inp.E_alpha, Z1=inp.Z1, ln_lambda=inp.ln_lambda)
    n_before = n_cl .* 1e19                          # 10^19 -> m^-3
    p_before = n_before .* (T_eq .* 1e3 .* IMAS.mks.e)  # keV -> Pa
    n_after = dfast(cp1d.ion[isep])
    p_after = pfast(cp1d.ion[isep])
    p_th = IMAS.pressure_thermal(cp1d)
    p_tot = p_th .+ sum(pfast(ion) for ion in cp1d.ion)

    layout := (1, 3)
    size --> (1300, 380)
    left_margin --> 6 * Plots.Measures.mm
    bottom_margin --> 8 * Plots.Measures.mm
    top_margin --> 5 * Plots.Measures.mm

    @series begin
        subplot := 1
        seriestype := :path
        linestyle := :dash
        linewidth := 2
        linecolor := :gray45
        label := " before (classical)"
        xguide := "rho_tor_norm"
        yguide := "n_fast [10¹⁹ m⁻³]"
        title := "$(lab) fast-ion density"
        rho, n_before ./ 1e19
    end
    @series begin
        subplot := 1
        linewidth := 2.5
        linecolor := :crimson
        label := " after (TJLFEP+ALPHA)"
        rho, n_after ./ 1e19
    end

    @series begin
        subplot := 2
        linestyle := :dash
        linewidth := 2
        linecolor := :gray45
        label := " before (classical)"
        xguide := "rho_tor_norm"
        yguide := "p_fast [kPa]"
        title := "$(lab) fast-ion pressure"
        rho, p_before ./ 1e3
    end
    @series begin
        subplot := 2
        linewidth := 2.5
        linecolor := :crimson
        label := " after"
        rho, p_after ./ 1e3
    end

    @series begin
        subplot := 3
        linestyle := :dash
        linewidth := 2
        linecolor := :black
        label := " thermal only"
        xguide := "rho_tor_norm"
        yguide := "pressure [kPa]"
        title := "total core pressure"
        rho, p_th ./ 1e3
    end
    @series begin
        subplot := 3
        linewidth := 2.5
        linecolor := :red
        label := " thermal + EP"
        rho, p_tot ./ 1e3
    end
end
