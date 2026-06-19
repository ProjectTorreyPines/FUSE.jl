# End-to-end smoke test for the NN pedestal predictor on a synthesized DIII-D
# shot. Exercises the live ZMQ → FPE channel mapping (all 32 inputs),
# the MSE history buffer round-trip (NPZ → dd._aux → encoder), and the final
# `predict_pedestal` call — the exact path ActorPedestal._step takes when
# `par.ne_from == :nn_predictor`.
#
# The scenario mirrors a DIII-D L-mode-to-H-mode transition snapshot with
# realistic actuator values; tweak the `SHOT_*` constants below to explore
# sensitivity. The smoke is intentionally *not* gated on exact prediction
# values — it asserts physical plausibility ranges only.
#
# Run from FUSE root:
#     julia --project=. test/smoke_pedestal_nn_shot.jl
#
# Requires `data/pedestal_nn/history_morning.npz` (checked into the repo).
# The ONNX weights are resolved through `FUSE_PEDESTAL_NN_DIR`; see
# `data/pedestal_nn/README.md` for container wiring.

using Printf
using Test
using Statistics: mean

using FUSE
import IMAS
import NPZ

# ── Test scenario ────────────────────────────────────────────────────────────
# Realistic DIII-D L→H actuator values (upper-H-mode-threshold regime).

const SHOT_ID        = 206804                   # last slot of history_morning.npz
const SHOT_T_NOW     = 3.500                    # s into the shot
const SHOT_IP        = 1.20e6                   # A  — plasma current
const SHOT_PR15V     = 1.18e6                   # A  — PR15V regulator current
const SHOT_BT        = -1.90                    # T  — vacuum toroidal field at R0
const SHOT_POHM      = 1.8e6                    # W  — ohmic power (zmq_Pohm passthrough)
const SHOT_PNBI      = 6.5e6                    # W  — total NBI power (zmq_Pnbi; /1e6 → MW at the FPE input)
const SHOT_PECH      = 2.3e6                    # W  — total ECH power (zmq_Pech; /1e6 → MW at the FPE input)
const SHOT_TNBI      = 5.2                      # N·m — NBI torque (zmq_Tnbi passthrough)
const SHOT_GAS_CALS  = (gasa=8.0, gasb=4.0, gasc=0.0, gasd=0.0, gase=0.0)  # Torr·L/s
const SHOT_I_COIL    = begin
    v = zeros(Float64, 24)
    v[1]  =  80.0;  v[4]  = -80.0                           # ecoila / ecoilb
    v[7:15] .= [-500.0,  -400.0,   300.0,    50.0, 200.0,
                -150.0,   -60.0,   -30.0,   -20.0]          # f1a..f9a
    v[16:24] .= [ 500.0,   400.0,  -300.0,   -50.0, -200.0,
                  150.0,    60.0,    30.0,    20.0]         # f1b..f9b
    v
end

const HISTORY_NPZ = joinpath(dirname(@__DIR__), "data", "pedestal_nn", "history_morning.npz")

# ── Helpers ──────────────────────────────────────────────────────────────────

_aux_push!(aux, key, t, value) = begin
    rec = get!(aux, key, (times=Float64[], values=Float64[]))
    push!(rec.times, t)
    push!(rec.values, value)
end

_aux_push_vec!(aux, key, t, vec) = begin
    rec = get!(aux, key, (times=Float64[], values=Vector{Float64}[]))
    push!(rec.times, t)
    push!(rec.values, collect(Float64, vec))
end

function seed_live_aux!(dd; t::Real=SHOT_T_NOW)
    aux = getfield(dd, :_aux)
    _aux_push!(aux, :zmq_Ip_avg,  t, SHOT_IP)
    _aux_push!(aux, :zmq_pr15v,   t, SHOT_PR15V)
    _aux_push!(aux, :zmq_Pohm,    t, SHOT_POHM)
    _aux_push!(aux, :zmq_Pnbi,    t, SHOT_PNBI)
    _aux_push!(aux, :zmq_Pech,    t, SHOT_PECH)
    _aux_push!(aux, :zmq_Tnbi,    t, SHOT_TNBI)
    _aux_push!(aux, :zmq_gasa_cal, t, SHOT_GAS_CALS.gasa)
    _aux_push!(aux, :zmq_gasb_cal, t, SHOT_GAS_CALS.gasb)
    _aux_push!(aux, :zmq_gasc_cal, t, SHOT_GAS_CALS.gasc)
    _aux_push!(aux, :zmq_gasd_cal, t, SHOT_GAS_CALS.gasd)
    _aux_push!(aux, :zmq_gase_cal, t, SHOT_GAS_CALS.gase)
    _aux_push_vec!(aux, :zmq_I_coil, t, SHOT_I_COIL)
    return aux
end

function seed_history_buffer!(dd; npz_path::AbstractString=HISTORY_NPZ,
                              shot_id::Integer=SHOT_ID)
    # Load the NPZ verbatim into dd._aux[:nn_history_buffer] matching
    # HistoryManager.save_state's schema.
    data = NPZ.npzread(npz_path)
    aux = getfield(dd, :_aux)
    aux[:nn_history_buffer] = (
        history_stats      = Float32.(data["history_stats"]),
        history_masks      = Float32.(data["history_masks"]),
        shot_ids           = Int64.(data["shot_ids"]),
        bzn_seconds        = Float64.(data["bzn_seconds"]),
        disrupt_seconds    = Float64.(data["disrupt_seconds"]),
        disrupt_coverage   = Float64.(data["disrupt_coverage"]),
        n_history_features = Int64.(data["n_history_features"]),
    )
    return aux[:nn_history_buffer]
end

function channel_report(nn::FUSE.PedestalNN, sequences::AbstractMatrix,
                        mask::AbstractVector; sigdigits::Int=4)
    # Per-channel provenance: live? value, training mean, |Δ| in training-σ units.
    μ = nn.signal_means
    σ = nn.signal_stds
    ncols = size(sequences, 2)
    @assert ncols == length(nn.signal_names) == length(μ) == length(σ) == length(mask)
    live_names  = String[]
    mean_names  = String[]
    rows = String[]
    for i in 1:ncols
        v   = sequences[1, i]
        m   = μ[i]
        s   = σ[i] > 0 ? σ[i] : one(eltype(σ))
        dz  = (v - m) / s
        is_live = mask[i] >= 0.5f0
        tag = is_live ? "LIVE " : "mean "
        push!(rows, @sprintf("  %s  %-12s  val=% .4g   mean=% .4g   Δ/σ=% .3g",
                             tag, nn.signal_names[i], v, m, dz))
        (is_live ? live_names : mean_names)
        push!(is_live ? live_names : mean_names, nn.signal_names[i])
    end
    return rows, live_names, mean_names
end

# ── Scenario ─────────────────────────────────────────────────────────────────

@testset "end-to-end NN pedestal smoke on shot $(SHOT_ID)" begin
    isfile(HISTORY_NPZ) || error("smoke_pedestal_nn_shot: history NPZ missing at $HISTORY_NPZ")

    dd = IMAS.dd()
    dd.global_time = SHOT_T_NOW
    # Seed Bt through IMAS so build_fpe_sequences_from_aux picks it up.
    # vacuum_toroidal_field has its own time array; stamp a single sample.
    dd.equilibrium.time = [SHOT_T_NOW]
    dd.equilibrium.vacuum_toroidal_field.b0 = [SHOT_BT]

    seed_live_aux!(dd)
    seed_history_buffer!(dd)

    nn = FUSE.load_pedestal_nn()
    @info "loaded PedestalNN" bundles=length(nn.bundles) onnx_dir=nn.onnx_dir

    sequences, mask = FUSE.build_fpe_sequences_from_aux(nn, dd)
    rows, live, mean_only = channel_report(nn, sequences, mask)

    println("\n── FPE input channel provenance ─────────────────────────────")
    println("LIVE  = live value from dd._aux[:zmq_*] / IMAS")
    println("mean  = channel left at training mean (z-scores to 0)\n")
    foreach(println, rows)
    println(@sprintf("\nsummary: %d/32 live, %d on training mean",
                     length(live), length(mean_only)))

    # Expected live set under this branch (post colleague's Pech + Tnbi wiring):
    # ip, ipspr15v, pohm, pinj, ech_total, tinj, gasa..gase_cal, ecoila, ecoilb,
    # f1a..f9b, bt = 3 + 4 + 5 + 2 + 18 + 1 = 33?  -> 32 channels total, all wired
    # except gasc/gasd/gase which are 0.0 in the scenario (still counted live).
    expected_live = Set(vcat(
        ["ip", "ipspr15v", "pohm", "pinj", "ech_total", "tinj"],
        ["gasa_cal", "gasb_cal", "gasc_cal", "gasd_cal", "gase_cal"],
        ["ecoila", "ecoilb"],
        ["f1a","f2a","f3a","f4a","f5a","f6a","f7a","f8a","f9a"],
        ["f1b","f2b","f3b","f4b","f5b","f6b","f7b","f8b","f9b"],
        ["bt"],
    ))
    live_set = Set(live)
    missing_live = setdiff(expected_live, live_set)
    @test isempty(missing_live)
    if !isempty(missing_live)
        @warn "Channels we expected live but got training-mean" channels=collect(missing_live)
    end
    @test length(live) >= 30  # all 32 if the scenario fills everything; 30 leaves a tiny margin

    # ── History buffer round-trip ─────────────────────────────────────────
    hist_buffer = FUSE.mse_history_from_aux(dd, nn)
    hist_hs, hist_hm, hist_aux = hist_buffer
    @test size(hist_hs)  == (1, 50, 458)
    @test size(hist_hm)  == (1, 50)
    @test size(hist_aux) == (1, 3)
    println("\n── MSE history ──────────────────────────────────────────────")
    println(@sprintf("  history_stats     : %s   range=[%.3g, %.3g]",
                     string(size(hist_hs)), minimum(hist_hs), maximum(hist_hs)))
    println(@sprintf("  history_masks     : sum=%.0f / %d",
                     sum(hist_hm), length(hist_hm)))
    println(@sprintf("  aux [bzn,dis,cov] : [% .3g, % .3g, % .3g]",
                     hist_aux[1], hist_aux[2], hist_aux[3]))

    # ── Predict ──────────────────────────────────────────────────────────
    out = FUSE.predict_pedestal(nn; sequences, signal_mask=mask, history=hist_buffer)

    # Mean across the FPE time window (same reduction ActorPedestal._step uses).
    ne   = mean(filter(isfinite, vec(out.edens_ped)))
    te   = mean(filter(isfinite, vec(out.te_ped)))
    ti   = mean(filter(isfinite, vec(out.ti_ped)))
    trot = mean(filter(isfinite, vec(out.t_rot_ped)))

    println("\n── Predictions (mean across FPE window) ─────────────────────")
    println(@sprintf("  ne_ped      = %7.3f × 10¹⁹ m⁻³", ne))
    println(@sprintf("  Te_ped      = %7.3f keV",        te))
    println(@sprintf("  Ti_ped      = %7.3f keV",        ti))
    println(@sprintf("  T_rot       = %7.3f krad/s",     trot))
    println(@sprintf("  hmode_prob  = %7.3f     (logit=% .3f)",
                     out.hmode_prob, out.hmode_logit))
    println(@sprintf("  is_h_mode   = %s", out.is_h_mode))

    # Physical plausibility envelopes (same ranges as the existing smoke test).
    @test 0.5  < ne   < 10.0
    @test 0.05 < te   < 3.0
    @test 0.05 < ti   < 5.0
    @test -50  < trot < 50
    @test 0.0  <= out.hmode_prob <= 1.0
end
