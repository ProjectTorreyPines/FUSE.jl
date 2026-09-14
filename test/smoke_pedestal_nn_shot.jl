# End-to-end smoke test for the fuse29 pedestal predictor on a synthesized
# DIII-D shot: exercises the live ZMQ → actuator mapping (all 29 channels),
# the per-shot recurrent-state tracker and the H-mode gate — the exact path
# ActorPedestal takes when `par.ne_from == :nn_predictor`, without needing an
# equilibrium or core profiles.
#
# The scenario is a beam-power ramp through the L→H threshold with realistic
# DIII-D actuator values; tweak the `SHOT_*` constants to explore sensitivity.
# The smoke is intentionally *not* gated on exact prediction values — it
# asserts physical plausibility ranges only.
#
# Run from FUSE root:
#     julia --project=. test/smoke_pedestal_nn_shot.jl
#
# The ONNX bundle is resolved through `FUSE_PEDESTAL_NN_DIR` (auto-fetched from
# HuggingFace if absent); see `data/pedestal_nn/README.md`.

using Printf
using Test
using Statistics: mean

using FUSE
import IMAS

# ── Test scenario ────────────────────────────────────────────────────────────

const SHOT_ID      = 206804
const SHOT_T0      = 0.5                        # s, first tick
const SHOT_TICKS   = 60                         # 3 s at 50 ms
const SHOT_BT      = -1.90                      # T  — signed vacuum toroidal field
const SHOT_PNBI_LO = 1.0e6                      # W  — NBI before the ramp
const SHOT_PNBI_HI = 6.5e6                      # W  — NBI after the ramp
const SHOT_PECH    = 0.0                        # W  — ECH
const SHOT_GAS     = (gasa=8.0, gasb=4.0, gasc=0.0, gasd=0.0, gase=0.0)  # Torr·L/s
const SHOT_I_COIL  = begin                      # A, PCS 24-vector (median-like DIII-D shape)
    v = zeros(Float64, 24)
    v[1] = -6900.0; v[4] = -6800.0                          # ecoila / ecoilb
    v[7:15] .= [-1500.0, -1440.0, -150.0, 1140.0, 1270.0, -2900.0, -2700.0, 1150.0, -580.0]   # f1a..f9a
    v[16:24] .= [-1460.0, -1450.0, -1170.0, 1390.0, 1570.0, -2840.0, -3040.0, 1740.0, -210.0] # f1b..f9b
    v
end

_aux_push!(aux, key, t, value) = begin
    rec = get!(aux, key, (times=Float64[], values=Float64[]))
    push!(rec.times, t); push!(rec.values, value)
end
_aux_push_vec!(aux, key, t, vec) = begin
    rec = get!(aux, key, (times=Float64[], values=Vector{Float64}[]))
    push!(rec.times, t); push!(rec.values, collect(Float64, vec))
end

pnbi_at(t) = t < SHOT_T0 + 1.0 ? SHOT_PNBI_LO : SHOT_PNBI_HI

function seed_live_aux!(dd, t)
    aux = getfield(dd, :_aux)
    _aux_push!(aux, :zmq_Pnbi, t, pnbi_at(t))
    _aux_push!(aux, :zmq_Pech, t, SHOT_PECH)
    for (k, v) in pairs(SHOT_GAS)
        _aux_push!(aux, Symbol("zmq_$(k)_cal"), t, v)
    end
    _aux_push_vec!(aux, :zmq_I_coil, t, SHOT_I_COIL)
end

@testset "end-to-end fuse29 smoke on synthetic shot $(SHOT_ID)" begin
    nn = FUSE.load_pedestal_nn()
    @info "loaded $(nn)"

    dd = IMAS.dd()
    dd.equilibrium.time = [SHOT_T0]   # coordinate of vacuum_toroidal_field.b0
    dd.equilibrium.vacuum_toroidal_field.b0 = [SHOT_BT]

    # channel provenance at the first tick
    dd.global_time = SHOT_T0
    seed_live_aux!(dd, SHOT_T0)
    u, live = FUSE.build_fuse29_actuators(nn, dd; source=:zmq, time=SHOT_T0)
    println("\n── fuse29 actuator provenance (t = $(SHOT_T0) s) ──────────────────")
    ch = nn.actuator_ranges["channels"]
    for (i, name) in enumerate(nn.actuator_names)
        println(@sprintf("  %s  %-10s val=% .4g   [p1=% .4g, p99=% .4g] %s",
                         live[i] ? "LIVE  " : "median", name, u[i],
                         ch[name]["pct"]["1"], ch[name]["pct"]["99"], ch[name]["units"]))
    end
    expected_live = Set(vcat(["pinj", "ech_total"], ["gas$(c)_cal" for c in "abcde"],
                             ["ecoila", "ecoilb"], ["f$(i)a" for i in 1:9], ["f$(i)b" for i in 1:9], ["bt"]))
    live_set = Set(nn.actuator_names[live])
    @test isempty(setdiff(expected_live, live_set))
    @test setdiff(Set(nn.actuator_names), live_set) == Set(["tinj"])   # no core_sources in this scenario
    @test isempty(FUSE.fuse29_check_units(nn, u))

    # step through the shot on the model's grid, exactly as fuse29_predict! does
    tr = FUSE.Fuse29Tracker(nn, SHOT_T0 - nn.period_s)
    trace = Vector{Float32}[]
    gate = Bool[]
    for k in 0:SHOT_TICKS-1
        t = SHOT_T0 + k * nn.period_s
        dd.global_time = t
        seed_live_aux!(dd, t)
        pred, h, n = FUSE.fuse29_advance!(tr, nn, t, τ -> FUSE.build_fuse29_actuators(nn, dd; source=:zmq, time=τ)[1])
        @test n == 1
        push!(trace, pred); push!(gate, h)
    end
    T = reduce(hcat, trace)'   # (ticks, 9)

    println("\n── fuse29 trace (every 10th tick) ─────────────────────────────────")
    println(@sprintf("  %6s %7s %8s %8s %8s %8s %8s %8s %8s %6s %4s", "t[s]", "ne", "te", "ti", "trot", "neped", "teped", "rhosym", "netop", "hmode", "gate"))
    for k in 1:10:SHOT_TICKS
        t = SHOT_T0 + (k - 1) * nn.period_s
        println(@sprintf("  %6.2f %7.3f %8.3f %8.3f %8.2f %8.3f %8.3f %8.3f %8.3f %6.3f %4s",
                         t, T[k, 1], T[k, 2], T[k, 3], T[k, 4], T[k, 5], T[k, 6], T[k, 7], T[k, 8], T[k, 9], gate[k] ? "H" : "L"))
    end

    # Physical plausibility envelopes (IO_CONTRACT typical values).
    settled = 21:SHOT_TICKS
    @test all(0.3 .< T[settled, 1] .< 10.0)          # ne, 1e19 m^-3
    @test all(0.05 .< T[settled, 2] .< 3.0)          # te at rho 0.85, keV
    @test all(0.05 .< T[settled, 3] .< 5.0)          # ti, keV
    @test all(-100 .< T[settled, 4] .< 100)          # trot, krad/s
    @test all(0 .<= T[:, 9] .<= 1)                   # hmode probability
    @test count(i -> gate[i] != gate[i-1], 2:SHOT_TICKS) <= 2   # the Schmitt gate does not chatter

    # gated record as the actor sees it
    pr = FUSE.fuse29_prediction(nn, trace[end], gate[end]; time=SHOT_T0 + (SHOT_TICKS - 1) * nn.period_s, n_ticks=1)
    println("\n  final: ", (; ne=pr.ne, te_ped=pr.te_ped, ti_ped=pr.ti_ped, hmode_prob=pr.hmode_prob, is_h_mode=pr.is_h_mode, rho_sym=pr.rho_sym))
    @test pr.is_h_mode == gate[end]
    @test pr.is_h_mode ? isfinite(pr.rho_sym) : isnan(pr.rho_sym)
end
