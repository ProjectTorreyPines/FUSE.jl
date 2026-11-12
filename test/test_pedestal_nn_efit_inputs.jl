#= ============================================================ =#
#  fuse_efit inputs on analytic equilibria: the 14 boundary-shape  #
#  quantities of fuse_efit_shape (minor radius, elongation,        #
#  triangularity, Shafranov shift, Luce squareness, drsep with     #
#  its sign, clip and saturation, strike points and their          #
#  presence flags, inner gap), and build_fuse_efit_inputs on a     #
#  synthetic dd (window mean and efit_valid of the :replay path,   #
#  Ip in MA, loop voltage from the boundary flux, actuators).      #
#  The shape tests need no ONNX bundle; the builder tests do       #
#  (FUSE_PEDESTAL_EFIT_NN_DIR) and are skipped otherwise.          #
#    cd FUSE && julia --project test/test_pedestal_nn_efit_inputs.jl
#= ============================================================ =#
using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
using Test
using FUSE
import IMAS

# ψ(R, Z) = x² + f(y), x = R − R0, y = Z, f(y) = y²/κ² − y⁴/(2 κ² yx²) + c·y
# Two X-points near (R0, ±yx). c = 0 is an exact double null; c > 0 lowers ψ at the bottom, so the lower
# X-point becomes the active one (lower single null, drsep < 0); c < 0 gives an upper single null.
const R0, KAP, YX = 1.7, sqrt(2.0), 1.2
f_of(y, c) = y^2 / KAP^2 - y^4 / (2 * KAP^2 * YX^2) + c * y
df_of(y, c) = 2y / KAP^2 - 2y^3 / (KAP^2 * YX^2) + c

function newton(g, dg, y)
    for _ in 1:100
        y -= g(y) / dg(y)
    end
    return y
end

function analytic(c)
    d2f(y) = 2 / KAP^2 - 6y^2 / (KAP^2 * YX^2)
    ya = newton(y -> df_of(y, c), d2f, 0.0)          # magnetic axis
    yl = newton(y -> df_of(y, c), d2f, -YX)          # lower X-point
    yu = newton(y -> df_of(y, c), d2f, YX)           # upper X-point
    psi_l, psi_u = f_of(yl, c), f_of(yu, c)
    return (ya=ya, yl=yl, yu=yu, psi_axis=f_of(ya, c), psi_l=psi_l, psi_u=psi_u, psi_x=min(psi_l, psi_u))
end

# equilibrium time slice with ψ on a rectangular grid and the boundary just inside the active separatrix
function synthetic_eqt!(eqt, c; n_boundary=180, ip=1.0e6, psi_offset=0.0)
    an = analytic(c)
    psi_b = an.psi_axis + (an.psi_x - an.psi_axis) * (1 - 1e-4)
    r = collect(range(R0 - 1.0, R0 + 1.0, 129))
    z = collect(range(-1.7, 1.7, 171))
    p2d = resize!(eqt.profiles_2d, 1)[1]
    p2d.grid_type.index = 1
    p2d.grid.dim1 = r
    p2d.grid.dim2 = z
    p2d.psi = [(R - R0)^2 + f_of(Z, c) + psi_offset for R in r, Z in z]
    eqt.profiles_1d.psi = [an.psi_axis + psi_offset, psi_b + psi_offset]
    eqt.global_quantities.magnetic_axis.r = R0
    eqt.global_quantities.magnetic_axis.z = an.ya
    eqt.global_quantities.psi_axis = an.psi_axis + psi_offset
    eqt.global_quantities.psi_boundary = psi_b + psi_offset
    eqt.global_quantities.ip = ip
    # boundary: x = ±sqrt(psi_b − f(y)) between the two zeros of psi_b − f(y) around the axis
    g(y) = psi_b - f_of(y, c)
    bisect(lo, hi) = begin
        for _ in 1:200
            mid = (lo + hi) / 2
            (g(lo) > 0) == (g(mid) > 0) ? (lo = mid) : (hi = mid)
        end
        (lo + hi) / 2
    end
    y_lo = bisect(an.ya, an.yl)
    y_hi = bisect(an.ya, an.yu)
    # an odd number of heights, so that the outermost point is a vertex (at the axis height when c = 0);
    # with an even number the two candidates tie and the discrete-vertex squareness of the upper and lower
    # quadrants differs by ~0.02, the resolution bias the model's INPUT_CONTRACT §4.3 describes
    ys = collect(range(y_lo, y_hi, n_boundary ÷ 2 + 1))
    xs = [sqrt(max(g(y), 0.0)) for y in ys]
    eqt.boundary.outline.r = vcat(R0 .+ xs, R0 .- reverse(xs)[2:end-1])
    eqt.boundary.outline.z = vcat(ys, reverse(ys)[2:end-1])
    return (an..., psi_b=psi_b, y_lo=y_lo, y_hi=y_hi, a=sqrt(psi_b - an.psi_axis))
end

function synthetic_eqt(c; kw...)
    dd = IMAS.dd()
    eqt = resize!(dd.equilibrium.time_slice, 1)[1]
    info = synthetic_eqt!(eqt, c; kw...)
    return eqt, info, dd
end

# rectangular wall around both X-points
const WALL_R = [R0 - 0.9, R0 + 0.9, R0 + 0.9, R0 - 0.9, R0 - 0.9]
const WALL_Z = [-1.5, -1.5, 1.5, 1.5, -1.5]

@testset "fuse_efit inputs" begin

@testset "fuse_efit_shape: double null" begin
    eqt, an, _ = synthetic_eqt(0.0)
    s = FUSE.fuse_efit_shape(eqt, WALL_R, WALL_Z)
    @test s !== nothing
    @test s.aminor ≈ an.a atol = 1e-3
    @test s.rsurf ≈ R0 atol = 1e-6
    @test s.zsurf ≈ 0.0 atol = 1e-6
    # the X-points are the top and the bottom of the boundary, as in EFIT: the outline stops ~6 mm short of them
    @test an.y_hi < an.yu - 3e-3
    @test s.kappa ≈ (an.yu - an.yl) / (2 * s.aminor) atol = 5e-4
    @test abs(s.tritop) < 1e-3 && abs(s.tribot) < 1e-3       # the top and bottom of this boundary sit at R0
    @test abs(s.shafranov) < 1e-6
    @test s.gapin ≈ (R0 - an.a) - (R0 - 0.9) atol = 1e-3
    # both X-points found, on the same flux surface: drsep = 0
    @test s.x_lower !== nothing && s.x_upper !== nothing
    @test s.x_lower[1] ≈ R0 atol = 1e-4
    @test s.x_lower[2] ≈ -YX atol = 1e-4
    @test s.x_upper[2] ≈ YX atol = 1e-4
    @test s.drsep ≈ 0.0 atol = 1e-3
    # all four legs reach the wall: x² = ψ_x − f(∓1.5) on the bottom and the top
    leg = sqrt(an.psi_x - f_of(-1.5, 0.0))
    @test s.rvsout_present && s.rvsou_present && s.rvsiu_present
    @test s.rvsout ≈ R0 + leg atol = 1e-2
    @test s.rvsou ≈ R0 + leg atol = 1e-2
    @test s.rvsiu ≈ R0 - leg atol = 1e-2
    # up-down symmetric and pointed at the X-points: equal, negative outer squareness
    @test s.zeta_ou ≈ s.zeta_ol atol = 1e-3
    @test -0.6 < s.zeta_ou < 0.0

    # an outline that rounds the X-points off, as a traced contour does: without the X-points its lowest and
    # highest vertices sit off R0 (a triangularity of ±0.05 or more) and 3 cm inside
    keep = findall(z -> an.y_lo + 0.03 < z < an.y_hi - 0.03, eqt.boundary.outline.z)
    eqt.boundary.outline.r = eqt.boundary.outline.r[keep]
    eqt.boundary.outline.z = eqt.boundary.outline.z[keep]
    k_low = argmin(eqt.boundary.outline.z)
    @test abs(eqt.boundary.outline.r[k_low] - R0) > 0.05 * an.a
    rounded = FUSE.fuse_efit_shape(eqt, WALL_R, WALL_Z)
    @test rounded.kappa ≈ s.kappa atol = 1e-6
    @test rounded.zsurf ≈ s.zsurf atol = 1e-6
    @test abs(rounded.tritop) < 1e-3 && abs(rounded.tribot) < 1e-3
end

@testset "fuse_efit_shape: squareness is IMAS's Luce squareness (0 for an ellipse)" begin
    eqt, _, _ = synthetic_eqt(0.0)
    θ = range(0, 2π, 361)[1:end-1]                            # 1° vertices: the extrema and the 45° diagonals fall on vertices
    eqt.boundary.outline.r = R0 .+ 0.6 .* cos.(θ)
    eqt.boundary.outline.z = 0.6 * 1.7 .* sin.(θ)
    s = FUSE.fuse_efit_shape(eqt)
    @test s.aminor ≈ 0.6 atol = 1e-9
    @test s.kappa ≈ 1.7 atol = 1e-9
    @test abs(s.zeta_ou) < 1e-6 && abs(s.zeta_ol) < 1e-6
    @test isnan(s.gapin) && !s.rvsout_present               # no wall: no gap, no strike points
    # a D shape: same numbers as IMAS.luce_squareness / miller_R_a_κ_δ_ζ on the closed polygon
    pr = R0 .+ 0.6 .* cos.(θ .+ 0.4 .* sin.(θ))
    pz = 0.6 * 1.7 .* sin.(θ .+ 0.1 .* sin.(2 .* θ))
    eqt.boundary.outline.r = pr
    eqt.boundary.outline.z = pz
    s = FUSE.fuse_efit_shape(eqt)
    m = IMAS.miller_R_a_κ_δ_ζ(vcat(pr, pr[1]), vcat(pz, pz[1]))
    @test s.rsurf ≈ m.R0 && s.aminor ≈ m.a && s.kappa ≈ m.κ
    @test s.tritop ≈ m.δu && s.tribot ≈ m.δl
    @test s.zeta_ou ≈ m.ζou atol = 1e-12
    @test s.zeta_ol ≈ m.ζol atol = 1e-12
    @test s.tritop > 0.3                                     # and a positive triangularity for a D pointing inward
end

@testset "fuse_efit_shape: drsep sign, clip and saturation; strike-point presence" begin
    # lower single null with the second separatrix 1-2 cm out: negative drsep, all strike points there
    c = 0.01
    eqt, an, _ = synthetic_eqt(c)
    @test an.psi_l < an.psi_u
    s = FUSE.fuse_efit_shape(eqt, WALL_R, WALL_Z)
    expected = sqrt(an.psi_l - an.psi_axis) - sqrt(an.psi_u - an.psi_axis)
    @test -0.1 < expected < -0.005
    @test s.drsep ≈ expected atol = 2e-3
    @test s.drsep_raw == s.drsep
    @test s.rvsout_present && s.rvsou_present && s.rvsiu_present
    @test s.rvsout ≈ R0 + sqrt(an.psi_l - f_of(-1.5, c)) atol = 1e-2
    @test s.rvsou ≈ R0 + sqrt(an.psi_u - f_of(1.5, c)) atol = 1e-2

    # upper single null: mirror image, positive drsep
    eqt, an, _ = synthetic_eqt(-c)
    s = FUSE.fuse_efit_shape(eqt, WALL_R, WALL_Z)
    @test s.drsep ≈ -expected atol = 2e-3

    # far lower single null: clipped to −0.1 m
    c = 0.12
    eqt, an, _ = synthetic_eqt(c)
    s = FUSE.fuse_efit_shape(eqt, WALL_R, WALL_Z)
    @test sqrt(an.psi_l - an.psi_axis) - sqrt(an.psi_u - an.psi_axis) < -0.1
    @test s.drsep == -0.1 && s.drsep_raw < -0.1

    # no upper X-point inside the wall: drsep saturates, the upper strike points are absent (NaN, never 0)
    eqt, an, _ = synthetic_eqt(0.01)
    wall_z = [-1.5, -1.5, 1.1, 1.1, -1.5]                   # above the plasma top (1.04 m), below the upper X-point (1.20 m)
    s = FUSE.fuse_efit_shape(eqt, WALL_R, wall_z)
    @test s.x_upper === nothing && s.x_lower !== nothing
    @test s.drsep_raw == -0.4 && s.drsep == -0.1
    @test s.rvsout_present && !s.rvsou_present && !s.rvsiu_present
    @test isnan(s.rvsou) && isnan(s.rvsiu) && isfinite(s.rvsout)

    # a slice without a usable boundary is not a plasma
    eqt, _, _ = synthetic_eqt(0.0)
    eqt.boundary.outline.r = eqt.boundary.outline.r[1:5]
    eqt.boundary.outline.z = eqt.boundary.outline.z[1:5]
    @test FUSE.fuse_efit_shape(eqt, WALL_R, WALL_Z) === nothing
end

have_bundle = withenv(FUSE.PEDESTAL_NN_AUTOFETCH_ENV => "0") do
    FUSE._fuse_efit_find_bundle_dir(abspath(FUSE.resolve_pedestal_efit_nn_dir()), :single) !== nothing
end
if !have_bundle
    @warn "fuse_efit bundle not on disk; skipping the build_fuse_efit_inputs tests"
else
    nn = FUSE.load_pedestal_efit_nn()
    col(name) = findfirst(==(name), nn.input_names)
    flag(name) = findfirst(==(name), nn.indicator_names)

    # experimental-like dd: equilibrium slices every 20 ms from 0.5 to 1.5 s, the boundary flux rising at 0.4 Wb/s
    function replay_dd(; dpsi_dt=0.4)
        rdd = IMAS.dd()
        times = collect(0.5:0.02:1.5)
        resize!(rdd.equilibrium.time_slice, length(times))
        for (k, t) in enumerate(times)
            eqt = rdd.equilibrium.time_slice[k]
            eqt.time = t
            synthetic_eqt!(eqt, 0.01; ip=1.0e6 + 1.0e5 * (t - 1.0), psi_offset=dpsi_dt * t)
        end
        rdd.equilibrium.time = times
        return rdd
    end

    function host_dd()
        dd = IMAS.dd()
        outline = resize!(resize!(resize!(dd.wall.description_2d, 1)[1].limiter.unit, 1), 1)[1].outline
        outline.r = WALL_R
        outline.z = WALL_Z
        t = collect(0.0:0.01:2.0)
        valve = resize!(dd.gas_injection.valve, 1)[1]
        valve.name = "GASA"
        valve.flow_rate.time = t
        valve.flow_rate.data = fill(10.0 * FUSE.PA_M3_PER_TORR_L, length(t))
        return dd
    end

    @testset "build_fuse_efit_inputs: :replay state" begin
        dd = host_dd()
        rdd = replay_dd()
        x, ind, live = FUSE.build_fuse_efit_inputs(nn, dd; source=:dd, state_source=:replay, time=1.0, replay_dd=rdd)
        @test ind[flag("efit_valid")] == 1 && ind[flag("vloop_valid")] == 1
        @test ind[flag("rvsout_present")] == 1 && ind[flag("rvsou_present")] == 1 && ind[flag("rvsiu_present")] == 1
        # window (0.95, 1.0]: slices at 0.96, 0.98, 1.00 -> mean Ip = 1e6 + 1e5 * (0.98 − 1.0), in MA
        @test x[col("ip")] ≈ (1.0e6 - 2.0e3) / 1e6 atol = 1e-6
        # dψ_boundary/dt (COCOS 11)
        @test x[col("vloop")] ≈ 0.4 atol = 1e-4
        one, _, _ = synthetic_eqt(0.01)
        s = FUSE.fuse_efit_shape(one, WALL_R, WALL_Z)
        for name in FUSE.FUSE_EFIT_SHAPE_NAMES
            @test x[col(name)] ≈ Float32(getproperty(s, Symbol(name))) atol = 1e-5
        end
        @test all(live[col(name)] for name in FUSE.FUSE_EFIT_SHAPE_NAMES)
        # actuators: the valve the dd carries, zero power for a dd without beams or gyrotrons, no other valves
        @test x[col("gasa_cal")] ≈ 10.0 atol = 1e-4
        @test x[col("pinj")] == 0 && x[col("tinj")] == 0 && x[col("ech_mw")] == 0
        @test !live[col("gasb_cal")] && isnan(x[col("gasb_cal")])
        @test isempty(FUSE.fuse_efit_check_units(nn, x, ind, Float32.(live)))

        # no slice in the tick window: no geometry, and no hold-last
        x, ind, live = FUSE.build_fuse_efit_inputs(nn, dd; source=:dd, state_source=:replay, time=1.7, replay_dd=rdd)
        @test ind[flag("efit_valid")] == 0 && ind[flag("vloop_valid")] == 0
        @test all(isnan(x[col(name)]) for name in FUSE.FUSE_EFIT_SHAPE_NAMES)
        @test !live[col("ip")] && !live[col("vloop")]
        # a reversed loop voltage keeps its sign
        x, _, _ = FUSE.build_fuse_efit_inputs(nn, dd; source=:dd, state_source=:replay, time=1.0, replay_dd=replay_dd(; dpsi_dt=-0.3))
        @test x[col("vloop")] ≈ -0.3 atol = 1e-4
        # :replay without a replay dd
        _, ind, _ = FUSE.build_fuse_efit_inputs(nn, dd; source=:dd, state_source=:replay, time=1.0)
        @test ind[flag("efit_valid")] == 0 && ind[flag("vloop_valid")] == 0
        @test_throws ErrorException FUSE.build_fuse_efit_inputs(nn, dd; source=:dd, state_source=:nope, time=1.0)
    end

    @testset "build_fuse_efit_inputs: :simulation state" begin
        dd = host_dd()
        for (k, t) in enumerate((0.95, 1.0))
            eqt = resize!(dd.equilibrium.time_slice, k)[k]
            eqt.time = t
            synthetic_eqt!(eqt, 0.01; ip=0.9e6, psi_offset=0.25 * t)
        end
        dd.equilibrium.time = [0.95, 1.0]
        x, ind, live = FUSE.build_fuse_efit_inputs(nn, dd; source=:dd, state_source=:simulation, time=1.0)
        @test ind[flag("efit_valid")] == 1
        # no core_profiles in this dd: the plasma current is the equilibrium's, and the loop voltage falls back on
        # the boundary-flux change of the last two slices
        @test x[col("ip")] ≈ 0.9 atol = 1e-6
        @test x[col("drsep")] < 0 && x[col("aminor")] > 0.5
        @test ind[flag("vloop_valid")] == 1
        @test x[col("vloop")] ≈ 0.25 atol = 1e-4
        # the equilibrium is too old for this tick
        _, ind, _ = FUSE.build_fuse_efit_inputs(nn, dd; source=:dd, state_source=:simulation, time=1.5)
        @test ind[flag("efit_valid")] == 0
        # ZMQ: the measured Ip streamed by GSLite wins over the equilibrium's
        aux = getfield(dd, :_aux)
        aux[:zmq_Ip_avg] = (times=[0.9, 1.0], values=[1.1e6, 1.2e6])
        aux[:zmq_Pnbi] = (times=[1.0], values=[2.5e6])
        x, _, live = FUSE.build_fuse_efit_inputs(nn, dd; source=:zmq, state_source=:simulation, time=1.0)
        @test x[col("ip")] ≈ 1.2 atol = 1e-6
        @test x[col("pinj")] ≈ 2.5 atol = 1e-6
        @test !live[col("gasa_cal")]                             # the :zmq path does not read dd.gas_injection
    end

    @testset "pedestal_nn_record! / pedestal_nn_history" begin
        dd = IMAS.dd()
        aux = getfield(dd, :_aux)
        @test FUSE.pedestal_nn_history(dd) === nothing
        x = fill(1.0f0, 25)
        x[3] = NaN32
        pred = Float32.(1:9)
        FUSE.pedestal_nn_record!(aux, :fuse_efit, nn.input_names, nn.indicator_names, nn.head_names, 1.0, x, ones(Float32, 5), pred, false, 1)
        FUSE.pedestal_nn_record!(aux, :fuse_efit, nn.input_names, nn.indicator_names, nn.head_names, 1.0, 2 .* x, ones(Float32, 5), 2 .* pred, true, 1) # re-run
        FUSE.pedestal_nn_record!(aux, :fuse_efit, nn.input_names, nn.indicator_names, nn.head_names, 1.05, x, ones(Float32, 5), pred, true, 1)
        h = FUSE.pedestal_nn_history(dd)
        @test h["model"] == "fuse_efit" && h["time"] == [1.0, 1.05]
        @test h["predictions"]["ne"] == [2.0, 1.0] && h["predictions"]["hmode"] == [18.0, 9.0]
        @test h["inputs"]["pinj"] == [2.0, 1.0] && h["inputs"]["ech_mw"] == [nothing, nothing]
        @test h["is_h_mode"] == [true, true] && h["indicators"]["efit_valid"] == [1.0, 1.0]
        # another model, or a rewind, starts a new record
        FUSE.pedestal_nn_record!(aux, :fuse29, ["a"], String[], nn.head_names, 2.0, [1.0f0], Float32[], pred, false, 1)
        @test FUSE.pedestal_nn_history(dd)["model"] == "fuse29" && FUSE.pedestal_nn_history(dd)["time"] == [2.0]
    end
end

end # fuse_efit inputs
