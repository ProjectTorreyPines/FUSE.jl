using FUSE
using IMAS
using Test
using ZMQ
using ProtoBuf
import FRESCO

# Pull in the generated protobuf types via the same path the actor uses,
# so we encode/decode against the exact stub the actor was compiled against.
const _ZMQ_PB = include(joinpath(dirname(pathof(FUSE)), "actors", "zmq_proto_generated", "zmq_messages_pb.jl"))
using ._ZMQ_PB: FUSERequest, WireDataForFUSE, WireDataFromFUSE, Ack

# Tiny encode/decode helpers — mirror _pb_send/_pb_recv inside the actor.
_zmq_encode(msg) = (io = IOBuffer(); ProtoBuf.encode(ProtoBuf.ProtoEncoder(io), msg); take!(io))
_zmq_decode(::Type{T}, raw) where {T} = ProtoBuf.decode(ProtoBuf.ProtoDecoder(IOBuffer(raw)), T)

# Keyword-style constructor for WireDataForFUSE that survives future field additions:
# uses ProtoBuf.default_values (generated alongside the struct) for any field the
# caller doesn't override.
function _make_WireDataForFUSE(; kwargs...)
    defaults = ProtoBuf.default_values(WireDataForFUSE)
    overrides = NamedTuple(kwargs)
    vals = (haskey(overrides, f) ? overrides[f] : defaults[f] for f in keys(defaults))
    return WireDataForFUSE(vals...)
end

@testset "ActorZMQ round-trip" begin
    # ipc:// + tempname() avoids TCP port collisions on shared CI; Linux-only.
    endpoint = "ipc://" * tempname()

    # --- Round-trip happy path: metadata-only WireDataForFUSE (no psizr) ---
    # This exercises the wire format, the FUSERequest/WireDataForFUSE/Ack handshake,
    # the version handshake (FUSE.SCHEMA_VERSION echoed on both sides), the try/catch
    # wrappers, and the Bt-storage branch — without dragging in a fully-initialized
    # equilibrium (which would require a full FUSE.init).
    ini, act = FUSE.case_parameters(:ITER; init_from=:scalars)
    dd = IMAS.dd()
    # send! reports equilibrium-derived quantities and fetches dd.equilibrium.time_slice[]
    # before its own try/catch, so give it a (bare) slice to land on. The slice is left
    # empty on purpose: send! must fall back to valid=false rather than throw.
    dd.global_time = 0.0
    dd.equilibrium.time = [0.0]
    resize!(dd.equilibrium.time_slice, 1)[1].time = 0.0
    act.ActorZMQ.enabled = true
    act.ActorZMQ.endpoint = endpoint
    act.ActorZMQ.timeout_ms = 5000  # short timeout so a bug fails the test instead of hanging on the default

    server = @async begin
        ctx = ZMQ.Context()
        sock = ZMQ.Socket(ctx, ZMQ.REP)
        ZMQ.bind(sock, endpoint)
        try
            req = _zmq_decode(FUSERequest, ZMQ.recv(sock))
            @test req.status == "ready"
            @test req.schema_version == FUSE.SCHEMA_VERSION
            resp = _make_WireDataForFUSE(;
                sim_time=0.1,
                Bt=2.1,
                has_Bt=true,
                schema_version=FUSE.SCHEMA_VERSION,
            )
            ZMQ.send(sock, _zmq_encode(resp))

            # Then a WireDataFromFUSE arrives; reply with an Ack with ok=true.
            _ = _zmq_decode(WireDataFromFUSE, ZMQ.recv(sock))
            ZMQ.send(sock, _zmq_encode(Ack("ack", true, "")))
        finally
            ZMQ.close(sock)
            ZMQ.close(ctx)
        end
    end

    actor = FUSE.ActorZMQ(dd, act)
    try
        FUSE.receive!(actor)
        @test !isempty(dd.equilibrium.vacuum_toroidal_field.b0)
        @test dd.equilibrium.vacuum_toroidal_field.b0[end] ≈ 2.1
        FUSE.send!(actor)
    finally
        FUSE.disconnect!(actor)
    end
    wait(server)
end

@testset "ActorZMQ rejects schema_version mismatch" begin
    # GSLite replies with schema_version=0 (pre-versioning binary). FUSE must
    # error loudly with both versions in the message and not proceed.
    endpoint = "ipc://" * tempname()
    ini, act = FUSE.case_parameters(:ITER; init_from=:scalars)
    dd = IMAS.dd()
    act.ActorZMQ.enabled = true
    act.ActorZMQ.endpoint = endpoint
    act.ActorZMQ.timeout_ms = 5000

    server = @async begin
        ctx = ZMQ.Context()
        sock = ZMQ.Socket(ctx, ZMQ.REP)
        ZMQ.bind(sock, endpoint)
        try
            _ = _zmq_decode(FUSERequest, ZMQ.recv(sock))
            resp = _make_WireDataForFUSE(; sim_time=0.0, schema_version=Int32(0))
            ZMQ.send(sock, _zmq_encode(resp))
        finally
            ZMQ.close(sock)
            ZMQ.close(ctx)
        end
    end

    actor = FUSE.ActorZMQ(dd, act)
    err = try
        FUSE.receive!(actor)
        nothing
    catch e
        e
    finally
        FUSE.disconnect!(actor)
    end
    wait(server)
    @test err !== nothing
    msg = sprint(showerror, err)
    @test occursin("schema mismatch", msg)
    @test occursin(string(FUSE.SCHEMA_VERSION), msg)
    @test occursin("version 0", msg)
end

@testset "ActorZMQ rejects ack.ok=false" begin
    # First exchange succeeds (so send! has somewhere to run); the second exchange
    # returns Ack(ok=false, error="solver diverged"). FUSE must surface the
    # rejection and tear down the coupled run.
    endpoint = "ipc://" * tempname()
    ini, act = FUSE.case_parameters(:ITER; init_from=:scalars)
    dd = IMAS.dd()
    # this testset calls send!, which fetches dd.equilibrium.time_slice[] outside its
    # own try/catch — without a slice it throws before sending, and the server below
    # blocks forever on recv. Same bare seed as the round-trip testset.
    dd.global_time = 0.0
    dd.equilibrium.time = [0.0]
    resize!(dd.equilibrium.time_slice, 1)[1].time = 0.0
    act.ActorZMQ.enabled = true
    act.ActorZMQ.endpoint = endpoint
    act.ActorZMQ.timeout_ms = 5000

    server = @async begin
        ctx = ZMQ.Context()
        sock = ZMQ.Socket(ctx, ZMQ.REP)
        ZMQ.bind(sock, endpoint)
        try
            _ = _zmq_decode(FUSERequest, ZMQ.recv(sock))
            ZMQ.send(sock, _zmq_encode(_make_WireDataForFUSE(;
                sim_time=0.1,
                Bt=2.1,
                has_Bt=true,
                schema_version=FUSE.SCHEMA_VERSION,
            )))
            _ = _zmq_decode(WireDataFromFUSE, ZMQ.recv(sock))
            ZMQ.send(sock, _zmq_encode(Ack("nack", false, "solver diverged")))
        finally
            ZMQ.close(sock)
            ZMQ.close(ctx)
        end
    end

    actor = FUSE.ActorZMQ(dd, act)
    err = try
        FUSE.receive!(actor)
        FUSE.send!(actor)
        nothing
    catch e
        e
    finally
        FUSE.disconnect!(actor)
    end
    wait(server)
    @test err !== nothing
    msg = sprint(showerror, err)
    @test occursin("GSLite rejected", msg)
    @test occursin("solver diverged", msg)
end

@testset "ActorZMQ stores Ip_avg=0.0 when has_Ip_avg=true" begin
    # Regression for the proto3 zero-elision trap: a legitimate 0 A average during
    # ramp-up / pre-breakdown must land in aux[:zmq_Ip_avg] when GSLite advertises
    # presence with has_Ip_avg=true. The earlier `!= 0.0` check (since removed)
    # would have silently dropped this; the has_* bool is what the GSLite team
    # MUST set on the C++ side, and this test pins that contract.
    endpoint = "ipc://" * tempname()
    ini, act = FUSE.case_parameters(:ITER; init_from=:scalars)
    dd = IMAS.dd()
    act.ActorZMQ.enabled = true
    act.ActorZMQ.endpoint = endpoint
    act.ActorZMQ.timeout_ms = 5000

    server = @async begin
        ctx = ZMQ.Context()
        sock = ZMQ.Socket(ctx, ZMQ.REP)
        ZMQ.bind(sock, endpoint)
        try
            _ = _zmq_decode(FUSERequest, ZMQ.recv(sock))
            ZMQ.send(sock, _zmq_encode(_make_WireDataForFUSE(;
                sim_time=0.0,
                Ip_avg=0.0,
                has_Ip_avg=true,
                schema_version=FUSE.SCHEMA_VERSION,
            )))
        finally
            ZMQ.close(sock)
            ZMQ.close(ctx)
        end
    end

    actor = FUSE.ActorZMQ(dd, act)
    try
        FUSE.receive!(actor)
    finally
        FUSE.disconnect!(actor)
    end
    wait(server)

    aux = getfield(dd, :_aux)
    @test haskey(aux, :zmq_Ip_avg)
    @test aux[:zmq_Ip_avg].values == [0.0]
    @test aux[:zmq_Ip_avg].times == [dd.global_time]
end

# ============================================================================
#  psizr commit: what receive! must leave behind for the rest of the substep
#
#  ActorZMQ writes into dd and then hands off to a substep sequence it does not
#  control (dynamic_plasma_actor.jl): run_sources -> run_sawteeth -> QED ->
#  run_equilibrium. Both tests below pin an invariant that sequence depends on,
#  and both correspond to a deterministic failure seen in coupled runs against
#  real GSLite (fresco_crash_logs, shot 207123).
# ============================================================================

# Analytic psi with a magnetic axis at (1.70, 0.0), increasing outward — enough
# for find_magnetic_axis / find_psi_boundary / flux_surfaces to succeed without
# dragging in a full FUSE.init.
_psi_analytic(r, z; R0=1.70, a=0.62, b=1.05) = ((r - R0) / a)^2 + (z / b)^2 - 1.0

# GSLite's fixed DIII-D 33x33 grid — the hardcoded fallback inside receive!,
# reached because the seeded dd carries a different (65x65) grid.
const _GS_R = collect(range(0.84, 2.54, length=33))
const _GS_Z = collect(range(-1.6, 1.6, length=33))
_gslite_psizr() = vec([_psi_analytic(r, z) for r in _GS_R, z in _GS_Z])

"""
A dd in the state receive! actually meets: a first wall, a 65x65 equilibrium
from a previous FRESCO solve, and NON-ZERO 1-D source profiles inherited from
the previous time slice (new_timeslice! deep-copies them forward).
"""
function _seeded_psizr_dd(; n=65, p0=1.0e4)
    dd = IMAS.dd()

    resize!(dd.wall.description_2d, 1)
    limiter = resize!(dd.wall.description_2d[1].limiter.unit, 1)[1]
    θ = range(0, 2π, length=97)
    limiter.outline.r = 1.70 .+ 0.72 .* cos.(θ)
    limiter.outline.z = 1.25 .* sin.(θ)

    dd.global_time = 0.0
    dd.equilibrium.time = [0.0]
    eqt = resize!(dd.equilibrium.time_slice, 1)[1]
    eqt.time = 0.0
    eqt.global_quantities.vacuum_toroidal_field.r0 = 1.6955
    eqt.global_quantities.vacuum_toroidal_field.b0 = 2.0

    R = collect(range(0.84, 2.54, length=n))
    Z = collect(range(-1.6, 1.6, length=n))
    p2d = resize!(eqt.profiles_2d, 1)[1]
    p2d.grid_type.index = 1
    p2d.grid.dim1 = R
    p2d.grid.dim2 = Z
    p2d.psi = [_psi_analytic(r, z) for r in R, z in Z]

    npsi = 101
    eqt1d = eqt.profiles_1d
    eqt1d.psi = collect(range(-1.0, 0.35, length=npsi))
    psin = (eqt1d.psi .- eqt1d.psi[1]) ./ (eqt1d.psi[end] - eqt1d.psi[1])
    eqt1d.pressure = p0 .* (1.0 .- psin) .^ 2
    eqt1d.dpressure_dpsi = -2.0 * p0 .* (1.0 .- psin) ./ (eqt1d.psi[end] - eqt1d.psi[1])
    eqt1d.f_df_dpsi = fill(-0.1, npsi)
    eqt1d.f = fill(2.0 * 1.6955, npsi)
    eqt.global_quantities.magnetic_axis.r = 1.70
    eqt.global_quantities.magnetic_axis.z = 0.0
    eqt.global_quantities.psi_axis = -1.0
    eqt.global_quantities.psi_boundary = 0.35

    return dd
end

# A GSLite stand-in that answers `n` successive receive! calls, each carrying psizr.
function _serve_psizr(endpoint, n)
    return @async begin
        ctx = ZMQ.Context()
        sock = ZMQ.Socket(ctx, ZMQ.REP)
        ZMQ.bind(sock, endpoint)
        try
            for _ in 1:n
                _ = _zmq_decode(FUSERequest, ZMQ.recv(sock))
                ZMQ.send(sock, _zmq_encode(_make_WireDataForFUSE(;
                    sim_time=0.0,
                    psizr=_gslite_psizr(),
                    schema_version=FUSE.SCHEMA_VERSION,
                )))
            end
        finally
            ZMQ.close(sock)
            ZMQ.close(ctx)
        end
    end
end

function _psizr_actor(dd, endpoint)
    _, act = FUSE.case_parameters(:ITER; init_from=:scalars)
    act.ActorZMQ.enabled = true
    act.ActorZMQ.endpoint = endpoint
    act.ActorZMQ.timeout_ms = 20000
    return FUSE.ActorZMQ(dd, act)
end

@testset "psizr commit keeps the equilibrium source profiles" begin
    # ActorFRESCO builds its source term with FRESCO.PressureJt(dd; ...), which
    # defaults to j_p_from=:equilibrium and therefore reads
    # equilibrium.profiles_1d.pressure / .j_tor. Zeroing those on commit — even
    # only on the first step — hands FRESCO p=0, Jt=0: it solves a vacuum field,
    # finds no closed boundary, and flux_bounds! assigns nothing into a Float64.
    # The inherited profiles are a perfectly good seed for flux_surfaces; keep them.
    endpoint = "ipc://" * tempname()
    dd = _seeded_psizr_dd()
    server = _serve_psizr(endpoint, 1)
    actor = _psizr_actor(dd, endpoint)
    try
        FUSE.receive!(actor)
    finally
        FUSE.disconnect!(actor)
    end
    wait(server)

    eqt1d = dd.equilibrium.time_slice[].profiles_1d
    @test maximum(eqt1d.pressure) > 0.0
    @test maximum(abs, eqt1d.f_df_dpsi) > 0.0
    @test maximum(abs, eqt1d.j_tor) > 0.0

    # the actual consumer: FRESCO must see a non-zero pressure profile
    profile = FRESCO.PressureJt(dd; grid=:psi_norm)
    @test maximum(abs, profile.pressure.(range(0, 1, 21))) > 0.0
end

@testset "successive psizr commits keep profiles_2d.phi" begin
    # Only the first-step branch rebuilds phi (via flux_surfaces). On later steps
    # the rebuild lives in ActorEquilibrium._finalize, which runs at the END of the
    # :run_equilibrium substep — but ActorSimpleNB reads profiles_2d.phi inside
    # :run_sources, which runs BEFORE it. So phi must survive the commit.
    endpoint = "ipc://" * tempname()
    dd = _seeded_psizr_dd()
    server = _serve_psizr(endpoint, 2)
    actor = _psizr_actor(dd, endpoint)
    try
        FUSE.receive!(actor)   # first step: full flux_surfaces
        @test IMAS.hasdata(dd.equilibrium.time_slice[].profiles_2d[1], :phi)
        FUSE.receive!(actor)   # subsequent step: boundary only
        @test IMAS.hasdata(dd.equilibrium.time_slice[].profiles_2d[1], :phi)
    finally
        FUSE.disconnect!(actor)
    end
    wait(server)
end

# TODO: additional cases worth covering once the round-trip is stable:
#   - done=true short-circuit: receive! disconnects without populating dd
#   - psizr length vs nR*nZ mismatch: error message mentions "ActorZMQ" and "psizr"
