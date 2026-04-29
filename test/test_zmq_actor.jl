using FUSE
using IMAS
using Test
using ZMQ
using ProtoBuf

# Pull in the generated protobuf types via the same path the actor uses,
# so we encode/decode against the exact stub the actor was compiled against.
const _ZMQ_PB = include(joinpath(dirname(pathof(FUSE)), "actors", "zmq_proto_generated", "zmq_messages_pb.jl"))
using ._ZMQ_PB: FUSERequest, DataForFUSE, DataFromFUSE, Ack

# Tiny encode/decode helpers — mirror _pb_send/_pb_recv inside the actor.
_zmq_encode(msg) = (io = IOBuffer(); ProtoBuf.encode(ProtoBuf.ProtoEncoder(io), msg); take!(io))
_zmq_decode(::Type{T}, raw) where {T} = ProtoBuf.decode(ProtoBuf.ProtoDecoder(IOBuffer(raw)), T)

# Keyword-style constructor for DataForFUSE that survives future field additions:
# uses ProtoBuf.default_values (generated alongside the struct) for any field the
# caller doesn't override.
function _make_DataForFUSE(; kwargs...)
    defaults = ProtoBuf.default_values(DataForFUSE)
    overrides = NamedTuple(kwargs)
    vals = (haskey(overrides, f) ? overrides[f] : defaults[f] for f in keys(defaults))
    return DataForFUSE(vals...)
end

@testset "ActorZMQ round-trip" begin
    # ipc:// + tempname() avoids TCP port collisions on shared CI; Linux-only.
    endpoint = "ipc://" * tempname()

    # --- Round-trip happy path: metadata-only DataForFUSE (no psizr) ---
    # This exercises the wire format, the FUSERequest/DataForFUSE/Ack handshake,
    # the version handshake (FUSE.SCHEMA_VERSION echoed on both sides), the try/catch
    # wrappers, and the Bt-storage branch — without dragging in a fully-initialized
    # equilibrium (which would require a full FUSE.init).
    ini, act = FUSE.case_parameters(:ITER; init_from=:scalars)
    dd = IMAS.dd()
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
            resp = _make_DataForFUSE(;
                sim_time=0.1,
                Bt=2.1,
                has_Bt=true,
                schema_version=FUSE.SCHEMA_VERSION,
            )
            ZMQ.send(sock, _zmq_encode(resp))

            # Then a DataFromFUSE arrives; reply with an Ack with ok=true.
            _ = _zmq_decode(DataFromFUSE, ZMQ.recv(sock))
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
            resp = _make_DataForFUSE(; sim_time=0.0, schema_version=Int32(0))
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
    act.ActorZMQ.enabled = true
    act.ActorZMQ.endpoint = endpoint
    act.ActorZMQ.timeout_ms = 5000

    server = @async begin
        ctx = ZMQ.Context()
        sock = ZMQ.Socket(ctx, ZMQ.REP)
        ZMQ.bind(sock, endpoint)
        try
            _ = _zmq_decode(FUSERequest, ZMQ.recv(sock))
            ZMQ.send(sock, _zmq_encode(_make_DataForFUSE(;
                sim_time=0.1,
                Bt=2.1,
                has_Bt=true,
                schema_version=FUSE.SCHEMA_VERSION,
            )))
            _ = _zmq_decode(DataFromFUSE, ZMQ.recv(sock))
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
            ZMQ.send(sock, _zmq_encode(_make_DataForFUSE(;
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

# TODO: additional cases worth covering once the round-trip is stable:
#   - done=true short-circuit: receive! disconnects without populating dd
#   - psizr length vs nR*nZ mismatch: error message mentions "ActorZMQ" and "psizr"
#     (will need a fully-initialized dd to reach the equilibrium branch)
#   - had_psizr semantics: two successive receive!s, only the second carries psizr;
#     full flux_surfaces path runs on the second call, not the first
