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

@testset "ActorZMQ" begin
    # ipc:// + tempname() avoids TCP port collisions on shared CI; Linux-only.
    endpoint = "ipc://" * tempname()

    # --- Round-trip happy path: metadata-only DataForFUSE (no psizr) ---
    # This exercises the wire format, the FUSERequest/DataForFUSE/Ack handshake,
    # the try/catch wrappers, and the Bt-storage branch — without dragging in a
    # fully-initialized equilibrium (which would require a full FUSE.init).
    ini, act = FUSE.case_parameters(:ITER; init_from=:scalars)
    dd = IMAS.dd()
    act.ActorZMQ.enabled = true
    act.ActorZMQ.endpoint = endpoint
    act.ActorZMQ.timeout_ms = 5000  # short timeout so a bug fails the test instead of hanging 30s

    server = @async begin
        ctx = ZMQ.Context()
        sock = ZMQ.Socket(ctx, ZMQ.REP)
        ZMQ.bind(sock, endpoint)
        try
            req = _zmq_decode(FUSERequest, ZMQ.recv(sock))
            @test req.status == "ready"
            resp = DataForFUSE(
                0.1,                  # sim_time
                false,                # done
                0.0,                  # Ip_latest (skip pulse_schedule write)
                0.0,                  # Ip_avg
                2.1,                  # Bt — the field we'll assert on
                0.0,                  # pr15v
                Float64[],            # I_coil
                Float64[],            # psizr (empty → skip equilibrium block)
                Float64[],            # r_grid
                Float64[],            # z_grid
                Float64[],            # pinj_per_beam
                Float64[],            # nbi_acc_voltage
                Float64[]             # gas_cal
            )
            ZMQ.send(sock, _zmq_encode(resp))

            # Then a DataFromFUSE arrives; reply with an Ack.
            _ = _zmq_decode(DataFromFUSE, ZMQ.recv(sock))
            ZMQ.send(sock, _zmq_encode(Ack("ack")))
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

    # TODO: additional cases worth covering once the round-trip is stable:
    #   - done=true short-circuit: receive! disconnects without populating dd
    #   - psizr length vs nR*nZ mismatch: error message mentions "ActorZMQ" and "psizr"
    #     (exercises the @assert→error fix; will need a fully-initialized dd to reach
    #      the equilibrium branch)
    #   - had_psizr semantics: two successive receive!s, only the second carries psizr;
    #     full flux_surfaces path runs on the second call, not the first
    #     (exercises the first_receive→had_psizr rename)
    #   - Ip_avg = 0.0 storage: blocked until proto3-default fix (issue #1) lands —
    #     a test today would freeze the existing wrong behavior.
end
