using FUSE
using Test
# Pkg.test's environment holds only FUSE and test/Project.toml: reach FUSE's own dependencies through FUSE
import FUSE: ZMQ, YAML
using FUSE.Serialization: deserialize
using FUSE.Logging: with_logger

# A stand-in for the digital twin's registry: it answers OK to `n` requests and keeps them.
function _stand_in_registry(n)
    ctx = ZMQ.Context()
    sock = ZMQ.Socket(ctx, ZMQ.REP)
    ZMQ.bind(sock, "tcp://127.0.0.1:*")
    port = split(sock.last_endpoint[1:end-1], ":")[end]
    asked = String[]
    served = @async begin
        for _ in 1:n
            push!(asked, String(ZMQ.recv(sock)))
            ZMQ.send(sock, "OK")
        end
        ZMQ.close(sock)
        ZMQ.close(ctx)
    end
    return "127.0.0.1:$port", asked, served
end

@testset "ActorZMQPublish" begin
    registry, asked, served = _stand_in_registry(2)
    act = FUSE.ParametersActors()
    act.ActorZMQPublish.enabled = true
    act.ActorZMQPublish.registry_host = registry
    dd = IMAS.dd()
    dd.global_time = 1.0
    dd.equilibrium.time = [1.0]
    eqt = resize!(dd.equilibrium.time_slice, 1)[1]
    eqt.time = 1.0
    eqt.global_quantities.ip = 1.0e6  # listed in zmq_publish_fields.yaml
    eqt.global_quantities.psi_axis = 0.1  # not listed

    actor = FUSE.ActorZMQPublish(dd, act)
    @test startswith(asked[1], "REGISTER|FUSE_CENTER|tcp://")

    sub = ZMQ.Socket(ZMQ.SUB)
    sub.rcvtimeo = 5000  # a bug fails the test instead of hanging it
    ZMQ.subscribe(sub, "")
    ZMQ.connect(sub, "tcp://127.0.0.1:" * split(asked[1], ":")[end])
    sleep(0.5)  # a PUB drops what it sends before the subscription reaches it

    # the dd has one of the listed fields: the slice holds that one, the rest are left out,
    # and the run logs them once, naming them
    logger = Test.TestLogger()
    with_logger(() -> FUSE.publish!(actor), logger)
    @test occursin("build.layer[:].outline.r", only(logger.logs).message)

    # once per run: a second slice adds no record
    logger = Test.TestLogger()
    with_logger(() -> FUSE.publish!(actor), logger)
    @test isempty(logger.logs)

    FUSE.disconnect!(actor)

    slice = deserialize(IOBuffer(ZMQ.recv(sub)))
    @test slice.global_time == 1.0
    @test slice.equilibrium.time_slice[1].time == 1.0
    @test slice.equilibrium.time_slice[1].global_quantities.ip == 1.0e6
    @test ismissing(slice.equilibrium.time_slice[1].global_quantities, :psi_axis)
    @test isempty(slice.build.layer)  # listed, but the dd has no build
    ZMQ.recv(sub)  # the second slice
    @test isempty(ZMQ.recv(sub))  # the end of the run
    wait(served)
    @test asked[2] == "DEREGISTER|FUSE_CENTER"
    ZMQ.close(sub)
end

@testset "ActorZMQPublish fields are in the dd" begin
    @testset "$path" for path in YAML.load_file(FUSE.PUBLISH_FIELDS_FILE)["sim"]
        @test IMAS.info(path) isa IMAS.IMASdd.Info  # throws for a path the dd doesn't have
    end
end

@testset "ActorZMQPublish is off by default" begin
    actor = FUSE.ActorZMQPublish(IMAS.dd(), FUSE.ParametersActors())
    @test actor.socket === nothing
    FUSE.publish!(actor)
    FUSE.disconnect!(actor)
end
