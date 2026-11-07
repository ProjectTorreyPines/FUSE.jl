#= ================ =#
#  ActorZMQPublish   #
#= ================ =#
# Publishes the dd at the end of each dynamic step, for viewers outside FUSE: the digital
# twin's FUSE center takes the slices on to FUSEGUI and Omniverse, and archives them.
#
# It binds a ZMQ PUB socket and registers it with the digital twin's registry as
# FUSE_CENTER, the way every DT runner registers. Each message is the dd at its global_time
# cut down to the fields in zmq_publish_fields.yaml, the ones the viewers read, and
# serialized; an empty message ends the run. A PUB never holds FUSE up: with nobody
# subscribed, it drops the slice.
#
# Julia Serialization only round-trips between the same Julia and IMAS/IMASdd versions,
# so a reader must run the ones FUSE runs.

import Serialization
import YAML

const PUBLISH_NAME = "FUSE_CENTER"
const PUBLISH_FIELDS_FILE = joinpath(@__DIR__, "zmq_publish_fields.yaml")

@actor_parameters_struct ActorZMQPublish{T} begin
    enabled::Entry{Bool} = Entry{Bool}("-", "Publish the dd at the end of each dynamic step over ZMQ, registered as $PUBLISH_NAME"; default=false)
    registry_host::Entry{String} = Entry{String}("-", "Digital twin registry the publisher registers with (host:port)"; default="")
end

mutable struct ActorZMQPublish{D,P} <: SingleAbstractActor{D,P}
    dd::IMAS.DD{D}
    par::OverrideParameters{P,FUSEparameters__ActorZMQPublish{P}}
    context::Union{Nothing,ZMQ.Context}
    socket::Union{Nothing,ZMQ.Socket}
    fields::Vector{Vector{String}}
end

function ActorZMQPublish(dd::IMAS.DD{D}, par::FUSEparameters__ActorZMQPublish{P}; kw...) where {D<:Real,P<:Real}
    logging_actor_init(ActorZMQPublish)
    par = OverrideParameters(par; kw...)
    return ActorZMQPublish{D,P}(dd, par, nothing, nothing, Vector{String}[])
end

"""
    ActorZMQPublish(dd::IMAS.DD, act::ParametersAllActors; kw...)

Publishes the dd at the end of each dynamic step over ZeroMQ, for viewers outside FUSE.

When enabled, it binds a PUB socket and registers it with the digital twin's registry as
`FUSE_CENTER`. In ActorDynamicPlasma:

  - Phase 2 end, after `send!(actor_zmq)`: `publish!(actor_publish)` sends the dd at its
    global_time, cut down to the fields in `zmq_publish_fields.yaml`, serialized
  - When the run ends, also by an error: `disconnect!(actor_publish)` sends an empty
    message and deregisters

A PUB never holds FUSE up: with nobody subscribed, it drops the slice. A slice that
fails to publish is logged and skipped, and the run goes on.
"""
function ActorZMQPublish(dd::IMAS.DD, act::ParametersAllActors; kw...)
    actor = ActorZMQPublish(dd, act.ActorZMQPublish; kw...)
    if actor.par.enabled
        connect!(actor)
    end
    return actor
end

"""
    connect!(actor::ActorZMQPublish)

Load the fields to publish, bind the PUB socket and register it with the registry as `FUSE_CENTER`.
"""
function connect!(actor::ActorZMQPublish)
    if actor.socket !== nothing
        return actor
    end
    actor.fields = [IMAS.i2p(path) for path in YAML.load_file(PUBLISH_FIELDS_FILE)["sim"]]
    ctx = ZMQ.Context()
    sock = ZMQ.Socket(ctx, ZMQ.PUB)
    sock.linger = 1000  # close waits up to 1 s, so the end-of-run message still goes out
    endpoint = "tcp://$(gethostname()):$(_bind_random_port(sock))"
    _registry_request(actor.par.registry_host, "REGISTER|$PUBLISH_NAME|$endpoint")
    @info "ActorZMQPublish: publishing on $endpoint, registered as $PUBLISH_NAME"
    actor.context = ctx
    actor.socket = sock
    return actor
end

"""
    publish!(actor::ActorZMQPublish)

Send the dd at its global_time, cut down to the fields in `zmq_publish_fields.yaml`.
A listed field the dd has no data for is left out: a DIII-D run has no LH or IC antennas.
"""
function publish!(actor::ActorZMQPublish)
    if actor.socket === nothing
        return actor
    end
    try
        slice = IMAS.dd()
        for path in actor.fields
            IMAS.selective_copy!(actor.dd, slice, path, actor.dd.global_time)
        end
        io = IOBuffer()
        Serialization.serialize(io, slice)
        ZMQ.send(actor.socket, take!(io))
    catch e
        # the viewers miss a step; the run itself is not worth stopping for them
        @warn "ActorZMQPublish: could not publish the slice at t=$(actor.dd.global_time) s" exception = e maxlog = 5
    end
    return actor
end

"""
    disconnect!(actor::ActorZMQPublish)

Send the empty end-of-run message, deregister, and close the socket.
"""
function disconnect!(actor::ActorZMQPublish)
    if actor.socket === nothing
        return actor
    end
    ZMQ.send(actor.socket, UInt8[])
    try
        _registry_request(actor.par.registry_host, "DEREGISTER|$PUBLISH_NAME")
    catch e
        @warn "ActorZMQPublish: could not deregister $PUBLISH_NAME" exception = e
    end
    ZMQ.close(actor.socket)
    ZMQ.close(actor.context)
    actor.socket = nothing
    actor.context = nothing
    return actor
end

# The digital twin's runners bind in this port range.
function _bind_random_port(sock::ZMQ.Socket)
    for _ in 1:100
        port = rand(5555:6000)
        try
            ZMQ.bind(sock, "tcp://*:$port")
            return port
        catch e
            e isa ZMQ.StateError || rethrow()
        end
    end
    return error("ActorZMQPublish: found no free port in 5555:6000")
end

# One request to the registry, which answers OK: REGISTER|name|endpoint or DEREGISTER|name.
function _registry_request(registry_host::AbstractString, message::String)
    ctx = ZMQ.Context()
    sock = ZMQ.Socket(ctx, ZMQ.REQ)
    sock.rcvtimeo = 5000
    sock.linger = 0
    try
        ZMQ.connect(sock, "tcp://$registry_host")
        ZMQ.send(sock, message)
        reply = String(ZMQ.recv(sock))
        reply == "OK" || error("ActorZMQPublish: the registry at $registry_host answered \"$reply\" to $message")
    finally
        ZMQ.close(sock)
        ZMQ.close(ctx)
    end
    return nothing
end
