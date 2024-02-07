import Jedis
import JSON

mutable struct Pubsub
    enabled::Bool
    identifier::String
    client::Jedis.Client
end

function pubsub(dd::IMAS.dd)
    aux = getfield(dd, :_aux)
    if :fuse_pubsub ∉ keys(aux)
        return false
    else
        return aux[:fuse_pubsub].enabled
    end
end

function pubsub(dd::IMAS.dd, enable::Bool; host::String="localhost", port::Int=55000, password::String="redispw")
    aux = getfield(dd, :_aux)
    if :fuse_pubsub ∈ keys(aux)
        aux[:fuse_pubsub].enabled = enable
    elseif enable
        client = Jedis.Client(; host, port, password, keepalive_enable=true)
        identifier = string(objectid(dd))
        aux[:fuse_pubsub] = Pubsub(enable, identifier, client)
    end
    return enable
end

function pubsub_client(dd)
    aux = getfield(dd, :_aux)
    if :fuse_pubsub ∈ keys(aux) # temporary fix
        identifier = aux[:fuse_pubsub].identifier
        delete!(aux, :fuse_pubsub)
        pubsub(dd, true)
        aux[:fuse_pubsub].identifier = identifier
    end
    if :fuse_pubsub ∈ keys(aux)
        return aux[:fuse_pubsub].client
    else
        return nothing
    end
end

function Jedis.publish(dd::IMAS.dd, channel_name::String, message::String)
    client = pubsub_client(dd)
    dd_channel = channel(dd, channel_name)
    return Jedis.publish(dd_channel, message; client)
end

function Jedis.publish(dd::IMAS.dd, channel_name::String; kw...)
    message = JSON.sprint(kw)
    return Jedis.publish(dd, channel_name, message)
end

function identifier(dd::IMAS.dd)
    aux = getfield(dd, :_aux)
    if :fuse_pubsub ∈ keys(aux)
        return aux[:fuse_pubsub].identifier
    else
        return string(objectid(dd))
    end
end

function channel(dd::IMAS.dd, channel_name::String)
    return "$(identifier(dd))__$(channel_name)"
end

function has_one_subscriber(dd::IMAS.dd, channel_name::String)
    subs = subscribers(dd, channel_name)
    @assert length(subs) <= 1 "Too many subscribers: $(subs)"
    return length(subs) == 1
end

function subscribers(dd::IMAS.dd, channel_name::String)
    aux = getfield(dd, :_aux)
    if :fuse_pubsub ∈ keys(aux)
        client = pubsub_client(dd)
        dd_channel_name = channel(dd, channel_name)
        return Jedis.execute("PUBSUB CHANNELS $(dd_channel_name)*", client)
    else
        return []
    end
end

function listen_and_wait(dd::IMAS.dd, channel_name::String; timeout::Float64)
    dd_channel = channel(dd, channel_name)
    client = pubsub_client(dd)

    # julia channel for blocking
    jchannel = Channel{String}(2)

    # async subscribe
    @async Jedis.subscribe(dd_channel; client) do msg
        return put!(jchannel, msg[end])
    end
    Jedis.wait_until_subscribed(client)

    # async timeout
    @async begin
        sleep(timeout)
        put!(jchannel, "")
    end

    # Wait for either the message to come or the timeout
    result = take!(jchannel)

    # unsubscribe
    Jedis.unsubscribe(dd_channel; client)

    if isempty(result)
        error("$(dd_channel) timeout")
    end

    return JSON.parse(result)
end