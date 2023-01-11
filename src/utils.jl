using InteractiveUtils: subtypes
import ForwardDiff

function unwrap(v, inplace=false)
    unwrapped = inplace ? v : copy(v)
    for i = 2:length(v)
        while (unwrapped[i] - unwrapped[i-1] >= pi)
            unwrapped[i] -= 2pi
        end
        while (unwrapped[i] - unwrapped[i-1] <= -pi)
            unwrapped[i] += 2pi
        end
    end
    return unwrapped
end

function IMAS.force_float(x::ForwardDiff.Dual)
    ## we purposly do not do it recursively since generally
    ## ForwardDiff.Dual of ForwardDiff.Dual is an indication of someghing going wrong
    # return force_float(x.value)
    return x.value
end

"""
    same_length_vectors(args...)

Returns scalars and vectors as vectors of the same lengths
For example:

    same_length_vectors(1, [2], [3,3,6], [4,4,4,4,4,4])
    
    4-element Vector{Vector{Int64}}:
    [1, 1, 1, 1, 1, 1]
    [2, 2, 2, 2, 2, 2]
    [3, 3, 6, 3, 3, 6]
    [4, 4, 4, 4, 4, 4]
"""
function same_length_vectors(args...)
    function length2(x)
        if ismissing(x)
            return 1
        else
            return length(x)
        end
    end
    n = maximum(map(length2, args))
    args = collect(map(x -> isa(x, Vector) ? x : [x], args))
    args = map(x -> vcat([x for k = 1:n]...)[1:n], args)
end

"""
    mirror_bound(x, l, u)

Return tuple with value of x bounded between l and u
The bounding is done by mirroring the value at the bound limits.
"""
function mirror_bound(x, l, u)
    d = (u - l) / 2.0
    c = (u + l) / 2.0
    x0 = (x .- c) / d
    while abs(x0) > 1.0
        if x0 < 1.0
            x0 = -2.0 - x0
        else
            x0 = 2.0 - x0
        end
    end
    return x0 * d + c
end

"""
    mirror_bound(x, l, u)

Return tuple with value of x bounded between l and u and error (cost) for going out of bounds
The bounding is done by mirroring the value at the bound limits.
"""
function mirror_bound_w_cost(x, l, u)
    y = mirror_bound.(x, l, u)
    return y, abs.((x .- l) .* (x .< l) .+ (x .- u) .* (x .> u))
end

# ******************************************
# Convex Hull
# ******************************************
import VacuumFields: convex_hull

# ******************************************
# TraceCAD
# ******************************************
import Images

"""
    dd_build_layers_to_ini(dd::IMAS.dd)

Utility function to convert layers in dd.build to layers in `ini.build.layers = layers = DataStructures.OrderedDict()`
"""
function dd_build_layers_to_ini(dd::IMAS.dd)
    for layer in dd.build.layer
        name = replace(layer.name, " " => "_")
        println("layers[:$name] = $(layer.thickness)")
    end
end

struct TraceCAD
    name::Symbol
    x_length::Real
    x_offset::Real
    y_offset::Real
end

function TraceCAD(device::Symbol)
    return TraceCAD(Val{device})
end

@recipe function plot_trace_cad(cad::TraceCAD)
    img = Images.load(joinpath(dirname(@__DIR__), "cases", "$(cad.name).jpg"))
    img .= img[end:-1:1, 1:end]
    x = LinRange(0, cad.x_length, size(img)[1]) .+ cad.x_offset
    y = LinRange(-0.5, 0.5, size(img)[2]) .* (size(img)[1] / size(img)[2]) * (x[end] - x[1]) .- cad.y_offset
    @series begin
        flip --> false
        x, y, img
    end
end

# ******************************************
# types
# ******************************************
"""
    returns enum from symbol
"""
function to_enum(smbl::Symbol)::Enum
    smbl = Symbol("_$(smbl)_")
    return @eval($smbl)
end

function to_enum(smbl::T where {T<:Enum})
    return smbl
end

"""
    concretetypes(type::Type)

List concrete subtypes of a given datatype
"""
function concretetypes(type::Type)
    concretetypes!(Any[], type)
end

function concretetypes!(out, type::Type)
    if !isabstracttype(type)
        push!(out, type)
    else
        foreach(T -> concretetypes!(out, T), subtypes(type))
    end
    out
end

# ******************************************
# FUSE
# ******************************************
function fuse()
    return """
███████╗██╗   ██╗███████╗███████╗
██╔════╝██║   ██║██╔════╝██╔════╝
█████╗  ██║   ██║███████╗█████╗  
██╔══╝  ██║   ██║╚════██║██╔══╝  
██║     ╚██████╔╝███████║███████╗
╚═╝      ╚═════╝ ╚══════╝╚══════╝
"""
end

"""
    warmup()

Function used to precompile the majority of FUSE
"""
function warmup()
    dd = IMAS.dd()
    return warmup(dd)
end

function warmup(dd::IMAS.dd)
    ini, act = case_parameters(:FPP; version=:v1_demount, init_from=:scalars)
    init(dd, ini, act)
    ActorWholeFacility(dd, act)
    IMAS.freeze(dd)
end