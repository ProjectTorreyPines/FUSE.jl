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

function atan_eq(r, z, r0, z0)
    if r[1] == r[end] && z[1] == z[end]
        r = r[1:end-1]
        z = z[1:end-1]
    end
    θ = unwrap(atan.(z .- z0, r .- r0))
    if θ[2] < θ[1]
        r = reverse(r)
        z = reverse(z)
        θ = reverse(θ)
    end
    return r, z, θ
end

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
    n = maximum(map(length, args))
    args = collect(map(x -> isa(x, Vector) ? x : [x], args))
    args = map(x -> vcat([x for k = 1:n]...)[1:n], args)
end

# ******************************************
# Convex Hull
# ******************************************
struct Point
    x::Float64
    y::Float64
end

function Base.isless(p::Point, q::Point)
    p.x < q.x || (p.x == q.x && p.y < q.y)
end

function isrightturn(p::Point, q::Point, r::Point)
    (q.x - p.x) * (r.y - p.y) - (q.y - p.y) * (r.x - p.x) < 0
end

function grahamscan(points::Vector{Point})
    sorted = sort(points)
    upperhull = halfhull(sorted)
    lowerhull = halfhull(reverse(sorted))
    [upperhull..., lowerhull[2:end-1]...]
end

function convex_hull(xy::Vector)
    return [(k.x, k.y) for k in grahamscan([Point(xx, yx) for (xx, yx) in xy])]
end

function halfhull(points::Vector{Point})
    halfhull = points[1:2]
    for p in points[3:end]
        push!(halfhull, p)
        while length(halfhull) > 2 && !isrightturn(halfhull[end-2:end]...)
            deleteat!(halfhull, length(halfhull) - 1)
        end
    end
    halfhull
end

function IMAS.force_float(x::ForwardDiff.Dual)
    ## we purposly do not do it recursively since generally
    ## ForwardDiff.Dual of ForwardDiff.Dual is an indication of someghing going wrong
    # return force_float(x.value)
    return x.value
end