import SpecialFunctions
import QuadGK
import Interpolations
import PolygonOps
import LibGEOS
import GeoInterface
import FuseUtils: mirror_bound

#= =============== =#
#  Shape functions  #
#= =============== =#
function layer_shape_message(shape_function_index)
    buf = IOBuffer()
    show(buf, MIME("text/plain"), IMAS.BuildLayerShape)
    valid_shapes = String(take!(buf))
    return """
        layer.shape=$(shape_function_index) is invalid. Valid options are:

        $(valid_shapes)

        shape + z_offset = +100"""
end

function initialize_shape_parameters(shape_function_index, r_obstruction, z_obstruction, r_start, r_end, clearance)
    shape_parameters = nothing
    if shape_function_index ∈ (Int(_offset_), Int(_negative_offset_), Int(_convex_hull_))
        return nothing
    else
        height = maximum(z_obstruction) - minimum(z_obstruction) + clearance * 2.0
        z_offset = (maximum(z_obstruction) + minimum(z_obstruction)) / 2.0
        r_center = (r_obstruction[argmax(z_obstruction)] + r_obstruction[argmin(z_obstruction)]) / 2.0

        shape_index_mod = shape_function_index
        is_z_offset = false
        if shape_index_mod > 100
            shape_index_mod = mod(shape_function_index, 100)
            is_z_offset = true
        end
        if shape_index_mod in (Int(_princeton_D_), Int(_mirror_princeton_D_))
            shape_parameters = Float64[]
        elseif shape_index_mod in (Int(_princeton_D_scaled_), Int(_mirror_princeton_D_scaled_))
            shape_parameters = [height]
        elseif shape_index_mod in (Int(_double_ellipse_), Int(_mirror_double_ellipse_))
            centerpost_height = (maximum(z_obstruction) - minimum(z_obstruction)) * 2.0 / 3.0
            shape_parameters = [r_center, centerpost_height, height]
        elseif shape_index_mod in (Int(_circle_ellipse_), Int(_mirror_circle_ellipse_))
            centerpost_height = (maximum(z_obstruction) - minimum(z_obstruction)) * 2.0 / 3.0
            shape_parameters = [centerpost_height, height]
        elseif shape_index_mod in (Int(_rectangle_ellipse_), Int(_mirror_rectangle_ellipse_))
            shape_parameters = [r_center, height]
        elseif shape_index_mod in (Int(_triple_arc_), Int(_mirror_triple_arc_))
            shape_parameters = [height, 0.5, 0.5, Float64(pi / 3), Float64(pi / 3)]
        elseif shape_index_mod == Int(_rectangle_)
            shape_parameters = [height]
        elseif shape_index_mod == Int(_racetrack_)
            shape_parameters = [height, 0.25]
        elseif shape_index_mod == Int(_silo_)
            shape_parameters = [height, height / 2.0]
        end
    end
    if shape_parameters === nothing
        error(layer_shape_message(shape_function_index))
    end
    if is_z_offset
        push!(shape_parameters, z_offset)
    end
    return shape_parameters
end

function shape_function(shape_function_index::Int; resolution::Float64)
    func = nothing
    if shape_function_index ∈ (Int(_offset_), Int(_negative_offset_), Int(_convex_hull_))
        return nothing
    else
        shape_index_mod = shape_function_index
        is_z_offset = false
        if shape_index_mod > 100
            shape_index_mod = mod(shape_function_index, 100)
            is_z_offset = true
        end
        if shape_index_mod in (Int(_princeton_D_), Int(_mirror_princeton_D_))
            func = princeton_D_approx
        elseif shape_index_mod in (Int(_princeton_D_scaled_), Int(_mirror_princeton_D_scaled_))
            func = princeton_D_scaled
        elseif shape_index_mod in (Int(_double_ellipse_), Int(_mirror_double_ellipse_))
            func = double_ellipse
        elseif shape_index_mod in (Int(_circle_ellipse_), Int(_mirror_circle_ellipse_))
            func = circle_ellipse
        elseif shape_index_mod in (Int(_rectangle_ellipse_), Int(_mirror_rectangle_ellipse_))
            func = rectangle_ellipse
        elseif shape_index_mod in (Int(_triple_arc_), Int(_mirror_triple_arc_))
            func = triple_arc
        elseif shape_index_mod == Int(_rectangle_)
            func = rectangle_shape
        elseif shape_index_mod == Int(_racetrack_)
            func = racetrack
        elseif shape_index_mod == Int(_silo_)
            func = silo
        end
    end
    if func === nothing
        error(layer_shape_message(shape_function_index))
    end

    # resolution
    rfunc(args...) = func(args...; resolution)

    # zoffset
    zfunc = rfunc
    if is_z_offset
        zfunc(args...) = begin
            R, Z = rfunc(args[1:end-1]...)
            Z .+= args[end]
            return R, Z
        end
    end

    # neg-D
    dfunc = zfunc
    if shape_index_mod in (
        Int(_mirror_princeton_D_),
        Int(_mirror_princeton_D_scaled_),
        Int(_mirror_double_ellipse_),
        Int(_mirror_circle_ellipse_),
        Int(_mirror_rectangle_ellipse_),
        Int(_mirror_triple_arc_)
    )
        dfunc(args...) = begin
            R, Z = zfunc(args...)
            R = -(R .- args[1]) .+ args[2]
            return reverse(R), reverse(Z)
        end
    end

    return dfunc
end

"""
    struveL(ν, z)

Modified Struve function.
Ported from https://github.com/gwater/Struve.jl
(Struve had some requirements conflicts with FUSE)
"""
function struveL(ν, z)
    integrate(f, a, b) = QuadGK.quadgk(f, a, b)[1]
    _M_integral(ν, z) = -2(0.5z)^ν / (sqrt(pi) * SpecialFunctions.gamma(ν + 0.5)) * integrate(t -> exp(-z * t) * (1 - t^2)^(ν - 0.5), 0, 1)
    return SpecialFunctions.besseli(ν, z) + _M_integral(ν, z)
end

"""
    ellipse(a::T, b::T, t0::T, t1::T, x0::T, z0::T; n_points::Integer=100) where {T<:Real}

Simple ellipse shape function
"""
function ellipse(a::T, b::T, t0::T, t1::T, x0::T, z0::T; n_points::Integer=100) where {T<:Real}
    t = range(t0, t1, n_points)
    x = a .* cos.(t) .+ x0
    z = b .* sin.(t) .+ z0
    return x, z
end

"""
    princeton_D_approx(r_start::T, r_end::T; n_points::Integer=100, resolution::Float64=1.0) where {T<:Real}

An pproximate version of the "Princeton-D" constant tension shape for a TF coil that is built with ellipses
References: Gralnick, S. L.; Tenney, F. H. Analytic Solutions for Constant‐tension Coil Shapes. J. Appl. Phys. 1976, 47, 7
"""
function princeton_D_approx(r_start::T, r_end::T; n_points::Integer=100, resolution::Float64=1.0) where {T<:Real}
    r1 = r_start
    r2 = r_end
    k = 0.5 * log(r2 / r1)
    r0 = sqrt(r1 * r2)

    # analytic equations for coil parameters
    centerpost_maxz = 2 * pi * r0 * k * SpecialFunctions.besseli(1, k) / 2 # Gralnick Eq. 28
    coil_maxz = pi * r0 * k * (SpecialFunctions.besseli(1, k) + struveL(1, k) + 2 / pi) / 2  # Gralnick Eq. 34

    # number of points
    ra_e1 = r0 - r1
    zb_e1 = coil_maxz - centerpost_maxz
    ra_e2 = r2 - r0
    zb_e2 = coil_maxz
    e1_e2_ratio = sqrt(ra_e1^2 + zb_e1^2) / sqrt(ra_e2^2 + zb_e2^2)
    n1_points = max(3, Int(round(n_points * resolution * e1_e2_ratio)))
    n2_points = Int(round(n_points * resolution))

    # make ellipses to approximate equal-tension arc
    r1, z1 = ellipse(r0 - r1, coil_maxz - centerpost_maxz, float(π), float(π / 2), r0, centerpost_maxz; n_points=n1_points)
    r2, z2 = ellipse(r2 - r0, coil_maxz, float(π / 2), float(0.0), r0, 0.0; n_points=n2_points)

    # stack
    R = [r1[1:end-1]; r2[1:end-1]; r2[end:-1:2]; r1[end:-1:1]; r1[1]]
    Z = [z1[1:end-1]; z2[1:end-1]; -z2[end:-1:2]; -z1[end:-1:1]; z1[1]]

    return R, Z
end

"""
    princeton_D_scaled(r_start::T, r_end::T, height::T; n_points::Integer=100, resolution::Float64=1.0) where {T<:Real}

This routine calculates a "shortened" TF coil shape that foregoes the equal-tension "Princeton-Dee" for
a squarer, more space-efficient shape. It replicates the inboard curve of the equal-tension arc, but
decreases the height of the coil to match a given value.
"""
function princeton_D_scaled(r_start::T, r_end::T, height::T; n_points::Integer=100, resolution::Float64=1.0) where {T<:Real}
    r1 = r_start
    r2 = r_end
    k = 0.5 * log(r2 / r1)
    r0 = sqrt(r1 * r2)

    centerpost_maxz = 2 * pi * r0 * k * SpecialFunctions.besseli(1, k) / 2 # Gralnick Eq. 28
    coil_maxz = pi * r0 * k * (SpecialFunctions.besseli(1, k) + struveL(1, k) + 2 / pi) / 2  # Gralnick Eq. 34

    inboard_curve_dz = coil_maxz - centerpost_maxz
    centerpost_maxz = height / 2
    coil_maxz = centerpost_maxz + inboard_curve_dz

    # number of points
    ra_e1 = r0 - r1
    zb_e1 = coil_maxz - centerpost_maxz
    ra_e2 = r2 - r0
    zb_e2 = coil_maxz
    e1_e2_ratio = sqrt(ra_e1^2 + zb_e1^2) / sqrt(ra_e2^2 + zb_e2^2)
    n1_points = max(3, Int(round(n_points * resolution * e1_e2_ratio)))
    n2_points = Int(round(n_points * resolution))

    # make ellipses to approximate equal-tension arc
    r1, z1 = ellipse(r0 - r1, coil_maxz - centerpost_maxz, float(π), float(π / 2), r0, centerpost_maxz; n_points=n1_points)
    r2, z2 = ellipse(r2 - r0, coil_maxz, float(π / 2), float(0.0), r0, 0.0; n_points=n2_points)

    # stack
    R = [r1[1:end-1]; r2[1:end-1]; r2[end:-1:2]; r1[end:-1:1]; r1[1]]
    Z = [z1[1:end-1]; z2[1:end-1]; -z2[end:-1:2]; -z1[end:-1:1]; z1[1]]

    return R, Z
end

"""
    double_ellipse(r_start::T, r_end::T, r_center::T, centerpost_height::T, height::T; n_points::Integer=100, resolution::Float64=1.0) where {T<:Real}

double ellipse shape
"""
function double_ellipse(r_start::T, r_end::T, r_center::T, centerpost_height::T, height::T; n_points::Integer=100, resolution::Float64=1.0) where {T<:Real}
    height = abs(height)
    centerpost_height = mirror_bound(abs(centerpost_height), 0.5 * height, height)
    r_center = mirror_bound(r_center, r_start, r_end)
    return double_ellipse(r_start, r_end, r_center, centerpost_height, 0.0, height; n_points, resolution)
end

function double_ellipse(r_start::T, r_end::T, height::T; n_points::Integer=100, resolution::Float64=1.0) where {T<:Real}
    return double_ellipse(r_start, r_end, (r_start + r_end) / 2, height, 0.0, height; n_points, resolution)
end

function double_ellipse(r_start::T, r_end::T, r_center::T, centerpost_height::T, outerpost_height::T, height::T; n_points::Integer=100, resolution::Float64=1.0) where {T<:Real}
    # inner
    r_e1 = r_center
    z_e1 = centerpost_height / 2.0
    ra_e1 = r_center - r_start
    zb_e1 = (height - centerpost_height) / 2.0

    # outer
    r_e2 = r_center
    z_e2 = outerpost_height / 2.0
    ra_e2 = r_end - r_center
    zb_e2 = (height - outerpost_height) / 2.0

    # number of points
    e1_e2_ratio = sqrt(ra_e1^2 + zb_e1^2) / sqrt(ra_e2^2 + zb_e2^2)
    n1_points = max(3, Int(round(n_points * resolution * e1_e2_ratio)))
    n2_points = Int(round(n_points * resolution))

    # ellipses
    r1, z1 = ellipse(ra_e1, zb_e1, float(π), float(π / 2), r_e1, z_e1; n_points=n1_points)
    r2, z2 = ellipse(ra_e2, zb_e2, float(π / 2), float(0.0), r_e2, z_e2; n_points=n2_points)

    # stack
    R = [r1[1:end-1]; r2[1:end-1]; r2[end:-1:2]; r1[end:-1:1]; r1[1]]
    Z = [z1[1:end-1]; z2[1:end-1]; -z2[end:-1:2]; -z1[end:-1:1]; z1[1]]

    return R, Z
end

"""
    circle_ellipse(r_start::T, r_end::T, centerpost_height::T, height::T; n_points::Integer=100, resolution::Float64=1.0) where {T<:Real}

circle ellipse shape (parametrization of TF coils used in GATM)

Special case of the double ellipse shape, where the inner ellipse is actually a circle
"""
function circle_ellipse(r_start::T, r_end::T, centerpost_height::T, height::T; n_points::Integer=100, resolution::Float64=1.0) where {T<:Real}
    height = abs(height)
    centerpost_height = mirror_bound(abs(centerpost_height), 0.0, height)
    r_circle = (height - centerpost_height) / 2.0
    r_center = r_start + r_circle
    return double_ellipse(r_start, r_end, r_center, centerpost_height, 0.0, height; n_points, resolution)
end

"""
    rectangle_ellipse(r_start::T, r_end::T, r_center::T, height::T; n_points::Integer=100, resolution::Float64=1.0) where {T<:Real}

Rectangle ellipse shape

The inner part of the TF is actually a rectangle
"""
function rectangle_ellipse(r_start::T, r_end::T, r_center::T, height::T; n_points::Integer=100, resolution::Float64=1.0) where {T<:Real}
    height = abs(height)
    r_center = mirror_bound(r_center, r_start, r_end)
    n_points = Int(round(n_points * resolution))
    r2, z2 = ellipse(r_end - r_center, height / 2.0, float(π / 2), float(0.0), r_center, 0.0; n_points)
    R = [r_start; r_start; r2[1:end-1]; r2[end:-1:1]; r_start]
    Z = [-z2[1]; z2[1]; z2[1:end-1]; -z2[end:-1:1]; -z2[1]]
    return R, Z
end

"""
    rectangle_shape(r_start::T, r_end::T, z_low::T, z_high::T; n_points::Int=400, resolution::Float64=0.0) where {T<:Real}

Asymmetric rectangular shape

NOTE: by default resolution=0.0, which returns 5 points
"""
function rectangle_shape(r_start::T, r_end::T, z_low::T, z_high::T; n_points::Int=400, resolution::Float64=0.0) where {T<:Real}
    if resolution == 0.0
        R = [r_start, r_start, r_end, r_end, r_start]
        Z = [z_low, z_high, z_high, z_low, z_low]
    else
        n_points = Int(round(n_points / 4 * resolution))
        R = vcat(
            range(r_start, r_start, n_points),
            range(r_start, r_end, n_points)[2:end],
            range(r_end, r_end, n_points)[2:end],
            range(r_end, r_start, n_points)[2:end])
        Z = vcat(
            range(z_low, z_high, n_points),
            range(z_high, z_high, n_points)[2:end],
            range(z_high, z_low, n_points)[2:end],
            range(z_low, z_low, n_points)[2:end])
    end
    return R, Z
end

"""
    rectangle_shape(r_start::T, r_end::T, height::T; n_points::Integer=400, resolution::Float64=0.0) where {T<:Real}

Symmetric rectangular shape

NOTE: by default resolution=0, which returns 5 points
"""
function rectangle_shape(r_start::T, r_end::T, height::T; n_points::Integer=400, resolution::Float64=0.0) where {T<:Real}
    Δ = height / 2.0
    return rectangle_shape(r_start, r_end, -Δ, Δ; n_points, resolution)
end

"""
    racetrack(r_start::T, r_end::T, z_low::T, z_high::T, radius_ratio::T; n_points::Int=400, resolution::Float64=1.0) where {T<:Real}

Asymmetric racetrack shape
"""
function racetrack(r_start::T, r_end::T, z_low::T, z_high::T, radius_ratio::T; n_points::Int=400, resolution::Float64=1.0) where {T<:Real}
    R, Z = rectangle_shape(r_start, r_end, z_low, z_high; n_points, resolution)
    if radius_ratio != 0.0
        radius_ratio = mirror_bound(abs(radius_ratio), 0.0, 1.0 - 1E-3)
        radius = min((z_high - z_low), (r_end - r_start)) / 2.0 * radius_ratio
        R, Z = buffer(R, Z, -radius)
        R, Z = buffer(R, Z, +radius)
    end
    return R, Z
end

"""
    racetrack(r_start::T, r_end::T, height::T, radius_ratio::T; n_points::Integer=400, resolution::Float64=1.0) where {T<:Real}

Symmetric racetrack shape
"""
function racetrack(r_start::T, r_end::T, height::T, radius_ratio::T; n_points::Integer=400, resolution::Float64=1.0) where {T<:Real}
    Δ = height / 2.0
    return racetrack(r_start, r_end, -Δ, Δ, radius_ratio; n_points, resolution)
end

"""
    triple_arc(
        r_start::T,
        r_end::T,
        height::T,
        small_radius::T,
        mid_radius::T,
        small_coverage::T,
        mid_coverage::T;
        min_small_radius_fraction::T=0.2,
        min_mid_radius_fraction::T=min_small_radius_fraction * 2.0,
        n_points::Integer=100,
        resolution::Float64=1.0) where {T<:Real}

TrippleArc shape. Angles are in radians.
"""
function triple_arc(
    r_start::T,
    r_end::T,
    height::T,
    small_radius::T,
    mid_radius::T,
    small_coverage::T,
    mid_coverage::T;
    n_points::Integer=400,
    resolution::Float64=1.0) where {T<:Real}

    n_points = Int(round(n_points / 4.0 * resolution))

    height = abs(height) / 2.0
    small_radius = mirror_bound(small_radius, 0.01, 1.0) * height
    mid_radius = mirror_bound(mid_radius, small_radius / height, 1.0) * height
    small_coverage = mirror_bound(small_coverage, 0.0, Float64(pi))
    mid_coverage = mirror_bound(mid_coverage, 0.0, Float64(pi))

    asum = small_coverage + mid_coverage

    # small arc
    theta = range(0, small_coverage, n_points)
    small_arc_R = r_start .+ small_radius .* (1 .- cos.(theta))
    small_arc_Z = height .+ small_radius .* sin.(theta)

    # mid arc
    theta = range(small_coverage, asum, n_points)
    mid_arc_R = small_arc_R[end] .+ mid_radius .* (cos.(small_coverage) .- cos.(theta))
    mid_arc_Z = small_arc_Z[end] .+ mid_radius .* (sin.(theta) .- sin.(small_coverage))

    # large arc
    theta = range(theta[end], pi, n_points)
    large_radius = mid_arc_Z[end] / sin(pi - asum)
    large_arc_R = mid_arc_R[end] .+ large_radius .* (cos.(pi .- theta) .- cos.(pi .- asum))
    large_arc_Z = mid_arc_Z[end] .- large_radius .* (sin(asum) .- sin.(pi .- theta))

    R = vcat(small_arc_R, mid_arc_R[2:end], large_arc_R[2:end])
    R = vcat(R, reverse(R)[2:end])
    Z = vcat(small_arc_Z, mid_arc_Z[2:end], large_arc_Z[2:end])
    Z = vcat(Z, -reverse(Z)[2:end])

    # Add vertical
    R = vcat(range(r_start, r_start, n_points), R)
    Z = vcat(range(-height, height, n_points), Z)

    # Resize to ensure r_start to r_end
    factor = (r_end - r_start) / (maximum(R) - minimum(R))
    Z = Z .* factor
    R = (R .- minimum(R)) .* factor .+ r_start

    return R, Z
end

"""
    xy_polygon(x::T, y::T) where {T<:AbstractVector{<:Real}}

Returns LibGEOS.Polygon from x and y arrays
"""
function xy_polygon(x::T, y::T) where {T<:AbstractVector{<:Real}}
    x = convert(Vector{Float64}, x)
    y = convert(Vector{Float64}, y)
    if (x[1] ≈ x[end]) && (y[1] ≈ y[end])
        coords = [[[x[i], y[i]] for i in 1:length(x)]]
        coords[1][end] .= coords[1][1]
    else
        coords = [[i > length(x) ? [x[1], y[1]] : [x[i], y[i]] for i in 1:length(x)+1]]
    end

    if (coords[1][1][1] != coords[1][end][1]) && (coords[1][1][1] ≈ coords[1][end][1])
        coords[1][end][1] = coords[1][1][1]
    end
    if (coords[1][1][2] != coords[1][end][2]) && (coords[1][1][2] ≈ coords[1][end][2])
        coords[1][end][2] = coords[1][1][2]
    end
    return LibGEOS.Polygon(coords)
end

"""
    xy_polygon(coords::T) where {T<:Vector{Tuple{Float64, Float64}}}

Returns LibGEOS.Polygon from vector of tuple with x, y points
"""
function xy_polygon(coords::T) where {T<:Vector{Tuple{Float64,Float64}}}
    return LibGEOS.Polygon([[[coord[1], coord[2]] for coord in coords]])
end

"""
    xy_polygon(layer::Union{IMAS.build__layer,IMAS.build__structure})

Returns LibGEOS.Polygon from IMAS build layer or structure
"""
function xy_polygon(layer::Union{IMAS.build__layer,IMAS.build__structure})
    return xy_polygon(layer.outline.r, layer.outline.z)
end

"""
    buffer(x::AbstractVector{T}, y::AbstractVector{T}, b::T)::Tuple{Vector{T},Vector{T}} where {T<:Real}

Buffer polygon defined by x,y arrays by a quantity b
"""
function buffer(x::AbstractVector{T}, y::AbstractVector{T}, b::T)::Tuple{Vector{Float64},Vector{Float64}} where {T<:Real}
    poly = xy_polygon(x, y)
    poly_b::LibGEOS.Polygon = LibGEOS.buffer(poly, b)
    @inline return get_xy(poly_b, Float64)
end

function get_xy(poly::LibGEOS.Polygon, T=Float64)
    # extract the LinearRing struct
    lr = first(GeoInterface.getgeom(poly))

    # preallocate a LibGEOS.Point struct and arrays to store (x,y) coordinates
    Npts = GeoInterface.ngeom(lr)
    x = zeros(T, Npts)
    y = zeros(T, Npts)
    pt = LibGEOS.Point(zero(T), zero(T))
    for k in 1:Npts
        pt = getPoint!(pt, lr, k)
        x[k] = GeoInterface.x(pt)::T
        y[k] = GeoInterface.y(pt)::T
    end
    return x, y
end

# An allocation free version of LibGEOS.getPoint
function getPoint!(
    pt::LibGEOS.Point,
    obj::LibGEOS.LinearRing,
    n::Integer,
    context::LibGEOS.GEOSContext=LibGEOS.get_context(obj)
)
    if !(0 < n <= LibGEOS.numPoints(obj, context))
        error(
            "LibGEOS: n=$n is out of bounds for LineString with numPoints=$(numPoints(obj, context))"
        )
    end
    result = LibGEOS.GEOSGeomGetPointN_r(context, obj, n - 1)
    if result == LibGEOS.C_NULL
        error("LibGEOS: Error in GEOSGeomGetPointN")
    end
    pt.ptr = result
    pt.context = context
    return pt
end

"""
    buffer(x::AbstractVector{T}, y::AbstractVector{T}, b_hfs::T, b_lfs::T)::Tuple{Vector{T},Vector{T}} where {T<:Real}

Buffer polygon defined by x,y arrays by a quantity b_hfs to the left and b_lfs to the right
"""
function buffer(x::AbstractVector{T}, y::AbstractVector{T}, b_hfs::T, b_lfs::T)::Tuple{Vector{T},Vector{T}} where {T<:Real}
    x_b, y_b = buffer(x, y, 0.5 * (b_lfs + b_hfs))
    x_offset = 0.5 * (b_lfs - b_hfs)
    x_b .+= x_offset
    return x_b, y_b
end

"""
    limit_curvature(x::AbstractVector{T}, y::AbstractVector{T}, max_curvature::Real)::Tuple{Vector{Float64},Vector{Float64}} where {T<:Real}

Limit maximum curvature of a polygon described by x,y arrays
"""
function limit_curvature(x::AbstractVector{T}, y::AbstractVector{T}, max_curvature::Real)::Tuple{Vector{Float64},Vector{Float64}} where {T<:Real}
    @assert max_curvature > 0.0
    x = convert(Vector{Float64}, x)
    y = convert(Vector{Float64}, y)
    max_curvature = convert(Float64, max_curvature)
    poly = xy_polygon(x, y)
    poly_b = LibGEOS.buffer(LibGEOS.buffer(poly, -max_curvature)::LibGEOS.Polygon, max_curvature)::LibGEOS.Polygon
    @inline return get_xy(poly_b, Float64)
end

"""
    volume_no_structures(layer::IMAS.build__layer, structures::IMAS.IDSvector{<:IMAS.build__structure})

Returns volume of the layer without structures
"""
function volume_no_structures(layer::IMAS.build__layer, structures::IMAS.IDSvector{<:IMAS.build__structure})
    vol = 0.0
    for structure in structures
        vol += layer_structure_intersect_volume(layer, structure)
    end
    return layer.volume - vol
end

# this expressions is added here because volume_no_structures is not a IMAS function
IMAS.dynamic_expressions["build.layer[:].volume_no_structures"] =
    (; build, layer, _...) -> volume_no_structures(layer, build.structure)

"""
    layer_structure_intersect_volume(layer::IMAS.build__layer, structure::IMAS.build__structure)

Returns volume of the intersection between build layer volume and structure volume
"""
function layer_structure_intersect_volume(layer::IMAS.build__layer, structure::IMAS.build__structure)
    if layer.type ∈ (Int(_in_), Int(_out_))
        return layer.volume
    elseif layer.type ∈ (Int(_tf_), Int(_plasma_))
        return layer.volume
    elseif layer.side ∈ (Int(_hfs_), Int(_lfs_))
        i = IMAS.index(layer)
        if layer.side == Int(_hfs_)
            layer_in = IMAS.parent(layer)[i+1]
        else
            layer_in = IMAS.parent(layer)[i-1]
        end
    end
    layer_poly = xy_polygon(layer)
    layer_in_poly = xy_polygon(layer_in)
    ring_poly = LibGEOS.difference(layer_poly, layer_in_poly)
    structure_poly = xy_polygon(structure)

    if structure.toroidal_extent == 2 * pi
        toroidal_angles = [0.0]
    else
        toroidal_angles = structure.toroidal_angles
    end

    vol = 0.0
    for poly in GeoInterface.coordinates(LibGEOS.intersection(ring_poly, structure_poly))
        pr = [v[1] for v in poly]
        pz = [v[2] for v in poly]
        vol += IMAS.area(pr, pz) * structure.toroidal_extent * length(toroidal_angles)
    end

    return vol
end

"""
    P_LH_threshold_from_scalars(Bt0::Real, nel::Real,surface_area::Real)

Calculates the L to H transition threshold according to 2008 scaling law
"""
function P_LH_threshold_from_scalars(Bt0::Real, nel::Real, surface_area::Real)
    return 0.049 * abs(Bt0)^0.8 * nel^0.72 * surface_area^0.94
end
"""
    approximate_surface_area(a::Real, R::Real ,κ::Real, δ::Real)

Approximation of the surface area of a miller geometry flux surface
"""
function approximate_surface_area(a::Real, R::Real, κ::Real, δ::Real)
    return 2pi^2 * R * a^2 * κ * (1.0 - 0.151 * δ * a / R) / (2pi * R) # m²
end

"""
    silo(r_start, r_end, height_start, height_end; n_points::Int=100, resolution::Float64=1.0)
"""
function silo(r_start::Real, r_end::Real, height_start::Real, height_end::Real; n_points::Int=100, resolution::Float64=1.0)
    n_points = Int(round(n_points * resolution))
    height_start = abs(height_start)
    height_end = abs(height_end)
    height_end = min(max(height_end, height_start * 0.0), height_start * 0.9)
    x, y = ellipse(r_end - r_start, height_start - height_end, 0.0, pi / 2, r_start, height_end; n_points)
    return vcat(r_start, r_start, r_end, x), vcat(height_start, 0.0, 0.0, y) .- (height_start / 2.0)
end

"""
    line_through_point(m, x0, y0, x)

return `y` values at `x` of line with gradient `m` going through point `(x0,y0)`
"""
function line_through_point(m::T, x0::T, y0::T, x::Vector{T}) where {T<:Real}
    return @. m * x + y0 - m * x0
end

mutable struct MXHboundary
    mxh::IMAS.MXH
    RX::Vector{<:Real}
    ZX::Vector{<:Real}
    r_boundary::Vector{<:Real}
    z_boundary::Vector{<:Real}
end

function Base.show(io::IO, ::MIME"text/plain", mxhb::MXHboundary)
    mxh = mxhb.mxh
    print(io, mxh)
    println(io, "RX: $(mxhb.RX)")
    println(io, "ZX: $(mxhb.ZX)")
    return nothing
end

@recipe function plot_mxhb(mxhb::MXHboundary)
    @series begin
        aspect_ratio --> :equal
        label --> ""
        lw --> 2
        color --> :black
        mxhb.r_boundary, mxhb.z_boundary
    end
    @series begin
        primary := false
        aspect_ratio --> :equal
        label --> ""
        mxhb.mxh
    end
    @series begin
        seriestype := :scatter
        primary := false
        aspect_ratio --> :equal
        label --> ""
        marker := :cross
        mxhb.RX, mxhb.ZX
    end
end

"""
    add_xpoint(
        mr::AbstractVector{T},
        mz::AbstractVector{T},
        R0::T,
        Z0::T;
        upper::Bool,
        α_multiplier::Float64) where {T<:Real}

Add a X-point to a boundary that does not have one.

The X-point is added at the point where there is the largest curvature in the boundary,
and is built in such a way to form 90° angle (which is what it should be for zero current at the LCFS).
The X-point is placed on a line that goes from (R0,Z0) through the point of maximum curvature.

Control of the X-point location can be achieved by modifying R0, Z0.
"""
function add_xpoint(
    mr::AbstractVector{T},
    mz::AbstractVector{T},
    R0::T,
    Z0::T;
    upper::Bool,
    α_multiplier::Float64) where {T<:Real}

    function cost(pr::Vector{Float64}, pz::Vector{Float64}, i::Integer, R0::T, Z0::T, α::Float64)
        return abs(add_xpoint(pr, pz, i, R0, Z0, α).θX - π / 2)
    end

    if upper
        i = argmax(abs.(IMAS.curvature(mr, mz)) .* (mz .> Z0))
    else
        i = argmax(abs.(IMAS.curvature(mr, mz)) .* (mz .< Z0))
    end

    res = Optim.optimize(α -> cost(mr, mz, i, R0, Z0, α), 1.0, 1.5, Optim.GoldenSection())
    return add_xpoint(mr, mz, i, R0, Z0, 1.0 + (res.minimizer[1] - 1.0) * α_multiplier)
end

function add_xpoint(pr::Vector{Float64}, pz::Vector{Float64}, i::Integer, R0::T, Z0::T, α::T) where {T<:Real}
    @assert IMAS.is_closed_polygon(pr, pz)
    RX = pr[i] .* α .+ R0 .* (1.0 .- α)
    ZX = pz[i] .* α .+ Z0 .* (1.0 .- α)

    pr = @views pr[1:end-1]
    pz = @views pz[1:end-1]
    n = length(pr)

    min_angle = Inf
    i_min = 0
    max_angle = -Inf
    i_max = 0
    for i in 1:n
        @inbounds angle = atan(pz[i] - ZX, pr[i] - RX)
        if angle < min_angle
            min_angle = angle
            i_min = i
        end
        if angle > max_angle
            max_angle = angle
            i_max = i
        end
    end

    index = [mod(idx - 1, n) + 1 for idx in i_max:i_max+n-(i_max-i_min)]
    R = @views [RX; pr[index]; RX]
    Z = @views [ZX; pz[index]; ZX]

    θX = max_angle - min_angle

    return (RX=RX, ZX=ZX, R=R, Z=Z, θX=θX)
end

"""
    MXHboundary(mxh::IMAS.MXH; upper_x_point::Bool, lower_x_point::Bool, n_points::Integer=0)

Return MXHboundary structure of boundary with x-points based on input MXH boundary parametrization
"""
function MXHboundary(mxh::IMAS.MXH; upper_x_point::Bool, lower_x_point::Bool, n_points::Integer=0)
    mr, mz = mxh()
    mxhb = MXHboundary(deepcopy(mxh), Float64[], Float64[], mr, mz)
    return MXHboundary!(mxhb; upper_x_point, lower_x_point, n_points)
end

function MXHboundary!(mxhb::MXHboundary; upper_x_point::Bool, lower_x_point::Bool, n_points::Integer=0, RU::Real=0.0, RL::Real=0.0)
    mr, mz = mxhb.mxh()

    if ~upper_x_point && ~lower_x_point && n_points === 0
        empty!(mxhb.RX)
        empty!(mxhb.ZX)
        mxhb.r_boundary = mr
        mxhb.z_boundary = mz
        return mxhb
    end

    R0 = mxhb.mxh.R0
    Z0 = mxhb.mxh.Z0
    if RU == 0.0
        RU = R0
    end
    if RL == 0.0
        RL = R0
    end

    R, Z = mr, mz
    RXU, ZXU, R1, Z1, _ = add_xpoint(R, Z, RU, Z0; upper=true, α_multiplier=(upper_x_point ? 1.0 : 2.0))
    if upper_x_point
        R = R1
        Z = Z1
    end

    RXL, ZXL, R2, Z2 = add_xpoint(R, Z, RL, Z0; upper=false, α_multiplier=(lower_x_point ? 1.0 : 2.0))
    if lower_x_point
        R = R2
        Z = Z2
    end

    RX = Float64[]
    ZX = Float64[]
    if upper_x_point
        push!(RX, RXU)
        push!(ZX, ZXU)
    end
    if lower_x_point
        push!(RX, RXL)
        push!(ZX, ZXL)
    end

    # resample boundary after convex_hull in such a way to preserve x-points and add proper curvature in the x-point region
    if upper_x_point + lower_x_point == 1
        if upper_x_point
            R, Z = IMAS.reorder_flux_surface!(R, Z, argmax(Z))
        else
            R, Z = IMAS.reorder_flux_surface!(R, Z, argmin(Z))
        end
        RR = @views [-(reverse(R[2:end-1]) .- R[1]) .+ R[1]; R[2:end-1]; -(reverse(R[2:end-1]) .- R[1]) .+ R[1]]
        ZZ = @views [-(reverse(Z[2:end-1]) .- Z[1]) .+ Z[1]; Z[2:end-1]; -(reverse(Z[2:end-1]) .- Z[1]) .+ Z[1]]
        RR, ZZ = IMAS.resample_plasma_boundary(RR, ZZ; n_points=length(R) * 2, method=:cubic)
        if upper_x_point
            I = ZZ .< Z[1]
        else
            I = ZZ .> Z[1]
        end
        R = [R[1]; RR[I]; R[1]]
        Z = [Z[1]; ZZ[I]; Z[1]]
    else
        R, Z = IMAS.reorder_flux_surface!(R, Z, argmax(Z))
        izmin = argmin(Z)
        R1 = @views R[1:izmin]
        Z1 = @views Z[1:izmin]
        R2 = @views R[izmin:end]
        Z2 = @views Z[izmin:end]
        ΔZ = (Z1[end] - Z[1])
        ZZ1 = @views [Z1[2:end-1] .- ΔZ; Z1[2:end-1]; Z1[2:end-1] .+ ΔZ]
        RR1 = @views [-(reverse(R1[2:end-1]) .- R1[1]) .+ R1[1]; R1[2:end-1]; -(reverse(R1[2:end-1]) .- R1[1]) .+ R1[1]]
        ZZ2 = @views [Z2[2:end-1] .+ ΔZ; Z2[2:end-1]; Z2[2:end-1] .- ΔZ]
        RR2 = @views [-(reverse(R2[2:end-1]) .- R2[1]) .+ R2[1]; R2[2:end-1]; -(reverse(R2[2:end-1]) .- R2[1]) .+ R2[1]]
        RR1, ZZ1 = IMAS.resample_2d_path(RR1, ZZ1; n_points=length(R) * 3, method=:cubic)
        RR2, ZZ2 = IMAS.resample_2d_path(RR2, ZZ2; n_points=length(R) * 3, method=:cubic)
        crossings = IMAS.intersection(RR1, ZZ1, RR2, ZZ2).crossings
        if crossings[1][2] > crossings[2][2]
            ((RXU, ZXU), (RXL, ZXL)) = crossings
        else
            ((RXU, ZXU), (RXL, ZXL)) = crossings
        end
        I1 = (ZZ1 .< ZXU) .&& (ZZ1 .> ZXL)
        I2 = (ZZ2 .< ZXU) .&& (ZZ2 .> ZXL)
        R = [RXU; RR1[I1]; RXL; RR2[I2]; RXU]
        Z = [ZXU; ZZ1[I1]; ZXL; ZZ2[I2]; ZXU]
    end

    IMAS.reorder_flux_surface!(R, Z, R0, Z0)

    mxhb.RX = [RXU, RXL]
    mxhb.ZX = [ZXU, ZXL]
    mxhb.r_boundary = R
    mxhb.z_boundary = Z

    return mxhb
end

function fitMXHboundary(mxh::IMAS.MXH, nx::Int; n_points::Int=0, target_area::Float64=0.0, target_volume::Float64=0.0, debug::Bool=false)
    return fitMXHboundary(mxh; upper_x_point=nx ∈ (1, 2), lower_x_point=nx ∈ (-1, 2), n_points, target_area, target_volume, debug)
end

"""
    fitMXHboundary(mxh::IMAS.MXH; upper_x_point::Bool, lower_x_point::Bool, n_points::Int=0, target_area::Float64=0.0, target_volume::Float64=0.0, debug::Bool=false)

Find boundary such that the output MXH parametrization (with x-points) matches the input MXH parametrization (without x-points)

Optimization can be further constrained to target specific plasma cross-sectional area and and volume
"""
function fitMXHboundary(mxh::IMAS.MXH; upper_x_point::Bool, lower_x_point::Bool, n_points::Int=0, target_area::Float64=0.0, target_volume::Float64=0.0, debug::Bool=false)
    if ~upper_x_point && ~lower_x_point
        return MXHboundary(mxh; upper_x_point, lower_x_point, n_points)
    end

    symmetric = IMAS.is_updown_symmetric(mxh) && !xor(upper_x_point, lower_x_point)

    pr, pz = mxh()
    i = argmax(pz)
    RXU = pr[i]
    ZXU = pz[i]
    i = argmin(pz)
    RXL = pr[i]
    ZXL = pz[i]
    pr, pz = IMAS.resample_plasma_boundary(pr, pz; n_points=100)

    M = 12
    N = min(M, length(mxh.s))

    # change mxhb0 parameters so that the boundary  with x-points fit has the desired elongation, triangularity, squareness when fit with mxh0
    mxhb0 = MXHboundary(mxh; upper_x_point, lower_x_point, n_points)
    mxh0 = IMAS.MXH(mxhb0.r_boundary, mxhb0.z_boundary, M)

    function mxhb_from_params!(mxhb0::MXHboundary, params::AbstractVector{<:Real}; upper_x_point::Bool, lower_x_point::Bool, n_points::Int)
        L = 1
        mxhb0.mxh.Z0 = params[L]
        L += 1
        mxhb0.mxh.κ = params[L]
        if symmetric
            mxhb0.mxh.c0 = 0.0
            mxhb0.mxh.s = params[L+1:end]
            mxhb0.mxh.c = zero(mxhb0.mxh.s)
        else
            L += 1
            mxhb0.mxh.c0 = params[L]
            mxhb0.mxh.s = params[L+1:L+Integer((end - L) / 2)]
            mxhb0.mxh.c = params[L+Integer((end - L) / 2)+1:end]
        end
        return MXHboundary!(mxhb0; upper_x_point, lower_x_point, n_points)
    end

    function cost(params::AbstractVector{<:Real}; mxhb0, mxh0, upper_x_point, lower_x_point, n_points, target_area::Float64, target_volume::Float64)
        # mxhb0: contains the boundary with the x-point
        # mxhb0.mxh: is the MXH parametrization that leads to mxhb0 once x-points are set
        # mxh0: is the MHX fit to the mxhb0 boundary

        mxhb_from_params!(mxhb0, params; upper_x_point, lower_x_point, n_points)
        IMAS.MXH!(mxh0, mxhb0.r_boundary, mxhb0.z_boundary)

        pr0, pz0 = IMAS.resample_plasma_boundary(mxhb0.r_boundary, mxhb0.z_boundary; n_points=length(pr))

        # X-points and non X-points halves
        c = 0.0
        if upper_x_point
            i = argmax(mxhb0.ZX)
            c += ((mxhb0.RX[i] - RXU)^2 + (mxhb0.ZX[i] - ZXU)^2)
        else
            i = pz .< mxhb0.mxh.Z0
            c += sum((pr0[i] .- pr[i]) .^ 2 .+ (pz0[i] .- pz[i]) .^ 2) / sum(i)
        end
        if lower_x_point
            i = argmin(mxhb0.ZX)
            c += ((mxhb0.RX[i] - RXL)^2 + (mxhb0.ZX[i] - ZXL)^2)
        else
            i = pz .> mxhb0.mxh.Z0
            c += sum((pr0[i] .- pr[i]) .^ 2 .+ (pz0[i] .- pz[i]) .^ 2) / sum(i)
        end

        # geometric center
        c += (mxh0.R0 - mxh.R0)^2
        c += (mxh0.Z0 - mxh.Z0)^2

        # normalize by R0
        c /= mxhb0.mxh.R0^2

        # other (already normalized) shape parameters
        if target_area == 0 && target_volume == 0
            c += (mxh0.ϵ - mxh.ϵ)^2
        end
        c += (mxh0.κ - mxh.κ)^2
        c += (mxh0.c0 - mxh.c0)^2
        c += sum((mxh0.s[1:N] .- mxh.s[1:N]) .^ 2)
        c += sum((mxh0.c[1:N] .- mxh.c[1:N]) .^ 2)

        # make MXH match boundary of MHXboundary
        x_point_area = IMAS.area(mxhb0.r_boundary, mxhb0.z_boundary)
        marea0 = IMAS.area(mxh0()...)
        c += (abs(x_point_area - marea0) / x_point_area)^2

        # To avoid MXH solutions with kinks force area and convex_hull area to match
        mr, mz = mxhb0.mxh()
        hull = convex_hull(mr, mz; closed_polygon=true)
        mrch = [r for (r, z) in hull]
        mzch = [z for (r, z) in hull]
        marea = IMAS.area(mr, mz)
        mareach = IMAS.area(mrch, mzch)
        c += (abs(mareach - marea) / mareach)^2

        # Optionally match a given target area and or volume
        if target_area > 0
            c += ((x_point_area - target_area) / target_area)^2
        end
        if target_volume > 0
            x_point_volume = IMAS.revolution_volume(mxhb0.r_boundary, mxhb0.z_boundary)
            c += ((x_point_volume - target_volume) / target_volume)^2
        end

        return sqrt(c)
    end

    if symmetric
        params = vcat(mxh.Z0, mxh.κ, mxh.s)
    else
        params = vcat(mxh.Z0, mxh.κ, mxh.c0, mxh.s, mxh.c)
    end
    res = Optim.optimize(
        x -> cost(x; mxhb0, mxh0, upper_x_point, lower_x_point, n_points, target_area, target_volume),
        params,
        Optim.NelderMead(),
        Optim.Options(; iterations=1000, g_tol=1E-5)
    )
    mxhb_from_params!(mxhb0, res.minimizer; upper_x_point, lower_x_point, n_points)
    IMAS.MXH!(mxh0, mxhb0.r_boundary, mxhb0.z_boundary)

    if debug
        println(res)

        println("mxh")
        println(mxh)

        println("mxh0")
        println(mxh0)

        println("mxhb0.mxh")
        println(mxhb0.mxh)

        display(plot(mxh0))
    end

    return mxhb0
end

"""
    free_boundary_private_flux_constraint(R::AbstractVector{T}, Z::AbstractVector{T}; upper_x_point::Bool, lower_x_point::Bool, fraction::Float64=0.25, n_points::Int=10) where {T<:Real}

Given a closed boundary R,Z returns Rp,Zp constraint points on plausible private flux regions
"""
function free_boundary_private_flux_constraint(
    R::AbstractVector{T},
    Z::AbstractVector{T};
    upper_x_point::Bool,
    lower_x_point::Bool,
    fraction::Float64=0.25,
    n_points::Int=10) where {T<:Real}

    Rp = T[]
    Zp = T[]

    if fraction <= 0.0
        return Rp, Zp
    end

    for x_point in [-1, 1]

        if x_point == -1 && !lower_x_point
            continue
        elseif x_point == 1 && !upper_x_point
            continue
        end

        pr = deepcopy(R)
        pz = deepcopy(Z)

        # up-down and left-right mirror of lcfs with respect to the X-point
        if x_point == 1
            pr, pz = IMAS.reorder_flux_surface!(pr, pz, argmin(pz))
            Rx = pr[argmax(pz)]
            Zx = pz[argmax(pz)]
        else
            pr, pz = IMAS.reorder_flux_surface!(pr, pz, argmax(pz))
            Rx = pr[argmin(pz)]
            Zx = pz[argmin(pz)]
        end
        pr = -(pr .- Rx) .+ Rx
        pz = -(pz .- Zx) .+ Zx

        # select points withing a given radius
        d = sqrt.((pr .- Rx) .^ 2 .+ (pz .- Zx) .^ 2)
        index = d .< abs(maximum(pz) - minimum(pz)) * fraction
        pr = pr[index]
        pz = pz[index]

        # resample to give requested number of points
        pr, pz = IMAS.resample_plasma_boundary(pr, pz; n_points)

        append!(Rp, pr)
        append!(Zp, pz)
    end

    return Rp, Zp
end

function boundary_shape(R0::Real; p=nothing)
    return boundary_shape(IMAS.MXH(R0, 3); p=p)
end

"""
    private_flux_regions_from_lcfs(mr::AbstractArray{T}, mz::AbstractArray{T}, len::T; n_points::Integer=2)

Traces private flux regions starting from LCFS boundary
"""
function private_flux_regions_from_lcfs(mr::AbstractArray{T}, mz::AbstractArray{T}, len::T; n_points::Integer=2) where {T<:Real}
    mr = deepcopy(mr)
    mz = deepcopy(mz)
    IMAS.reorder_flux_surface!(mr, mz, sum(mr) / length(mr), sum(mz) / length(mz))

    iu = argmax(mz)
    il = argmin(mz)

    r_hfs = reverse(mr[iu:il])
    z_hfs = reverse(mz[iu:il])

    r_lfs = vcat(mr[il:end-1], mr[1:iu])
    z_lfs = vcat(mz[il:end-1], mz[1:iu])

    zz_u_hfs = range(mz[iu] + len, mz[iu], n_points)
    zz_u_lfs = range(mz[iu], mz[iu] + len, n_points)
    zz_l_hfs = range(mz[il] - len, mz[il], n_points)
    zz_l_lfs = range(mz[il], mz[il] - len, n_points)
    rr_u_hfs = IMAS.interp1d(z_lfs, r_lfs).(zz_u_hfs)
    rr_u_lfs = IMAS.interp1d(z_hfs, r_hfs).(zz_u_lfs)
    rr_l_hfs = IMAS.interp1d(z_lfs, r_lfs).(zz_l_hfs)
    rr_l_lfs = IMAS.interp1d(z_hfs, r_hfs).(zz_l_lfs)

    rr_u = vcat(rr_u_hfs[1:end-1], rr_u_lfs[2:end])
    zz_u = vcat(zz_u_hfs[1:end-1], zz_u_lfs[2:end])
    rr_l = vcat(rr_l_hfs[1:end-1], rr_l_lfs[2:end])
    zz_l = vcat(zz_l_hfs[1:end-1], zz_l_lfs[2:end])

    return rr_u, zz_u, rr_l, zz_l
end

"""
    check_ped_finder(profile::AbstractVector{<:Real},psi_norm::AbstractVector{<:Real})

Plots the profile and fitted profile for pedestal finding (region outside pedestal not important for fit)
"""
function check_ped_finder(profile::AbstractVector{<:Real}, psi_norm::AbstractVector{<:Real})
    ped_height, ped_width = IMAS.pedestal_finder(profile, psi_norm)
    plot(psi_norm, profile; label="original profile")
    return plot!(
        psi_norm,
        IMAS.Hmode_profiles(profile[end], ped_height, profile[1], length(profile), 2.0, 2.0, ped_width);
        label="fitted profile (pedestal region is important only)"
    )
end

"""
    ip_from_q_star(a::Real, B0::Real, q_star::Real, R0::Real, κ::Real)

Calculate ip from the cylindrical equivalent safety factor
"""
function ip_from_q_cyl(a::Real, B0::Real, q_cyl::Real, R0::Real, κ::Real)
    return 1e6 * abs((5 * a^2 * B0) / (R0 * (q_cyl)) * (1 + κ^2) / 2)
end

"""
    B0_from_Bmax(Bmax::Real, distance_tf_to_plasma::Real, a::Real, R0::Real)

Calculate the on axis magnetic field from the maximum magnetic field at the tf coils
"""
function B0_from_b_max(b_max::Real, distance_tf_to_plasma::Real, a::Real, R0::Real)
    return b_max * (1 - (a + distance_tf_to_plasma) / R0)
end