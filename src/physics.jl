import SpecialFunctions
import QuadGK
import Interpolations
import PolygonOps
import Interact

#= =============== =#
#  Shape functions  #
#= =============== =#

function layer_shape_message(shape_function_index)
    return "layer.shape=$(shape_function_index) is invalid. Valid options are:
        -2: Offset & convex-hull
        -1: Offset
         1: Priceton D exact  (shape_parameters = [])
         2: Priceton D approx (shape_parameters = [])
         3: Priceton D scaled (shape_parameters = [height])
         4: rectangle         (shape_parameters = [height])
         5: tripple-arc       (shape_parameters = [height, small_radius, mid_radius, small_coverage, mid_coverage])
         6: miller            (shape_parameters = [elongation, triangularity])
         7: spline            (shape_parameters = [hfact, rz...)
         8: silo              (shape_parameters = [h_start, h_end)
       10x: shape + z_offset  (shape_parameters = [..., z_offset]) 
10x: shape + z_offset  (shape_parameters = [..., z_offset]) 
       10x: shape + z_offset  (shape_parameters = [..., z_offset]) 
       "
end

function initialize_shape_parameters(shape_function_index, r_obstruction, z_obstruction, r_start, r_end, target_clearance)
    height = maximum(z_obstruction) - minimum(z_obstruction) + target_clearance * 2.0
    z_offset = (maximum(z_obstruction) + minimum(z_obstruction)) / 2.0
    shape_parameters = nothing
    if shape_function_index in [Int(_offset_), Int(_negative_offset_), Int(_convex_hull_)]
        return nothing
    else
        shape_index_mod = shape_function_index
        is_negative_D = false
        if shape_index_mod > 1000
            shape_index_mod = mod(shape_function_index, 1000)
            is_negative_D = true
        end
        is_z_offset = false
        if shape_index_mod > 100
            shape_index_mod = mod(shape_function_index, 100)
            is_z_offset = true
        end
        if shape_index_mod == Int(_princeton_D_exact_)
            shape_parameters = Real[]
        elseif shape_index_mod == Int(_princeton_D_)
            shape_parameters = Real[]
        elseif shape_index_mod == Int(_princeton_D_scaled_)
            shape_parameters = [height]
        elseif shape_index_mod == Int(_rectangle_)
            shape_parameters = [height]
        elseif shape_index_mod == Int(_triple_arc_)
            shape_parameters = [log10(height), log10(1E-3), log10(1E-3), log10(45), log10(45)]
        elseif shape_index_mod == Int(_miller_)
            _, imaxr = findmax(r_obstruction)
            _, iminr = findmin(r_obstruction)
            _, imaxz = findmax(z_obstruction)
            _, iminz = findmin(z_obstruction)
            r_at_max_z, max_z = r_obstruction[imaxz], z_obstruction[imaxz]
            r_at_min_z, min_z = r_obstruction[iminz], z_obstruction[iminz]
            z_at_max_r, max_r = z_obstruction[imaxr], r_obstruction[imaxr]
            z_at_min_r, min_r = z_obstruction[iminr], r_obstruction[iminr]
            a = 0.5 * (max_r - min_r)
            b = 0.5 * (max_z - min_z)
            R = 0.5 * (max_r + min_r)
            elongation = b / a
            triup = (R - r_at_max_z) / a
            tridown = (R - r_at_min_z) / a
            shape_parameters = [elongation, (triup + tridown) / 2.0]
        elseif shape_index_mod == Int(_spline_)
            n = 1
            R = range(r_start, r_end, length=2 + n)[2:end-1]
            Z = range(height / 2.0, height / 2.0, length=2 + n)[2:end-1]
            shape_parameters = Float64[0.8]
            for (r, z) in zip(R, Z)
                append!(shape_parameters, [r, z])
            end
        elseif shape_index_mod == Int(_silo_)
            shape_parameters = [height, height / 2.0]
        end
    end
    if shape_parameters === nothing
        error(shape_function_index)
    end
    if is_z_offset
        push!(shape_parameters, z_offset)
    end
    return shape_parameters
end

function shape_function(shape_function_index)
    func = nothing
    if shape_function_index in [Int(_offset_), Int(_negative_offset_), Int(_convex_hull_)]
        return nothing
    else
        shape_index_mod = shape_function_index
        is_negative_D = false
        if shape_index_mod > 1000
            shape_index_mod = mod(shape_function_index, 1000)
            is_negative_D = true
        end
        is_z_offset = false
        if shape_index_mod > 100
            shape_index_mod = mod(shape_function_index, 100)
            is_z_offset = true
        end
        if shape_index_mod in Int(_princeton_D_exact_)
            func = princeton_D_exact
        elseif shape_index_mod == Int(_princeton_D_)
            func = princeton_D_approx
        elseif shape_index_mod == Int(_princeton_D_scaled_)
            func = princeton_D_scaled
        elseif shape_index_mod == Int(_rectangle_)
            func = (r_start, r_end, height) -> rectangle_shape(r_start, r_end, height; n_points=100)
        elseif shape_index_mod == Int(_triple_arc_)
            func = triple_arc
        elseif shape_index_mod == Int(_miller_)
            func = miller_Rstart_Rend
        elseif shape_index_mod == Int(_spline_)
            func = spline_shape
        elseif shape_index_mod == Int(_silo_)
            func = silo
        end
    end
    if func === nothing
        error(layer_shape_message(shape_index_mod))
    end

    # zoffset
    zfunc = func
    if is_z_offset
        zfunc(args...) = begin
            R, Z = func(args[1:end-1]...)
            Z .+= args[end]
            return R, Z
        end
    end

    # neg-D
    dfunc = zfunc
    if is_negative_D
        dfunc(args...) = begin
            R, Z = zfunc(args...)
            R = -(R .- args[1]) .+ args[2]
            return reverse(R), reverse(Z)
        end
    end

    # uniform resampling
    function resampled_zfunc(args...; resample=true)
        if resample
            return IMAS.resample_2d_line(dfunc(args...)...)
        else
            return dfunc(args...)
        end
    end

    return resampled_zfunc
end
"""
    optimize_shape(r_obstruction, z_obstruction, target_clearance, func, r_start, r_end, shape_parameters; verbose=false, time_limit=60)

Find shape parameters that generate smallest shape and target clearance from an obstruction
"""
function optimize_shape(r_obstruction, z_obstruction, target_clearance, func, r_start, r_end, shape_parameters; verbose=false, time_limit=60)
    if length(shape_parameters) in [0, 1]
        func(r_start, r_end, shape_parameters...)
        return shape_parameters
    end

    function cost_TF_shape(r_obstruction, z_obstruction, rz_obstruction, target_clearance, func, r_start, r_end, shape_parameters; verbose=false)
        R, Z = func(r_start, r_end, shape_parameters...)

        # disregard near r_start and r_end where optimizer has no control and shape is allowed to go over obstruction
        index = abs.((R .- r_start) .* (R .- r_end)) .> (target_clearance / 1000)
        R = R[index]
        Z = Z[index]

        # no polygon crossings  O(N)
        inpoly = [PolygonOps.inpolygon((r, z), rz_obstruction) for (r, z) in zip(R, Z)]
        cost_inside = sum(inpoly)

        # target clearance  O(1)
        minimum_distance = IMAS.minimum_distance_two_shapes(R, Z, r_obstruction, z_obstruction)
        cost_min_clearance = (minimum_distance - target_clearance) / target_clearance
        mean_distance_error = IMAS.mean_distance_error_two_shapes(R, Z, r_obstruction, z_obstruction, target_clearance)
        cost_mean_distance = mean_distance_error / target_clearance

        # favor up/down symmetric solutions
        cost_up_down_symmetry = abs(maximum(Z) + minimum(Z)) / maximum(abs.(Z))

        if verbose
            @show minimum_distance
            @show mean_distance_error
            @show target_clearance
            @show cost_min_clearance^2
            @show cost_mean_distance^2
            @show cost_inside^2
            @show cost_up_down_symmetry^2
        end

        # return cost
        return cost_min_clearance^2 + cost_mean_distance^2 + cost_inside^2 + 0.1 * cost_up_down_symmetry^2
    end

    rz_obstruction = collect(zip(r_obstruction, z_obstruction))
    initial_guess = copy(shape_parameters)
    # res = optimize(shape_parameters-> cost_TF_shape(r_obstruction, z_obstruction, rz_obstruction, target_clearance, func, r_start, r_end, shape_parameters),
    #                initial_guess, Newton(), Optim.Options(time_limit=time_limit); autodiff=:forward)
    res = Optim.optimize(shape_parameters -> cost_TF_shape(r_obstruction, z_obstruction, rz_obstruction, target_clearance, func, r_start, r_end, shape_parameters),
        initial_guess, length(shape_parameters) == 1 ? Optim.BFGS() : Optim.NelderMead(), Optim.Options(time_limit=time_limit); autodiff=:forward)
    if verbose
        println(res)
    end
    shape_parameters = Optim.minimizer(res)
    # R, Z = func(r_start, r_end, shape_parameters...; resample=false)
    # plot(func(r_start, r_end, initial_guess...); markershape=:x)
    # plot!(r_obstruction, z_obstruction, ; markershape=:x)
    # display(plot!(R, Z; markershape=:x, aspect_ratio=:equal))
    # cost_TF_shape(r_obstruction, z_obstruction, rz_obstruction, obstruction_area, target_clearance, func, r_start, r_end, shape_parameters; verbose=true)
    return shape_parameters
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
    ellipse(a, b, t0, t1, x0, z0; n_points=100)

Simple ellipse shape function
"""
function ellipse(a, b, t0, t1, x0, z0; n_points=100)
    t = LinRange(t0, t1, n_points)
    x = a .* cos.(t) .+ x0
    z = b .* sin.(t) .+ z0
    return x, z
end

"""
    princeton_D_approx(r_start, r_end; n_points=100)

This retuns the an approximate version of the "Princeton-D" constant tension shape for a TF coil that is built with ellipses
References: Gralnick, S. L.; Tenney, F. H. Analytic Solutions for Constant‐tension Coil Shapes. J. Appl. Phys. 1976, 47, 7
"""
function princeton_D_approx(r_start, r_end; n_points=100)
    r1 = r_start
    r2 = r_end
    k = 0.5 * log(r2 / r1)
    r0 = sqrt(r1 * r2)

    # analytic equations for coil parameters
    centerpost_maxz = 2 * pi * r0 * k * SpecialFunctions.besseli(1, k) / 2 # Gralnick Eq. 28
    coil_maxz = pi * r0 * k * (SpecialFunctions.besseli(1, k) + struveL(1, k) + 2 / pi) / 2  # Gralnick Eq. 34

    # make ellipses to approximate equal-tension arc
    x1, z1 = ellipse(r0 - r1, coil_maxz - centerpost_maxz, pi, pi / 2, r0, centerpost_maxz; n_points)
    x2, z2 = ellipse(r2 - r0, coil_maxz, pi / 2, 0, r0, 0; n_points)

    x = vcat(x1[1:end-1], x2)
    z = vcat(z1[1:end-1], z2)

    x = vcat(x, x[end-1:-1:1], x[1])
    z = vcat(z, -z[end-1:-1:1], z[1])

    return x, z
end

"""
    function princeton_D_scaled(r_start, r_end, height; n_points=100)

This routine calculates a "shortened" TF coil shape that foregoes the equal-tension "Princeton-Dee" for 
a squater, more space-efficient shape. It replicates the inboard curve of the equal-tension arc, but
decreases the height of the coil to match a given value.
"""
function princeton_D_scaled(r_start, r_end, height; n_points=100)
    r1 = r_start
    r2 = r_end
    k = 0.5 * log(r2 / r1)
    r0 = sqrt(r1 * r2)

    centerpost_maxz = 2 * pi * r0 * k * SpecialFunctions.besseli(1, k) / 2 # Gralnick Eq. 28
    coil_maxz = pi * r0 * k * (SpecialFunctions.besseli(1, k) + struveL(1, k) + 2 / pi) / 2  # Gralnick Eq. 34

    centerpost_height = height - (coil_maxz - centerpost_maxz) * 2

    inboard_curve_dz = coil_maxz - centerpost_maxz
    centerpost_maxz = centerpost_height / 2
    coil_maxz = centerpost_maxz + inboard_curve_dz

    # make ellipses to connect points
    x1, z1 = ellipse(r0 - r1, coil_maxz - centerpost_maxz, pi, pi / 2, r0, centerpost_maxz; n_points)
    x2, z2 = ellipse(r2 - r0, coil_maxz, pi / 2, 0, r0, 0; n_points)

    x = vcat(x1[1:end-1], x2)
    z = vcat(z1[1:end-1], z2)

    x = vcat(x, x[end-1:-1:1], x[1])
    z = vcat(z, -z[end-1:-1:1], z[1])

    return x, z
end

"""
    rectangle_shape(r_start::Real, r_end::Real, z_low::Real, z_high::Real; n_points = 5)

Asymmetric rectangular contour
layer[:].shape = 2
"""
function rectangle_shape(r_start::Real, r_end::Real, z_low::Real, z_high::Real; n_points::Int=5)
    if n_points == 5
        R = [r_start, r_end, r_end, r_start, r_start]
        Z = [z_low, z_low, z_high, z_high, z_low]
    else
        R = vcat(
            range(r_start, r_end; length=n_points),
            range(r_end, r_end; length=n_points)[2:end],
            range(r_end, r_start; length=n_points)[2:end],
            range(r_start, r_start; length=n_points)[2:end],
            r_start)
        Z = vcat(
            range(z_low, z_low; length=n_points),
            range(z_low, z_high; length=n_points)[2:end],
            range(z_high, z_high; length=n_points)[2:end],
            range(z_high, z_low; length=n_points)[2:end],
            z_low)
    end
    return R, Z
end

"""
    rectangle_shape(r_start::Real, r_end::Real, height::Real; n_points = 5)

Symmetric rectangular contour
"""
function rectangle_shape(r_start::Real, r_end::Real, height::Real; n_points::Int=5)
    Δ = height / 2.0
    return rectangle_shape(r_start, r_end, -Δ, Δ; n_points)
end

"""
    function triple_arc(r_start::Real,
        r_end::Real,
        height::Real,
        small_radius::Real,
        mid_radius::Real,
        small_coverage::Real,
        mid_coverage::Real;
        min_small_radius_fraction::Real=0.5,
        min_mid_radius_fraction::Real=min_small_radius_fraction*2.0,
        n_points::Int=400)

TrippleArc contour
Angles are in degrees
height, small_radius, mid_radius, small_coverage, mid_coverage are 10^exponent (to ensure positiveness)
"""
function triple_arc(
    r_start::Real,
    r_end::Real,
    height::Real,
    small_radius::Real,
    mid_radius::Real,
    small_coverage::Real,
    mid_coverage::Real;
    min_small_radius_fraction::Real=0.2,
    min_mid_radius_fraction::Real=min_small_radius_fraction * 2.0,
    n_points::Int=400)

    height = 10^height / 2.0
    small_radius = 10^small_radius + height * min_small_radius_fraction
    mid_radius = 10^mid_radius + height * min_mid_radius_fraction
    small_coverage = 10^small_coverage * pi / 180
    mid_coverage = 10^mid_coverage * pi / 180

    asum = small_coverage + mid_coverage
    n_points = floor(Int, n_points / 4)

    # small arc
    theta = LinRange(0, small_coverage, n_points)
    small_arc_R = r_start .+ small_radius .* (1 .- cos.(theta))
    small_arc_Z = height .+ small_radius .* sin.(theta)

    # mid arc
    theta = LinRange(small_coverage, asum, n_points)
    mid_arc_R = small_arc_R[end] .+ mid_radius .* (cos.(small_coverage) .- cos.(theta))
    mid_arc_Z = small_arc_Z[end] .+ mid_radius .* (sin.(theta) .- sin.(small_coverage))

    # large arc
    theta = LinRange(theta[end], pi, n_points)
    large_radius = mid_arc_Z[end] / sin(pi - asum)
    large_arc_R = mid_arc_R[end] .+ large_radius .* (cos.(pi .- theta) .- cos.(pi .- asum))
    large_arc_Z = mid_arc_Z[end] .- large_radius .* (sin(asum) .- sin.(pi .- theta))

    R = vcat(small_arc_R, mid_arc_R[2:end], large_arc_R[2:end])
    R = vcat(R, reverse(R)[2:end])
    Z = vcat(small_arc_Z, mid_arc_Z[2:end], large_arc_Z[2:end])
    Z = vcat(Z, -reverse(Z)[2:end])

    # Add vertical
    R = vcat(LinRange(r_start, r_start, n_points), R)
    Z = vcat(LinRange(-height, height, n_points), Z)

    # Resize to ensure r_start to r_end
    factor = (r_end - r_start) / (maximum(R) - minimum(R))
    Z = Z .* factor
    R = (R .- minimum(R)) .* factor .+ r_start

    return R, Z
end

"""
    miller(R0, inverse_aspect_ratio, elongation, triangularity, n_points)

Miller contour
layer[:].shape = 4
"""
function miller(R0::Real, rmin_over_R0::Real, elongation::Real, triangularity::Real; n_points::Int=401)
    θ = range(0, 2 * pi; length=n_points)
    triangularity = mirror_bound(triangularity, -1.0, 1.0)
    δ₀ = asin(triangularity)
    R = R0 * (1 .+ rmin_over_R0 .* cos.(θ .+ δ₀ * sin.(θ)))
    Z = R0 * (rmin_over_R0 * elongation * sin.(θ))
    R[end] = R[1]
    Z[end] = Z[1]
    return R, Z
end

"""
    miller_Rstart_Rend(r_start, r_end, elongation, triangularity, n_points)

Miller contour
"""
function miller_Rstart_Rend(r_start::Real, r_end::Real, elongation::Real, triangularity::Real; n_points::Int=401)
    return miller((r_end + r_start) / 2.0, (r_end - r_start) / (r_end + r_start), elongation, triangularity; n_points)
end

"""
    spline_shape(r::Real, z::Real; n_points::Int=101)

Spline contour
"""
function spline_shape(r::Vector{T}, z::Vector{T}; n_points::Int=101) where {T<:Real}
    r = vcat(r[1], r[1], r, r[end], r[end])
    z = vcat(0, z[1] / 2, z, z[end] / 2, 0)
    d = cumsum(sqrt.(vcat(0, diff(r)) .^ 2.0 .+ vcat(0, diff(z)) .^ 2.0))

    itp_r = Interpolations.interpolate(d, r, Interpolations.FritschButlandMonotonicInterpolation())
    itp_z = Interpolations.interpolate(d, z, Interpolations.FritschButlandMonotonicInterpolation())

    D = LinRange(d[1], d[end], n_points)
    R, Z = itp_r.(D), itp_z.(D)
    R[end] = R[1]
    Z[end] = Z[1]
    return R, Z
end

function spline_shape(r_start::Real, r_end::Real, hfact::Real, rz...; n_points::Int=101)
    rz = collect(rz)
    R = rz[1:2:end]
    Z = rz[2:2:end]
    hfact_max = 0.5 + (minimum(R) - r_start) / (r_end - r_start) / 2.0
    hfact = min(abs(hfact), hfact_max)
    h = maximum(Z) * hfact
    r = vcat(r_start, R, r_end, reverse(R), r_start)
    z = vcat(h, Z, 0, -reverse(Z), -h)
    return spline_shape(r, z; n_points=n_points)
end

function xy_polygon(x::Vector{<:Real}, y::Vector{<:Real})
    x = deepcopy(x)
    y = deepcopy(y)
    if (x[1] != x[end]) && (x[1] ≈ x[end])
        x[end] = x[1]
    end
    if (y[1] != y[end]) && (y[1] ≈ y[end])
        y[end] = y[1]
    end
    if (x[1] != x[end]) || (y[1] != y[end])
        push!(x, x[1])
        push!(y, y[1])
    end
    coords = [collect(map(collect, zip(x, y)))]
    return LibGEOS.Polygon(coords)
end

function xy_polygon(layer::Union{IMAS.build__layer,IMAS.build__structure})
    return xy_polygon(layer.outline.r, layer.outline.z)
end

"""
    volume_no_structures(layer::IMAS.build__layer, structures::IMAS.IDSvector{IMAS.build__structure})

Returns volume of the layer without structures
"""
function volume_no_structures(layer::IMAS.build__layer, structures::IMAS.IDSvector{IMAS.build__structure})
    vol = 0.0
    for structure in structures
        vol += layer_structure_intersect_volume(layer, structure)
    end
    return layer.volume - vol
end

IMAS.expressions["build.layer[:].volume_no_structures"] =
    (; build, layer, _...) -> volume_no_structures(layer, build.structure)

"""
    layer_structure_intersect_volume(layer::IMAS.build__layer, structure::IMAS.build__structure)

Returns volume of the intersection between build layer volume and structure volume
"""
function layer_structure_intersect_volume(layer::IMAS.build__layer, structure::IMAS.build__structure)
    if layer.type in [Int(_in_), Int(_out_)]
        return layer.volume
    elseif layer.type in [Int(_tf_), Int(_plasma_)]
        return layer.volume
    elseif layer.fs ∈ [Int(_hfs_), Int(_lfs_)]
        i = IMAS.index(layer)
        if layer.fs == Int(_hfs_)
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
    for poly in LibGEOS.coordinates(LibGEOS.intersection(ring_poly, structure_poly))
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
    silo(r_start, r_end, height_start, height_end)  
"""
function silo(r_start, r_end, height_start, height_end)
    height_start = abs(height_start)
    height_end = abs(height_end)
    height_end = min(max(height_end, height_start * 0.0), height_start * 0.9)
    x, y = ellipse(r_end - r_start, height_start - height_end, 0, pi / 2, r_start, height_end)
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
    rX::Vector{<:Real}
    zX::Vector{<:Real}
    r_boundary::Vector{<:Real}
    z_boundary::Vector{<:Real}
end

function add_xpoint(mr::Vector{T}, mz::Vector{T}, i::Integer, R0::T, Z0::T, α::T) where {T<:Real}
    rX = mr[i] .* α .+ R0 .* (1.0 .- α)
    zX = mz[i] .* α .+ Z0 .* (1.0 .- α)
    rrrzzz = FUSE.convex_hull(collect(zip(vcat(mr, rX), vcat(mz, zX))); closed_polygon=true)
    rrr = [r for (r, z) in rrrzzz]
    zzz = [z for (r, z) in rrrzzz]
    return rrr, zzz
end

function add_xpoint(mr::Vector{T}, mz::Vector{T}, R0::T, Z0::T; upper::Bool) where {T<:Real}

    function cost(mr, mz, i, R0, Z0, α)
        rrr, zzz = add_xpoint(mr, mz, i, R0, Z0, α[1])
        tmp = 1.0 - maximum(abs.(IMAS.curvature(rrr, zzz)))
        return tmp
    end

    if upper
        i = argmax(abs.(IMAS.curvature(mr, mz)) .* (mz .> Z0))
    else
        i = argmax(abs.(IMAS.curvature(mr, mz)) .* (mz .< Z0))
    end
    res = FUSE.Optim.optimize(α -> cost(mr, mz, i, R0, Z0, α), 1.0, 1.5, FUSE.Optim.GoldenSection())
    α = res.minimizer[1]
    rX = mr[i] .* α .+ R0 .* (1.0 .- α)
    zX = mz[i] .* α .+ Z0 .* (1.0 .- α)

    return rX, zX
end

function MXH_boundary(mxh::IMAS.MXH; upper_x_point::Bool, lower_x_point::Bool)
    mr, mz = mxh()
    R0 = mxh.R0
    Z0 = mxh.Z0

    rX = Float64[]
    zX = Float64[]
    if upper_x_point
        rXU, zXU = add_xpoint(mr, mz, R0, Z0; upper=true)
        push!(rX, rXU)
        push!(zX, zXU)
    end
    if lower_x_point
        rXL, zXL = add_xpoint(mr, mz, R0, Z0; upper=false)
        push!(rX, rXL)
        push!(zX, zXL)
    end

    rrrzzz = FUSE.convex_hull(collect(zip(vcat(mr, rX), vcat(mz, zX))); closed_polygon=true)
    rrr = [r for (r, z) in rrrzzz]
    zzz = [z for (r, z) in rrrzzz]
    rrr, zzz = IMAS.resample_2d_line(rrr, zzz; npoints=length(mr))
    IMAS.reorder_flux_surface!(rrr, zzz, R0, Z0)

    return MXHboundary(mxh, rX, zX, rrr, zzz)
end

"""
    boundary_shape(mxh::IMAS.MXH; upper_x_point::Bool, lower_x_point::Bool, p=nothing)

Plot and manipulate Miller Extended Harmonic (MXH) boundary
"""
function boundary_shape(mxh::IMAS.MXH; upper_x_point::Bool, lower_x_point::Bool, p::Union{Nothing,Plots.Plot})
    n = 101
    if length(mxh.c) != 3
        mxh = IMAS.MXH(mxh()..., 3)
    end
    Interact.@manipulate for ϵ in Interact.slider(LinRange(0.0, 1.0, n), value=mxh.ϵ, label="ϵ"),
        κ in Interact.slider(LinRange(1, 3, n), value=mxh.κ, label="κ"),
        #tilt in Interact.slider(LinRange(-1,1,n);value=mxh.c0,label="tilt"),
        #ovality in Interact.slider(LinRange(-1,1,n);value=mxh.c[1],label="ovality"),
        triangularity in Interact.slider(LinRange(-1, 1, n), value=mxh.s[1], label="δ"),
        #s_shape in Interact.slider(LinRange(-1,1,n);value=mxh.c[2],label="s1"),
        #squareness in Interact.slider(LinRange(-1, 1, n), value=-(mxh.s[2] + mxh.s[3]) / 2.0, label="ζ"),
        squareness in Interact.slider(LinRange(-1, 1, n), value=-mxh.s[2], label="ζ"),
        #c3 in Interact.slider(LinRange(-1,1,n);value=mxh.c[3],label="s2"),
        pentagonnes in Interact.slider(LinRange(-1, 1, n); value=mxh.s[3], label="⬠")

        mxh.ϵ = ϵ
        mxh.κ = κ
        #mxh.c0=tilt
        #mxh.c[1]=ovality
        mxh.s[1] = triangularity
        #mxh.c[2]=s_shape
        mxh.s[2] = -squareness
        #mxh.c[3]=c3
        mxh.s[3] = -pentagonnes

        if p === nothing
            q = plot()
        else
            q = deepcopy(p)
            plot(q)
        end

        if upper_x_point || lower_x_point
            mxhb = MXH_boundary(mxh; upper_x_point, lower_x_point)
            plot!(q, mxh, color=:gray, linewidth=1.5)
            plot!(q, mxhb.r_boundary, mxhb.z_boundary, color=:black, linewidth=2, label="")
        else
            plot!(q, mxh, color=:black, linewidth=2, label="")
        end
    end
end

function boundary_shape(R0::Real; p=nothing)
    return boundary_shape(IMAS.MXH(R0, 3); p=p)
end
