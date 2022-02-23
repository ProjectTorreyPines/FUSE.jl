import LibGEOS
import Interpolations
import Contour
import LazySets
import DataStructures

#= ================== =#
#  init core_profiles  #
#= ================== =#
function init_core_profiles(cp::IMAS.core_profiles; kw...)
    for item in keys(kw)
        IMAS.set_time_array(cp.global_quantities, item, kw[item])
    end
    return cp
end

#= ========== =#
#  init build  #
#= ========== =#
"""
    init_build(bd::IMAS.build; verbose = false, layers...)

Initialize build IDS based on center stack layers (thicknesses)

NOTE: layer[:].type and layer[:].material follows from naming of layers
*   0 ...gap... : vacuum
*   1 OH: ohmic coil
*   2 TF: toroidal field coil
*   3 shield...: neutron shield
*   4 blanket...: neutron blanket
*   5 wall....: 
*  -1 ...plasma...: 

layer[:].hfs is set depending on if "hfs" or "lfs" appear in the name

layer[:].identifier is created as a hash of then name removing "hfs" or "lfs"
"""
function init_build(bd::IMAS.build; verbose = false, layers...)
    # empty build IDS
    empty!(bd)
    # assign layers
    resize!(bd.layer, length([layer_name for (layer_name, layer_thickness) in layers if layer_thickness >= 0.0]))
    k = 0
    radius_start = 0.0
    radius_end = 0.0
    for (layer_name, layer_thickness) in layers
        if layer_thickness < 0.0
            continue
        end
        k += 1
        layer = bd.layer[k]
        layer.thickness = layer_thickness
        layer.name = replace(String(layer_name), "_" => " ")
        if occursin("gap", lowercase(layer.name))
            layer.type = 0
            layer.material = "vacuum"
        elseif uppercase(layer.name) == "OH"
            layer.type = 1
        elseif occursin("TF", uppercase(layer.name))
            layer.type = 2
        elseif occursin("shield", lowercase(layer.name))
            layer.type = 3
        elseif occursin("blanket", lowercase(layer.name))
            layer.type = 4
        elseif occursin("wall", lowercase(layer.name))
            layer.type = 5
        end
        if occursin("hfs", lowercase(layer.name))
            layer.hfs = 1
        elseif occursin("lfs", lowercase(layer.name))
            layer.hfs = -1
        else
            layer.hfs = 0
        end
        if occursin("plasma", lowercase(layer.name))
            layer.type = -1
            layer.material = "vacuum"
        end
        layer.identifier = UInt(hash(replace(replace(lowercase(layer.name), "hfs" => ""), "lfs" => "")))
        radius_end += layer.thickness
        if verbose
            println("$(rpad(string(k),4)) $(rpad(layer.name, 17)) Δr=$(@sprintf("%3.3f",layer.thickness))    R=[$(@sprintf("%3.3f",radius_start)) <=> $(@sprintf("%3.3f",radius_end))]")
        end
        radius_start = radius_end
    end
    if ismissing(bd.layer[end], :material) || bd.layer[end].material != "vacuum"
        error("Material of last layer ($(bd.layer[end].name)) must be `vacuum`")
    end

    return bd
end

function init_build(bd::IMAS.build, layers::AbstractDict; verbose = false)
    nt = (; zip([Symbol(k) for k in keys(layers)], values(layers))...)
    init_build(bd; verbose, nt...)
end

"""
    init_build(
        bd::IMAS.build,
        eq::IMAS.equilibrium;
        tf_shape_index::Int = 3,
        is_nuclear_facility::Bool = true,
        pf_inside_tf::Bool = false,
        pf_outside_tf::Bool = true)

Initialization of build IDS based on equilibrium time_slice
"""
function init_build(dd::IMAS.dd;
    tf_shape_index::Int = 3,
    is_nuclear_facility::Bool = true,
    pf_inside_tf::Bool = false,
    pf_outside_tf::Bool = true,
    verbose::Bool = false)

    if !(ismissing(dd.wall.description_2d, [1, :limiter, :unit, 1, :outline, :r]))
        rmin = minimum(dd.wall.description_2d[1].limiter.unit[1].outline.r)
        rmax = maximum(dd.wall.description_2d[1].limiter.unit[1].outline.r)
    else
        eqt = dd.equilibrium.time_slice[]
        rmin = eqt.boundary.geometric_axis.r - eqt.boundary.minor_radius
        rmax = eqt.boundary.geometric_axis.r + eqt.boundary.minor_radius
        gap = (rmax - rmin) / 20.0 # plasma-wall gap
        gap = (rmax - rmin) / 20.0 # plasma-wall gap
        rmin -= gap
        rmax += gap
    end

    if is_nuclear_facility
        n_hfs_layers = 6
        dr = rmin / n_hfs_layers
        init_build(dd.build;
            verbose,
            gap_OH = dr * 2.0,
            OH = dr,
            hfs_TF = dr,
            gap_hfs_TF_shield = pf_inside_tf ? 0 : -1,
            hfs_shield = dr / 2.0,
            hfs_blanket = dr,
            hfs_wall = dr / 2.0,
            plasma = rmax - rmin,
            lfs_wall = dr / 2.0,
            lfs_blanket = dr * 2,
            lfs_shield = dr / 2.0,
            gap_lfs_TF_shield = dr * (pf_inside_tf ? 4 : -1),
            lfs_TF = dr,
            gap_cryostat = dr * (pf_outside_tf ? 5 : 1))

    else
        n_hfs_layers = 4.5
        dr = rmin / n_hfs_layers
        init_build(dd.build;
            verbose,
            gap_OH = dr * 2.0,
            OH = dr,
            hfs_TF = dr,
            gap_hfs_TF_wall = pf_inside_tf ? 0 : -1,
            hfs_wall = dr / 2.0,
            plasma = rmax - rmin,
            lfs_wall = dr / 2.0,
            gap_lfs_TF_wall = dr * (pf_inside_tf ? 2.25 : -1),
            lfs_TF = dr,
            gap_cryostat = dr * (pf_outside_tf ? 3 : 1))
    end

    # TF coils
    dd.build.tf.coils_n = 16

    # cross-section outlines
    if !(ismissing(dd.wall.description_2d, [1, :limiter, :unit, 1, :outline, :r]))
        build_cx(dd.build, dd.wall.description_2d[1].limiter.unit[1].outline.r, dd.wall.description_2d[1].limiter.unit[1].outline.z, tf_shape_index)
    else
        build_cx(dd.build, eqt, tf_shape_index)
    end

    return dd.build
end

function init_build(dd::IMAS.dd, gasc::GASC;
    no_small_gaps::Bool = true,
    tf_shape_index::Int = 3,
    verbose::Bool = false)
    gasc = gasc.solution

    # build
    norm = gasc["OUTPUTS"]["radial build"]["innerPlasmaRadius"]

    radial_build = DataStructures.OrderedDict()
    radial_build["gap_OH"] = gasc["OUTPUTS"]["radial build"]["innerSolenoidRadius"]
    radial_build["OH"] = gasc["INPUTS"]["radial build"]["rbOH"] * norm

    radial_build["hfs_gap_TF"] = gasc["INPUTS"]["radial build"]["gapTFOH"] * norm
    radial_build["hfs_TF"] = gasc["INPUTS"]["radial build"]["rbTF"] * norm
    if no_small_gaps
        radial_build["hfs_TF"] += radial_build["hfs_gap_TF"]
        pop!(radial_build, "hfs_gap_TF")
    end

    radial_build["hfs_gap_shield"] = gasc["INPUTS"]["radial build"]["gapBlanketCoil"] * norm
    radial_build["hfs_shield"] = gasc["INPUTS"]["radial build"]["rbInnerShield"] * norm
    if no_small_gaps
        radial_build["hfs_shield"] += radial_build["hfs_gap_shield"]
        pop!(radial_build, "hfs_gap_shield")
    end
    radial_build["hfs_blanket"] = gasc["INPUTS"]["radial build"]["rbInnerBlanket"] * norm

    radial_build["hfs_wall"] = gasc["INPUTS"]["radial build"]["gapInnerBlanketWall"] * norm
    radial_build["plasma"] = (gasc["INPUTS"]["radial build"]["majorRadius"] - sum(values(radial_build))) * 2
    radial_build["lfs_wall"] = gasc["INPUTS"]["radial build"]["gapOuterBlanketWall"] * norm

    radial_build["lfs_blanket"] = gasc["INPUTS"]["radial build"]["rbOuterBlanket"] * norm
    radial_build["lfs_shield"] = gasc["INPUTS"]["radial build"]["rbOuterShield"] * norm
    radial_build["lfs_gap_shield"] = gasc["INPUTS"]["radial build"]["gapBlanketCoil"] * norm
    if no_small_gaps
        radial_build["lfs_shield"] += radial_build["lfs_gap_shield"]
        pop!(radial_build, "lfs_gap_shield")
    end

    radial_build["lfs_TF"] = radial_build["hfs_TF"]
    radial_build["lfs_gap_TF"] = gasc["INPUTS"]["radial build"]["gapTFOH"] * norm
    if no_small_gaps
        radial_build["lfs_TF"] += radial_build["lfs_gap_TF"]
        pop!(radial_build, "lfs_gap_TF")
    end

    radial_build["gap_cryostat"] = radial_build["gap_OH"] * 3

    # thin layers can cause LibGEOS to crash
    min_fraction_thin_wall = 0.02
    if no_small_gaps && (radial_build["hfs_wall"] < min_fraction_thin_wall * norm)
        radial_build["hfs_blanket"] -= (min_fraction_thin_wall * norm - radial_build["hfs_wall"])
        radial_build["hfs_wall"] = min_fraction_thin_wall * norm
    end
    if no_small_gaps && (radial_build["lfs_wall"] < min_fraction_thin_wall * norm)
        radial_build["lfs_blanket"] -= (min_fraction_thin_wall * norm - radial_build["lfs_wall"])
        radial_build["lfs_wall"] = min_fraction_thin_wall * norm
    end

    init_build(dd.build, radial_build; verbose)

    # TF coils
    dd.build.tf.coils_n = 16

    # cross-section outlines
    build_cx(dd.build, dd.equilibrium.time_slice[], tf_shape_index)

    return dd
end


function init_build(dd::IMAS.dd, par::Parameters)
    init_from = par.general.init_from

    if init_from == :gasc
        gasc = GASC(par.gasc.filename, par.gasc.case)
        init_build(dd, gasc; par.gasc.no_small_gaps, tf_shape_index = 3)

    elseif init_from == :ods
        dd1 = IMAS.json2imas(par.ods.filename)
        if length(keys(dd1.build)) > 0
            dd.build = dd1.build
        else
            if length(keys(dd1.wall)) > 0
                dd.wall = dd1.wall
            end
            init_from = :scalars
        end
    end

    if init_from == :scalars
        init_build(dd;
            tf_shape_index = 3,
            is_nuclear_facility = par.build.is_nuclear_facility,
            pf_inside_tf = (par.pf_active.n_pf_coils_inside > 0),
            pf_outside_tf = (par.pf_active.n_pf_coils_outside > 0))
    end
    return dd
end

function wall_miller_conformal(bd, layer_type, elongation, triangularity; n_points = 101)
    if layer_type == -1
        Rstart = IMAS.get_build(bd, type = layer_type).start_radius
        Rend = IMAS.get_build(bd, type = layer_type).end_radius
        line = miller_Rstart_Rend(Rstart, Rend, elongation, triangularity; n_points)
        return line, line
    else
        Rstart_lfs = IMAS.get_build(bd, type = layer_type, hfs = -1).start_radius
        Rend_lfs = IMAS.get_build(bd, type = layer_type, hfs = -1).end_radius
        Rstart_hfs = IMAS.get_build(bd, type = layer_type, hfs = 1).start_radius
        Rend_hfs = IMAS.get_build(bd, type = layer_type, hfs = 1).end_radius
        inner_line = miller_Rstart_Rend(Rend_hfs, Rstart_lfs, elongation, triangularity; n_points)
        outer_line = miller_Rstart_Rend(Rstart_hfs, Rend_lfs, elongation, triangularity; n_points)
        return inner_line, outer_line
    end
end

"""
    build_cx(bd::IMAS.build, eqt::IMAS.equilibrium__time_slice, tf_shape_index::Int)

Translates 1D build to 2D cross-sections starting from equilibrium
"""
function build_cx(bd::IMAS.build, eqt::IMAS.equilibrium__time_slice, tf_shape_index::Int)
    # Inner radii of the plasma
    R_hfs_plasma = IMAS.get_build(bd, type = -1).start_radius
    R_lfs_plasma = IMAS.get_build(bd, type = -1).end_radius

    # Plasma as buffered convex-hull polygon of LCFS and strike points
    ψb = IMAS.find_psi_boundary(eqt)
    ψa = eqt.profiles_1d.psi[1]
    δψ = 0.10 # this sets the length of the strike divertor legs
    r_in, z_in, _ = IMAS.flux_surface(eqt, ψb * (1 - δψ) + ψa * δψ, true)
    Z0 = eqt.global_quantities.magnetic_axis.z
    rlcfs, zlcfs, _ = IMAS.flux_surface(eqt, ψb, true)
    theta = range(0.0, 2 * pi, length = 101)
    private_extrema = []
    private = IMAS.flux_surface(eqt, ψb, false)
    a = 0
    for (pr, pz) in private
        if sign(pz[1] - Z0) != sign(pz[end] - Z0)
            # open flux surface does not encicle the plasma
            continue
        elseif minimum_distance_two_shapes(pr, pz, rlcfs, zlcfs) > (maximum(zlcfs) - minimum(zlcfs)) / 20
            # secondary Xpoint far away
            continue
        elseif (sum(pz) - Z0) < 0
            # lower private region
            index = argmax(pz)
            a = minimum(z_in) - minimum(zlcfs)
            a = min(a, pz[index] - minimum(pz))
        else
            # upper private region
            index = argmin(pz)
            a = maximum(zlcfs) - maximum(z_in)
            a = min(a, maximum(pz) - pz[index])
        end
        Rx = pr[index]
        Zx = pz[index]
        append!(private_extrema, IMAS.intersection(a .* cos.(theta) .+ Rx, a .* sin.(theta) .+ Zx, pr, pz))
    end
    h = [[r, z] for (r, z) in vcat(collect(zip(rlcfs, zlcfs)), private_extrema)]
    hull = LazySets.convex_hull(h)
    R = [r for (r, z) in hull]
    R[R.<R_hfs_plasma] .= R_hfs_plasma
    R[R.>R_lfs_plasma] .= R_lfs_plasma
    Z = [z for (r, z) in hull]
    hull_poly = xy_polygon(R, Z)
    plasma_poly = LibGEOS.buffer(hull_poly, ((R_lfs_plasma - R_hfs_plasma) - (maximum(rlcfs) - minimum(rlcfs))) / 2.0)

    # make the divertor domes in the plasma
    δψ = 0.05 # how close to the LCFS shoudl the divertor plates be
    for (pr, pz) in IMAS.flux_surface(eqt, ψb * (1 - δψ) + ψa * δψ, false)
        if pr[1] != pr[end]
            pz[1] = pz[1] * 2
            pz[end] = pz[end] * 2
            plasma_poly = LibGEOS.difference(plasma_poly, xy_polygon(pr, pz))
        end
    end

    # plasma first wall
    pr = [v[1] for v in LibGEOS.coordinates(plasma_poly)[1]]
    pz = [v[2] for v in LibGEOS.coordinates(plasma_poly)[1]]

    # build cx
    return build_cx(bd, pr, pz, tf_shape_index)
end

"""
    build_cx(bd::IMAS.build, pr, pz, tf_shape_index::Int)

Translates 1D build to 2D cross-sections starting from R and Z coordinates of plasma first wall
"""
function build_cx(bd::IMAS.build, pr::Vector{Float64}, pz::Vector{Float64}, tf_shape_index::Int)
    # plasma
    IMAS.get_build(bd, type = -1).outline.r = pr
    IMAS.get_build(bd, type = -1).outline.z = pz

    # all layers between plasma and OH
    plasma_to_oh = []
    valid = false
    for (k, layer) in reverse(collect(enumerate(bd.layer)))
        # stop once you see the OH
        if layer.type == 1
            valid = false
            break
        end
        if valid
            push!(plasma_to_oh, k)
        end
        # valid starting from the plasma
        if layer.type == -1
            valid = true
        end
    end
    shape_set = false
    for (n, k) in enumerate(plasma_to_oh)
        # layer that preceeds the TF (or shield) sets the TF (and shield) shape
        if (!shape_set) && (n < length(plasma_to_oh)) && (bd.layer[plasma_to_oh[n+1]].type in [2, 3])
            FUSE.optimize_shape(bd, k, tf_shape_index)
            shape_set = true
            # everything else is conformal convex hull
        else
            FUSE.optimize_shape(bd, k, -2)
        end
    end

    # plug
    L = 0
    R = IMAS.get_build(bd, type = 1).start_radius
    D = minimum(IMAS.get_build(bd, type = 5, hfs = 1).outline.z)
    U = maximum(IMAS.get_build(bd, type = 5, hfs = 1).outline.z)
    bd.layer[1].outline.r, bd.layer[1].outline.z = rectangle_shape(L, R, D, U)

    # oh
    L = IMAS.get_build(bd, type = 1).start_radius
    R = IMAS.get_build(bd, type = 1).end_radius
    bd.layer[2].outline.r, bd.layer[2].outline.z = rectangle_shape(L, R, D, U)

    # cryostat
    L = 0
    R = bd.layer[end].end_radius
    D = minimum(IMAS.get_build(bd, type = 2, hfs = 1).outline.z) - bd.layer[end].thickness
    U = maximum(IMAS.get_build(bd, type = 2, hfs = 1).outline.z) + bd.layer[end].thickness
    bd.layer[end].outline.r, bd.layer[end].outline.z = rectangle_shape(L, R, D, U)

    # set the toroidal thickness of the TF coils based on the innermost radius and the number of coils
    bd.tf.thickness = 2 * π * IMAS.get_build(bd, type = 2, hfs = 1).start_radius / bd.tf.coils_n
    return bd
end

"""
    build_cx(bd::IMAS.build, pr, pz, tf_shape_index::Int)

Translates 1D build to 2D cross-sections starting from R and Z coordinates of plasma first wall
"""
function build_cx(bd::IMAS.build, wl::IMAS.wall__description_2d___limiter__unit___outline, tf_shape_index::Int)
    return build_cx(bd, wl.r, wl.z, tf_shape_index)
end

function optimize_shape(bd::IMAS.build, layer_index::Int, tf_shape_index::Int)
    # properties of current layer
    layer = bd.layer[layer_index]
    id = bd.layer[layer_index].identifier
    r_start = layer.start_radius
    r_end = IMAS.get_build(bd, identifier = layer.identifier, hfs = -1).end_radius
    hfs_thickness = layer.thickness
    lfs_thickness = IMAS.get_build(bd, identifier = id, hfs = -1).thickness
    target_minimum_distance = (hfs_thickness + lfs_thickness) / 2.0

    # obstruction
    oR = bd.layer[layer_index+1].outline.r
    oZ = bd.layer[layer_index+1].outline.z

    # only update shape if that's not been set before
    # this is to allow external overriding of default shape setting
    if ismissing(layer, :shape)
        layer.shape = tf_shape_index
    end

    # handle offset and offset & convex-hull
    if layer.shape in [-1, -2]
        poly = LibGEOS.buffer(xy_polygon(oR, oZ), (hfs_thickness + lfs_thickness) / 2.0)
        layer.outline.r = [v[1] .+ (lfs_thickness .- hfs_thickness) / 2.0 for v in LibGEOS.coordinates(poly)[1]]
        layer.outline.z = [v[2] for v in LibGEOS.coordinates(poly)[1]]
        if layer.shape == -2
            h = [[r, z] for (r, z) in collect(zip(layer.outline.r, layer.outline.z))]
            hull = LazySets.convex_hull(h)
            layer.outline.r = vcat([r for (r, z) in hull], hull[1][1])
            layer.outline.z = vcat([z for (r, z) in hull], hull[1][2])
        end
        # handle shapes
    else
        up_down_symmetric = false
        if abs(sum(oZ) / sum(abs.(oZ))) < 1E-2
            up_down_symmetric = true
        end

        if up_down_symmetric
            layer.shape = mod(layer.shape, 100)
        else
            layer.shape = mod(layer.shape, 100) + 100
        end

        func = shape_function(layer.shape)
        if ismissing(layer, :shape_parameters)
            layer.shape_parameters = init_shape_parameters(layer.shape, oR, oZ, r_start, r_end, target_minimum_distance)
        end
        layer.outline.r, layer.outline.z = func(r_start, r_end, layer.shape_parameters...)
        layer.shape_parameters = optimize_shape(oR, oZ, target_minimum_distance, func, r_start, r_end, layer.shape_parameters)
        layer.outline.r, layer.outline.z = func(r_start, r_end, layer.shape_parameters...)
    end
end

