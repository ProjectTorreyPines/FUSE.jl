import LibGEOS
import Interpolations
import Contour
import DataStructures

import IMAS: BuildLayerType, _plasma_, _gap_, _oh_, _tf_, _shield_, _blanket_, _wall_, _vessel_
import IMAS: BuildLayerSide, _lfs_, _lhfs_, _hfs_
import IMAS: BuildLayerShape, _offset_, _convex_hull_, _princeton_D_, _rectangle_, _triple_arc_, _miller_, _spline_

#= ========== =#
#  init build  #
#= ========== =#
function init_build(dd::IMAS.dd, par::Parameters)
    init_from = par.general.init_from

    if init_from == :gasc
        gasc = GASC(par.gasc.filename, par.gasc.case)
        init_radial_build(dd.build, gasc; no_small_gaps = par.gasc.no_small_gaps)

    elseif init_from == :ods
        dd1 = IMAS.json2imas(par.ods.filename)
        if length(keys(dd1.wall)) > 0
            dd.wall = dd1.wall
        end
        if length(keys(dd1.build)) > 0
            dd.build = dd1.build
        else
            init_from = :scalars
        end
    end

    if init_from == :scalars
        if ismissing(par.build, :layers)
            init_radial_build(
                dd.build,
                dd.equilibrium.time_slice[],
                first_wall(dd.wall);
                is_nuclear_facility = par.build.is_nuclear_facility,
                pf_inside_tf = (par.pf_active.n_pf_coils_inside > 0),
                pf_outside_tf = (par.pf_active.n_pf_coils_outside > 0))
        else
            init_radial_build(dd.build, par.build.layers; verbose = true)
        end
    end

    # cross-section outlines
    sh = Symbol("_$(par.tf.shape)_")
    build_cx(dd; tf_shape = @eval($sh))

    # TF coils
    dd.build.tf.coils_n = par.tf.n_coils
    # set the toroidal thickness of the TF coils based on the innermost radius and the number of coils
    dd.build.tf.thickness = 2 * π * IMAS.get_build(dd.build, type = _tf_, fs = _hfs_).start_radius / dd.build.tf.coils_n

    # assign materials
    assign_build_layers_materials(dd, par)

    return dd
end

"""
    init_build(bd::IMAS.build; verbose = false, layers...)

Initialize build IDS based on center stack layers (thicknesses)

NOTE: layer[:].type and layer[:].material follows from naming of layers
*   0 ...gap... : vacuum
*   1 OH: ohmic coil
*   2 TF: toroidal field coil
*   6 vessel...: vacuum vessel
*   3 shield...: neutron shield
*   4 blanket...: neutron blanket
*   5 wall....: first wall
*  -1 ...plasma...: 

layer[:].fs is set depending on if "hfs" or "lfs" appear in the name

layer[:].identifier is created as a hash of then name removing "hfs" or "lfs"
"""
function init_radial_build(bd::IMAS.build; verbose::Bool = false, layers...)
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
            layer.type = Int(_gap_)
        end
        if occursin("plasma", lowercase(layer.name))
            layer.type = Int(_plasma_)
        elseif uppercase(layer.name) == "OH"
            layer.type = Int(_oh_)
        elseif occursin("TF", uppercase(layer.name))
            layer.type = Int(_tf_)
        elseif occursin("shield", lowercase(layer.name))
            layer.type = Int(_shield_)
        elseif occursin("blanket", lowercase(layer.name))
            layer.type = Int(_blanket_)
        elseif occursin("wall", lowercase(layer.name))
            layer.type = Int(_wall_)
        elseif occursin("vessel", lowercase(layer.name))
            layer.type = Int(_vessel_)
        end
        if occursin("hfs", lowercase(layer.name))
            layer.fs = Int(_hfs_)
        elseif occursin("lfs", lowercase(layer.name))
            layer.fs = Int(_lfs_)
        else
            layer.fs = Int(_lhfs_)
        end

        layer.identifier = UInt(hash(replace(replace(lowercase(layer.name), "hfs" => ""), "lfs" => "")))
        radius_end += layer.thickness
        if verbose
            println("$(rpad(string(k),4)) $(rpad(layer.name, 17)) Δr=$(@sprintf("%3.3f",layer.thickness))    R=[$(@sprintf("%3.3f",radius_start)) <=> $(@sprintf("%3.3f",radius_end))]")
        end
        radius_start = radius_end
    end

    return bd
end

function init_radial_build(bd::IMAS.build, layers::AbstractDict; verbose = false)
    nt = (; zip([Symbol(k) for k in keys(layers)], values(layers))...)
    init_radial_build(bd; verbose, nt...)
end

"""
    init_radial_build(
        bd::IMAS.build,
        eqt=IMAS.equilibrium__time_slice;
        is_nuclear_facility::Bool = true,
        pf_inside_tf::Bool = false,
        pf_outside_tf::Bool = true,
        verbose::Bool = false)

Initialization of build IDS based on equilibrium time_slice
"""
function init_radial_build(
    bd::IMAS.build,
    eqt::IMAS.equilibrium__time_slice,
    wall::T where {T<:Union{IMAS.wall__description_2d___limiter__unit___outline,Missing}};
    is_nuclear_facility::Bool = true,
    pf_inside_tf::Bool = false,
    pf_outside_tf::Bool = true,
    verbose::Bool = false)

    if wall !== missing
        rmin = minimum(wall.r)
        rmax = maximum(wall.r)
    else
        rmin = eqt.boundary.geometric_axis.r - eqt.boundary.minor_radius
        rmax = eqt.boundary.geometric_axis.r + eqt.boundary.minor_radius
        gap = (rmax - rmin) / 20.0 # plasma-wall gap
        gap = (rmax - rmin) / 20.0 # plasma-wall gap
        rmin -= gap
        rmax += gap
    end

    if is_nuclear_facility
        n_hfs_layers = 6.125
        dr = rmin / n_hfs_layers
        init_radial_build(bd;
            verbose,
            gap_OH = dr * 2.0,
            OH = dr,
            hfs_TF = dr,
            gap_hfs_vacuum_vessel = dr / 8.0,
            gap_hfs_coils = pf_inside_tf ? 0 : -1,
            hfs_shield = dr / 2.0,
            hfs_blanket = dr,
            hfs_wall = dr / 2.0,
            plasma = rmax - rmin,
            lfs_wall = dr / 2.0,
            lfs_blanket = dr * 2,
            lfs_shield = dr / 2.0,
            gap_lfs_coils = dr * (pf_inside_tf ? 4 : -1),
            gap_lfs_vacuum_vessel = dr / 8.0,
            lfs_TF = dr,
            gap_cryostat = dr * (pf_outside_tf ? 5 : 1))

    else
        n_hfs_layers = 4.5
        dr = rmin / n_hfs_layers
        init_radial_build(bd;
            verbose,
            gap_OH = dr * 2.0,
            OH = dr,
            hfs_TF = dr,
            gap_hfs_coils = pf_inside_tf ? 0 : -1,
            hfs_wall = dr / 2.0,
            plasma = rmax - rmin,
            lfs_wall = dr / 2.0,
            gap_lfs_coils = dr * (pf_inside_tf ? 2.25 : -1),
            lfs_TF = dr,
            gap_cryostat = dr * (pf_outside_tf ? 3 : 1))
    end

    return bd
end

"""
    init_radial_build(
        bd::IMAS.build,
        gasc::GASC;
        no_small_gaps::Bool = true,
        verbose::Bool = false)

Initialization of radial build based on equilibrium GASC output
"""
function init_radial_build(
    bd::IMAS.build,
    gasc::GASC;
    no_small_gaps::Bool = true,
    vacuum_vessel::Float64 = 0.1,
    verbose::Bool = false)

    layers = gascrb2layers(gasc.solution["INPUTS"]["radial build"]; no_small_gaps, vacuum_vessel)

    init_radial_build(bd, layers; verbose)

    return bd
end

"""
    function gascrb2layers(
        gascrb::Dict;
        no_small_gaps::Bool = true,
        vacuum_vessel::Float64 = 0.1)

Convert GASC ["INPUTS"]["radial build"] to FUSE build layers dictionary
"""
function gascrb2layers(
    gascrb::Dict;
    no_small_gaps::Bool = true,
    vacuum_vessel::Float64 = 0.1)

    majorRadius = gascrb["majorRadius"]
    aspectRatio = gascrb["aspectRatio"]
    minorRadius = majorRadius / aspectRatio
    innerPlasmaRadius = majorRadius - minorRadius
    norm = innerPlasmaRadius

    layers = DataStructures.OrderedDict()
    for run in 1:2
        layers["OH"] = gascrb["rbOH"] * norm

        layers["hfs_gap_TF"] = gascrb["gapTFOH"] * norm
        layers["hfs_TF"] = gascrb["rbTF"] * norm
        if no_small_gaps
            layers["hfs_TF"] += layers["hfs_gap_TF"]
            pop!(layers, "hfs_gap_TF")
        end

        if vacuum_vessel > 0.0
            layers["gap_hfs_vacuum_vessel"] = gascrb["rbInnerBlanket"] * norm * vacuum_vessel
        end

        layers["hfs_gap_shield"] = gascrb["gapBlanketCoil"] * norm
        layers["hfs_shield"] = gascrb["rbInnerShield"] * norm
        if no_small_gaps
            layers["hfs_shield"] += layers["hfs_gap_shield"]
            pop!(layers, "hfs_gap_shield")
        end
        layers["hfs_blanket"] = gascrb["rbInnerBlanket"] * norm * (1 - vacuum_vessel)

        layers["hfs_wall"] = gascrb["gapInnerBlanketWall"] * norm

        if run == 1
            between_gapOH_and_plasma = sum(values(layers))
            empty!(layers)
            layers["gap_OH"] = innerPlasmaRadius - between_gapOH_and_plasma
        end
    end

    layers["plasma"] = (gascrb["majorRadius"] - sum(values(layers))) * 2
    layers["lfs_wall"] = gascrb["gapOuterBlanketWall"] * norm

    layers["lfs_blanket"] = gascrb["rbOuterBlanket"] * norm * (1 - vacuum_vessel)
    layers["lfs_shield"] = gascrb["rbOuterShield"] * norm
    layers["lfs_gap_shield"] = gascrb["gapBlanketCoil"] * norm
    if no_small_gaps
        layers["lfs_shield"] += layers["lfs_gap_shield"]
        pop!(layers, "lfs_gap_shield")
    end

    if vacuum_vessel > 0.0
        layers["gap_lfs_vacuum_vessel"] = gascrb["rbOuterBlanket"] * norm * vacuum_vessel
    end

    layers["lfs_TF"] = layers["hfs_TF"]
    layers["lfs_gap_TF"] = gascrb["gapTFOH"] * norm
    if no_small_gaps
        layers["lfs_TF"] += layers["lfs_gap_TF"]
        pop!(layers, "lfs_gap_TF")
    end

    layers["gap_cryostat"] = layers["gap_OH"] * 3

    # thin layers can cause LibGEOS to crash
    # if wall is too thin, then thicken it at the expense of the blanket
    min_fraction_thin_wall = 0.02
    if no_small_gaps && (layers["hfs_wall"] < min_fraction_thin_wall * norm)
        layers["hfs_blanket"] -= (min_fraction_thin_wall * norm - layers["hfs_wall"])
        layers["hfs_wall"] = min_fraction_thin_wall * norm
    end
    if no_small_gaps && (layers["lfs_wall"] < min_fraction_thin_wall * norm)
        layers["lfs_blanket"] -= (min_fraction_thin_wall * norm - layers["lfs_wall"])
        layers["lfs_wall"] = min_fraction_thin_wall * norm
    end

    return layers
end

"""
    first_wall(wall::IMAS.wall)

return outline of first wall
"""
function first_wall(wall::IMAS.wall)
    if (!ismissing(wall.description_2d, [1, :limiter, :unit, 1, :outline, :r])) && (length(wall.description_2d[1].limiter.unit[1].outline.r) > 5)
        return wall.description_2d[1].limiter.unit[1].outline
    else
        return missing
    end
end

"""
    wall_from_eq(dd)

Generate first wall outline starting from an equilibrium
"""
function wall_from_eq(bd::IMAS.build, eqt::IMAS.equilibrium__time_slice; divertor_length_length_multiplier::Real = 1.0)
    # Inner radii of the plasma
    R_hfs_plasma = IMAS.get_build(bd, type = _plasma_).start_radius
    R_lfs_plasma = IMAS.get_build(bd, type = _plasma_).end_radius

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
        a *= divertor_length_length_multiplier
        cr = a .* cos.(theta) .+ Rx
        cz = a .* sin.(theta) .+ Zx
        append!(private_extrema, IMAS.intersection(cr, cz, pr, pz))
    end
    h = [[r, z] for (r, z) in vcat(collect(zip(rlcfs, zlcfs)), private_extrema)]
    hull = convex_hull(h)
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

    return pr, pz
end

"""
    build_cx(dd::IMAS.dd; tf_shape_index::Int)

Translates 1D build to 2D cross-sections starting either wall information
If wall information is missing, then the first wall information is generated starting from equilibrium time_slice
"""
function build_cx(dd::IMAS.dd; tf_shape::BuildLayerShape)
    wall = first_wall(dd.wall)
    if wall === missing
        pr, pz = wall_from_eq(dd.build, dd.equilibrium.time_slice[])
        resize!(dd.wall.description_2d, 1)
        resize!(dd.wall.description_2d[1].limiter.unit, 1)
        dd.wall.description_2d[1].limiter.unit[1].outline.r = pr
        dd.wall.description_2d[1].limiter.unit[1].outline.z = pz
        wall = first_wall(dd.wall)
    end
    build_cx(dd.build, wall.r, wall.z, tf_shape)
end

"""
    build_cx(bd::IMAS.build, pr::Vector{Float64}, pz::Vector{Float64}, tf_shape_index::Int)

Translates 1D build to 2D cross-sections starting from R and Z coordinates of plasma first wall
"""
function build_cx(bd::IMAS.build, pr::Vector{Float64}, pz::Vector{Float64}, tf_shape::BuildLayerShape)
    # plasma
    IMAS.get_build(bd, type = _plasma_).outline.r = pr
    IMAS.get_build(bd, type = _plasma_).outline.z = pz

    # all layers between plasma and OH
    plasma_to_oh = []
    valid = false
    for (k, layer) in reverse(collect(enumerate(bd.layer)))
        # stop once you see the OH
        if layer.type == Int(_oh_)
            valid = false
            break
        end
        if valid
            push!(plasma_to_oh, k)
        end
        # valid starting from the plasma
        if layer.type == Int(_plasma_)
            valid = true
        end
    end
    shape_set = false
    for (n, k) in enumerate(plasma_to_oh)
        # layer that preceeds the TF (or shield) sets the TF (and shield) shape
        if (!shape_set) && (n < length(plasma_to_oh)) && (bd.layer[plasma_to_oh[n+1]].type in [Int(_tf_), Int(_shield_)])
            FUSE.optimize_shape(bd, k, tf_shape)
            shape_set = true
            # everything else is conformal convex hull
        else
            FUSE.optimize_shape(bd, k, _convex_hull_)
        end
    end

    # plug
    L = 0
    R = IMAS.get_build(bd, type = _oh_).start_radius
    D = minimum(IMAS.get_build(bd, type = _wall_, fs = _hfs_).outline.z)
    U = maximum(IMAS.get_build(bd, type = _wall_, fs = _hfs_).outline.z)
    bd.layer[1].outline.r, bd.layer[1].outline.z = rectangle_shape(L, R, D, U)

    # oh
    L = IMAS.get_build(bd, type = _oh_).start_radius
    R = IMAS.get_build(bd, type = _oh_).end_radius
    bd.layer[2].outline.r, bd.layer[2].outline.z = rectangle_shape(L, R, D, U)

    # cryostat
    L = 0
    R = bd.layer[end].end_radius
    D = minimum(IMAS.get_build(bd, type = _tf_, fs = _hfs_).outline.z) - bd.layer[end].thickness
    U = maximum(IMAS.get_build(bd, type = _tf_, fs = _hfs_).outline.z) + bd.layer[end].thickness
    bd.layer[end].outline.r, bd.layer[end].outline.z = rectangle_shape(L, R, D, U)

    return bd
end

"""
    optimize_shape(bd::IMAS.build, layer_index::Int, tf_shape::BuildLayerShape)

Generates outline of layer in such a way to maintain minimum distance from inner layer
"""
function optimize_shape(bd::IMAS.build, layer_index::Int, tf_shape::BuildLayerShape)
    # properties of current layer
    layer = bd.layer[layer_index]
    id = bd.layer[layer_index].identifier
    r_start = layer.start_radius
    r_end = IMAS.get_build(bd, identifier = layer.identifier, fs = _lfs_).end_radius
    hfs_thickness = layer.thickness
    lfs_thickness = IMAS.get_build(bd, identifier = id, fs = _lfs_).thickness
    target_minimum_distance = (hfs_thickness + lfs_thickness) / 2.0

    # obstruction
    oR = bd.layer[layer_index+1].outline.r
    oZ = bd.layer[layer_index+1].outline.z

    # only update shape if that is not been set before
    # this is to allow external overriding of default shape setting
    if ismissing(layer, :shape)
        layer.shape = Int(tf_shape)
    end

    # handle offset and offset & convex-hull
    if layer.shape in [-1, -2]
        poly = LibGEOS.buffer(xy_polygon(oR, oZ), (hfs_thickness + lfs_thickness) / 2.0)
        layer.outline.r = [v[1] .+ (lfs_thickness .- hfs_thickness) / 2.0 for v in LibGEOS.coordinates(poly)[1]]
        layer.outline.z = [v[2] for v in LibGEOS.coordinates(poly)[1]]
        if layer.shape == -2
            h = [[r, z] for (r, z) in collect(zip(layer.outline.r, layer.outline.z))]
            hull = convex_hull(h)
            layer.outline.r, layer.outline.z = IMAS.resample_2d_line(vcat([r for (r, z) in hull], hull[1][1]), vcat([z for (r, z) in hull], hull[1][2]))
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

function assign_build_layers_materials(dd::IMAS.dd, par::Parameters)
    bd = dd.build
    for layer in bd.layer
        if layer.type == Int(_plasma_)
            layer.material = par.build.is_nuclear_facility ? "DT_plasma" : "DD_plasma"
        elseif layer.type == Int(_gap_)
            layer.material = "Vacuum"
        elseif layer.type == Int(_oh_)
            layer.material = par.material.wall
        elseif layer.type == Int(_tf_)
            layer.material = par.tf.technology.material
        elseif layer.type == Int(_shield_)
            layer.material = par.material.shield
        elseif layer.type == Int(_blanket_)
            layer.material = par.material.blanket
        elseif layer.type == Int(_wall_)
            layer.material = par.material.wall
        elseif layer.type == Int(_vessel_)
            layer.material = layer.material = "Water, Liquid"
        end
    end
end
