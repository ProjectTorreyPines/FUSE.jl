import LibGEOS
import GeoInterface
import Interpolations
import DataStructures

import IMAS: BuildLayerType, _plasma_, _gap_, _oh_, _tf_, _shield_, _blanket_, _wall_, _vessel_, _cryostat_, _divertor_
import IMAS: BuildLayerSide, _lfs_, _lhfs_, _hfs_, _in_, _out_
import IMAS: BuildLayerShape, _offset_, _negative_offset_, _convex_hull_, _princeton_D_exact_, _princeton_D_, _princeton_D_scaled_, _rectangle_, _triple_arc_, _miller_, _square_miller_, _spline_, _silo_

#= ========================================== =#
#  Visualization of IMAS.build.layer as table  #
#= ========================================== =#
function DataFrames.DataFrame(layers::IMAS.IDSvector{<:IMAS.build__layer})
    df = DataFrames.DataFrame(group=String[], details=String[], type=String[], ΔR=Float64[], R_start=Float64[], R_end=Float64[], material=String[], area=Float64[], volume=Float64[])
    for layer in layers
        group = replace(string(BuildLayerSide(layer.fs)), "_" => "")
        type = replace(string(BuildLayerType(layer.type)), "_" => "")
        type = replace(type, r"^gap" => "")
        details = replace(lowercase(layer.name), r"^[hl]fs " => "")
        details = replace(details, r"^gap .*" => "")
        details = replace(details, r"\b" * type * r"\b" => "")
        material = getproperty(layer, :material, "?")
        material = split(material, ",")[1]
        material = replace(material, "Vacuum" => "")
        area = getproperty(layer, :area, NaN)
        volume = getproperty(layer, :volume, NaN)
        push!(df, [group, details, type, layer.thickness, layer.start_radius, layer.end_radius, material, area, volume])
    end
    return df
end

function Base.show(io::IO, ::MIME"text/plain", layers::IMAS.IDSvector{<:IMAS.build__layer})
    old_lines = get(ENV, "LINES", missing)
    old_columns = get(ENV, "COLUMNS", missing)
    df = DataFrames.DataFrame(layers)
    try
        ENV["LINES"] = 1000
        ENV["COLUMNS"] = 1000
        return show(io::IO, df)
    finally
        if old_lines === missing
            delete!(ENV, "LINES")
        else
            ENV["LINES"] = old_lines
        end
        if old_columns === missing
            delete!(ENV, "COLUMNS")
        else
            ENV["COLUMNS"] = old_columns
        end
    end
end

#= ========== =#
#  init build  #
#= ========== =#
"""
    init_build(dd::IMAS.dd, ini::ParametersAllInits, act::ParametersAllActors)

Initialize `dd.build` starting from `ini` and `act` parameters
"""
function init_build(dd::IMAS.dd, ini::ParametersAllInits, act::ParametersAllActors)
    init_from = ini.general.init_from

    if init_from == :ods
        dd1 = IMAS.json2imas(ini.ods.filename)
        if length(keys(dd1.wall)) > 0
            dd.wall = dd1.wall
        end
        if length(dd1.build.layer) > 0
            dd.build = dd1.build
        else
            init_from = :scalars
        end
    end

    # if layers are not filled explicitly, then generate them from fractions in ini.build.
    if ismissing(ini.build, :layers)
        layers = layers_meters_from_fractions(
            dd.equilibrium.time_slice[],
            IMAS.first_wall(dd.wall);
            shield=ini.build.shield,
            blanket=ini.build.blanket,
            vessel=ini.build.vessel,
            plasma_gap=ini.build.plasma_gap,
            pf_inside_tf=(ini.pf_active.n_pf_coils_inside > 0),
            pf_outside_tf=(ini.pf_active.n_pf_coils_outside > 0))
    else
        layers = deepcopy(ini.build.layers)
    end

    # scale radial build layers based on equilibrium R0 and ϵ
    iplasma = findfirst(key -> key in ["plasma", :plasma], collect(keys(layers)))
    R_hfs_build = sum([d for (k, d) in enumerate(values(layers)) if k < iplasma])
    if ismissing(ini.equilibrium, :R0)
        R0 = R_hfs_build + collect(values(layers))[iplasma] / 2.0
    else
        R0 = ini.equilibrium.R0
    end
    if ismissing(ini.equilibrium, :ϵ)
        a_gap = collect(values(layers))[iplasma] / 2.0
    else
        a_gap = ini.equilibrium.ϵ * R0 * (1.0 + ini.build.plasma_gap)
    end
    R_hfs = R0 - a_gap
    for key in keys(layers)
        if key in ["plasma", :plasma]
            layers[key] = a_gap * 2.0
        else
            layers[key] = layers[key] * R_hfs / R_hfs_build
        end
    end

    # populate dd.build with radial build layers
    init_build(dd.build, layers)

    # set the TF shape
    tf_to_plasma = IMAS.get_build(dd.build, fs=_hfs_, return_only_one=false, return_index=true)
    dd.build.layer[tf_to_plasma[1]].shape = Int(_offset_)
    dd.build.layer[tf_to_plasma[2]].shape = Int(to_enum(ini.tf.shape))
    # set all other shapes
    for k in tf_to_plasma[2:end]
        dd.build.layer[k+1].shape = Int(_convex_hull_)
    end
    k = tf_to_plasma[end]
    if (dd.build.layer[k].type == Int(_wall_)) && ((dd.build.layer[k-1].type == Int(_blanket_)) || (dd.build.layer[k-1].type == Int(_shield_)))
        dd.build.layer[k].shape = Int(_offset_)
    end
    if ini.build.n_first_wall_conformal_layers >= 0
        for k in tf_to_plasma[2:end-ini.build.n_first_wall_conformal_layers]
            dd.build.layer[k+1].shape = Int(_negative_offset_)
        end
    end

    # 2D build cross-section
    ActorCXbuild(dd, act)

    # TF coils
    dd.build.tf.coils_n = ini.tf.n_coils
    # set the toroidal thickness of the TF coils based on the innermost radius and the number of coils
    dd.build.tf.wedge_thickness = 2 * π * IMAS.get_build(dd.build, type=_tf_, fs=_hfs_).start_radius / dd.build.tf.coils_n
    # ripple
    dd.build.tf.ripple = ini.tf.ripple

    # center stack solid mechanics
    dd.solid_mechanics.center_stack.bucked = Int(ini.center_stack.bucked)
    dd.solid_mechanics.center_stack.noslip = Int(ini.center_stack.noslip)
    dd.solid_mechanics.center_stack.plug = Int(ini.center_stack.plug)

    # assign materials
    assign_build_layers_materials(dd, ini)

    return dd
end

"""
    init_build(bd::IMAS.build; layers...)

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
function init_build(bd::IMAS.build; layers...)
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
        if occursin("gap ", lowercase(layer.name))
            layer.type = Int(_gap_)
        elseif lowercase(layer.name) == "plasma"
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
        elseif occursin("cryostat", lowercase(layer.name))
            layer.type = Int(_cryostat_)
        end
        if occursin("hfs", lowercase(layer.name))
            layer.fs = Int(_hfs_)
        elseif occursin("lfs", lowercase(layer.name))
            layer.fs = Int(_lfs_)
        else
            if layer.type == Int(_plasma_)
                layer.fs = Int(_lhfs_)
            elseif k < length(layers) / 2
                layer.fs = Int(_in_)
            elseif k > length(layers) / 2
                layer.fs = Int(_out_)
            end
        end
        layer.identifier = UInt(hash(replace(replace(lowercase(layer.name), "hfs" => ""), "lfs" => ""))) % Int
        radius_end += layer.thickness
        radius_start = radius_end
    end

    return bd
end

function init_build(bd::IMAS.build, layers::AbstractDict)
    nt = (; zip([Symbol(k) for k in keys(layers)], values(layers))...)
    init_build(bd; nt...)
end

"""
    function init_build(
        bd::IMAS.build,
        eqt::IMAS.equilibrium__time_slice,
        wall::T where {T<:Union{IMAS.wall__description_2d___limiter__unit___outline,Missing}};
        blanket::Float64 = 1.0,
        shield::Float64 = 0.5,
        vessel::Float64 = .125,
        pf_inside_tf::Bool = false,
        pf_outside_tf::Bool = true)

Initialization of build IDS based on equilibrium time_slice
"""
function layers_meters_from_fractions(
    eqt::IMAS.equilibrium__time_slice,
    wall::T where {T<:Union{IMAS.wall__description_2d___limiter__unit___outline,Missing}};
    blanket::Float64=1.0,
    shield::Float64=0.5,
    vessel::Float64=0.125,
    plasma_gap::Float64=0.1,
    pf_inside_tf::Bool=false,
    pf_outside_tf::Bool=true)

    if wall !== missing
        rmin = minimum(wall.r)
        rmax = maximum(wall.r)
    else
        rmin = eqt.boundary.geometric_axis.r - eqt.boundary.minor_radius
        rmax = eqt.boundary.geometric_axis.r + eqt.boundary.minor_radius
        gap = (rmax - rmin) / 2.0 * plasma_gap
        rmin -= gap
        rmax += gap
    end

    # express layer thicknesses as fractions
    layers = DataStructures.OrderedDict()
    layers[:gap_OH] = 2.0
    layers[:OH] = 1.0
    layers[:hfs_TF] = 1.0
    if vessel > 0.0
        layers[:gap_hfs_vacuum_vessel] = vessel
    end
    layers[:gap_hfs_coils] = pf_inside_tf ? 0 : -1
    if shield > 0.0
        layers[:hfs_shield] = shield
    end
    if blanket > 0.0
        layers[:hfs_blanket] = blanket
    end
    layers[:hfs_wall] = 0.5
    layers[:plasma] = rmax - rmin
    layers[:lfs_wall] = 0.5
    if blanket > 0.0
        layers[:lfs_blanket] = blanket * 2.0
    end
    if shield > 0.0
        layers[:lfs_shield] = shield
    end
    layers[:gap_lfs_coils] = 1.0 * (pf_inside_tf ? 2.25 : -1)
    if vessel > 0.0
        layers[:gap_lfs_vacuum_vessel] = vessel
    end
    layers[:lfs_TF] = 1.0
    layers[:gap_cryostat] = 1.0 * (pf_outside_tf ? 3 : 1)
    if blanket > 0.0
        layers[:cryostat] = 0.5
    end

    # from fractions to meters
    n_hfs_layers = 0.0
    for layer in keys(layers)
        if layer == :plasma
            break
        end
        if layers[layer] > 0
            n_hfs_layers += layers[layer]
        end
    end
    dr = rmin / n_hfs_layers
    for layer in keys(layers)
        if layer != :plasma
            layers[layer] = layers[layer] * dr
        end
    end

    return layers
end
