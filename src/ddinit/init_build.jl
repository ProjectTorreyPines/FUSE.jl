import LibGEOS
import GeoInterface
import Interpolations
import OrderedCollections

import IMAS: BuildLayerType, _plasma_, _gap_, _oh_, _tf_, _shield_, _blanket_, _wall_, _vessel_, _cryostat_, _divertor_
import IMAS: BuildLayerSide, _lfs_, _lhfs_, _hfs_, _in_, _out_
import IMAS: BuildLayerShape, _offset_, _negative_offset_, _convex_hull_, _princeton_D_exact_, _princeton_D_, _princeton_D_scaled_, _rectangle_, _double_ellipse_, _triple_arc_,
    _miller_, _square_miller_, _spline_, _silo_

#= ========================================== =#
#  Visualization of IMAS.build.layer as table  #
#= ========================================== =#
function DataFrames.DataFrame(layers::IMAS.IDSvector{<:IMAS.build__layer})

    df = DataFrames.DataFrame(;
        group=String[],
        details=String[],
        type=String[],
        ΔR=Float64[],
        R_start=Float64[],
        R_end=Float64[],
        material=String[],
        area=Float64[],
        volume=Float64[]
    )

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
    init_build!(dd::IMAS.dd, ini::ParametersAllInits, act::ParametersAllActors, dd1::IMAS.dd)

Initialize `dd.build` starting from `ini` and `act` parameters
"""
function init_build!(dd::IMAS.dd, ini::ParametersAllInits, act::ParametersAllActors, dd1::IMAS.dd)
    TimerOutputs.reset_timer!("init_build")
    TimerOutputs.@timeit timer "init_build" begin
        init_from = ini.general.init_from

        if init_from == :ods
            if length(keys(dd1.wall)) > 0
                dd.wall = dd1.wall
            end
            if length(dd1.build.layer) > 0
                dd.build = dd1.build
            else
                init_from = :scalars
            end
        end

        eqt = dd.equilibrium.time_slice[]

        # if layers are not filled explicitly, then generate them from fractions in ini.build
        if ismissing(ini.build, :layers)
            layers = layers_meters_from_fractions(
                eqt,
                IMAS.first_wall(dd.wall);
                shield=ini.build.shield,
                blanket=ini.build.blanket,
                vessel=ini.build.vessel,
                plasma_gap=ini.build.plasma_gap,
                pf_inside_tf=(ini.pf_active.n_coils_inside > 0),
                pf_outside_tf=(ini.pf_active.n_coils_outside > 0))
        else
            layers = ini.build.layers
            # scale radial build layers based on equilibrium R0, a, and the requested plasma_gap
            scale_build_layers!(layers, dd.equilibrium.vacuum_toroidal_field.r0, eqt.boundary.minor_radius, ini.build.plasma_gap)
        end

        # populate dd.build with radial build layers
        init_build!(dd.build, layers)

        # divertors
        if ini.build.divertors == :from_x_points
            if ismissing(ini.equilibrium, :xpoints)
                # if number of divertors is not set explicitly, get it from the ODS
                n_divertors = length(eqt.boundary.x_point)
                if n_divertors == 0
                    ini.equilibrium.xpoints = :none
                elseif n_divertors == 1
                    if eqt.boundary.x_point[1].z > eqt.boundary.geometric_axis.z
                        ini.equilibrium.xpoints = :upper
                    else
                        ini.equilibrium.xpoints = :lower
                    end
                elseif n_divertors == 2
                    ini.equilibrium.xpoints = :double
                end
            end
            ini.build.divertors = ini.equilibrium.xpoints
        end
        if ini.build.divertors ∈ (:upper, :double)
            dd.build.divertors.upper.installed = 1
        else
            dd.build.divertors.upper.installed = 0
        end
        if ini.build.divertors ∈ (:lower, :double)
            dd.build.divertors.lower.installed = 1
        else
            dd.build.divertors.lower.installed = 0
        end

        # set the TF shape
        tf_to_plasma = IMAS.get_build_indexes(dd.build.layer; fs=_hfs_)
        plama_to_tf = collect(reverse(tf_to_plasma))
        # set all shapes to convex hull by default
        for k in tf_to_plasma
            dd.build.layer[k].shape = Int(_convex_hull_)
        end
        # set TF (shape is set by inner)
        dd.build.layer[tf_to_plasma[2]].shape = Int(ini.tf.shape)
        # first layer is a offset
        k = plama_to_tf[1]
        if (dd.build.layer[k].type == Int(_wall_)) && ((dd.build.layer[k-1].type == Int(_blanket_)) || (dd.build.layer[k-1].type == Int(_shield_)))
            dd.build.layer[k].shape = Int(_offset_)
        end

        if ini.build.n_first_wall_conformal_layers >= 0
            # for k in plama_to_tf[ini.build.n_first_wall_conformal_layers:end-1]
            #     dd.build.layer[k+1].shape = Int(_offset_)
            # end
            dd.build.layer[plama_to_tf[ini.build.n_first_wall_conformal_layers]].shape = Int(ini.tf.shape)
        end

        # if ini.build.n_first_wall_conformal_layers >= 0
        #     dd.build.layer[tf_to_plasma[1]].shape = Int(_offset_)
        #     dd.build.layer[tf_to_plasma[2]].shape = Int(_offset_)
        #     dd.build.layer[plama_to_tf[ini.build.n_first_wall_conformal_layers]].shape = Int(ini.tf.shape)
        # end

        # 2D build cross-section
        ActorCXbuild(dd, act)

        # number of TF coils
        dd.build.tf.coils_n = ini.tf.n_coils
        # target TF ripple
        dd.build.tf.ripple = ini.tf.ripple

        # center stack solid mechanics
        dd.solid_mechanics.center_stack.bucked = Int(ini.center_stack.bucked)
        dd.solid_mechanics.center_stack.noslip = Int(ini.center_stack.noslip)
        dd.solid_mechanics.center_stack.plug = Int(ini.center_stack.plug)

        # assign coils and CS technologies
        assign_technologies(dd, ini)

        # assign radial build materials
        assign_build_layers_materials(dd, ini)

        return dd
    end
end

"""
    init_build!(bd::IMAS.build; layers...)

Initialize build IDS based on center stack layers (thicknesses)

NOTE: layer[:].type and layer[:].material follows from naming of layers

  - 0 ...gap... : vacuum
  - 1 OH: ohmic coil
  - 2 TF: toroidal field coil
  - 6 vessel...: vacuum vessel
  - 3 shield...: neutron shield
  - 4 blanket...: neutron blanket
  - 5 wall....: first wall
  - -1 ...plasma...:

layer[:].fs is set depending on if "hfs" or "lfs" appear in the name

layer[:].identifier is handled via IMAS.jl expressions
"""
function init_build!(bd::IMAS.build; layers...)
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
        layer.name = replace(string(layer_name), "_" => " ")
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
        radius_end += layer.thickness
        radius_start = radius_end
    end

    return bd
end

function init_build!(bd::IMAS.build, layers::AbstractDict)
    nt = (; zip([Symbol(k) for k in keys(layers)], values(layers))...)
    return init_build!(bd; nt...)
end

"""
    function init_build!(
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
    layers = OrderedCollections.OrderedDict{Symbol,Float64}()
    layers[:gap_OH] = 2.0
    layers[:OH] = 1.0
    layers[:hfs_TF] = 1.0
    if vessel > 0.0
        layers[:hfs_vacuum_vessel_wall_outer] = vessel * 0.1
        layers[:gap_hfs_vacuum_vessel] = vessel * 0.8
        layers[:hfs_vacuum_vessel_wall_inner] = vessel * 0.1
    end
    layers[:gap_hfs_coils] = pf_inside_tf ? 0 : -1
    if shield > 0.0
        layers[:hfs_shield] = shield
    end
    if blanket > 0.0
        layers[:hfs_blanket] = blanket
    end
    layers[:hfs_wall] = 0.5
    layers[:plasma] = 0.0 # this number does not matter
    layers[:lfs_wall] = 0.5
    if blanket > 0.0
        layers[:lfs_blanket] = blanket * 2.0
    end
    if shield > 0.0
        layers[:lfs_shield] = shield
    end
    layers[:gap_lfs_coils] = 1.0 * (pf_inside_tf ? 2.25 : -1)
    if vessel > 0.0
        layers[:lfs_vacuum_vessel_wall_inner] = vessel * 0.1
        layers[:gap_lfs_vacuum_vessel] = vessel * 0.8
        layers[:lfs_vacuum_vessel_wall_outer] = vessel * 0.1
    end
    layers[:lfs_TF] = 1.0
    layers[:gap_cryostat] = 1.0 * (pf_outside_tf ? 3 : 1)
    if blanket > 0.0
        layers[:cryostat] = 0.5
    end

    # from fractions to meters
    scale_build_layers!(layers, (rmax + rmin) / 2.0, (rmax - rmin) / 2.0, 0.0)

    return layers
end

function scale_build_layers!(layers::OrderedCollections.OrderedDict{Symbol,Float64}, R0::Float64, a::Float64, gap_fraction::Float64)
    gap = a * gap_fraction
    plasma_start = R0 - a - gap
    layer_plasma_start = 0.0
    for (layer, thickness) in layers
        if layer == :plasma
            break
        end
        if thickness > 0.0
            layer_plasma_start += thickness
        end
    end
    factor = plasma_start / layer_plasma_start
    for (layer, thickness) in layers
        if layer == :plasma
            layers[layer] = 2.0 * (a + gap)
        else
            layers[layer] = thickness * factor
        end
    end
end

function assign_build_layers_materials(dd::IMAS.dd, ini::ParametersAllInits)
    bd = dd.build
    for (k, layer) in enumerate(bd.layer)
        if k == 1 && ini.center_stack.plug
            layer.material = ini.material.wall
        elseif layer.type == Int(_plasma_)
            layer.material = any((layer.type in (Int(_blanket_), Int(_shield_)) for layer in bd.layer)) ? "DT_plasma" : "DD_plasma"
        elseif layer.type == Int(_gap_)
            layer.material = "Vacuum"
        elseif layer.type == Int(_oh_)
            layer.material = bd.oh.technology.material
        elseif layer.type == Int(_tf_)
            layer.material = bd.tf.technology.material
        elseif layer.type == Int(_shield_)
            layer.material = ini.material.shield
        elseif layer.type == Int(_blanket_)
            layer.material = ini.material.blanket
        elseif layer.type == Int(_wall_)
            layer.material = ini.material.wall
        elseif layer.type == Int(_vessel_)
            layer.material = "Water, Liquid"
        elseif layer.type == Int(_cryostat_)
            layer.material = ini.material.wall
        end
    end
end

function assign_technologies(dd::IMAS.dd, ini::ParametersAllInits)
    # plug
    if ini.center_stack.plug
        mechanical_technology(dd, :pl)
    end

    # oh
    coil_technology(dd.build.oh.technology, ini.oh.technology, :oh)
    mechanical_technology(dd, :oh)

    # tf
    coil_technology(dd.build.tf.technology, ini.tf.technology, :tf)
    mechanical_technology(dd, :tf)

    #pf
    return coil_technology(dd.build.pf_active.technology, ini.pf_active.technology, :pf_active)
end

function mechanical_technology(dd::IMAS.dd, what::Symbol)
    if what != :pl && getproperty(dd.build, what).technology.material == "Copper"
        material = pure_copper
    else
        material = stainless_steel
    end
    setproperty!(dd.solid_mechanics.center_stack.properties.yield_strength, what, material.yield_strength)
    setproperty!(dd.solid_mechanics.center_stack.properties.poisson_ratio, what, material.poisson_ratio)
    return setproperty!(dd.solid_mechanics.center_stack.properties.young_modulus, what, material.young_modulus)
end