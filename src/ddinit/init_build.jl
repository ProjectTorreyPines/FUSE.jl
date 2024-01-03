import LibGEOS
import GeoInterface
import Interpolations
import OrderedCollections

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
            if !isempty(dd1.wall.description_2d)
                dd.wall = dd1.wall
            end
            if length(dd1.build.layer) > 0
                dd.build = dd1.build
            else
                init_from = :scalars
            end
        end

        eqt = dd.equilibrium.time_slice[]

        # scale radial build layers based on equilibrium R0, a, and the requested plasma_gap
        scale_build_layers!(ini.build.layers, dd.equilibrium.vacuum_toroidal_field.r0, eqt.boundary.minor_radius, ini.build.plasma_gap)

        # populate dd.build with radial build layers
        init_build!(dd.build, ini.build.layers)

        # divertors
        if ini.build.divertors == :from_x_points
            if ismissing(ini.equilibrium, :xpoints)
                # if number of x-points is not set explicitly, get it from the pulse_schedule
                upper = any(x_point.z .> eqt.boundary.geometric_axis.z for x_point in dd.pulse_schedule.position_control.x_point)
                lower = any(x_point.z .< eqt.boundary.geometric_axis.z for x_point in dd.pulse_schedule.position_control.x_point)
                if !upper && !lower
                    ini.equilibrium.xpoints = :none
                elseif upper && !lower
                    ini.equilibrium.xpoints = :upper
                elseif lower && upper
                    ini.equilibrium.xpoints = :lower
                else
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
    init_build!(bd::IMAS.build, layers::ParametersVector{<:FUSEparameters__build_layer})

Initialize dd.build.layers from ini.build.layers
"""
function init_build!(bd::IMAS.build, layers::ParametersVector{<:FUSEparameters__build_layer})
    empty!(bd.layer)

    k = 0
    for ini_layer in layers
        if ini_layer.thickness < 0.0
            continue
        end
        k += 1
        layer = resize!(bd.layer, k)[k]

        layer.name = ini_layer.name
        layer.thickness = ini_layer.thickness
        layer.type = Int(ini_layer.type)
        layer.side = Int(ini_layer.side)
        if !ismissing(ini_layer, :material)
            layer.material = ini_layer.material
        end
        if !ismissing(ini_layer, :shape)
            layer.shape = Int(ini_layer.shape)
        end
    end

    return bd
end

"""
    setproperty!(parameters_build::FUSEparameters__build{T}, field::Symbol, layers::AbstractDict{Symbol,<:Real}) where {T<:Real}

Allows users to initialize layers from a dictionary
"""
function Base.setproperty!(parameters_build::FUSEparameters__build{T}, field::Symbol, layers::AbstractDict{Symbol,<:Real}) where {T<:Real}
    @assert field == :layers
    for (k, (name, thickness)) in enumerate(layers)
        layer = FUSEparameters__build_layer{T}()
        push!(parameters_build.layers, layer)

        # name
        layer.name = replace(string(name), "_" => " ")

        # thickness
        layer.thickness = thickness

        # type
        @show layer.name
        if occursin("OH", uppercase(layer.name)) && occursin("hfs", lowercase(layer.name))
            layer.type = :oh
        elseif occursin("gap ", lowercase(layer.name))
            layer.type = :gap
        elseif lowercase(layer.name) == "plasma"
            layer.type = :plasma
        elseif occursin("OH", uppercase(layer.name))
            layer.type = :oh
        elseif occursin("TF", uppercase(layer.name))
            layer.type = :tf
        elseif occursin("shield", lowercase(layer.name))
            layer.type = :shield
        elseif occursin("blanket", lowercase(layer.name))
            layer.type = :blanket
        elseif occursin("wall", lowercase(layer.name))
            layer.type = :wall
        elseif occursin("vessel", lowercase(layer.name))
            layer.type = :vessel
        elseif occursin("cryostat", lowercase(layer.name))
            layer.type = :cryostat
            layer.shape = :silo
        end

        # side
        if occursin("hfs", lowercase(layer.name))
            layer.side = :hfs
        elseif occursin("lfs", lowercase(layer.name))
            layer.side = :lfs
        else
            if layer.type == _plasma_
                layer.side = :lhfs
            elseif k < length(layers) / 2
                layer.side = :in
            elseif k > length(layers) / 2
                layer.side = :out
            end
        end
    end
end

"""
    dict2par!(dct::AbstractDict, par::ParametersVector{<:FUSEparameters__build_layer})

Custom reading from file of FUSEparameters__build_layer
"""
function SimulationParameters.dict2par!(dct::AbstractDict, par::ParametersVector{<:FUSEparameters__build_layer})
    return parent(par).layers = dct
end


"""
    par2ystr(par::ParametersVector{<:FUSEparameters__build_layer}, txt::Vector{String})

Custom writing to file for FUSEparameters__build_layer
"""
function SimulationParameters.par2ystr(par::ParametersVector{<:FUSEparameters__build_layer}, txt::Vector{String}; show_info::Bool=true, skip_defaults::Bool=false)
    for parameter in par
        p = SimulationParameters.path(parameter)
        sp = SimulationParameters.spath(p)
        depth = (count(".", sp) + count("[", sp) - 1) * 2
        pre = " "^depth
        push!(txt, string(pre, replace(getproperty(parameter, :name, "_each_layer_name_and_thickness_"), " " => "_"), ": ", repr(getproperty(parameter, :thickness, 0.0))))
    end
    return txt
end

"""
    to_index(layers::Vector{FUSEparameters__build_layer{T}}, name::Symbol) where {T<:Real}

Allows accesing parameters layers by their Symbol
"""
function Base.to_index(layers::Vector{FUSEparameters__build_layer{T}}, name::Symbol) where {T<:Real}
    tmp = findfirst(x -> x.name == replace(string(name), "_" => " "), layers)
    if tmp === nothing
        error("Valid ini.build.layers are: $([Symbol(replace(layer.name," " => "_")) for layer in layers])")
    end
    return tmp
end

"""
    scale_build_layers!(layers::ParametersVector{<:FUSEparameters__build_layer}, R0::Float64, a::Float64, plasma_gap_fraction::Float64)

Scale build layers thicknesses so that the plasma equilibrium is in the middle of the plasma layer
"""
function scale_build_layers!(layers::ParametersVector{<:FUSEparameters__build_layer}, R0::Float64, a::Float64, plasma_gap_fraction::Float64)
    gap = a * plasma_gap_fraction
    plasma_start = R0 - a - gap
    layer_plasma_start = 0.0
    for layer in layers
        if layer.type == _plasma_
            break
        end
        if layer.thickness > 0.0
            layer_plasma_start += layer.thickness
        end
    end
    factor = plasma_start / layer_plasma_start
    for (k, layer) in enumerate(layers)
        if layer.type == _plasma_
            layers[k].thickness = 2.0 * (a + gap)
        else
            layers[k].thickness = layer.thickness * factor
        end
    end
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
        if !ismissing(layer, :material)
            continue
        end
        if k == 1 && ini.center_stack.plug
            layer.material = "Steel, Stainless 316"
        elseif layer.type == Int(_plasma_)
            layer.material = any((layer.type in (Int(_blanket_), Int(_shield_)) for layer in bd.layer)) ? "DT_plasma" : "DD_plasma"
        elseif layer.type == Int(_gap_)
            layer.material = "Vacuum"
        elseif layer.type == Int(_oh_)
            layer.material = bd.oh.technology.material
        elseif layer.type == Int(_tf_)
            layer.material = bd.tf.technology.material
        elseif layer.type == Int(_shield_)
            layer.material = "Steel, Stainless 316"
        elseif layer.type == Int(_blanket_)
            layer.material = "lithium-lead"
        elseif layer.type == Int(_wall_)
            layer.material = "Tungsten"
        elseif layer.type == Int(_vessel_)
            layer.material = "Water, Liquid"
        elseif layer.type == Int(_cryostat_)
            layer.material = "Steel, Stainless 316"
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

"""
    layers_meters_from_fractions(; blanket::Float64, shield::Float64, vessel::Float64, pf_inside_tf::Bool, pf_outside_tf::Bool, thin_vessel_walls::Bool=false)

Handy function for initializing layers based on few scalars
"""
function layers_meters_from_fractions(; blanket::Float64, shield::Float64, vessel::Float64, pf_inside_tf::Bool, pf_outside_tf::Bool, thin_vessel_walls::Bool=false)

    # express layer thicknesses as fractions
    layers = OrderedCollections.OrderedDict{Symbol,Float64}()
    layers[:gap_OH] = 2.0
    layers[:OH] = 1.0
    layers[:hfs_TF] = 1.0
    if vessel > 0.0
        if thin_vessel_walls
            layers[:hfs_vacuum_vessel_wall_outer] = 0.1
            layers[:hfs_vacuum_vessel] = vessel
            layers[:hfs_vacuum_vessel_wall_inner] = 0.1
        else
            layers[:hfs_vacuum_vessel] = vessel
        end
    end
    layers[:gap_hfs_coils] = pf_inside_tf ? 0 : -1
    if shield > 0.0
        layers[:hfs_shield] = shield
    end
    if blanket > 0.0
        layers[:hfs_blanket] = blanket
    end
    layers[:hfs_wall] = 0.1
    layers[:plasma] = 0.0 # this number does not matter
    layers[:lfs_wall] = 0.1
    if blanket > 0.0
        layers[:lfs_blanket] = blanket * 2.0
    end
    if shield > 0.0
        layers[:lfs_shield] = shield
    end
    layers[:gap_lfs_coils] = 1.0 * (pf_inside_tf ? 2.25 : -1)
    if vessel > 0.0
        if thin_vessel_walls
            layers[:lfs_vacuum_vessel_wall_inner] = 0.1
            layers[:lfs_vacuum_vessel] = vessel
            layers[:lfs_vacuum_vessel_wall_outer] = 0.1
        else
            layers[:lfs_vacuum_vessel] = vessel
        end
    end
    layers[:lfs_TF] = 1.0
    if blanket > 0.0
        layers[:gap_cryostat] = 1.0 * (pf_outside_tf ? 3 : 1)
        layers[:cryostat] = 0.5
    else
        layers[:gap_world] = 1.0 * (pf_outside_tf ? 3 : 1)
    end

    return layers
end
