import Interpolations
import OrderedCollections

#= ========== =#
#  init build  #
#= ========== =#
"""
    init_build!(dd::IMAS.dd, ini::ParametersAllInits, act::ParametersAllActors, dd1::IMAS.dd=IMAS.dd())

Initialize `dd.build` starting from `ini` and `act` parameters
"""
function init_build!(dd::IMAS.dd, ini::ParametersAllInits, act::ParametersAllActors, dd1::IMAS.dd=IMAS.dd())
    TimerOutputs.reset_timer!("init_build")
    TimerOutputs.@timeit timer "init_build" begin
        init_from = ini.general.init_from

        if init_from == :ods
            if length(dd1.build.layer) > 0
                dd.build = deepcopy(dd1.build)
            else
                init_from = :scalars
            end
        end

        eqt = dd.equilibrium.time_slice[]

        # scale radial build layers based on equilibrium R0, a, and the requested plasma_gap
        scale_build_layers!(ini.build.layers, eqt.boundary.geometric_axis.r, eqt.boundary.minor_radius, ini.build.plasma_gap)

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

        # set build layer shapes
        plama_to_tf = IMAS.get_build_indexes(dd.build.layer; fs=_lfs_)

        # everything inside of the TF is _negative_offset_ by default
        for k in plama_to_tf
            dd.build.layer[k].shape = Int(_negative_offset_)
        end

        # Allow first few layers to have shapes conformal with the first wall
        for k in 1:ini.build.n_first_wall_conformal_layers
            dd.build.layer[plama_to_tf[k]].shape = Int(_convex_hull_)
        end

        # first wall is of type offset instead of convex hull, to allow for concave shape
        k = plama_to_tf[1]
        if (dd.build.layer[k].type == Int(_wall_)) &&
           ((dd.build.layer[k+1].type == Int(_blanket_)) || (dd.build.layer[k+1].type == Int(_shield_)) || (dd.build.layer[k+1].type == Int(_wall_)))
            dd.build.layer[k].shape = Int(_offset_)
        end

        # set TF shape (shape is set by inner layer)
        dd.build.layer[plama_to_tf[end-1]].shape = Int(ini.tf.shape)
        dd.build.layer[plama_to_tf[end]].shape = Int(_convex_hull_)

        # number of TF coils
        dd.build.tf.coils_n = ini.tf.n_coils
        # target TF ripple
        dd.build.tf.ripple = ini.tf.ripple

        # 2D build cross-section
        ActorCXbuild(dd, act)

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
        @assert ini_layer.thickness >= 0.0
        k += 1
        layer = resize!(bd.layer, k)[k]

        layer.name = ini_layer.name
        layer.thickness = ini_layer.thickness
        layer.type = Int(ini_layer.type)
        layer.side = Int(ini_layer.side)
        if !ismissing(ini_layer, :material)
            layer.material = String(ini_layer.material)
        end
        if !ismissing(ini_layer, :shape)
            layer.shape = Int(ini_layer.shape)
        end
    end

    return bd
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

"""
    wall_radii(R0::Float64, a::Float64, plasma_gap_fraction::Float64)

Returns the hfs and lfs radii of the wall
"""
function wall_radii(R0::Float64, a::Float64, plasma_gap_fraction::Float64)
    gap = a * plasma_gap_fraction
    r_hfs = R0 - a - gap
    r_lfs = R0 + a + gap
    return (r_hfs=r_hfs, r_lfs=r_lfs)
end

function assign_build_layers_materials(dd::IMAS.dd, ini::ParametersAllInits)
    bd = dd.build
    for (k, layer) in enumerate(bd.layer)
        if !ismissing(layer, :material)
            continue
        end
        if k == 1 && ini.center_stack.plug
            layer.material = "steel"
        elseif layer.type == Int(_plasma_)
            layer.material = "plasma"
        elseif layer.type == Int(_gap_)
            layer.material = "vacuum"
        elseif layer.type == Int(_oh_)
            layer.material = bd.oh.technology.material
        elseif layer.type == Int(_tf_)
            layer.material = bd.tf.technology.material
        elseif layer.type == Int(_shield_)
            layer.material = "steel"
        elseif layer.type == Int(_blanket_)
            layer.material = "lithium_lead"
        elseif layer.type == Int(_wall_)
            if bd.layer[k-1].type == Int(_plasma_) || bd.layer[k+1].type == Int(_plasma_)
                layer.material = "tungsten"
            else
                layer.material = "steel"
            end
        elseif layer.type == Int(_vessel_)
            layer.material = "water"
        elseif layer.type == Int(_cryostat_)
            layer.material = "steel"
        end
    end
end

function assign_technologies(dd::IMAS.dd, ini::ParametersAllInits)
    # plug
    if ini.center_stack.plug
        IMAS.mechanical_technology(dd, :pl)
    end

    # oh
    IMAS.coil_technology(dd.build.oh.technology, ini.oh.technology, :oh)
    IMAS.mechanical_technology(dd, :oh)

    # tf
    IMAS.coil_technology(dd.build.tf.technology, ini.tf.technology, :tf)
    IMAS.mechanical_technology(dd, :tf)

    #pf
    return IMAS.coil_technology(dd.build.pf_active.technology, ini.pf_active.technology, :pf_active)
end

"""
    layers_meters_from_fractions(;
        lfs_multiplier::Float64=0.0,
        wall::Float64=0.1,
        blanket::Float64,
        shield::Float64,
        vessel::Float64,
        pf_inside_tf::Bool,
        pf_outside_tf::Bool,
        thin_vessel_walls::Bool=false
    )

Handy function for initializing layers based on few scalars
"""
function layers_meters_from_fractions(;
    lfs_multiplier::Float64,
    wall::Float64,
    blanket::Float64,
    shield::Float64,
    vessel::Float64,
    pf_inside_tf::Bool,
    pf_outside_tf::Bool,
    thin_vessel_walls::Bool=false
)

    lfs_asymmetry(x) = x * lfs_multiplier

    # express layer thicknesses as fractions
    layers = OrderedCollections.OrderedDict{Symbol,Float64}()
    layers[:gap_OH] = 2.0
    layers[:OH] = 1.0
    layers[:hfs_TF] = 1.0
    @assert vessel >= 0.0
    if thin_vessel_walls
        layers[:hfs_vacuum_vessel_wall_outer] = 0.1 * vessel
        layers[:hfs_vacuum_vessel] = 0.8 * vessel
        layers[:hfs_vacuum_vessel_wall_inner] = 0.1 * vessel
    else
        layers[:hfs_vacuum_vessel] = vessel
    end
    if pf_inside_tf
        layers[:gap_hfs_coils] = 0.0
    end
    if shield > 0.0
        layers[:hfs_shield] = shield
    end
    if blanket > 0.0
        layers[:hfs_blanket] = blanket
    end
    layers[:hfs_wall] = wall
    layers[:plasma] = 0.0 # this number does not matter
    layers[:lfs_wall] = lfs_asymmetry(wall)
    if blanket > 0.0
        layers[:lfs_blanket] = lfs_asymmetry(blanket * 2.0)
    end
    if shield > 0.0
        layers[:lfs_shield] = lfs_asymmetry(shield)
    end
    if pf_inside_tf
        layers[:gap_lfs_coils] = lfs_asymmetry(2.25)
    end
    if vessel > 0.0
        if thin_vessel_walls
            layers[:lfs_vacuum_vessel_wall_inner] = lfs_asymmetry(0.1 * vessel)
            layers[:lfs_vacuum_vessel] = lfs_asymmetry(0.8 * vessel)
            layers[:lfs_vacuum_vessel_wall_outer] = lfs_asymmetry(0.1 * vessel)
        else
            layers[:lfs_vacuum_vessel] = lfs_asymmetry(vessel)
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
