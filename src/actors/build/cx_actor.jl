#= ============= =#
#  cross-section  #
#= ============= =#
Base.@kwdef mutable struct FUSEparameters__ActorCXbuild{T} <: ParametersActor where {T<:Real}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    rebuild_wall::Entry{Bool} = Entry(Bool, "-", "Rebuild wall based on equilibrium"; default=false)
    do_plot::Entry{Bool} = Entry(Bool, "-", "plot"; default=false)
end

mutable struct ActorCXbuild <: ReactorAbstractActor
    dd::IMAS.dd
    par::FUSEparameters__ActorCXbuild
    act::ParametersAllActors
    function ActorCXbuild(dd::IMAS.dd, par::FUSEparameters__ActorCXbuild, act::ParametersAllActors; kw...)
        logging_actor_init(ActorCXbuild)
        par = par(kw...)
        return new(dd, par, act)
    end
end

"""
    ActorCXbuild(dd::IMAS.dd, act::ParametersAllActors; kw...)

Generates the 2D cross section of the tokamak build

!!! note 
    Manipulates data in `dd.build`
"""
function ActorCXbuild(dd::IMAS.dd, act::ParametersAllActors; kw...)
    par = act.ActorCXbuild
    actor = ActorCXbuild(dd, par, act; kw...)
    step(actor)
    finalize(actor)
    if par.do_plot
        plot(dd.build)
        display(plot!(dd.build; cx=false))
    end
    return actor
end

function _step(actor::ActorCXbuild)
    dd = actor.dd
    par = actor.par
    act = actor.act

    # empty layer outlines and structures
    for layer in dd.build.layer
        empty!(layer.outline)
    end
    empty!(dd.build.structure)
    actor_pfs = nothing
    wall = PlasmaFacingSurfaces.get_merged_wall_outline(dd.wall)
    if wall === missing || par.rebuild_wall
        actor_pfs = ActorPlasmaFacingSurfaces(dd, act)
        wall = PlasmaFacingSurfaces.get_merged_wall_outline(dd.wall)
    end
    pr = wall.r
    pz = wall.z

    build_cx!(dd.build, pr, pz)

    #divertor_regions!(dd.build, dd.equilibrium.time_slice[])

    blanket_regions!(dd.build, dd.equilibrium.time_slice[])
    
    #JG: redo to adjust the pfs to current build. Not sure where to set a switch to call adjust2build procedure. probabbly in cx_build?
    if actor_pfs !==nothing 
        PlasmaFacingSurfaces.adjust2build(actor_pfs.pfs,dd)
        actor_pfs.pfs()
        finalize(actor_pfs)
        wall = PlasmaFacingSurfaces.get_merged_wall_outline(dd.wall)
        pr = wall.r
        pz = wall.z
        build_cx!(dd.build, pr, pz)

    #divertor_regions!(dd.build, dd.equilibrium.time_slice[])

        blanket_regions!(dd.build, dd.equilibrium.time_slice[])
    end

    return actor
end

"""
    build_cx!(bd::IMAS.build, pr::Vector{Float64}, pz::Vector{Float64})

Translates 1D build to 2D cross-sections starting from R and Z coordinates of plasma first wall
"""
function build_cx!(bd::IMAS.build, pr::Vector{Float64}, pz::Vector{Float64})
    ipl = IMAS.get_build(bd, type=_plasma_, return_index=true)
    itf = IMAS.get_build(bd, type=_tf_, fs=_hfs_, return_index=true)

    # rescale _plasma_ outline to match 1D radial build
    # start_radius = bd.layer[ipl].start_radius
    # end_radius = bd.layer[ipl].end_radius
    # pr1 = minimum(pr)
    # pr2 = maximum(pr)
    # fact = (end_radius - start_radius) / (pr2 - pr1)
    # pz .= pz .* fact
    # pr .= (pr .- pr1) .* fact .+ start_radius
    bd.layer[ipl].outline.r = pr
    bd.layer[ipl].outline.z = pz

    coils_inside = any([contains(lowercase(l.name), "coils") for l in bd.layer])

    # all layers between plasma and OH
    # k-1 means the layer outside (ie. towards the tf)
    # k   is the current layer
    # k+1 means the layer inside (ie. towards the plasma)
    # 
    # forward pass: from plasma to TF _convex_hull_ and then desired TF shape
    tf_to_plasma = IMAS.get_build(bd, fs=_hfs_, return_only_one=false, return_index=true)
    plasma_to_tf = reverse(tf_to_plasma)
    for k in plasma_to_tf
        original_shape = bd.layer[k].shape
        if k == itf + 1
            # layer that is inside of the TF sets TF shape
            layer_shape = BuildLayerShape(mod(mod(bd.layer[k].shape, 1000), 100))
            @debug "forward $(bd.layer[k].name) $(layer_shape)"
            optimize_shape(bd, k + 1, k, layer_shape; tight=!coils_inside)
        else
            # everything else is conformal convex hull
            @debug "forward $(bd.layer[k].name) CONVEX HULL"
            optimize_shape(bd, k + 1, k, _convex_hull_)
            bd.layer[k].shape = original_shape
        end
    end
    # NOTE:: To debug it helps to comment out the subsequent reverse and forward passes
    # reverse pass: from TF to plasma only with negative offset
    for k in tf_to_plasma[2:end]
        if bd.layer[k+1].shape == Int(_negative_offset_)
            @debug "reverse $(bd.layer[k].name) NEGATIVE OFFSET"
            optimize_shape(bd, k, k + 1, _negative_offset_)
        end
    end
    # forward pass: from plasma to TF with desired shapes
    for k in plasma_to_tf[1:end-1]
        if bd.layer[k].shape == Int(_negative_offset_)
            break
        else
            layer_shape = BuildLayerShape(mod(mod(bd.layer[k].shape, 1000), 100))
            @debug "reverse $(bd.layer[k].name) $(layer_shape)"
            optimize_shape(bd, k + 1, k, layer_shape)
        end
    end

    # _in_
    D = minimum(IMAS.get_build(bd, type=_tf_, fs=_hfs_).outline.z)
    U = maximum(IMAS.get_build(bd, type=_tf_, fs=_hfs_).outline.z)
    for k in IMAS.get_build(bd, fs=_in_, return_index=true, return_only_one=false)
        L = bd.layer[k].start_radius
        R = bd.layer[k].end_radius
        bd.layer[k].outline.r, bd.layer[k].outline.z = rectangle_shape(L, R, D, U)
    end

    # _out_
    iout = IMAS.get_build(bd, fs=_out_, return_index=true, return_only_one=false)
    if lowercase(bd.layer[iout[end]].name) == "cryostat"
        olfs = IMAS.get_build(bd, fs=_lfs_, return_index=true, return_only_one=false)[end]
        optimize_shape(bd, olfs, iout[end], _silo_)
        for k in reverse(iout[2:end])
            optimize_shape(bd, k, k - 1, _negative_offset_)
        end
    else
        for k in iout
            L = 0.0
            R = bd.layer[k].end_radius
            D = minimum(bd.layer[k-1].outline.z) - bd.layer[k].thickness
            U = maximum(bd.layer[k-1].outline.z) + bd.layer[k].thickness
            bd.layer[k].outline.r, bd.layer[k].outline.z = rectangle_shape(L, R, D, U)
        end
    end

    return bd
end

function divertor_regions!(bd::IMAS.build, eqt::IMAS.equilibrium__time_slice)
    R0 = eqt.global_quantities.magnetic_axis.r
    Z0 = eqt.global_quantities.magnetic_axis.z

    ipl = IMAS.get_build(bd, type=_plasma_, return_index=true)
    plasma_poly = xy_polygon(bd.layer[ipl])
    ψb = IMAS.find_psi_boundary(eqt)
    rlcfs, zlcfs, _ = IMAS.flux_surface(eqt, ψb, true)
    linear_plasma_size = maximum(zlcfs) - minimum(zlcfs)

    wall_poly = xy_polygon(bd.layer[ipl-1])
    for ltype in [_blanket_, _shield_, _wall_,]
        iwl = IMAS.get_build(bd, type=ltype, fs=_hfs_, return_index=true, raise_error_on_missing=false)
        if iwl !== missing
            wall_poly = xy_polygon(bd.layer[iwl])
            break
        end
    end

    divertors = IMAS.IDSvectorElement[]
    for x_point in eqt.boundary.x_point
        Rx = x_point.r
        Zx = x_point.z
        if IMAS.minimum_distance_two_shapes(rlcfs, zlcfs, [Rx], [Zx]) > (linear_plasma_size / 20)
            # secondary Xpoint far away
            continue
        end

        m = (Zx - Z0) / (Rx - R0)
        xx = [0, R0 * 2.0]
        yy = line_through_point(-1.0 ./ m, Rx, Zx, xx)
        pr = vcat(xx, reverse(xx), xx[1])
        pz = vcat(yy, [Zx * 5, Zx * 5], yy[1])

        domain = xy_polygon(pr, pz)
        divertor_poly = LibGEOS.intersection(wall_poly, domain)
        divertor_poly = LibGEOS.difference(divertor_poly, plasma_poly)

        pr = [v[1] for v in GeoInterface.coordinates(divertor_poly)[1]]
        pz = [v[2] for v in GeoInterface.coordinates(divertor_poly)[1]]

        # assign to build structure
        if Zx > Z0
            name = "Upper divertor"
        else
            name = "Lower divertor"
        end
        structure = resize!(bd.structure, "type" => Int(_divertor_), "name" => name)
        structure.material = "Tungsten"
        structure.outline.r = pr
        structure.outline.z = pz
        structure.toroidal_extent = 2pi

        push!(divertors, structure)
    end

    return divertors
end

function blanket_regions!(bd::IMAS.build, eqt::IMAS.equilibrium__time_slice)
    R0 = eqt.global_quantities.magnetic_axis.r

    layers = bd.layer
    iblanket = IMAS.get_build(bd; type=_blanket_, fs=_lfs_, return_index=true, raise_error_on_missing=false)
    if iblanket === missing
        return IMAS.IDSvectorElement[]
    end
    layer = layers[iblanket]
    layer_in = layers[iblanket-1]

    layer_poly = xy_polygon(layer)
    layer_in_poly = xy_polygon(layer_in)
    ring_poly = LibGEOS.difference(layer_poly, layer_in_poly)
    for structure in [structure for structure in bd.structure if structure.type == Int(_divertor_)]
        structure_poly = xy_polygon(structure)
        ring_poly = LibGEOS.difference(ring_poly, structure_poly)
    end

    geometries = LibGEOS.getGeometries(ring_poly)
    blankets = IMAS.IDSvectorElement[]
    for poly in geometries
        coords = GeoInterface.coordinates(poly)
        pr = [v[1] for v in coords[1]]
        pz = [v[2] for v in coords[1]]

        # assign to build structure
        if length(geometries) == 2
            if sum(pr) / length(pr) > R0
                name = "LFS blanket"
            else
                name = "HFS blanket"
            end
        else
            name = "blanket"
        end

        structure = resize!(bd.structure, "type" => Int(_blanket_), "name" => name)
        structure.outline.r = pr
        structure.outline.z = pz
        structure.toroidal_extent = 2pi

        push!(blankets, structure)
    end

    return blankets
end

"""
    optimize_shape(bd::IMAS.build, obstr_index::Int, layer_index::Int, shape::BuildLayerShape)

Generates outline of layer in such a way to maintain minimum distance from inner layer
"""
function optimize_shape(bd::IMAS.build, obstr_index::Int, layer_index::Int, shape::BuildLayerShape; tight::Bool=true)
    layer = bd.layer[layer_index]
    obstr = bd.layer[obstr_index]
    # display("Layer $layer_index = $(layer.name)")
    # display("Obstr $obstr_index = $(obstr.name)")
    if layer.fs == Int(_out_)
        l_start = 0.0
        l_end = layer.end_radius
        o_start = 0.0
        o_end = obstr.end_radius
    else
        if obstr.fs in [Int(_lhfs_), Int(_out_)]
            o_start = obstr.start_radius
            o_end = obstr.end_radius
        else
            o_start = obstr.start_radius
            o_end = IMAS.get_build(bd, identifier=obstr.identifier, fs=_lfs_).end_radius
        end
        l_start = layer.start_radius
        if layer.type == Int(_plasma_)
            l_end = layer.end_radius
        else
            l_end = IMAS.get_build(bd, identifier=layer.identifier, fs=_lfs_).end_radius
        end
    end
    hfs_thickness = o_start - l_start
    lfs_thickness = l_end - o_end
    oR = obstr.outline.r
    oZ = obstr.outline.z

    if layer.fs == Int(_out_)
        target_minimum_distance = lfs_thickness
    else
        if tight
            target_minimum_distance = min(hfs_thickness, lfs_thickness)
        else
            target_minimum_distance = sqrt(hfs_thickness^2 + lfs_thickness^2) / 2.0
        end
    end
    r_offset = (lfs_thickness .- hfs_thickness) / 2.0

    # update shape
    layer.shape = Int(shape)

    # handle offset, negative offset, negative offset, and convex-hull
    if layer.shape in [Int(_offset_), Int(_negative_offset_), Int(_convex_hull_)]
        poly = LibGEOS.buffer(xy_polygon(oR, oZ), (hfs_thickness + lfs_thickness) / 2.0)
        R = [v[1] .+ r_offset for v in GeoInterface.coordinates(poly)[1]]
        Z = [v[2] for v in GeoInterface.coordinates(poly)[1]]
        if layer.shape == Int(_convex_hull_)
            hull = convex_hull(R, Z; closed_polygon=true)
            R = [r for (r, z) in hull]
            Z = [z for (r, z) in hull]
            # resample disabled because this can lead to outlines of different layers to be crossing
            # R, Z = IMAS.resample_2d_path(R, Z)
        end
        layer.outline.r, layer.outline.z = R, Z

    else # handle shapes
        if layer.shape > 1000
            layer.shape = mod(layer.shape, 1000)
        end
        if layer.shape > 100
            layer.shape = mod(layer.shape, 100)
        end

        if layer.shape == Int(_silo_)
            is_up_down_symmetric = false
        elseif abs(sum(oZ) / sum(abs.(oZ))) < 1E-2
            is_up_down_symmetric = true
        else
            is_up_down_symmetric = false
        end

        is_negative_D = false
        if layer.shape != Int(_silo_)
            _, imaxr = findmax(oR)
            _, iminr = findmin(oR)
            _, imaxz = findmax(oZ)
            _, iminz = findmin(oZ)
            r_at_max_z, max_z = oR[imaxz], oZ[imaxz]
            r_at_min_z, min_z = oR[iminz], oZ[iminz]
            z_at_max_r, max_r = oZ[imaxr], oR[imaxr]
            z_at_min_r, min_r = oZ[iminr], oR[iminr]
            a = 0.5 * (max_r - min_r)
            R = 0.5 * (max_r + min_r)
            δu = (R - r_at_max_z) / a
            δl = (R - r_at_min_z) / a
            if δu + δl < -0.1
                is_negative_D = true
            end
        end

        if is_negative_D
            layer.shape = layer.shape + 1000
        end

        if !is_up_down_symmetric
            layer.shape = layer.shape + 100
        end

        func = shape_function(layer.shape)
        layer.shape_parameters = initialize_shape_parameters(layer.shape, oR, oZ, l_start, l_end, target_minimum_distance)

        layer.outline.r, layer.outline.z = func(l_start, l_end, layer.shape_parameters...)
        layer.shape_parameters = optimize_shape(oR, oZ, target_minimum_distance, func, l_start, l_end, layer.shape_parameters)
        layer.outline.r, layer.outline.z = func(l_start, l_end, layer.shape_parameters...; resample=false)
    end

    IMAS.reorder_flux_surface!(layer.outline.r, layer.outline.z)
    # display(plot(layer.outline.r, layer.outline.z; aspect_ratio=:equal))
end

function assign_build_layers_materials(dd::IMAS.dd, ini::ParametersAllInits)
    bd = dd.build
    for (k, layer) in enumerate(bd.layer)
        if k == 1 && ini.center_stack.plug
            layer.material = ini.material.wall
        elseif layer.type == Int(_plasma_)
            layer.material = any([layer.type in [Int(_blanket_), Int(_shield_)] for layer in dd.build.layer]) ? "DT_plasma" : "DD_plasma"
        elseif layer.type == Int(_gap_)
            layer.material = "Vacuum"
        elseif layer.type == Int(_oh_)
            layer.material = dd.build.oh.technology.material
        elseif layer.type == Int(_tf_)
            layer.material = dd.build.tf.technology.material
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

