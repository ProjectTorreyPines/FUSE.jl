#= ============= =#
#  cross-section  #
#= ============= =#
Base.@kwdef mutable struct FUSEparameters__ActorCXbuild{T} <: ParametersActor where {T<:Real}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    rebuild_wall::Entry{Bool} = Entry{Bool}("-", "Rebuild wall based on equilibrium"; default=true)
    do_plot::Entry{Bool} = Entry{Bool}("-", "Plot"; default=false)
end

mutable struct ActorCXbuild{D,P} <: ReactorAbstractActor
    dd::IMAS.dd{D}
    par::FUSEparameters__ActorCXbuild{P}
    function ActorCXbuild(dd::IMAS.dd{D}, par::FUSEparameters__ActorCXbuild{P}; kw...) where {D<:Real,P<:Real}
        logging_actor_init(ActorCXbuild)
        par = par(kw...)
        return new{D,P}(dd, par)
    end
end

"""
    ActorCXbuild(dd::IMAS.dd, act::ParametersAllActors; kw...)

Generates the 2D cross section of the tokamak build

!!! note 
    Manipulates data in `dd.build`
"""
function ActorCXbuild(dd::IMAS.dd, act::ParametersAllActors; kw...)
    actor = ActorCXbuild(dd, act.ActorCXbuild; kw...)
    step(actor)
    finalize(actor)
    if actor.par.do_plot
        plot(dd.build)
        display(plot!(dd.build; cx=false))
    end
    return actor
end

function _step(actor::ActorCXbuild)
    dd = actor.dd
    par = actor.par

    bd = dd.build
    eqt = dd.equilibrium.time_slice[]

    # If wall information is missing, then the first wall information is generated starting from equilibrium time_slice
    wall = IMAS.first_wall(dd.wall)
    if wall === missing || par.rebuild_wall
        pr, pz = wall_from_eq(bd, eqt)
    else
        pr = wall.r
        pz = wall.z
    end

    # empty layer outlines and structures
    for layer in bd.layer
        empty!(layer.outline)
    end
    empty!(bd.structure)

    build_cx!(bd, pr, pz)

    divertor_regions!(bd, eqt, dd.divertors)

    blanket_regions!(bd, eqt)

    if wall === missing || par.rebuild_wall
        plasma = IMAS.get_build_layer(bd.layer, type=_plasma_)
        resize!(dd.wall.description_2d, 1)
        resize!(dd.wall.description_2d[1].limiter.unit, 1)
        dd.wall.description_2d[1].limiter.unit[1].outline.r = plasma.outline.r
        dd.wall.description_2d[1].limiter.unit[1].outline.z = plasma.outline.z
    end

    IMAS.find_strike_points!(eqt, dd.divertors)

    return actor
end

"""
    wall_from_eq(bd::IMAS.build, eqt::IMAS.equilibrium__time_slice; divertor_length_fraction::Real=0.2)

Generate first wall outline starting from an equilibrium
"""
function wall_from_eq(bd::IMAS.build, eqt::IMAS.equilibrium__time_slice; divertor_length_fraction::Real=0.2)
    R0 = eqt.global_quantities.magnetic_axis.r
    Z0 = eqt.global_quantities.magnetic_axis.z

    # lcfs
    ψb = IMAS.find_psi_boundary(eqt)
    rlcfs, zlcfs, _ = IMAS.flux_surface(eqt, ψb, true)

    # Set the radial build thickness of the plasma
    plasma = IMAS.get_build_layer(bd.layer, type=_plasma_)
    a = (minimum(rlcfs) - plasma.start_radius)
    plasma.thickness = maximum(rlcfs) - minimum(rlcfs) + 2.0 * a
    R_hfs_plasma = plasma.start_radius
    R_lfs_plasma = plasma.end_radius

    # main chamber (clip elements that go beyond plasma radial build thickness)
    R, Z = buffer(rlcfs, zlcfs, a)
    R[R.<R_hfs_plasma] .= R_hfs_plasma
    R[R.>R_lfs_plasma] .= R_lfs_plasma
    Z = (Z .- Z0) .* 1.05 .+ Z0
    wall_poly = xy_polygon(R, Z)

    t = LinRange(0, 2π, 31)

    # divertor lengths
    linear_plasma_size = maximum(zlcfs) - minimum(zlcfs)
    max_divertor_length = linear_plasma_size * divertor_length_fraction

    # private flux regions
    private = IMAS.flux_surface(eqt, ψb, false)
    for (pr, pz) in private
        if sign(pz[1] - Z0) != sign(pz[end] - Z0)
            # open flux surface does not encicle the plasma
            continue
        end

        # xpoint between lcfs and private region
        index = IMAS.minimum_distance_two_shapes(pr, pz, rlcfs, zlcfs; return_index=true)
        Rx = (pr[index[1]] + rlcfs[index[2]]) / 2.0
        Zx = (pz[index[1]] + zlcfs[index[2]]) / 2.0
        d = sqrt((pr[index[1]] - rlcfs[index[2]])^2 + (pz[index[1]] - zlcfs[index[2]])^2)
        if d > linear_plasma_size / 5
            continue
        end

        # check that this divertor is in fact installed
        if Zx > Z0 && (ismissing(bd.divertors.upper, :installed) || bd.divertors.upper.installed != 1)
            continue
        elseif Zx < Z0 && (ismissing(bd.divertors.lower, :installed) || bd.divertors.lower.installed != 1)
            continue
        end

        max_d = maximum(sqrt.((Rx .- pr) .^ 2.0 .+ (Zx .- pz) .^ 2.0))
        divertor_length = min(max_d * 0.8, max_divertor_length)

        # limit extent of private flux regions
        circle = collect(zip(divertor_length .* cos.(t) .+ Rx, sign(Zx) .* divertor_length .* sin.(t) .+ Zx))
        circle[1] = circle[end]
        slot = [(rr, zz) for (rr, zz) in zip(pr, pz) if PolygonOps.inpolygon((rr, zz), circle) == 1 && rr >= R_hfs_plasma && rr <= R_lfs_plasma]
        pr = [rr for (rr, zz) in slot]
        pz = [zz for (rr, zz) in slot]

        # remove private flux region from wall (necessary because of Z expansion)
        # this may fail if the private region is small or weirdly shaped
        try
            wall_poly = LibGEOS.difference(wall_poly, xy_polygon(pr, pz))
        catch e
            if !(typeof(e) <: LibGEOS.GEOSError)
                rethrow(e)
            end
        end

        # add the divertor slots
        α = 0.2
        pr = vcat(pr, R0 * α + Rx * (1 - α))
        pz = vcat(pz, Z0 * α + Zx * (1 - α))
        slot = LibGEOS.buffer(xy_polygon(pr, pz), a)
        wall_poly = LibGEOS.union(wall_poly, slot)
    end

    # vertical clip
    wall_poly = LibGEOS.difference(wall_poly, xy_polygon(rectangle_shape(0.0, R_hfs_plasma, 100.0)...))
    wall_poly = LibGEOS.difference(wall_poly, xy_polygon(rectangle_shape(R_lfs_plasma, 10 * R_lfs_plasma, 100.0)...))

    # round corners
    wall_poly = LibGEOS.buffer(wall_poly, -a / 4)
    wall_poly = LibGEOS.buffer(wall_poly, a / 4)

    pr = [v[1] for v in GeoInterface.coordinates(wall_poly)[1]]
    pz = [v[2] for v in GeoInterface.coordinates(wall_poly)[1]]

    pr, pz = IMAS.resample_2d_path(pr, pz; step=0.1)

    return pr, pz
end

function divertor_regions!(bd::IMAS.build, eqt::IMAS.equilibrium__time_slice, divertors::IMAS.divertors)
    R0 = eqt.global_quantities.magnetic_axis.r
    Z0 = eqt.global_quantities.magnetic_axis.z

    # plasma poly (this sets the divertor's pfs)
    ipl = IMAS.get_build_index(bd.layer, type=_plasma_)
    plasma_poly = xy_polygon(bd.layer[ipl])
    pl_r = bd.layer[ipl].outline.r
    pl_z = bd.layer[ipl].outline.z
    pl_r[1] = pl_r[end]
    pl_z[1] = pl_z[end]
    IMAS.reorder_flux_surface!(pl_r, pl_z, R0, Z0)

    # wall poly (this sets how back the divertor structure goes)
    wall_poly = xy_polygon(bd.layer[ipl-1])
    for ltype in (_blanket_, _shield_, _wall_)
        iwls = IMAS.get_build_indexes(bd.layer, type=ltype, fs=_hfs_)
        if !isempty(iwls)
            wall_poly = xy_polygon(bd.layer[iwls[1]])
            break
        end
    end

    ψb = IMAS.find_psi_boundary(eqt)
    rlcfs, zlcfs, _ = IMAS.flux_surface(eqt, ψb, true)
    linear_plasma_size = maximum(zlcfs) - minimum(zlcfs)

    empty!(divertors)

    private = IMAS.flux_surface(eqt, ψb, false)
    for (pr, pz) in private
        if sign(pz[1] - Z0) != sign(pz[end] - Z0)
            # open flux surface does not encicle the plasma
            continue
        end

        # xpoint between lcfs and private region
        index = IMAS.minimum_distance_two_shapes(pr, pz, rlcfs, zlcfs; return_index=true)
        Rx = (pr[index[1]] + rlcfs[index[2]]) / 2.0
        Zx = (pz[index[1]] + zlcfs[index[2]]) / 2.0
        d = sqrt((pr[index[1]] - rlcfs[index[2]])^2 + (pz[index[1]] - zlcfs[index[2]])^2)
        if d > linear_plasma_size / 5
            continue
        end

        # check that this divertor is in fact installed
        if Zx > Z0 && (ismissing(bd.divertors.upper, :installed) || bd.divertors.upper.installed != 1)
            continue
        elseif Zx < Z0 && (ismissing(bd.divertors.lower, :installed) || bd.divertors.lower.installed != 1)
            continue
        end

        if Zx > Z0
            ul_name = "upper"
        else
            ul_name = "lower"
        end

        # Define divertor region by tracing a line going through X-point
        # and perpendicular to line connecting X-point to magnetic axis
        m = (Zx - Z0) / (Rx - R0)
        xx = [0.0, R0 * 2.0]
        yy = line_through_point(-1.0 ./ m, Rx, Zx, xx)
        domain_r = vcat(xx, reverse(xx), xx[1])
        domain_z = vcat(yy, [Zx * 5.0, Zx * 5.0], yy[1])
        domain_poly = xy_polygon(domain_r, domain_z)
        divertor_poly = LibGEOS.intersection(wall_poly, domain_poly)
        divertor_poly = LibGEOS.difference(divertor_poly, plasma_poly)

        # Assign to build structure
        coords = GeoInterface.coordinates(divertor_poly)
        structure = resize!(bd.structure, "type" => Int(_divertor_), "name" => "$ul_name divertor")
        structure.material = "Tungsten"
        structure.outline.r = [v[1] for v in coords[1]]
        structure.outline.z = [v[2] for v in coords[1]]
        structure.toroidal_extent = 2pi

        # now find divertor plasma facing surfaces
        indexes, crossings = IMAS.intersection(xx, yy, pl_r, pl_z)
        divertor_r = [crossings[1][1]; pl_r[indexes[1][2]+1:indexes[2][2]+1]; crossings[2][1]]
        divertor_z = [crossings[1][2]; pl_z[indexes[1][2]+1:indexes[2][2]+1]; crossings[2][2]]

        # split inner/outer based on line connecting X-point to magnetic axis
        m = (Zx - Z0) / (Rx - R0)
        xx = [0.0, R0 * 2.0]
        yy = line_through_point(m, Rx, Zx, xx)
        indexes, crossings = IMAS.intersection(xx, yy, divertor_r, divertor_z)

        # add target info
        divertor = resize!(divertors.divertor, length(divertors.divertor) + 1)[end]
        for io_name in ("inner", "outer")
            target = resize!(divertor.target, "name" => "$ul_name $io_name target")
            tile = resize!(target.tile, "name" => "$ul_name $io_name tile")
            if ul_name == "upper"
                if io_name == "outer"
                    divertor_pfs_r = divertor_r[indexes[1][2]:end]
                    divertor_pfs_z = divertor_z[indexes[1][2]:end]
                else
                    divertor_pfs_r = divertor_r[1:indexes[1][2]]
                    divertor_pfs_z = divertor_z[1:indexes[1][2]]
                end
            else
                if io_name == "outer"
                    divertor_pfs_r = divertor_r[1:indexes[1][2]]
                    divertor_pfs_z = divertor_z[1:indexes[1][2]]
                else
                    divertor_pfs_r = divertor_r[indexes[1][2]:end]
                    divertor_pfs_z = divertor_z[indexes[1][2]:end]
                end
            end
            tile.surface_outline.r = divertor_pfs_r
            tile.surface_outline.z = divertor_pfs_z
        end
    end

    return nothing
end

function blanket_regions!(bd::IMAS.build, eqt::IMAS.equilibrium__time_slice)
    R0 = eqt.global_quantities.magnetic_axis.r

    layers = bd.layer
    iblankets = IMAS.get_build_indexes(bd.layer; type=_blanket_, fs=_lfs_)
    if isempty(iblankets)
        return nothing
    end
    iblanket = iblankets[1]
    layer = layers[iblanket]
    layer_in = layers[iblanket-1]

    layer_poly = xy_polygon(layer)
    layer_in_poly = xy_polygon(layer_in)
    ring_poly = LibGEOS.difference(layer_poly, layer_in_poly)
    for structure in (structure for structure in bd.structure if structure.type == Int(_divertor_))
        structure_poly = xy_polygon(structure)
        ring_poly = LibGEOS.difference(ring_poly, structure_poly)
    end

    geometries = LibGEOS.getGeometries(ring_poly)
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
    end

    return nothing
end

"""
    build_cx!(bd::IMAS.build, pr::Vector{Float64}, pz::Vector{Float64})

Translates 1D build to 2D cross-sections starting from R and Z coordinates of plasma first wall
"""
function build_cx!(bd::IMAS.build, pr::Vector{Float64}, pz::Vector{Float64})
    plasma = IMAS.get_build_layer(bd.layer, type=_plasma_)

    # _plasma_ outline scaled to match 1D radial build
    start_radius = plasma.start_radius
    end_radius = plasma.end_radius
    pr1 = minimum(pr)
    pr2 = maximum(pr)
    fact = (end_radius - start_radius) / (pr2 - pr1)
    pz .= pz .* fact
    pr .= (pr .- pr1) .* fact .+ start_radius
    plasma.outline.r = pr
    plasma.outline.z = pz

    coils_inside = any([contains(lowercase(l.name), "coils") for l in bd.layer])

    # all layers between plasma and TF
    # k-1 means the layer outside (ie. towards the tf)
    # k   is the current layer
    # k+1 means the layer inside (ie. towards the plasma)
    tf_to_plasma = IMAS.get_build_indexes(bd.layer, fs=_hfs_)
    plasma_to_tf = reverse(tf_to_plasma)
    for k in plasma_to_tf
        layer_shape = BuildLayerShape(mod(mod(bd.layer[k].shape, 1000), 100))
        @debug "$(bd.layer[k].name) $(layer_shape)"
        optimize_shape(bd, k + 1, k, layer_shape; tight=!coils_inside)
    end

    # _in_
    TF = IMAS.get_build_layer(bd.layer, type=_tf_, fs=_hfs_)
    D = minimum(TF.outline.z)
    U = maximum(TF.outline.z)
    if coils_inside
        # generally the OH does not go higher than the PF coils
        D += 2.0 * TF.thickness
        U -= 2.0 * TF.thickness
    end
    for k in IMAS.get_build_indexes(bd.layer, fs=_in_)
        L = bd.layer[k].start_radius
        R = bd.layer[k].end_radius
        bd.layer[k].outline.r, bd.layer[k].outline.z = rectangle_shape(L, R, D, U)
    end

    # _out_
    iout = IMAS.get_build_indexes(bd.layer, fs=_out_)
    if lowercase(bd.layer[iout[end]].name) == "cryostat"
        olfs = IMAS.get_build_indexes(bd.layer, fs=_lfs_)[end]
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
        if obstr.fs in (Int(_lhfs_), Int(_out_))
            o_start = obstr.start_radius
            o_end = obstr.end_radius
        else
            o_start = obstr.start_radius
            o_end = IMAS.get_build_layer(bd.layer, identifier=obstr.identifier, fs=_lfs_).end_radius
        end
        l_start = layer.start_radius
        if layer.type == Int(_plasma_)
            l_end = layer.end_radius
        else
            l_end = IMAS.get_build_layer(bd.layer, identifier=layer.identifier, fs=_lfs_).end_radius
        end
    end
    hfs_thickness = o_start - l_start
    lfs_thickness = l_end - o_end
    oR = obstr.outline.r
    oZ = obstr.outline.z
    if layer.fs == Int(_out_)
        target_clearance = lfs_thickness * 1.2
        use_curvature = false
    else
        use_curvature = true
        if tight
            target_clearance = min(hfs_thickness, lfs_thickness)
        else
            target_clearance = sqrt(hfs_thickness^2 + lfs_thickness^2) / 2.0
        end
    end
    r_offset = (lfs_thickness .- hfs_thickness) / 2.0

    # update shape
    layer.shape = Int(shape)

    # handle offset, negative offset, and convex-hull
    if layer.shape in (Int(_offset_), Int(_negative_offset_), Int(_convex_hull_))
        R, Z = buffer(oR, oZ, (hfs_thickness + lfs_thickness) / 2.0)
        R .+= r_offset
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
        layer.shape_parameters = initialize_shape_parameters(layer.shape, oR, oZ, l_start, l_end, target_clearance)

        layer.outline.r, layer.outline.z = func(l_start, l_end, layer.shape_parameters...)
        layer.shape_parameters = optimize_shape(oR, oZ, target_clearance, func, l_start, l_end, layer.shape_parameters; use_curvature)
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
            layer.material = any((layer.type in (Int(_blanket_), Int(_shield_)) for layer in dd.build.layer)) ? "DT_plasma" : "DD_plasma"
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
