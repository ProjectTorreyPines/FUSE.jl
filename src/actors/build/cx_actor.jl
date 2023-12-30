#= ============= =#
#  cross-section  #
#= ============= =#
Base.@kwdef mutable struct FUSEparameters__ActorCXbuild{T} <: ParametersActor where {T<:Real}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    rebuild_wall::Entry{Bool} = Entry{Bool}("-", "Rebuild wall based on equilibrium"; default=true)
    n_points::Entry{Int} = Entry{Int}("-", "Number of points used for cross-sectional outlines"; default=101)
    do_plot::Entry{Bool} = Entry{Bool}("-", "Plot"; default=false)
end

mutable struct ActorCXbuild{D,P} <: ReactorAbstractActor{D,P}
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
    if isempty(wall.r) || par.rebuild_wall
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

    build_cx!(bd, pr, pz; par.n_points)

    divertor_regions!(bd, eqt, dd.divertors)

    blanket_regions!(bd, eqt)

    if isempty(wall.r) || par.rebuild_wall
        plasma = IMAS.get_build_layer(bd.layer; type=_plasma_)
        resize!(dd.wall.description_2d, 1)
        resize!(dd.wall.description_2d[1].limiter.unit, 1)
        dd.wall.description_2d[1].limiter.unit[1].outline.r = plasma.outline.r
        dd.wall.description_2d[1].limiter.unit[1].outline.z = plasma.outline.z
    end

    IMAS.find_strike_points!(eqt, dd.divertors)

    return actor
end

"""
    segmented_wall(eq_r::AbstractVector{T}, eq_z::AbstractVector{T}, gap::T, max_segment_relative_error::T; symmetric::Union{Bool,Nothing}=nothing) where {T<:Real}

Generate a segmented first wall outline starting from an equilibrium boundary

If `symmetric === nothing` then symmetry of the plasma is automatically determined
"""
function segmented_wall(eq_r::AbstractVector{T}, eq_z::AbstractVector{T}, gap::T, max_segment_relative_error::T; symmetric::Union{Bool,Nothing}=nothing) where {T<:Real}
    mxh = IMAS.MXH(eq_r, eq_z, 4)

    # automatic symmetry detection
    if symmetric === nothing
        symmetric = sum(abs.(mxh.c)) < 1E-3
    end

    # R,Z from theta
    Θ = LinRange(-π / 2 * 3 / 2, π / 2 * 3 / 2, 100) # overshoot
    R = [r for (r, z) in mxh.(Θ)]
    Z = [z for (r, z) in mxh.(Θ)]
    Z = (Z .- mxh.Z0) .* 1.1 .+ mxh.Z0

    # correct hfs to have a flat center stack wall
    R += IMAS.interp1d([Θ[1], Θ[argmax(Z)], Θ[argmin(Z)], Θ[end]], [minimum(eq_r) - R[1], 0.0, 0.0, minimum(eq_r) - R[end]]).(Θ)

    # flat midplane wall
    R, Z = buffer(R, Z, gap * (1.0 + max_segment_relative_error / 2.0))
    R[R.>(maximum(eq_r)+gap)] .= maximum(eq_r) + gap
    R, Z = buffer(R, Z, -gap * (1.0 + max_segment_relative_error / 2.0))
    circshift!(Z, argmin(R))
    circshift!(R, argmin(R))

    # segments
    R, Z = IMAS.rdp_simplify_2d_path(R, Z, gap * max_segment_relative_error)
    if symmetric
        R = (R .+ reverse(R)) / 2.0
        Z = ((Z .- mxh.Z0) .- reverse(Z .- mxh.Z0)) / 2.0 .+ mxh.Z0
    end

    # rounded joints
    R, Z = buffer(R, Z, gap * (1.0 + max_segment_relative_error))

    return R, Z
end

"""
    wall_from_eq(bd::IMAS.build, eqt::IMAS.equilibrium__time_slice; max_divertor_length_fraction_z_plasma::Real=0.2)

Generate first wall and divertors outline starting from an equilibrium
"""
function wall_from_eq(bd::IMAS.build, eqt::IMAS.equilibrium__time_slice; max_divertor_length_fraction_z_plasma::Real=0.2)
    R0 = eqt.global_quantities.magnetic_axis.r
    Z0 = eqt.global_quantities.magnetic_axis.z

    # lcfs
    ψb = IMAS.find_psi_boundary(eqt)
    ((rlcfs, zlcfs),), _ = IMAS.flux_surface(eqt, ψb, :closed)

    # Set the radial build thickness of the plasma
    plasma = IMAS.get_build_layer(bd.layer; type=_plasma_)
    gap = (minimum(rlcfs) - plasma.start_radius)
    plasma.thickness = maximum(rlcfs) - minimum(rlcfs) + 2.0 * gap
    R_hfs_plasma = plasma.start_radius
    R_lfs_plasma = plasma.end_radius

    # main chamber (clip elements that go beyond plasma radial build thickness)
    R, Z = segmented_wall(rlcfs, zlcfs, gap, 0.75)
    wall_poly = xy_polygon(R, Z)

    t = LinRange(0, 2π, 31)

    # divertor lengths
    linear_z_plasma_size = maximum(zlcfs) - minimum(zlcfs)
    max_divertor_length = linear_z_plasma_size * max_divertor_length_fraction_z_plasma

    detected_upper = bd.divertors.upper.installed
    detected_lower = bd.divertors.lower.installed

    # private flux regions sorted by distance from lcfs
    private, _ = IMAS.flux_surface(eqt, ψb, :open)
    private = collect(private)
    sort!(private; by=p -> IMAS.minimum_distance_two_shapes(p..., rlcfs, zlcfs))

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
        if d > linear_z_plasma_size / 5
            continue
        end

        # check that this divertor is in fact installed
        if Zx > Z0 && (ismissing(bd.divertors.upper, :installed) || bd.divertors.upper.installed != 1)
            continue
        elseif Zx < Z0 && (ismissing(bd.divertors.lower, :installed) || bd.divertors.lower.installed != 1)
            continue
        end

        # distance from points in private flux region to X-point
        Dx = sqrt.((Rx .- pr) .^ 2.0 .+ (Zx .- pz) .^ 2.0)
        divertor_length = minimum((Dx[1], Dx[end], max_divertor_length))

        # limit extent of private flux regions
        circle_r = divertor_length .* cos.(t) .+ Rx
        circle_z = sign(Zx) .* divertor_length .* sin.(t) .+ Zx
        circle = collect(zip(circle_r, circle_z))
        circle[1] = circle[end]
        inner_slot = [(rr, zz) for (rr, zz) in zip(pr, pz) if PolygonOps.inpolygon((rr, zz), circle) == 1]
        pr1 = [rr for (rr, zz) in inner_slot]
        pz1 = [zz for (rr, zz) in inner_slot]
        if isempty(pr1)
            continue
        end
        inner_slot_poly = xy_polygon(convex_hull(pr, pz; closed_polygon=true))

        # do not add more than one private flux region for each of the x-points
        if Zx > Z0
            if detected_upper == 0
                continue
            end
            detected_upper -= 1
        else
            if detected_lower == 0
                continue
            end
            detected_lower -= 1
        end

        # remove private flux region from wall (necessary because of Z expansion)
        # this may fail if the private region is small or weirdly shaped
        try
            wall_poly = LibGEOS.difference(wall_poly, inner_slot_poly)
        catch e
            if typeof(e) <: LibGEOS.GEOSError
                continue
            else
                rethrow(e)
            end
        end

        # add the divertor slots
        α = 0.25
        pr2 = vcat(pr1, R0 * α + Rx * (1 - α))
        pz2 = vcat(pz1, Z0 * α + Zx * (1 - α))

        slot_convhull = xy_polygon(convex_hull(pr2, pz2; closed_polygon=true))
        slot = LibGEOS.difference(slot_convhull, inner_slot_poly)
        slot = LibGEOS.buffer(slot, gap)

        scale = 1.001
        Rc1, Zc1 = IMAS.centroid(pr1, pz1)
        pr3 = vcat((pr1 .- Rc1) .* scale .+ Rc1, reverse(pr1))
        pz3 = vcat((pz1 .- Zc1) .* scale .+ Zc1, reverse(pz1))
        slot = LibGEOS.union(slot, LibGEOS.buffer(xy_polygon(pr3, pz3), gap))

        wall_poly = LibGEOS.union(wall_poly, slot)
    end

    # detect if equilibrium has x-points to define build of divertors
    if detected_upper != 0 || detected_lower != 0
        # plot(wall_poly)
        # display(plot!(eqt; cx=true, show_x_points=true))
        @warn(
            "Equilibrium does not allow building the right number of upper ($(bd.divertors.upper.installed)→$(-detected_upper+bd.divertors.upper.installed)) and lower ($(bd.divertors.lower.installed)→$(-detected_lower+bd.divertors.lower.installed)) divertors."
        )
    end

    # round corners
    corner_radius = gap / 4
    wall_poly = LibGEOS.buffer(wall_poly, -corner_radius)
    wall_poly = LibGEOS.buffer(wall_poly, corner_radius)

    # vertical clip
    wall_poly = LibGEOS.difference(wall_poly, xy_polygon(rectangle_shape(0.0, R_hfs_plasma, 100.0)...))
    wall_poly = LibGEOS.difference(wall_poly, xy_polygon(rectangle_shape(R_lfs_plasma, 10 * R_lfs_plasma, 100.0)...))

    pr = [v[1] for v in GeoInterface.coordinates(wall_poly)[1]]
    pz = [v[2] for v in GeoInterface.coordinates(wall_poly)[1]]

    try
        pr, pz = IMAS.resample_2d_path(pr, pz; step=0.1, method=:linear)
    catch e
        pp = plot(wall_poly; aspect_ratio=:equal)
        for (pr, pz) in private
            plot!(pp, pr, pz; label="")
        end
        display(pp)
        rethrow(e)
    end

    return pr, pz
end

function divertor_regions!(bd::IMAS.build, eqt::IMAS.equilibrium__time_slice, divertors::IMAS.divertors)
    R0 = eqt.global_quantities.magnetic_axis.r
    Z0 = eqt.global_quantities.magnetic_axis.z

    # plasma poly (this sets the divertor's pfs)
    ipl = IMAS.get_build_index(bd.layer; type=_plasma_)
    plasma_poly = xy_polygon(bd.layer[ipl])
    pl_r = bd.layer[ipl].outline.r
    pl_z = bd.layer[ipl].outline.z
    pl_r[1] = pl_r[end]
    pl_z[1] = pl_z[end]
    IMAS.reorder_flux_surface!(pl_r, pl_z, R0, Z0)

    # wall poly (this sets how back the divertor structure goes)
    wall_poly = xy_polygon(bd.layer[ipl-1])
    for ltype in (_blanket_, _shield_, _wall_)
        iwls = IMAS.get_build_indexes(bd.layer; type=ltype, fs=_hfs_)
        if !isempty(iwls)
            wall_poly = xy_polygon(bd.layer[iwls[1]])
            break
        end
    end
    wall_rz = [v for v in GeoInterface.coordinates(wall_poly)[1]]

    ψb = IMAS.find_psi_boundary(eqt)
    ((rlcfs, zlcfs),), _ = IMAS.flux_surface(eqt, ψb, :closed)
    linear_plasma_size = maximum(zlcfs) - minimum(zlcfs)

    empty!(divertors)

    detected_upper = bd.divertors.upper.installed
    detected_lower = bd.divertors.lower.installed

    private, _ = IMAS.flux_surface(eqt, ψb, :open)
    private = collect(private)
    sort!(private; by=p -> IMAS.minimum_distance_two_shapes(p..., rlcfs, zlcfs))
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

        # private flux region must be in the wall
        prz = [(rr, zz) for (rr, zz) in zip(pr, pz) if PolygonOps.inpolygon((rr, zz), wall_rz) == 1]
        pr = [rr for (rr, zz) in prz]
        pz = [zz for (rr, zz) in prz]
        if length(pr) < 5
            continue
        end

        if Zx > Z0
            if detected_upper == 0
                continue
            end
            detected_upper -= 1
        else
            if detected_lower == 0
                continue
            end
            detected_lower -= 1
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
        wall_domain_poly = LibGEOS.intersection(wall_poly, domain_poly)
        divertor_poly = LibGEOS.difference(wall_domain_poly, plasma_poly)

        # Assign to build structure
        coords = GeoInterface.coordinates(divertor_poly)
        structure = resize!(bd.structure, "type" => Int(_divertor_), "name" => "$ul_name divertor")
        structure.material = "Tungsten"
        try
            structure.outline.r = [v[1] for v in coords[1]]
            structure.outline.z = [v[2] for v in coords[1]]
            structure.toroidal_extent = 2π
        catch e
            p = plot(wall_poly; alpha=0.3, aspect_ratio=:equal)
            plot!(domain_poly; alpha=0.3)
            plot!(plasma_poly; alpha=0.3)
            display(p)
            rethrow(e)
        end

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

function IMAS.resample_2d_path(layer::IMAS.build__layer; method::Symbol=:linear, kw...)
    layer.outline.r, layer.outline.z = IMAS.resample_2d_path(layer.outline.r, layer.outline.z; method, kw...)
    return layer
end

"""
    build_cx!(bd::IMAS.build, pr::Vector{Float64}, pz::Vector{Float64}; n_points::Int)

Translates 1D build to 2D cross-sections starting from R and Z coordinates of plasma first wall
"""
function build_cx!(bd::IMAS.build, pr::Vector{Float64}, pz::Vector{Float64}; n_points::Int)
    plasma = IMAS.get_build_layer(bd.layer; type=_plasma_)

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
    tf_to_plasma = IMAS.get_build_indexes(bd.layer; fs=_hfs_)
    plasma_to_tf = reverse(tf_to_plasma)
    for k in plasma_to_tf
        layer = bd.layer[k]
        layer_shape = BuildLayerShape(mod(mod(layer.shape, 1000), 100))
        @debug "$(layer.name) $(layer_shape)"
        layer.shape, layer.shape_parameters = optimize_shape(bd, k + 1, k, layer_shape; tight=!coils_inside, resolution=n_points / 201.0)
    end

    # resample
    for k in tf_to_plasma[1:end-1]
        layer = bd.layer[k]
        IMAS.resample_2d_path(layer; n_points, method=:linear)
    end

    # _in_
    TF = IMAS.get_build_layer(bd.layer; type=_tf_, fs=_hfs_)
    D = (minimum(plasma.outline.z) * 2 + minimum(TF.outline.z)) / 3.0
    U = (maximum(plasma.outline.z) * 2 + maximum(TF.outline.z)) / 3.0
    if coils_inside
        # generally the OH does not go higher than the PF coils
        D += 2.0 * TF.thickness
        U -= 2.0 * TF.thickness
    end
    for k in IMAS.get_build_indexes(bd.layer; fs=_in_)
        layer = bd.layer[k]
        L = layer.start_radius
        R = layer.end_radius
        layer.outline.r, layer.outline.z = rectangle_shape(L, R, D, U)
    end

    # _out_
    iout = IMAS.get_build_indexes(bd.layer; fs=_out_)
    if lowercase(bd.layer[iout[end]].name) == "cryostat"
        olfs = IMAS.get_build_indexes(bd.layer; fs=_lfs_)[end]
        optimize_shape(bd, olfs, iout[end], BuildLayerShape(mod(mod(bd.layer[iout[end]].shape, 1000), 100)))
        for k in reverse(iout[2:end])
            optimize_shape(bd, k, k - 1, _negative_offset_)
        end
    else
        for k in iout
            layer = bd.layer[k]
            L = 0.0
            R = layer.end_radius
            D = minimum(bd.layer[k-1].outline.z) - layer.thickness
            U = maximum(bd.layer[k-1].outline.z) + layer.thickness
            layer.outline.r, layer.outline.z = rectangle_shape(L, R, D, U)
        end
    end

    return bd
end

"""
    optimize_shape(bd::IMAS.build, obstr_index::Int, layer_index::Int, shape::BuildLayerShape; tight::Bool=true, resolution::Float64=1.0)

Generates outline of layer in such a way to maintain minimum distance from inner layer
"""
function optimize_shape(bd::IMAS.build, obstr_index::Int, layer_index::Int, shape_enum::BuildLayerShape; tight::Bool=true, resolution::Float64=1.0)
    shape = Int(shape_enum)

    layer = bd.layer[layer_index]
    obstr = bd.layer[obstr_index]
    # display("Layer $layer_index = $(layer.name)")
    # display("Obstr $obstr_index = $(obstr.name)")

    if layer.side == Int(_out_)
        l_start = 0.0
        l_end = layer.end_radius
        o_start = 0.0
        o_end = obstr.end_radius
    else
        if obstr.side in (Int(_lhfs_), Int(_out_))
            o_start = obstr.start_radius
            o_end = obstr.end_radius
        else
            o_start = obstr.start_radius
            o_end = IMAS.get_build_layer(bd.layer; identifier=obstr.identifier, fs=_lfs_).end_radius
        end
        l_start = layer.start_radius
        if layer.type == Int(_plasma_)
            l_end = layer.end_radius
        else
            l_end = IMAS.get_build_layer(bd.layer; identifier=layer.identifier, fs=_lfs_).end_radius
        end
    end

    hfs_thickness = o_start - l_start
    lfs_thickness = l_end - o_end

    oR = obstr.outline.r
    oZ = obstr.outline.z

    if layer.side == Int(_out_)
        target_clearance = lfs_thickness * 1.2
        use_curvature = false
    else
        use_curvature = shape_enum == IMAS._rectangle_ ? false : true
        if tight
            target_clearance = min(hfs_thickness, lfs_thickness)
        else
            target_clearance = sqrt(hfs_thickness^2 + lfs_thickness^2) / 2.0
        end
    end
    r_offset = (lfs_thickness .- hfs_thickness) / 2.0

    # handle offset, negative offset, and convex-hull
    if shape in (Int(_offset_), Int(_negative_offset_), Int(_convex_hull_))
        R, Z = buffer(oR, oZ, (hfs_thickness + lfs_thickness) / 2.0)
        R .+= r_offset
        if shape == Int(_convex_hull_)
            hull = convex_hull(R, Z; closed_polygon=true)
            R = [r for (r, z) in hull]
            Z = [z for (r, z) in hull]
        end
        layer.outline.r, layer.outline.z = R, Z
        shape_parameters = Float64[]

    else # handle shapes
        if shape > 1000
            shape = mod(shape, 1000)
        end
        if shape > 100
            shape = mod(shape, 100)
        end

        if shape == Int(_silo_)
            is_up_down_symmetric = false
        elseif abs(sum(oZ) / sum(abs, oZ)) < 1E-2
            is_up_down_symmetric = true
        else
            is_up_down_symmetric = false
        end

        is_negative_D = false
        if shape != Int(_silo_)
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
            shape = shape + 1000
        end

        if !is_up_down_symmetric
            shape = shape + 100
        end

        func = shape_function(shape; resolution)
        shape_parameters = initialize_shape_parameters(shape, oR, oZ, l_start, l_end, target_clearance)

        layer.outline.r, layer.outline.z = func(l_start, l_end, shape_parameters...)
        shape_parameters = optimize_shape(oR, oZ, target_clearance, func, l_start, l_end, shape_parameters; use_curvature)
        layer.outline.r, layer.outline.z = func(l_start, l_end, shape_parameters...; resample=false)
    end

    IMAS.reorder_flux_surface!(layer.outline.r, layer.outline.z)

    return Int(shape_enum), shape_parameters
end

