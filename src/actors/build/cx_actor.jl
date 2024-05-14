import LibGEOS
import GeoInterface

#= ============= =#
#  cross-section  #
#= ============= =#
Base.@kwdef mutable struct FUSEparameters__ActorCXbuild{T<:Real} <: ParametersActorBuild{T}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    _time::Float64 = NaN
    rebuild_wall::Entry{Bool} = Entry{Bool}("-", "Rebuild wall based on equilibrium"; default=true)
    n_points::Entry{Int} = Entry{Int}("-", "Number of points used for cross-sectional outlines"; default=101)
    do_plot::Entry{Bool} = act_common_parameters(; do_plot=false)
end

mutable struct ActorCXbuild{D,P} <: CompoundAbstractActor{D,P}
    dd::IMAS.dd{D}
    par::FUSEparameters__ActorCXbuild{P}
    act::ParametersAllActors
    function ActorCXbuild(dd::IMAS.dd{D}, par::FUSEparameters__ActorCXbuild{P}, act::ParametersAllActors{P}; kw...) where {D<:Real,P<:Real}
        logging_actor_init(ActorCXbuild)
        par = par(kw...)
        return new{D,P}(dd, par, act)
    end
end

"""
    ActorCXbuild(dd::IMAS.dd, act::ParametersAllActors; kw...)

Generates the 2D cross section of the tokamak build

!!! note

    Manipulates data in `dd.build`
"""
function ActorCXbuild(dd::IMAS.dd, act::ParametersAllActors; kw...)
    actor = ActorCXbuild(dd, act.ActorCXbuild, act; kw...)
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
    act = actor.act

    bd = dd.build
    eqt = dd.equilibrium.time_slice[]
    wall = dd.wall

    # If wall information is missing, then the first wall information is generated starting from equilibrium time_slice
    wall_outline = IMAS.first_wall(wall)
    if isempty(wall_outline.r) || par.rebuild_wall
        wall_from_eq!(wall, eqt, bd)
    end

    # empty layer outlines and structures
    for layer in bd.layer
        empty!(layer.outline)
    end
    empty!(bd.structure)

    build_cx!(bd, wall, dd.pf_active; par.n_points)

    divertor_regions!(bd, eqt, dd.divertors)

    blanket_regions!(bd, eqt)

    if act.ActorLFSsizing.maintenance != :none
        port_regions!(bd; act.ActorLFSsizing.maintenance, act.ActorLFSsizing.tor_modularity, act.ActorLFSsizing.pol_modularity)
    end

    IMAS.find_strike_points!(eqt, dd.divertors)

    return actor
end

"""
    segmented_wall(eq_r::AbstractVector{T}, eq_z::AbstractVector{T}, gap::T, max_segment_relative_error::T; symmetric::Union{Bool,Nothing}=nothing) where {T<:Real}

Generate a segmented first wall outline starting from an equilibrium boundary.

If `max_segment_relative_error == 0.0` then no segmentation is applied.

If `symmetric === nothing` then symmetry of the plasma is automatically determined.
"""
function segmented_wall(eq_r::AbstractVector{T}, eq_z::AbstractVector{T}, gap::T, max_segment_relative_error::T; symmetric::Union{Bool,Nothing}=nothing) where {T<:Real}
    if max_segment_relative_error < 0.0
        R, Z = buffer(eq_r, eq_z, gap)

    else
        mxh = IMAS.MXH(eq_r, eq_z, 4)

        # automatic symmetry detection
        if symmetric === nothing
            symmetric = sum(abs.(mxh.c)) / length(mxh.c) < 1E-3
        end

        # negative δ
        δ = sin(mxh.s[1])
        if δ < 0.0
            Θ = LinRange(π / 4, 2 * π - π / 4, 100) # overshoot
        else
            Θ = LinRange(-π * 3 / 4, π * 3 / 4, 100) # overshoot
        end

        # R,Z from theta
        R = [r for (r, z) in mxh.(Θ)]
        Z = [z for (r, z) in mxh.(Θ)]
        Z = (Z .- mxh.Z0) .* 1.1 .+ mxh.Z0

        # correct hfs to have a flat center stack wall
        if δ < 0.0
            R += IMAS.interp1d([Θ[1], Θ[argmax(Z)], Θ[argmin(Z)], Θ[end]], [maximum(eq_r) - R[1], 0.0, 0.0, maximum(eq_r) - R[end]]).(Θ)
        else
            R += IMAS.interp1d([Θ[1], Θ[argmax(Z)], Θ[argmin(Z)], Θ[end]], [minimum(eq_r) - R[1], 0.0, 0.0, minimum(eq_r) - R[end]]).(Θ)
        end

        if max_segment_relative_error == 0.0
            R, Z = buffer(R, Z, gap)

        else
            if symmetric # this works because points are distributed like Θ
                R = (R .+ reverse(R)) / 2.0
                Z = ((Z .- mxh.Z0) .- reverse(Z .- mxh.Z0)) / 2.0 .+ mxh.Z0
            end

            # segments
            R, Z = IMAS.rdp_simplify_2d_path(R, Z, gap * max_segment_relative_error)

            # rounded joints
            R, Z = buffer(R, Z, gap * (1.0 + max_segment_relative_error))
        end
    end

    return R, Z
end

"""
    wall_from_eq!(wall::IMAS.wall, eqt::IMAS.equilibrium__time_slice, bd::IMAS.build; max_divertor_length_fraction_z_plasma::Real=0.15)

Generate first wall and divertors outline starting from an equilibrium and radial build
"""
function wall_from_eq!(wall::IMAS.wall, eqt::IMAS.equilibrium__time_slice, bd::IMAS.build; max_divertor_length_fraction_z_plasma::Real=0.15)
    # Set the radial build thickness of the plasma vacuum chamber
    plasma = IMAS.get_build_layer(bd.layer; type=_plasma_)
    rlcfs, zlcfs = eqt.boundary.outline.r, eqt.boundary.outline.z
    gap = (minimum(rlcfs) - plasma.start_radius)
    plasma.thickness = maximum(rlcfs) - minimum(rlcfs) + 2.0 * gap

    upper_divertor = ismissing(bd.divertors.upper, :installed) ? false : Bool(bd.divertors.upper.installed)
    lower_divertor = ismissing(bd.divertors.lower, :installed) ? false : Bool(bd.divertors.lower.installed)

    return wall_from_eq!(wall, eqt, gap; upper_divertor, lower_divertor, max_divertor_length_fraction_z_plasma)
end

function wall_from_eq!(
    wall::IMAS.wall,
    eqt::IMAS.equilibrium__time_slice,
    gap::Float64;
    upper_divertor::Bool,
    lower_divertor::Bool,
    max_divertor_length_fraction_z_plasma::Real)

    upper_divertor = Int(upper_divertor)
    lower_divertor = Int(lower_divertor)

    # domain of the equilibrium includes some buffer for divertor slots
    div_gap = gap
    eqt2d = findfirst(:rectangular, eqt.profiles_2d)
    req = eqt2d.grid.dim1
    req = [req[1] + div_gap, req[end] - div_gap, req[end] - div_gap, req[1] + div_gap, req[1] + div_gap]
    zeq = eqt2d.grid.dim2
    zeq = [zeq[1] + div_gap, zeq[1] + div_gap, zeq[end] - div_gap, zeq[end] - div_gap, zeq[1] + div_gap]
    eq_domain = collect(zip(req, zeq))

    # lcfs and magnetic axis
    ((rlcfs, zlcfs),), ψb = IMAS.flux_surface(eqt, eqt.global_quantities.psi_boundary, :closed)
    RA = eqt.global_quantities.magnetic_axis.r
    ZA = eqt.global_quantities.magnetic_axis.z

    # Set the radial build thickness of the plasma
    R_hfs_plasma = minimum(rlcfs) - gap
    R_lfs_plasma = maximum(rlcfs) + gap

    # main chamber (clip elements that go beyond plasma radial build thickness)
    R, Z = segmented_wall(rlcfs, zlcfs, gap, 0.)
    wall_poly = xy_polygon(R, Z)

    # divertor lengths
    linear_z_plasma_size = maximum(zlcfs) - minimum(zlcfs)
    max_divertor_length = linear_z_plasma_size * max_divertor_length_fraction_z_plasma

    # private flux regions sorted by distance from lcfs
    private, _ = IMAS.flux_surface(eqt, ψb, :open)
    private = collect(private)
    sort!(private; by=p -> IMAS.minimum_distance_two_shapes(p..., rlcfs, zlcfs))

    t = LinRange(0, 2π, 31)
    detected_upper = upper_divertor
    detected_lower = lower_divertor
    for (pr, pz) in private
        if sign(pz[1] - ZA) != sign(pz[end] - ZA)
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
        if Zx > ZA && upper_divertor == 0
            continue
        elseif Zx < ZA && lower_divertor == 0
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
        inner_slot = [(rr, zz) for (rr, zz) in zip(pr, pz) if PolygonOps.inpolygon((rr, zz), circle) == 1 && PolygonOps.inpolygon((rr, zz), eq_domain) == 1]
        pr1 = [rr for (rr, zz) in inner_slot]
        pz1 = [zz for (rr, zz) in inner_slot]
        if isempty(pr1)
            continue
        end
        inner_slot_poly = xy_polygon(convex_hull(pr, pz; closed_polygon=true))

        # do not add more than one private flux region for each of the x-points
        if Zx > ZA
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
        pr2 = vcat(pr1, RA * α + Rx * (1 - α))
        pz2 = vcat(pz1, ZA * α + Zx * (1 - α))

        slot_convhull = xy_polygon(convex_hull(pr2, pz2; closed_polygon=true))
        slot = LibGEOS.difference(slot_convhull, inner_slot_poly)
        slot = LibGEOS.buffer(slot, div_gap)

        scale = 1.00
        Rc1, Zc1 = IMAS.centroid(pr1, pz1)
        pr3 = vcat((pr1 .- Rc1) .* scale .+ Rc1, reverse(pr1))
        pz3 = vcat((pz1 .- Zc1) .* scale .+ Zc1, reverse(pz1))
        slot = LibGEOS.union(slot, LibGEOS.buffer(xy_polygon(pr3, pz3), div_gap))

        wall_poly = LibGEOS.union(wall_poly, slot)
    end

    # detect if equilibrium has x-points to define build of divertors
    if detected_upper != 0 || detected_lower != 0
        # plot(wall_poly)
        # display(plot!(eqt; cx=true, show_x_points=true))
        @warn(
            "Equilibrium does not allow building the right number of upper ($(upper_divertor)→$(-detected_upper+upper_divertor)) and lower ($(lower_divertor)→$(-detected_lower+lower_divertor)) divertors."
        )
    end

    # round corners
    corner_radius = div_gap / 4
    wall_poly = LibGEOS.buffer(wall_poly, -corner_radius)
    wall_poly = LibGEOS.buffer(wall_poly, corner_radius)

    # vertical clip
    wall_poly = LibGEOS.difference(wall_poly, xy_polygon(rectangle_shape(0.0, R_hfs_plasma, 100.0)...))
    wall_poly = LibGEOS.difference(wall_poly, xy_polygon(rectangle_shape(R_lfs_plasma, 10 * R_lfs_plasma, 100.0)...))

    pr = [v[1] for v in GeoInterface.coordinates(wall_poly)[1]]
    pz = [v[2] for v in GeoInterface.coordinates(wall_poly)[1]]

    try
        pr, pz = IMAS.resample_2d_path(pr, pz; n_points=101, method=:linear)
    catch e
        pp = plot(wall_poly; aspect_ratio=:equal)
        for (pr, pz) in private
            plot!(pp, pr, pz; label="")
        end
        display(pp)
        rethrow(e)
    end

    # update the wall IDS
    resize!(wall.description_2d, 1)
    resize!(wall.description_2d[1].limiter.unit, 1)
    wall.description_2d[1].limiter.unit[1].outline.r = pr
    wall.description_2d[1].limiter.unit[1].outline.z = pz

    return wall
end

function divertor_regions!(bd::IMAS.build, eqt::IMAS.equilibrium__time_slice, divertors::IMAS.divertors)
    RA = eqt.global_quantities.magnetic_axis.r
    ZA = eqt.global_quantities.magnetic_axis.z

    # plasma poly (this sets the divertor's pfs)
    ipl = IMAS.get_build_index(bd.layer; type=_plasma_)
    plasma_poly = xy_polygon(bd.layer[ipl])
    pl_r = bd.layer[ipl].outline.r
    pl_z = bd.layer[ipl].outline.z
    pl_r[1] = pl_r[end]
    pl_z[1] = pl_z[end]
    IMAS.reorder_flux_surface!(pl_r, pl_z, RA, ZA)

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

    ψb, _ = IMAS.find_psi_boundary(eqt)
    ((rlcfs, zlcfs),), _ = IMAS.flux_surface(eqt, ψb, :closed)
    linear_plasma_size = maximum(zlcfs) - minimum(zlcfs)

    empty!(divertors)

    detected_upper = bd.divertors.upper.installed
    detected_lower = bd.divertors.lower.installed

    private, _ = IMAS.flux_surface(eqt, ψb, :open)
    private = collect(private)
    sort!(private; by=p -> IMAS.minimum_distance_two_shapes(p..., rlcfs, zlcfs))
    for (pr, pz) in private
        if sign(pz[1] - ZA) != sign(pz[end] - ZA)
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
        if Zx > ZA && (ismissing(bd.divertors.upper, :installed) || bd.divertors.upper.installed != 1)
            continue
        elseif Zx < ZA && (ismissing(bd.divertors.lower, :installed) || bd.divertors.lower.installed != 1)
            continue
        end

        # private flux region must be in the wall
        prz = [(rr, zz) for (rr, zz) in zip(pr, pz) if PolygonOps.inpolygon((rr, zz), wall_rz) == 1]
        pr = [rr for (rr, zz) in prz]
        pz = [zz for (rr, zz) in prz]
        if length(pr) < 5
            continue
        end

        if Zx > ZA
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

        if Zx > ZA
            ul_name = "upper"
        else
            ul_name = "lower"
        end

        # Define divertor region by tracing a line going through X-point
        # and perpendicular to line connecting X-point to magnetic axis
        m = (Zx - ZA) / (Rx - RA)
        xx = [0.0, RA * 2.0]
        yy = line_through_point(-1.0 ./ m, Rx, Zx, xx)
        domain_r = vcat(xx, reverse(xx), xx[1])
        domain_z = vcat(yy, [Zx * 5.0, Zx * 5.0], yy[1])
        domain_poly = xy_polygon(domain_r, domain_z)
        wall_domain_poly = LibGEOS.intersection(wall_poly, domain_poly)
        divertor_poly = LibGEOS.difference(wall_domain_poly, plasma_poly)

        # Assign to build structure
        coords = GeoInterface.coordinates(divertor_poly)
        structure = resize!(bd.structure, "type" => Int(_divertor_), "name" => "$ul_name divertor")
        structure.material = "tungsten"
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
        m = (Zx - ZA) / (Rx - RA)
        xx = [0.0, RA * 2.0]
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
    RA = eqt.global_quantities.magnetic_axis.r

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
            if sum(pr) / length(pr) > RA
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
    port_regions!(bd::IMAS.build; maintenance::Symbol=:none, tor_modularity::Int=2, pol_modularity::Int=1)

Defines maintenance port structures
"""
function port_regions!(bd::IMAS.build; maintenance::Symbol=:none, tor_modularity::Int=2, pol_modularity::Int=1)
    vessel = IMAS.get_build_layer(bd.layer; type=_vessel_, fs=_lfs_)
    cryostat = IMAS.get_build_layer(bd.layer; type=_cryostat_)

    # calculate radial port locations
    vessel_maxz = maximum(vessel.outline.z)
    vessel_minz = minimum(vessel.outline.z)
    horz_port_maxr = maximum(cryostat.outline.r) * 1.1
    horz_port_maxz = 0.25 * (vessel_maxz - vessel_minz) + vessel_minz
    vert_port_maxz = maximum(cryostat.outline.z)

    # outer segment of the vacuum vessel
    istart = argmin(vessel.outline.z)
    pr = circshift(vessel.outline.r, 1 - istart)
    pz = circshift(vessel.outline.z, 1 - istart)
    index = argmax(pz)
    vR = [pr[index:end]; pr[1]]
    vZ = [pz[index:end]; pz[1]]

    # build port structure based on maintenance type
    if maintenance == :none
        return nothing

    elseif maintenance == :vertical
        rVP_hfs_ib, rVP_hfs_ob, rVP_lfs_ib, rVP_lfs_ob = IMAS.vertical_maintenance(bd; tor_modularity, pol_modularity)
        pr = vcat(rVP_hfs_ib, rVP_hfs_ib, vR, vR[end], horz_port_maxr, horz_port_maxr, rVP_lfs_ob, rVP_lfs_ob, rVP_hfs_ib)
        pz = vcat(vert_port_maxz, vZ[1], vZ, vessel_minz, vessel_minz, horz_port_maxz, horz_port_maxz, vert_port_maxz, vert_port_maxz)
        name = "vertical maintenance port"

    elseif maintenance == :horizontal
        pr = vcat(vR, vR[end], horz_port_maxr, horz_port_maxr, vR[1])
        pz = vcat(vZ, vessel_minz, vessel_minz, vessel_maxz, vessel_maxz)
        name = "horizontal maintenance port"
    end

    structure = resize!(bd.structure, "type" => Int(_port_), "name" => name)
    structure.outline.r = pr
    structure.outline.z = pz

    return nothing
end

function IMAS.resample_2d_path(layer::IMAS.build__layer; method::Symbol=:linear, kw...)
    layer.outline.r, layer.outline.z = IMAS.resample_2d_path(layer.outline.r, layer.outline.z; method, kw...)
    return layer
end

"""
    build_cx!(bd::IMAS.build, wall::IMAS.wall, pfa::IMAS.pf_active; n_points::Int)

Translates 1D build to 2D cross-sections starting from R and Z coordinates of plasma first wall
"""
function build_cx!(bd::IMAS.build, wall::IMAS.wall, pfa::IMAS.pf_active; n_points::Int)
    plasma = IMAS.get_build_layer(bd.layer; type=_plasma_)

    wall_outline = IMAS.first_wall(wall)
    pr = wall_outline.r
    pz = wall_outline.z

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

    # up-down symmetric plasma
    is_z_offset = false
    if abs(sum(pz) / sum(abs, pz)) > 1E-2
        is_z_offset = true
    end

    # up-down symmetric plasma
    is_negative_D = false
    _, imaxr = findmax(pr)
    _, iminr = findmin(pr)
    _, imaxz = findmax(pz)
    _, iminz = findmin(pz)
    r_at_max_z, _ = pr[imaxz], pz[imaxz]
    r_at_min_z, _ = pr[iminz], pz[iminz]
    _, max_r = pz[imaxr], pr[imaxr]
    _, min_r = pz[iminr], pr[iminr]
    a = 0.5 * (max_r - min_r)
    R = 0.5 * (max_r + min_r)
    δu = (R - r_at_max_z) / a
    δl = (R - r_at_min_z) / a
    if δu + δl < -0.1
        is_negative_D = true
    end

    #plot()
    vertical_clearance = 1.0
    obstruction_outline = nothing
    plasma_to_tf = reverse(IMAS.get_build_indexes(bd.layer; fs=IMAS._hfs_))
    pushfirst!(plasma_to_tf, plasma_to_tf[1] + 1)

    # @show plasma_to_tf
    # plot(bd.layer[plasma_to_tf[1]])

    k = 2
    while k <= length(plasma_to_tf)
        layer = bd.layer[plasma_to_tf[k]]
        layer_shape = IMAS.BuildLayerShape(mod(mod(layer.shape, 1000), 100))
        if layer_shape == IMAS._negative_offset_
            # go forward until find a normal layer
            neg_off_layers = []
            for kk in k:length(plasma_to_tf)
                layer = bd.layer[plasma_to_tf[kk]]
                layer_shape = IMAS.BuildLayerShape(mod(mod(layer.shape, 1000), 100))
                if layer_shape == IMAS._negative_offset_
                    push!(neg_off_layers, kk)
                else
                    #@show "A", layer.name, layer_shape
                    layer.shape, layer.shape_parameters = FUSE.optimize_layer_outline(
                        bd,
                        plasma_to_tf[k-1],
                        plasma_to_tf[kk],
                        layer_shape;
                        vertical_clearance,
                        resolution=n_points / 201.0,
                        obstruction_outline,
                        is_negative_D,
                        is_z_offset
                    )
                    k = kk
                    #display(plot!(bd.layer[plasma_to_tf[k]]))
                    break
                end
            end
            # go back with negative offset layers
            for kk in reverse(neg_off_layers)
                layer = bd.layer[plasma_to_tf[kk]]
                layer_shape = IMAS.BuildLayerShape(mod(mod(layer.shape, 1000), 100))
                #                @show "C", layer.name, layer_shape
                layer.shape, layer.shape_parameters = FUSE.optimize_layer_outline(
                    bd,
                    plasma_to_tf[kk+1],
                    plasma_to_tf[kk],
                    layer_shape;
                    vertical_clearance,
                    resolution=n_points / 201.0,
                    obstruction_outline,
                    is_negative_D,
                    is_z_offset
                )
                #                display(plot!(bd.layer[plasma_to_tf[kk]]))
            end

        else
            #            @show "B", layer.name, layer_shape
            layer.shape, layer.shape_parameters = FUSE.optimize_layer_outline(
                bd,
                plasma_to_tf[k-1],
                plasma_to_tf[k],
                layer_shape;
                vertical_clearance,
                resolution=n_points / 201.0,
                obstruction_outline,
                is_negative_D,
                is_z_offset
            )
            #            display(plot!(bd.layer[plasma_to_tf[k]]))
        end

        k += 1
    end

    # _in_
    TF = IMAS.get_build_layer(bd.layer; type=_tf_, fs=_hfs_)
    D = minimum(TF.outline.z) + TF.thickness / 2.0
    U = maximum(TF.outline.z) - TF.thickness / 2.0
    for layer in IMAS.get_build_layers(bd.layer; fs=_in_)
        L = layer.start_radius
        R = layer.end_radius
        layer.outline.r, layer.outline.z = rectangle_shape(L, R, D, U)
    end

    # _out_
    iout = IMAS.get_build_indexes(bd.layer; fs=_out_)
    if lowercase(bd.layer[iout[end]].name) == "cryostat"
        olfs = IMAS.get_build_indexes(bd.layer; fs=_lfs_)[end]
        # if !isempty(pfa.coil)
        #     obstruction_outline = convex_outline(pfa.coil)
        # end
        optimize_layer_outline(bd, olfs, iout[end], BuildLayerShape(mod(mod(bd.layer[iout[end]].shape, 1000), 100)); vertical_clearance=1.3)
        for k in reverse(iout[2:end])
            optimize_layer_outline(bd, k, k - 1, _negative_offset_)
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

    # resample
    for layer in bd.layer
        if layer.type != Int(_plasma_) || getproperty(layer, :shape, _undefined_) != Int(_offset_)
            layer.outline.r, layer.outline.z = IMAS.rdp_simplify_2d_path(layer.outline.r, layer.outline.z, 1E-3)
            layer.outline.r, layer.outline.z = IMAS.split_long_segments(layer.outline.r, layer.outline.z, n_points)
        end
    end

    return bd
end

"""
    optimize_layer_outline(bd::IMAS.build, obstr_index::Int, layer_index::Int, shape_enum::BuildLayerShape; vertical_clearance::Float64=1.0, resolution::Float64=1.0, obstruction_outline=nothing)

Generates outline of layer in such a way to maintain minimum distance from inner layer
"""
function optimize_layer_outline(
    bd::IMAS.build,
    obstr_index::Int,
    layer_index::Int,
    shape_enum::BuildLayerShape;
    vertical_clearance::Float64=1.0,
    resolution::Float64=1.0,
    obstruction_outline=nothing,
    is_z_offset::Bool=false,
    is_negative_D::Bool=false
)
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
    if obstruction_outline !== nothing
        hull = convex_hull(vcat(oR, obstruction_outline.r), vcat(oZ, obstruction_outline.z); closed_polygon=true)
        oR = [r for (r, z) in hull]
        oZ = [z for (r, z) in hull]
    end

    # handle offset, negative offset, and convex-hull
    if shape in (Int(_offset_), Int(_negative_offset_), Int(_convex_hull_))
        R, Z = buffer(oR, oZ, hfs_thickness, lfs_thickness)
        if shape == Int(_convex_hull_)
            hull = convex_hull(R, Z; closed_polygon=true)
            R = [r for (r, z) in hull]
            Z = [z for (r, z) in hull]
        end
        layer.outline.r, layer.outline.z = R, Z
        shape_parameters = Float64[]

    else # handle shapes
        use_curvature = true
        if layer.side == Int(_out_) || shape_enum ∈ (IMAS._rectangle_, IMAS._silo_)
            use_curvature = false
            target_clearance = lfs_thickness
        else
            target_clearance = min(hfs_thickness, lfs_thickness)
        end

        target_clearance *= vertical_clearance

        shape = mod(mod(shape, 1000), 100)
        if is_negative_D
            shape = shape + 1000
        end
        if is_z_offset
            shape = shape + 100
        end

        func = shape_function(shape; resolution)
        shape_parameters = initialize_shape_parameters(shape, oR, oZ, l_start, l_end, target_clearance)

        layer.outline.r, layer.outline.z = func(l_start, l_end, shape_parameters...)
        shape_parameters = optimize_outline(oR, oZ, target_clearance, func, l_start, l_end, shape_parameters; use_curvature)
        layer.outline.r, layer.outline.z = func(l_start, l_end, shape_parameters...; resample=false)
    end

    IMAS.reorder_flux_surface!(layer.outline.r, layer.outline.z)

    return Int(shape_enum), shape_parameters
end

"""
    convex_outline(coils::AbstractVector{IMAS.pf_active__coil{T}})::IMAS.pf_active__coil___element___geometry__outline{T} where {T<:Real}

returns convex-hull outline of all the pf coils elements in the input array
"""
function convex_outline(coils::AbstractVector{IMAS.pf_active__coil{T}})::IMAS.pf_active__coil___element___geometry__outline{T} where {T<:Real}
    out = IMAS.pf_active__coil___element___geometry__outline()
    out.r = T[]
    out.z = T[]

    for coil in coils
        for element in coil.element
            oute = IMAS.outline(element)
            append!(out.r, oute.r)
            append!(out.z, oute.z)
        end
    end

    points = convex_hull(out.r, out.z; closed_polygon=true)
    empty!(out.z)
    empty!(out.r)
    out.r = [p[1] for p in points]
    out.z = [p[2] for p in points]

    return out
end

function convex_outline(rail::IMAS.build__pf_active__rail{T})::IMAS.pf_active__coil___element___geometry__outline{T} where {T<:Real}
    rails = parent(rail)
    irail = IMAS.index(rail)
    coil_start = sum([rails[k].coils_number for k in 1:irail-1])
    coil_end = coil_start + rail.coils_number
    dd = IMAS.top_dd(rail)
    coils = [dd.pf_active.coil[k] for k in coil_start:coil_end]
    return convex_outline(coils)
end

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
        volume=Float64[],
        shape=String[])

    for layer in layers
        group = replace(string(BuildLayerSide(layer.side)), "_" => "")
        type = replace(string(BuildLayerType(layer.type)), "_" => "")
        type = replace(type, r"^gap" => "")
        details = replace(lowercase(layer.name), r"^[hl]fs " => "", r"^gap .*" => "", r"\b" * type * r"\b" => "")
        material = getproperty(layer, :material, "?")
        material = split(material, ",")[1]
        material = replace(material, "Vacuum" => "")
        area = getproperty(layer, :area, NaN)
        volume = getproperty(layer, :volume, NaN)
        shape = string(IMAS.BuildLayerShape(mod(mod(getproperty(layer, :shape, Int(IMAS._undefined_)), 1000), 100)))
        shape = replace(shape, "_undefined_" => "", r"^_" => "", r"_$" => "", "_" => " ")
        push!(df, [group, details, type, layer.thickness, layer.start_radius, layer.end_radius, material, area, volume, shape])
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