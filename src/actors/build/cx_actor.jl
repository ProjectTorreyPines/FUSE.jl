import LibGEOS
import GeoInterface

#= ============= =#
#  cross-section  #
#= ============= =#
Base.@kwdef mutable struct FUSEparameters__ActorCXbuild{T<:Real} <: ParametersActor{T}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    _time::Float64 = NaN
    rebuild_wall::Entry{Bool} = Entry{Bool}("-", "Rebuild wall based on equilibrium, even if dd.wall is already filled"; default=true)
    layers_aware_of_pf_coils::Entry{Bool} = Entry{Bool}("-", "Build layers are aware of pf_active coils"; default=true)
    divertor_hfs_size_fraction::Entry{T} = Entry{T}(
        "-",
        "Divertor size on the high-field-side as fraction of plasma minor radius";
        default=0.50,
        check=x -> @assert x >= 0 "divertor_hfs_size_fraction must be >= 0.0")
    divertor_lfs_size_fraction::Entry{T} = Entry{T}(
        "-",
        "Divertor size on the low-field-side as fraction of plasma minor radius";
        default=0.50,
        check=x -> @assert x >= 0 "divertor_lfs_size_fraction must be >= 0.0")
    n_points::Entry{Int} = Entry{Int}("-", "Number of points used for cross-sectional outlines"; default=101, check=x -> @assert x >= 0 "n_points must be > 0")
    do_plot::Entry{Bool} = act_common_parameters(; do_plot=false)
end

mutable struct ActorCXbuild{D,P} <: CompoundAbstractActor{D,P}
    dd::IMAS.dd{D}
    par::FUSEparameters__ActorCXbuild{P}
    act::ParametersAllActors{P}
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

    # If wall information is missing, then the first wall information is generated starting from equilibrium time_slice
    fw = IMAS.first_wall(dd.wall)
    if isempty(fw.r) || par.rebuild_wall
        wall_from_eq!(dd.wall, eqt, bd, dd.pulse_schedule.position_control, dd.pf_active; par.divertor_hfs_size_fraction, par.divertor_lfs_size_fraction)
        fw = IMAS.first_wall(dd.wall)
        IMAS.flux_surfaces(eqt, fw.r, fw.z) # output of flux_surfaces depends on the wall
    end

    # empty layer outlines and structures
    for layer in bd.layer
        empty!(layer.outline)
    end
    empty!(bd.structure)

    # layers
    if par.layers_aware_of_pf_coils
        pfa = dd.pf_active
    else
        pfa = IMAS.pf_active()
    end
    build_cx!(bd, dd.wall, pfa; par.n_points)

    # divertors + find strike points on divertors
    divertor_regions!(bd, eqt, dd.divertors, fw.r, fw.z)
    psi_first_open = IMAS.find_psi_boundary(eqt, fw.r, fw.z; raise_error_on_not_open=true).first_open
    IMAS.find_strike_points!(eqt, fw.r, fw.z, psi_first_open, dd.divertors)

    # blankets
    blanket_regions!(bd, eqt)

    # maintenance ports
    if act.ActorLFSsizing.maintenance != :none
        port_regions!(bd; act.ActorLFSsizing.maintenance, act.ActorLFSsizing.tor_modularity, act.ActorLFSsizing.pol_modularity)
    end

    return actor
end

"""
wall_from_eq!(wall::IMAS.wall, eqt::IMAS.equilibrium__time_slice, bd::IMAS.build, pc::IMAS.pulse_schedule__position_control, pfa::IMAS.pf_active; divertor_hfs_size_fraction::Float64, divertor_lfs_size_fraction::Float64)

Generate first wall and divertors outline starting from an equilibrium and radial build
"""
function wall_from_eq!(wall::IMAS.wall, eqt::IMAS.equilibrium__time_slice, bd::IMAS.build, pc::IMAS.pulse_schedule__position_control, pfa::IMAS.pf_active; divertor_hfs_size_fraction::Float64, divertor_lfs_size_fraction::Float64)
    # Set the radial build thickness of the plasma vacuum chamber
    plasma = IMAS.get_build_layer(bd.layer; type=_plasma_)
    rlcfs = eqt.boundary.outline.r
    gap = (minimum(rlcfs) - plasma.start_radius)
    #plasma.thickness = maximum(rlcfs) - minimum(rlcfs) + 2.0 * gap

    upper_divertor = ismissing(bd.divertors.upper, :installed) ? false : Bool(bd.divertors.upper.installed)
    lower_divertor = ismissing(bd.divertors.lower, :installed) ? false : Bool(bd.divertors.lower.installed)

    return wall_from_eq!(wall, eqt, pc, pfa, gap; upper_divertor, lower_divertor, divertor_hfs_size_fraction, divertor_lfs_size_fraction)
end

function wall_from_eq!(
    wall::IMAS.wall{T},
    eqt::IMAS.equilibrium__time_slice{T},
    pc::IMAS.pulse_schedule__position_control,
    pfa::IMAS.pf_active{T},
    gap::Float64;
    upper_divertor::Bool,
    lower_divertor::Bool,
    divertor_hfs_size_fraction::Float64,
    divertor_lfs_size_fraction::Float64) where {T<:Real}

    upper_divertor = Int(upper_divertor)
    lower_divertor = Int(lower_divertor)

    div_gap = gap / 2

    # domain of the equilibrium includes some buffer for divertor slots
    eqt2d = findfirst(:rectangular, eqt.profiles_2d)
    req = eqt2d.grid.dim1
    req = [req[1] + div_gap, req[end] - div_gap, req[end] - div_gap, req[1] + div_gap, req[1] + div_gap]
    zeq = eqt2d.grid.dim2
    zeq = [zeq[1] + div_gap, zeq[1] + div_gap, zeq[end] - div_gap, zeq[end] - div_gap, zeq[1] + div_gap]
    eq_domain = collect(zip(req, zeq))

    # lcfs and magnetic axis
    ψb = eqt.global_quantities.psi_boundary
    pf_wall_r, pf_wall_z = IMAS.first_wall(eqt, pfa)
    ((rlcfs, zlcfs),) = IMAS.flux_surface(eqt, ψb, :closed, pf_wall_r, pf_wall_z)

    # convex hull of equilibrium lcfs with ideal boundary from pulse_schedule
    pc_r, pc_z = IMAS.boundary(pc)
    hull = IMAS.convex_hull([rlcfs; pc_r], [zlcfs; pc_z]; closed_polygon=true)
    rlcfs = [r for (r, z) in hull]
    zlcfs = [z for (r, z) in hull]

    # uniform resample of lcfs
    rlcfs, zlcfs = IMAS.resample_2d_path(rlcfs, zlcfs; n_points=201, method=:linear)
    IMAS.reorder_flux_surface!(rlcfs, zlcfs, argmax(rlcfs))

    # Set the radial build thickness of the plasma
    R_hfs_plasma = minimum(rlcfs) - gap
    R_lfs_plasma = maximum(rlcfs) + gap
    R0 = (R_hfs_plasma + R_lfs_plasma) / 2.0
    Z0 = (minimum(zlcfs) + maximum(zlcfs)) / 2.0

    # main chamber (clip elements that go beyond plasma radial build thickness)
    R, Z = IMAS.rdp_simplify_2d_path(rlcfs, zlcfs, gap / 2)
    wall_poly = LibGEOS.buffer(xy_polygon(R, Z), gap * 3 / 2)

    # divertor lengths
    minor_radius = (maximum(rlcfs) - minimum(rlcfs)) / 2.0

    # private flux regions sorted by distance from lcfs
    private = IMAS.flux_surface(eqt, ψb, :open, pf_wall_r, pf_wall_z)
    sort!(private; by=p -> IMAS.minimum_distance_polygons_vertices(p..., rlcfs, zlcfs))

    t = LinRange(0, 2π, 31)
    detected_upper = upper_divertor
    detected_lower = lower_divertor
    for (pr, pz) in private
        if sign(pz[1] - Z0) != sign(pz[end] - Z0)
            # open flux surface does not encicle the plasma
            continue
        end

        # xpoint between lcfs and private region
        index = IMAS.minimum_distance_polygons_vertices(pr, pz, rlcfs, zlcfs; return_index=true)
        Rx = (pr[index[1]] + rlcfs[index[2]]) / 2.0
        Zx = (pz[index[1]] + zlcfs[index[2]]) / 2.0
        d = sqrt((pr[index[1]] - rlcfs[index[2]])^2 + (pz[index[1]] - zlcfs[index[2]])^2)
        if d > minor_radius
            continue
        end

        # check that this divertor is in fact installed
        if Zx > Z0 && upper_divertor == 0
            continue
        elseif Zx < Z0 && lower_divertor == 0
            continue
        end

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

        # distance from points in private flux region to X-point
        Dx = sqrt.((Rx .- pr) .^ 2.0 .+ (Zx .- pz) .^ 2.0)
        L = cumsum(sqrt.(IMAS.gradient(pr) .^ 2 .+ IMAS.gradient(pz) .^ 2))
        L .-= L[argmin(Dx)]
        if pr[argmin(Dx)+1] < pr[argmin(Dx)-1]
            L .*= -1.0
        end

        # limit extent of private flux regions
        divertor_hfs_size = min(maximum(-L), divertor_hfs_size_fraction * minor_radius)
        divertor_lfs_size = min(maximum(L), divertor_lfs_size_fraction * minor_radius)
        index = (-divertor_hfs_size .< L) .&& (L .< divertor_lfs_size)
        pr = pr[index]
        pz = pz[index]
        L = L[index]

        # slot carving
        L = ((1.0 .- abs.(L) ./ maximum(abs.(L))) .* (1.0 + max(divertor_hfs_size, divertor_lfs_size)) * 3.0 .+ 1.0) .* div_gap / 2.0
        POLYs = [xy_polygon(IMAS.thick_line_polygon(pr[i], pz[i], pr[i+1], pz[i+1], L[i], L[i+1])) for i in 1:length(pr)-1]
        POLYs[1] = LibGEOS.buffer(POLYs[1], div_gap)
        for k in 2:length(POLYs)
            POLYs[1] = LibGEOS.union(POLYs[1], LibGEOS.buffer(POLYs[k], div_gap))
        end
        slot_poly = POLYs[1]

        wall_poly = LibGEOS.union(wall_poly, slot_poly)
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
    corner_radius = gap / 4
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
        pp = plot(wall_poly; aspect_ratio=:equal, alpha=0.5)
        for (pr, pz) in private
            plot!(pp, pr, pz; label="")
        end
        display(pp)
        rethrow(e)
    end

    # update the wall IDS
    IMAS.first_wall!(wall, pr, pz)

    return wall
end

function divertor_regions!(
    bd::IMAS.build{T},
    eqt::IMAS.equilibrium__time_slice{T},
    divertors::IMAS.divertors{T},
    wall_r::AbstractVector{T},
    wall_z::AbstractVector{T}) where {T<:Real}

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
    backwall_layer = bd.layer[ipl-1]
    for ltype in (_blanket_, _shield_, _wall_)
        iwls = IMAS.get_build_indexes(bd.layer; type=ltype, fs=_hfs_)
        if !isempty(iwls)
            backwall_layer = bd.layer[iwls[end]]
            break
        end
    end
    backwall_poly = xy_polygon(backwall_layer)
    backwall_rz = [v for v in GeoInterface.coordinates(backwall_poly)[1]]

    ψb = eqt.profiles_1d.psi[end]
    ((rlcfs, zlcfs),) = IMAS.flux_surface(eqt, ψb, :closed, wall_r, wall_z)
    minor_radius = (maximum(rlcfs) - minimum(rlcfs)) / 2.0

    empty!(divertors)

    detected_upper = bd.divertors.upper.installed
    detected_lower = bd.divertors.lower.installed

    private = IMAS.flux_surface(eqt, ψb, :open, wall_r, wall_z)
    sort!(private; by=p -> IMAS.minimum_distance_polygons_vertices(p..., rlcfs, zlcfs))
    for (pr, pz) in private
        if sign(pz[1] - ZA) != sign(pz[end] - ZA)
            # open flux surface does not encicle the plasma
            continue
        end

        # xpoint between lcfs and private region
        index = IMAS.minimum_distance_polygons_vertices(pr, pz, rlcfs, zlcfs; return_index=true)
        Rx = (pr[index[1]] + rlcfs[index[2]]) / 2.0
        Zx = (pz[index[1]] + zlcfs[index[2]]) / 2.0
        d = sqrt((pr[index[1]] - rlcfs[index[2]])^2 + (pz[index[1]] - zlcfs[index[2]])^2)
        if d > minor_radius
            continue
        end

        # check that this divertor is in fact installed
        if Zx > ZA && (ismissing(bd.divertors.upper, :installed) || bd.divertors.upper.installed != 1)
            continue
        elseif Zx < ZA && (ismissing(bd.divertors.lower, :installed) || bd.divertors.lower.installed != 1)
            continue
        end

        # private flux region must be in the wall
        prz = [(rr, zz) for (rr, zz) in zip(pr, pz) if PolygonOps.inpolygon((rr, zz), backwall_rz) == 1]
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

        α = 5.0
        domain_r = vcat(xx, reverse(xx), xx[1])
        domain_z = vcat(yy, [(Zx - ZA) * α + ZA, (Zx - ZA) * α + ZA], yy[1])
        domain_poly = xy_polygon(domain_r, domain_z)
        backwall_domain_poly = try
            LibGEOS.intersection(backwall_poly, domain_poly)
        catch e
            display(plot(backwall_poly; aspect_ratio=:equal))
            display(plot!(eqt; cx=true))
            display(scatter!([RA, Rx], [ZA, Zx]))
            display(plot!(xx, yy))
            display(plot!(domain_poly; alpha=0.5))
            rethrow(e)
        end
        divertor_poly = LibGEOS.difference(backwall_domain_poly, plasma_poly)

        # Assign to build structure
        coords = GeoInterface.coordinates(divertor_poly)
        structure = resize!(bd.structure, "type" => Int(_divertor_), "name" => "$ul_name divertor")
        structure.material = "tungsten"
        structure.toroidal_extent = 2π

        try
            if typeof(coords[1][1]) <: Vector{Vector{Float64}}
                # if poly intersection finds multiple domains
                # pick the domain with the most points
                coords = coords[argmax(map(length, coords[1]))][1]
            else
                # if poly intersection finds single domain
                coords = coords[1]
            end
            structure.outline.r = [v[1] for v in coords]
            structure.outline.z = [v[2] for v in coords]
        catch e
            p = plot(eqt; cx=true)
            plot!(backwall_poly; alpha=0.3)
            plot!(domain_poly; alpha=0.3)
            plot!(plasma_poly; alpha=0.3)
            display(p)
            rethrow(e)
        end

        # find divertor plasma facing surfaces starting from divertor structure
        divertor_r = Float64[]
        divertor_z = Float64[]
        for (r, z) in zip(structure.outline.r, structure.outline.z)
            if any(sqrt.((r .- pl_r) .^ 2 .+ (z .- pl_z) .^ 2) .< 1E-3)
                push!(divertor_r, r)
                push!(divertor_z, z)
            end
        end

        # split inner/outer strike plates based on line connecting X-point to magnetic axis
        m = (Zx - ZA) / (Rx - RA)
        xx = [0.0, RA * 2.0]
        yy = line_through_point(m, Rx, Zx, xx)
        indexes, crossings = IMAS.intersection(xx, yy, divertor_r, divertor_z)

        # add target info
        if !isempty(indexes)
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
    vessel = IMAS.get_build_layers(bd.layer; type=_vessel_, fs=_lfs_)[1]
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

function obstructing_coils!(
    ocoils::Vector{Vector{IMAS.pf_active__coil{T}}},
    coils::IMAS.IDSvector{IMAS.pf_active__coil{T}},
    layer::IMAS.build__layer{T},
    direction::Int) where {T<:Real}

    if !isempty(coils)
        opposite_layer = IMAS.opposite_side_layer(layer)
        for clayer in (layer, opposite_layer)
            if !ismissing(clayer, :coils_inside)
                if direction > 0
                    push!(ocoils, coils[clayer.coils_inside])
                else
                    pop!(ocoils)
                end
            end
        end
    end
    ocoils = IMAS.pf_active__coil{T}[coil for coil in vcat(ocoils...)]
    obstruction_outline = convex_outline(ocoils)
    return obstruction_outline
end

"""
    build_cx!(bd::IMAS.build{T}, wall::IMAS.wall{T}, pfa::IMAS.pf_active{T}; n_points::Int, verbose::Bool=false) where {T<:Real}

Translates 1D build to 2D cross-sections starting from R and Z coordinates of plasma first wall
"""
function build_cx!(bd::IMAS.build{T}, wall::IMAS.wall{T}, pfa::IMAS.pf_active{T}; n_points::Int, verbose::Bool=false) where {T<:Real}
    # NOTE: In this function all operations are done on HFS, leaving LFS layers as expressions
    plasma = IMAS.get_build_layer(bd.layer; type=_plasma_)

    # "plasma" layer is the first wall
    fw = IMAS.first_wall(wall)
    plasma.outline.r = fw.r
    plasma.outline.z = fw.z

    # detect offset in z
    if IMAS.is_z_offset(plasma.outline.r, plasma.outline.z)
        z_offset = true
    else
        z_offset = false
    end

    plasma_to_tf = reverse(IMAS.get_build_indexes(bd.layer; fs=IMAS._hfs_))
    pushfirst!(plasma_to_tf, plasma_to_tf[1] + 1) # this is the plasma
    # @show plasma_to_tf
    # plot(bd.layer[plasma_to_tf[1]])

    vertical_clearance = 1.0
    kl = 2
    ocoils = Vector{Vector{IMAS.pf_active__coil{T}}}()
    obstruction_outline = nothing
    while kl <= length(plasma_to_tf)
        layer = bd.layer[plasma_to_tf[kl]]
        layer_shape = IMAS.BuildLayerShape(mod(layer.shape, 100))

        if layer_shape !== IMAS._negative_offset_
            verbose && @show "N", IMAS.index(layer), layer.name, layer_shape
            obstruction_outline = obstructing_coils!(ocoils, pfa.coil, layer, +1)
            layer.shape, layer.shape_parameters = FUSE.optimize_layer_outline(
                bd,
                plasma_to_tf[kl-1],
                plasma_to_tf[kl],
                layer_shape;
                vertical_clearance,
                resolution=n_points / 201.0,
                obstruction_outline,
                z_offset)
            if verbose
                plot(ocoils)
                display(plot!(bd.layer[plasma_to_tf[kl]]))
            end

        else
            # go forward until find a normal layer
            neg_off_layers = []
            for kk in kl:length(plasma_to_tf)
                layer = bd.layer[plasma_to_tf[kk]]
                obstruction_outline = obstructing_coils!(ocoils, pfa.coil, layer, +1)
                layer_shape = IMAS.BuildLayerShape(mod(layer.shape, 100))
                if layer_shape != IMAS._negative_offset_
                    verbose && @show "F", IMAS.index(layer), layer.name, layer_shape
                    layer.shape, layer.shape_parameters = FUSE.optimize_layer_outline(
                        bd,
                        plasma_to_tf[kl-1],
                        plasma_to_tf[kk],
                        layer_shape;
                        vertical_clearance,
                        resolution=n_points / 201.0,
                        obstruction_outline,
                        z_offset)
                    kl = kk
                    if verbose
                        #plot(ocoils)
                        display(plot!(bd.layer[plasma_to_tf[kl]]))
                    end
                    break
                else
                    # deal with negative offsets later
                    push!(neg_off_layers, kk)
                end
            end

            # go back with negative offset layers
            for kk in reverse(neg_off_layers)
                layer = bd.layer[plasma_to_tf[kk]]
                obstruction_outline = obstructing_coils!(ocoils, pfa.coil, layer, -1)
                layer_shape = IMAS.BuildLayerShape(mod(layer.shape, 100))
                verbose && @show "B", IMAS.index(layer), layer.name, layer_shape
                layer.shape, layer.shape_parameters = optimize_layer_outline(
                    bd,
                    plasma_to_tf[kk+1],
                    plasma_to_tf[kk],
                    layer_shape;
                    vertical_clearance,
                    resolution=n_points / 201.0,
                    obstruction_outline=nothing,
                    z_offset)
                if verbose
                    plot(ocoils)
                    display(plot!(bd.layer[plasma_to_tf[kk]]))
                end
            end
        end

        kl += 1
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
        layer_shape = BuildLayerShape(mod(bd.layer[iout[end]].shape, 100))
        optimize_layer_outline(bd, olfs, iout[end], layer_shape; vertical_clearance=1.3)
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
    optimize_layer_outline(
        bd::IMAS.build,
        obstr_index::Int,
        layer_index::Int,
        shape_enum::BuildLayerShape;
        vertical_clearance::Float64=1.0,
        resolution::Float64=1.0,
        obstruction_outline=nothing,
        z_offset::Bool=false)

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
    z_offset::Bool=false)

    shape = Int(shape_enum)

    layer = bd.layer[layer_index]
    obstr = bd.layer[obstr_index]
    # display("Layer $layer_index = $(layer.name)")
    # display("Obstr $obstr_index = $(obstr.name)")

    o_start = obstr.start_radius
    l_start = layer.start_radius
    if obstr.side in (Int(_lhfs_), Int(_out_))
        o_end = obstr.end_radius
    else
        o_end = IMAS.opposite_side_layer(obstr).end_radius
    end
    if layer.side in (Int(_lhfs_), Int(_out_))
        l_end = layer.end_radius
    else
        l_end = IMAS.opposite_side_layer(layer).end_radius
    end
    if layer.side == Int(_out_) || obstr.side == Int(_out_)
        l_start = 0.0
        o_start = 0.0
    end

    # handle internal obstructions (in addition to inner layer)
    oR = obstr.outline.r
    oZ = obstr.outline.z
    if obstruction_outline !== nothing
        hull = IMAS.convex_hull(vcat(oR, obstruction_outline.r), vcat(oZ, obstruction_outline.z); closed_polygon=true)
        oR = [r for (r, z) in hull]
        oZ = [z for (r, z) in hull]
        o_start = minimum(oR)
        o_end = maximum(oR)
    end

    hfs_thickness = o_start - l_start
    lfs_thickness = l_end - o_end

    # handle offset, negative offset, and convex-hull
    if shape in (Int(_offset_), Int(_negative_offset_), Int(_convex_hull_), Int(_miller_))
        # vertical clearance logic is to make things as tight as possible
        # both when having a positive and negative offset
        # NOTE: special logic when obstruction starts on axis (eg. cryostat)
        if o_start == 0.0
            vert_thickness = lfs_thickness
        elseif (hfs_thickness + lfs_thickness) > 0
            vert_thickness = min(hfs_thickness, lfs_thickness)
        else
            vert_thickness = max(hfs_thickness, lfs_thickness)
        end
        R, Z = buffer(oR, oZ, hfs_thickness, lfs_thickness, vert_thickness)

        if shape in (Int(_convex_hull_), Int(_miller_))
            hull = IMAS.convex_hull(R, Z; closed_polygon=true)
            R = [r for (r, z) in hull]
            Z = [z for (r, z) in hull]
        end
        if shape == Int(_miller_)
            R, Z = IMAS.resample_2d_path(R, Z; method=:linear)
            R, Z = IMAS.MXH(R, Z, 4)()
        end
        layer.outline.r, layer.outline.z = R, Z
        shape_parameters = Float64[]

    else # handle shapes
        shape = mod(shape, 100)
        if z_offset
            shape = shape + 100
        end

        hfs_thickness = abs(hfs_thickness)
        lfs_thickness = abs(lfs_thickness)
        oZ = oZ .* vertical_clearance

        func = shape_function(shape; resolution)
        initial_clerance = max(hfs_thickness, lfs_thickness) * vertical_clearance
        shape_parameters0 = initialize_shape_parameters(shape, oR, oZ, l_start, l_end, initial_clerance)
        shape_parameters = optimize_outline(oR, oZ, hfs_thickness, lfs_thickness, func, l_start, l_end, shape_parameters0)
        layer.outline.r, layer.outline.z = func(l_start, l_end, shape_parameters...)
    end

    IMAS.reorder_flux_surface!(layer.outline.r, layer.outline.z)

    return Int(shape_enum), shape_parameters
end

"""
    optimize_outline(r_obstruction, z_obstruction, hfs_thickness, lfs_thickness, func, r_start, r_end, shape_parameters; verbose=false)

Find shape parameters that generate smallest shape and target clearance from an obstruction
"""
function optimize_outline(
    r_obstruction::Vector{Float64},
    z_obstruction::Vector{Float64},
    hfs_thickness::Float64,
    lfs_thickness::Float64,
    func::Function,
    r_start::Float64,
    r_end::Float64,
    shape_parameters::Vector{Float64};
    verbose::Bool=false)

    initial_guess = deepcopy(shape_parameters)

    if length(shape_parameters) in (0, 1)
        func(r_start, r_end, shape_parameters...)

    else
        function cost_shape(
            obstruction_buffered_poly,
            obstruction_buffered_area::Float64,
            func::Function,
            r_start::Float64,
            r_end::Float64,
            shape_parameters::Vector{Float64};
            verbose::Bool=false)

            R, Z = func(r_start, r_end, shape_parameters...)

            RZ_poly = xy_polygon(R, Z)

            try
                intersection_poly = LibGEOS.intersection(RZ_poly, obstruction_buffered_poly)
                intersection_area = LibGEOS.area(intersection_poly)

                difference_poly = LibGEOS.difference(RZ_poly, obstruction_buffered_poly)
                difference_area = LibGEOS.area(difference_poly)

                cost_overlap = (intersection_area - obstruction_buffered_area) * 100.0
                cost_overspill = difference_area

                if verbose
                    @show cost_overlap^2
                    @show cost_overspill^2
                    println()
                end

                return norm((cost_overlap, cost_overspill))
            catch e
                if typeof(e) <: LibGEOS.GEOSError
                    if verbose
                        plot(RZ_poly; alpha=0.5)
                        plot!(obstruction_buffered_poly; alpha=0.5)
                        display(plot!())
                    end
                    return Inf
                end
            end
        end

        # buffer the obstruction, this is used for minimizing mean distance error
        if (hfs_thickness + lfs_thickness) > 0
            vert_thickness = min(hfs_thickness, lfs_thickness)
        else
            vert_thickness = max(hfs_thickness, lfs_thickness)
        end
        r_obstruction_buffered, z_obstruction_buffered = buffer(r_obstruction, z_obstruction, hfs_thickness, lfs_thickness, vert_thickness)
        obstruction_buffered_poly = xy_polygon(r_obstruction_buffered, z_obstruction_buffered)
        obstruction_buffered_area = LibGEOS.area(obstruction_buffered_poly)

        res = Optim.optimize(
            shape_parameters -> cost_shape(
                obstruction_buffered_poly,
                obstruction_buffered_area,
                func,
                r_start,
                r_end,
                shape_parameters
            ),
            initial_guess, length(shape_parameters) == 1 ? Optim.BFGS() : Optim.NelderMead(), Optim.Options(; iterations=1000))
        shape_parameters = Optim.minimizer(res)

        if verbose
            cost_shape(
                obstruction_buffered_poly,
                obstruction_buffered_area,
                func,
                r_start,
                r_end,
                shape_parameters;
                verbose=true
            )
            println(res)
        end
    end

    # R, Z = func(r_start, r_end, shape_parameters...)
    # plot(r_obstruction, z_obstruction; label="obstruction")
    # plot!(r_obstruction_buffered, z_obstruction_buffered; label="obstruction-buffered")
    # plot!(buffer(R, Z, -target_clearance)...; label="final-clearance")
    # plot!(func(r_start, r_end, shape_parameters...); label="final", lw=2)
    # display(plot!(; aspect_ratio=:equal))

    return shape_parameters
end

function convex_outline(coil::Union{IMAS.pf_active__coil{T},IMAS.pf_passive__loop{T}}) where {T<:Real}
    r = T[]
    z = T[]

    for element in coil.element
        oute = IMAS.outline(element)
        append!(r, oute.r)
        append!(z, oute.z)
    end

    points = IMAS.convex_hull(r, z; closed_polygon=true)
    r = [p[1] for p in points]
    z = [p[2] for p in points]

    return (r=r, z=z)
end

"""
    convex_outline(coils::AbstractVector{IMAS.pf_active__coil{T}}) where {T<:Real}

returns convex-hull outline of all the pf coils elements in the input array
"""
function convex_outline(coils::AbstractVector{IMAS.pf_active__coil{T}}) where {T<:Real}
    r = T[]
    z = T[]

    for coil in coils
        for element in coil.element
            oute = IMAS.outline(element)
            append!(r, oute.r)
            append!(z, oute.z)
        end
    end

    points = IMAS.convex_hull(r, z; closed_polygon=true)
    r = [p[1] for p in points]
    z = [p[2] for p in points]

    return (r=r, z=z)
end

function convex_outline(rail::IMAS.build__pf_active__rail{T}) where {T<:Real}
    rails = parent(rail)
    irail = IMAS.index(rail)
    coil_start = sum([rails[k].coils_number for k in 1:irail-1]) + 1
    coil_end = coil_start + rail.coils_number - 1
    dd = IMAS.top_dd(rail)
    coils = IMAS.pf_active__coil{T}[dd.pf_active.coil[k] for k in coil_start:coil_end if :shaping ∈ dd.pf_active.coil[k].function && :flux ∉ dd.pf_active.coil[k].function]
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
        shape = string(IMAS.BuildLayerShape(mod(getproperty(layer, :shape, Int(IMAS._undefined_)), 100)))
        shape = replace(shape, "_undefined_" => "", r"^_" => "", r"_$" => "", "_" => " ")
        push!(df, [group, details, type, layer.thickness, layer.start_radius, layer.end_radius, material, area, volume, shape])
    end

    return df
end

function Base.show(io::IO, mime::MIME"text/plain", layers::IMAS.IDSvector{<:IMAS.build__layer})
    old_lines = get(ENV, "LINES", missing)
    old_columns = get(ENV, "COLUMNS", missing)
    df = DataFrames.DataFrame(layers)
    try
        ENV["LINES"] = 1000
        ENV["COLUMNS"] = 1000
        return show(io, mime, df)
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