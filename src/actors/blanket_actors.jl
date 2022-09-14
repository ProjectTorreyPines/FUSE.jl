#= ============ =#
#  ActorBlanket  #
#= ============ =#

import NNeutronics

mutable struct ActorBlanket <: ReactorAbstractActor
    dd::IMAS.dd
    par::ParametersActor
    blanket_multiplier::Real
    thermal_power_extraction_efficiency::Real
end

function ParametersActor(::Type{Val{:ActorBlanket}})
    par = ParametersActor(nothing)
    par.blanket_multiplier = Entry(Real, "", "Neutron thermal power multiplier in blanket"; default=1.2)
    par.model = Switch([:midplane, :sectors], "", "Blanket model"; default=:sectors)
    par.thermal_power_extraction_efficiency = Entry(
        Real,
        "",
        "Fraction of thermal power that is carried out by the coolant at the blanket interface, rather than being lost in the surrounding strutures.";
        default=1.0
    )
    return par
end

"""
    ActorBlanket(dd::IMAS.dd, act::ParametersAllActors; kw...)

Blanket actor

!!! note 
    Stores data in `dd.blanket`
"""
function ActorBlanket(dd::IMAS.dd, act::ParametersAllActors; kw...)
    par = act.ActorBlanket(kw...)
    actor = ActorBlanket(dd, par)
    step(actor)
    finalize(actor)
    return actor
end

function ActorBlanket(dd::IMAS.dd, par::ParametersActor; kw...)
    logging_actor_init(ActorBlanket)
    par = par(kw...)
    return ActorBlanket(dd, par, par.blanket_multiplier, par.thermal_power_extraction_efficiency)
end

function _step(actor::ActorBlanket)
    @warn "currently only using GAMBL materials in blanket and first wall."
    dd = actor.dd
    empty!(dd.blanket)
    blanket = IMAS.get_build(dd.build, type=_blanket_, fs=_hfs_, raise_error_on_missing=false)
    if blanket === missing
        @warn "No blanket present for ActorBlanket to do anything"
        return actor
    end

    if actor.par.model == :midplane
        midplane_blanket_model(actor)
    elseif actor.par.model == :sectors
        sectors_blanket_model(actor)
    else
        error("ActorBlanket model `$(actor.par.model)` is not supported")
    end

    # Optimize Li6/Li7 ratio to obtain target TBR
    function target_TBR(blanket_model, Li6, dd, target=nothing)
        total_tritium_breeding_ratio = 0.0
        total_power_neutrons = sum([bm.time_slice[].power_incident_neutrons for bm in dd.blanket.module])
        for bm in dd.blanket.module
            bmt = bm.time_slice[]
            bm.layer[2].material = @sprintf("lithium-lead: Li6/7=%3.3f", Li6)
            bmt.tritium_breeding_ratio = NNeutronics.TBR(blanket_model, [dl.thickness for dl in bm.layer]..., Li6)
            total_tritium_breeding_ratio += bmt.tritium_breeding_ratio * bmt.power_incident_neutrons / total_power_neutrons
        end
        if target === nothing
            return total_tritium_breeding_ratio
        else
            return (total_tritium_breeding_ratio - target)^2
        end
    end

    blanket_model_1d = NNeutronics.Blanket()
    res = Optim.optimize(Li6 -> target_TBR(blanket_model_1d, Li6, dd, dd.target.tritium_breeding_ratio), 0.0, 100.0, Optim.GoldenSection(), rel_tol=1E-6)
    total_tritium_breeding_ratio = target_TBR(blanket_model_1d, res.minimizer, dd)

    @ddtime(dd.blanket.tritium_breeding_ratio = total_tritium_breeding_ratio)

    return actor
end

"""
    midplane_blanket_model(actor::ActorBlanket)

Uses wall,blanket,shield thicknesses at midplane.
Calculates wall neutron/radiation loading for hfs and lfs blankets.
"""
function midplane_blanket_model(actor::ActorBlanket)
    dd = actor.dd
    eqt = dd.equilibrium.time_slice[]
    nnt = dd.neutronics.time_slice[]

    total_power_neutrons = IMAS.fusion_power(dd.core_profiles.profiles_1d[]) .* 4 / 5
    total_power_neutrons = sum(nnt.wall_loading.power)
    total_power_radiated = 0.0 # IMAS.radiative_power(dd.core_profiles.profiles_1d[])

    blankets = [structure for structure in dd.build.structure if structure.type == Int(_blanket_)]
    resize!(dd.blanket.module, length(blankets))

    for (k, structure) in enumerate(blankets)
        bm = dd.blanket.module[k]
        bm.name = structure.name

        # evaluate neutron_capture_fraction
        tmp = convex_hull(vcat(eqt.boundary.outline.r, structure.outline.r), vcat(eqt.boundary.outline.z, structure.outline.z); closed_polygon=true)
        index = findall(x -> x == 1, [IMAS.PolygonOps.inpolygon((r, z), tmp) for (r, z) in zip(dd.neutronics.first_wall.r, dd.neutronics.first_wall.z)])
        neutron_capture_fraction = sum(nnt.wall_loading.power[index]) / total_power_neutrons

        # evaluate radiative_capture_fraction (could be done better, since ratiation source may not be coming from core)
        radiative_capture_fraction = neutron_capture_fraction

        resize!(bm.time_slice)
        bmt = bm.time_slice[]

        bmt.power_incident_neutrons = total_power_neutrons .* neutron_capture_fraction
        bmt.power_incident_radiated = total_power_radiated .* radiative_capture_fraction

        bmt.power_thermal_neutrons = bmt.power_incident_neutrons * actor.blanket_multiplier
        bmt.power_thermal_radiated = bmt.power_incident_radiated

        bmt.power_thermal_extracted = actor.thermal_power_extraction_efficiency * (bmt.power_thermal_neutrons + bmt.power_thermal_radiated)

        # blanket layer structure (designed to handle missing wall and/or missing shield)
        resize!(bm.layer, 3)
        if sum(structure.outline.r) / length(structure.outline.r) < eqt.boundary.geometric_axis.r # HFS
            d1 = IMAS.get_build(dd.build, type=_wall_, fs=_hfs_, raise_error_on_missing=false)
            d2 = IMAS.get_build(dd.build, type=_blanket_, fs=_hfs_)
            d3 = IMAS.get_build(dd.build, type=_shield_, fs=_hfs_, return_only_one=false, raise_error_on_missing=false)
            if d3 !== missing
                d3 = d3[end]
            end
        else  # LFS
            d1 = IMAS.get_build(dd.build, type=_wall_, fs=_lfs_)
            d2 = IMAS.get_build(dd.build, type=_blanket_, fs=_lfs_)
            d3 = IMAS.get_build(dd.build, type=_shield_, fs=_lfs_, return_only_one=false, raise_error_on_missing=false)
            if d3 !== missing
                d3 = d3[1]
            end
        end
        for (kl, dl) in enumerate(bm.layer)
            dl.name = "dummy layer $kl"
            dl.thickness = 0.0
            dl.material = "vacuum"
        end
        for (kl, dl) in enumerate([d1, d2, d3])
            if dl !== missing
                for field in [:name, :thickness, :material]
                    setproperty!(bm.layer[kl], field, getproperty(dl, field))
                end
            end
        end
    end

    return actor
end

"""
    sectors_blanket_model(actor::ActorBlanket)

Splits poloidal space into sectors and calculates wall,blanket,shield thicknesses for each segment.
alculates wall neutron/radiation loading for each sector.
"""
function sectors_blanket_model(actor::ActorBlanket)

    num_of_sections = 25
    debug_plot = true

    dd = actor.dd

    dd = actor.dd
    eqt = dd.equilibrium.time_slice[]
    nnt = dd.neutronics.time_slice[]

    total_power_neutrons = IMAS.fusion_power(dd.core_profiles.profiles_1d[]) .* 4 / 5
    total_power_neutrons = sum(nnt.wall_loading.power)
    total_power_radiated = 0.0 # IMAS.radiative_power(dd.core_profiles.profiles_1d[])

    resize!(dd.blanket.module, num_of_sections)
    for k in 1:num_of_sections
        bm = dd.blanket.module[k]
        resize!(bm.layer, 3)
    end

    # magnetic axis (where most of fusion is likely to occur)
    R0 = eqt.global_quantities.magnetic_axis.r
    Z0 = eqt.global_quantities.magnetic_axis.z

    # layers
    wall_layers = []
    for layer in dd.build.layer
        if layer.type == Int(_plasma_)
            push!(wall_layers, layer)
        elseif layer.type == Int(_wall_) && layer.fs == Int(_hfs_)
            push!(wall_layers, layer)
        end
    end

    # structures
    divertor_structs = []
    blanket_structs = []
    for layer in dd.build.structure
        if layer.type == Int(_blanket_)
            push!(blanket_structs, layer)
        elseif layer.type == Int(_divertor_)
            push!(divertor_structs, layer)
        end
    end

    # angles
    angle_bounds = range(0, -2π, num_of_sections + 1)
    angles = angle_bounds[1:end-1] .+ diff(angle_bounds) / 2.0

    wall_intersections_per_angle = Vector{Tuple{Float64,Float64}}[]
    blanket_intersections_per_angle = Vector{Tuple{Float64,Float64}}[]
    divertor_intersections_per_angle = Vector{Tuple{Float64,Float64}}[]
    current_coords = Tuple{Float64,Float64}[]

    # find intersections
    for angle in angles
        r_star = [cos(angle) * 1000 + R0, R0]
        z_star = [sin(angle) * 1000 + Z0, Z0]

        current_coords = Tuple{Float64,Float64}[]
        for layer in wall_layers
            append!(current_coords, IMAS.intersection(
                r_star, z_star,
                layer.outline.r, layer.outline.z)
            )
        end
        push!(wall_intersections_per_angle, current_coords)

        current_coords = Tuple{Float64,Float64}[]
        for layer in blanket_structs
            append!(current_coords, IMAS.intersection(
                r_star, z_star,
                layer.outline.r, layer.outline.z)
            )
        end
        push!(blanket_intersections_per_angle, current_coords)

        current_coords = Tuple{Float64,Float64}[]
        for layer in divertor_structs
            append!(current_coords, IMAS.intersection(
                r_star, z_star,
                layer.outline.r, layer.outline.z)
            )
        end
        push!(divertor_intersections_per_angle, current_coords)
    end

    for (k, intersection_points) in enumerate(wall_intersections_per_angle)
        if length(intersection_points) == 0
            wall_thickness = 0.0
        elseif length(intersection_points) == 2
            wall_thickness = norm((intersection_points[1][1] - intersection_points[2][1], intersection_points[1][2] - intersection_points[2][2]))
        elseif length(intersection_points) == 4
            wall_thickness = norm((intersection_points[1][1] - intersection_points[2][1], intersection_points[1][2] - intersection_points[2][2])) + norm((intersection_points[3][1] - intersection_points[4][1], intersection_points[3][2] - intersection_points[4][2]))
        else
            error("number of intersections $(length(intersection_points))")
        end
        bm = dd.blanket.module[k]
        bm.layer[1].name = "wall"
        bm.layer[1].thickness = wall_thickness
    end

    for (k, intersection_points) in enumerate(blanket_intersections_per_angle)
        if length(intersection_points) == 0
            blanket_thickness = 0.0
        elseif length(intersection_points) == 2
            blanket_thickness = norm((intersection_points[1][1] - intersection_points[2][1], intersection_points[1][2] - intersection_points[2][2]))
        elseif length(intersection_points) == 4
            blanket_thickness = norm((intersection_points[1][1] - intersection_points[2][1], intersection_points[1][2] - intersection_points[2][2])) + norm((intersection_points[3][1] - intersection_points[4][1], intersection_points[3][2] - intersection_points[4][2]))
        else
            error("number of intersections $(length(intersection_points))")
        end
        bm = dd.blanket.module[k]
        bm.layer[2].name = "blanket"
        bm.layer[2].thickness = blanket_thickness
    end

    for (k, intersection_points) in enumerate(divertor_intersections_per_angle)
        if length(intersection_points) == 0
            divertor_thickness = 0.0
        elseif length(intersection_points) == 2
            divertor_thickness = norm((intersection_points[1][1] - intersection_points[2][1], intersection_points[1][2] - intersection_points[2][2]))
        elseif length(intersection_points) == 4
            divertor_thickness = norm((intersection_points[1][1] - intersection_points[2][1], intersection_points[1][2] - intersection_points[2][2])) + norm((intersection_points[3][1] - intersection_points[4][1], intersection_points[3][2] - intersection_points[4][2]))
        else
            error("number of intersections $(length(intersection_points))")
        end
        bm = dd.blanket.module[k]
        bm.layer[3].name = "divertor"
        bm.layer[3].thickness = divertor_thickness
    end

    # neutron wall loading as a function of the poloidal angle
    wall_r = dd.neutronics.first_wall.r .- R0
    wall_z = dd.neutronics.first_wall.z .- Z0
    IMAS.reorder_flux_surface!(wall_r, wall_z, 0.0, 0.0)
    wall_angles = atan.(wall_z, wall_r)
    wall_angles .-= wall_angles[1]
    wall_angles = unwrap(wall_angles)

    # fractional wall loading for each section
    for (k, (angle_end, angle_start)) in enumerate(zip(angle_bounds[1:end-1], angle_bounds[2:end]))
        bm = dd.blanket.module[k]
        index = findall(x -> (x > angle_start) && (x <= angle_end), wall_angles)

        neutron_capture_fraction = sum(nnt.wall_loading.power[index]) / total_power_neutrons
        # evaluate radiative_capture_fraction (could be done better, since ratiation source may not be coming from
        radiative_capture_fraction = neutron_capture_fraction

        resize!(bm.time_slice)
        bmt = bm.time_slice[]

        bmt.power_incident_neutrons = total_power_neutrons .* neutron_capture_fraction
        bmt.power_incident_radiated = total_power_radiated .* radiative_capture_fraction

        bmt.power_thermal_neutrons = bmt.power_incident_neutrons * actor.blanket_multiplier
        bmt.power_thermal_radiated = bmt.power_incident_radiated

        bmt.power_thermal_extracted = actor.thermal_power_extraction_efficiency * (bmt.power_thermal_neutrons + bmt.power_thermal_radiated)
    end

    if debug_plot
        p = plot(legend=false, xlabel="R (m)", ylabel="Z (m)", size=(600, 600), aspect_ratio=:equal)
        plot!(p, [R0], [Z0], marker=:cross)
        xvals = []
        yvals = []
        for intersections in divertor_intersections_per_angle
            for intersection in intersections
                append!(xvals, intersection[1])
                append!(yvals, intersection[2])
            end
        end
        for intersections in blanket_intersections_per_angle
            for intersection in intersections
                append!(xvals, intersection[1])
                append!(yvals, intersection[2])
            end
        end
        for intersections in wall_intersections_per_angle
            for intersection in intersections
                append!(xvals, intersection[1])
                append!(yvals, intersection[2])
            end
        end

        plot!(p, xvals, yvals, seriestype=:scatter)
        for layer in vcat(wall_layers, blanket_structs, divertor_structs)
            rcoords = layer.outline.r
            zcoords = layer.outline.z
            plot!(p, rcoords, zcoords)
        end
        display(p)
    end

    return actor
end