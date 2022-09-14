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
    par.model = Switch([:TBR_1D, :TBR_2D], "", "Blanket model"; default=:TBR_2D)
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
    blanket = IMAS.get_build(actor.dd.build, type=_blanket_, fs=_hfs_, raise_error_on_missing=false)
    if blanket === missing
        @warn "No blanket present for ActorBlanket to do anything"
        return actor
    end
    if actor.par.model == :TBR_1D
        TBR_1D(actor)
    elseif actor.par.model == :TBR_2D
        TBR_2D(actor)
    end
    return actor
end

function TBR_1D(actor::ActorBlanket)
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

        # evaluate radiative_capture_fraction (to be fixed)
        radiative_capture_fraction = 1.0 / length(dd.blanket.module)

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

    # Optimize Li6/Li7 ratio to obtain target TBR
    function target_TBR(blanket_model, Li6, dd, target=nothing)
        total_tritium_breeding_ratio = 0.0
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
    TBR_2D(actor::ActorBlanket)

Calculates 2-dimensional TBR from dd using NNeutronics.TBR(). Calculates layer
thicknesses and wall loading at poloidal angles, a local TBR at each segment, 
then sums for a 2D TBR calculation. Blanket Li6 enrichment given to function in %.
"""
function TBR_2D(actor::ActorBlanket)
    Li6 = 90.0
    num_of_sections = 25
    debug_plot = true

    dd = actor.dd
    resize!(dd.blanket.module, num_of_sections)
    bm = dd.blanket.module

    TBRdd = Dict()

    # magnetic axis (where most of fusion is likely to occur)
    eqt = dd.equilibrium.time_slice[]
    R0 = eqt.global_quantities.magnetic_axis.r
    Z0 = eqt.global_quantities.magnetic_axis.z

    # layers
    wall_layer_names = String[]
    for layer in dd.build.layer
        if layer.type == Int(_plasma_)
            TBRdd[layer.name] = Dict()
            TBRdd[layer.name]["r_coordinates"] = push!(layer.outline.r)
            TBRdd[layer.name]["z_coordinates"] = push!(layer.outline.z)
            push!(wall_layer_names, layer.name)
        elseif layer.type == Int(_wall_) && layer.fs == Int(_hfs_)
            TBRdd[layer.name] = Dict()
            TBRdd[layer.name]["r_coordinates"] = push!(layer.outline.r)
            TBRdd[layer.name]["z_coordinates"] = push!(layer.outline.z)
            push!(wall_layer_names, layer.name)
        end
    end

    # structures
    divertor_layer_names = String[]
    blanket_layer_names = String[]
    for layer in dd.build.structure
        TBRdd[layer.name] = Dict()
        TBRdd[layer.name]["r_coordinates"] = push!(layer.outline.r)
        TBRdd[layer.name]["z_coordinates"] = push!(layer.outline.z)
        if layer.type == Int(_blanket_)
            push!(blanket_layer_names, layer.name)
        elseif layer.type == Int(_divertor_)
            push!(divertor_layer_names, layer.name)
        end
    end

    # angles
    angles = range(0, 2π - 2π / num_of_sections, num_of_sections)
    angle_bounds = angles .- π / num_of_sections

    TBRdd["divertor_intersections_per_angle"] = Vector{Tuple{Float64,Float64}}[]
    TBRdd["blanket_intersections_per_angle"] = Vector{Tuple{Float64,Float64}}[]
    TBRdd["wall_intersections_per_angle"] = Vector{Tuple{Float64,Float64}}[]
    current_coords = Tuple{Float64,Float64}[]

    # find intersections
    for angle in angles
        r_star = [cos(angle) * 1000 + R0, R0]
        z_star = [sin(angle) * 1000 + Z0, Z0]

        current_coords = Tuple{Float64,Float64}[]
        for layer_name in blanket_layer_names
            append!(current_coords, IMAS.intersection(
                r_star, z_star,
                TBRdd[layer_name]["r_coordinates"], TBRdd[layer_name]["z_coordinates"])
            )
        end
        push!(TBRdd["blanket_intersections_per_angle"], current_coords)

        current_coords = Tuple{Float64,Float64}[]
        for layer_name in wall_layer_names
            append!(current_coords, IMAS.intersection(
                r_star, z_star,
                TBRdd[layer_name]["r_coordinates"], TBRdd[layer_name]["z_coordinates"])
            )
        end
        push!(TBRdd["wall_intersections_per_angle"], current_coords)

        current_coords = Tuple{Float64,Float64}[]
        for layer_name in divertor_layer_names
            append!(current_coords, IMAS.intersection(
                r_star, z_star,
                TBRdd[layer_name]["r_coordinates"], TBRdd[layer_name]["z_coordinates"])
            )
        end
        push!(TBRdd["divertor_intersections_per_angle"], current_coords)
    end

    for (k, intersection_points) in enumerate(TBRdd["wall_intersections_per_angle"])
        if length(intersection_points) == 0
            wall_thickness = 0.0
        elseif length(intersection_points) == 2
            wall_thickness = norm(vcat(intersection_points[1][1] - intersection_points[2][1], intersection_points[1][2] - intersection_points[2][2]))
        elseif length(intersection_points) == 4
            wall_thickness = norm(vcat(intersection_points[1][1] - intersection_points[2][1], intersection_points[1][2] - intersection_points[2][2])) + norm(vcat(intersection_points[3][1] - intersection_points[4][1], intersection_points[3][2] - intersection_points[4][2]))
        else
            error("number of intersections $(length(intersection_points))")
        end
        resize!(bm[k].layer, 3)
        bm[k].layer[1].thickness = wall_thickness
    end
    for (k, intersection_points) in enumerate(TBRdd["blanket_intersections_per_angle"])
        if length(intersection_points) == 0
            blanket_thickness = 0.0
        elseif length(intersection_points) == 2
            blanket_thickness = norm(vcat(intersection_points[1][1] - intersection_points[2][1], intersection_points[1][2] - intersection_points[2][2]))
        elseif length(intersection_points) == 4
            blanket_thickness = norm(vcat(intersection_points[1][1] - intersection_points[2][1], intersection_points[1][2] - intersection_points[2][2])) + norm(vcat(intersection_points[3][1] - intersection_points[4][1], intersection_points[3][2] - intersection_points[4][2]))
        else
            error("number of intersections $(length(intersection_points))")
        end
        bm[k].layer[2].thickness = blanket_thickness
    end
    for (k, intersection_points) in enumerate(TBRdd["divertor_intersections_per_angle"])
        if length(intersection_points) == 0
            divertor_thickness = 0.0
        elseif length(intersection_points) == 2
            divertor_thickness = norm(vcat(intersection_points[1][1] - intersection_points[2][1], intersection_points[1][2] - intersection_points[2][2]))
        elseif length(intersection_points) == 4
            divertor_thickness = norm(vcat(intersection_points[1][1] - intersection_points[2][1], intersection_points[1][2] - intersection_points[2][2])) + norm(vcat(intersection_points[3][1] - intersection_points[4][1], intersection_points[3][2] - intersection_points[4][2]))
        else
            error("number of intersections $(length(intersection_points))")
        end
        bm[k].layer[3].thickness = divertor_thickness
    end

    # neutron wall loading as a function of the poloidal angle
    fw_power_list = dd.neutronics.time_slice[].wall_loading.power
    wall_angles = atan.(dd.neutronics.time_slice[].wall_loading.flux_z.-Z0, dd.neutronics.time_slice[].wall_loading.flux_r.-R0)
    wall_angles = unwrap(wall_angles)
    wall_angles .-= wall_angles[end]

    fw_total_power = sum(dd.neutronics.time_slice[].wall_loading.power)
    fw_fractional_power_list = [0.0 for x in angles]
    for (wall_angle_index, wall_angle) in enumerate(wall_angles)
        @assert wall_angle <= 2 * pi "Angle exceeds 2π, problem somewhere"
        for (angle_bin_index, angle_bin_lower_bound) in enumerate(angle_bounds)
            if angle_bin_lower_bound == last(angle_bounds)
                if wall_angle >= angle_bin_lower_bound && wall_angle < first(angle_bounds) + 2 * pi
                    fw_fractional_power_list[angle_bin_index] += fw_power_list[wall_angle_index] / fw_total_power
                end
            elseif angle_bin_lower_bound < 0.0
                angle_bin_upper_bound = angle_bounds[angle_bin_index+1]
                if wall_angle >= 2 * pi + angle_bin_lower_bound || wall_angle < angle_bin_upper_bound
                    fw_fractional_power_list[angle_bin_index] += fw_power_list[wall_angle_index] / fw_total_power
                end
            else
                angle_bin_upper_bound = angle_bounds[angle_bin_index+1]
                if wall_angle >= angle_bin_lower_bound && wall_angle < angle_bin_upper_bound
                    fw_fractional_power_list[angle_bin_index] += fw_power_list[wall_angle_index] / fw_total_power
                end
            end
        end
    end

    tbr = 0.0
    blanket_model_section = NNeutronics.Blanket()
    for (index, power_frac) in enumerate(fw_fractional_power_list)
        if bm[index].layer[2].thickness == 0
            local_tbr = 0.0
        else
            local_tbr = NNeutronics.TBR(blanket_model_section, bm[index].layer[1].thickness, bm[index].layer[2].thickness, bm[index].layer[3].thickness, Li6)
        end
        tbr += local_tbr * power_frac
        # println("wall:",TBRdd["wall_thickness_per_angle"][index]/100,"  blanket:",TBRdd["blanket_thickness_per_angle"][index]/100," divertor:", TBRdd["divertor_thickness_per_angle"][index]/100," TBR:",local_tbr)
    end
    println("Total TBR: ", tbr)

    if debug_plot
        p = plot(legend=false, xlabel="R (m)", ylabel="Z (m)", size=(600, 600))
        xvals = []
        yvals = []
        for intersections in TBRdd["divertor_intersections_per_angle"]
            for intersection in intersections
                append!(xvals, intersection[1])
                append!(yvals, intersection[2])
            end
        end
        for intersections in TBRdd["blanket_intersections_per_angle"]
            for intersection in intersections
                append!(xvals, intersection[1])
                append!(yvals, intersection[2])
            end
        end
        for intersections in TBRdd["wall_intersections_per_angle"]
            for intersection in intersections
                append!(xvals, intersection[1])
                append!(yvals, intersection[2])
            end
        end

        plot!(p, xvals, yvals, seriestype=:scatter)
        for layer in vcat(wall_layer_names, vcat(blanket_layer_names, divertor_layer_names))
            rcoords = TBRdd[layer]["r_coordinates"]
            zcoords = TBRdd[layer]["z_coordinates"]
            plot!(p, rcoords, zcoords)
        end
        display(p)
    end

    return actor
end