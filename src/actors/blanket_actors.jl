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
    if actor.par.model==:TBR_1D
        return TBR_1D(actor)
    elseif actor.par.model==:TBR_2D
        return TBR_2D(actor)
    end
end



function TBR_1D( actor::ActorBlanket)
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
    TBR_2D(dd::IMAS.dd, Li6::Float64; do_plot::Bool=false)

Calculates 2-dimensional TBR from dd using NNeutronics.TBR(). Calculates layer
thicknesses and wall loading at poloidal angles, a local TBR at each segment, 
then sums for a 2D TBR calculation. Blanket Li6 enrichment given to function
in %.
"""
function TBR_2D(actor::ActorBlanket)
    Li6=90.0
    do_plot=true
    dd = actor.dd
    println("WARNING: currently only using GAMBL materials in blanket and first wall.")
    
    TBRdd=Dict()
    plasma_midpoint=()
    blanket_present=false
    wall_layer_names=[]
    blanket_layer_names=[]
    divertor_layer_names=[]

    for layer in dd.build.layer
        if layer.type==-1
            plasma_midpoint=[(layer.end_radius*100+layer.start_radius*100)/2]
            plasma_midpoint=hcat(plasma_midpoint,0)
            TBRdd[layer.name]=Dict()
            TBRdd[layer.name]["r_coordinates"]=push!(layer.outline.r.*100)
            TBRdd[layer.name]["z_coordinates"]=push!(layer.outline.z.*100)
            push!(wall_layer_names,layer.name)
        elseif layer.type == 5 && occursin("hfs",layer.name)
            TBRdd[layer.name]=Dict()
            TBRdd[layer.name]["r_coordinates"]=push!(layer.outline.r.*100)
            TBRdd[layer.name]["z_coordinates"]=push!(layer.outline.z.*100)
            push!(wall_layer_names,layer.name)
        elseif layer.type==4
            blanket_present=true
        end
    end
    if !blanket_present
        error("No blanket present, TBR not evaluated.")
    end

    for layer in dd.build.structure
        TBRdd[layer.name]=Dict()
        TBRdd[layer.name]["r_coordinates"]=push!(layer.outline.r.*100)
        TBRdd[layer.name]["z_coordinates"]=push!(layer.outline.z.*100)
        if layer.type==4
            push!(blanket_layer_names, layer.name)
        elseif layer.type==8
            push!(divertor_layer_names, layer.name)
        end
    end

    num_of_sections=25
    TBRdd["angles"]=range(0,2π-2π/num_of_sections,num_of_sections)
    TBRdd["angle_boundaries"]=TBRdd["angles"].-π/num_of_sections
    fw_power_list=dd.neutronics.time_slice[1].wall_loading.power

    TBRdd["r_values_by_angle"]=Dict()
    TBRdd["z_values_by_angle"]=Dict()

    for num in range(1,num_of_sections)
        TBRdd["r_values_by_angle"][num]=zeros(0)
        TBRdd["z_values_by_angle"][num]=zeros(0)
        push!(TBRdd["r_values_by_angle"][num], (cos(TBRdd["angles"][num])*1000+plasma_midpoint[1]))
        push!(TBRdd["r_values_by_angle"][num], plasma_midpoint[1])
        push!(TBRdd["z_values_by_angle"][num], (sin(TBRdd["angles"][num])*1000+plasma_midpoint[2]))
        push!(TBRdd["z_values_by_angle"][num], plasma_midpoint[2])
    end

    TBRdd["divertor_intersections_per_angle"]=[]
    TBRdd["blanket_intersections_per_angle"]=[]
    TBRdd["wall_intersections_per_angle"]=[]

    for (angle_idx, r_values) in enumerate(TBRdd["r_values_by_angle"])
        current_coords=[]
        for layer_name in blanket_layer_names
            append!(current_coords,IMAS.intersection(
                TBRdd["r_values_by_angle"][angle_idx],TBRdd["z_values_by_angle"][angle_idx],
                TBRdd[layer_name]["r_coordinates"], TBRdd[layer_name]["z_coordinates"])
                )
        end
        push!(TBRdd["blanket_intersections_per_angle"], current_coords)
        current_coords=[]
        for layer_name in wall_layer_names
            append!(current_coords,IMAS.intersection(
                TBRdd["r_values_by_angle"][angle_idx],TBRdd["z_values_by_angle"][angle_idx],
                TBRdd[layer_name]["r_coordinates"], TBRdd[layer_name]["z_coordinates"])
                )
        end
        push!(TBRdd["wall_intersections_per_angle"], current_coords)
        current_coords=[]
        for layer_name in divertor_layer_names
            append!(current_coords,IMAS.intersection(
                TBRdd["r_values_by_angle"][angle_idx],TBRdd["z_values_by_angle"][angle_idx],
                TBRdd[layer_name]["r_coordinates"], TBRdd[layer_name]["z_coordinates"])
                )
        end
        push!(TBRdd["divertor_intersections_per_angle"], current_coords)
    end

    TBRdd["divertor_thickness_per_angle"]=[]
    TBRdd["blanket_thickness_per_angle"]=[]
    TBRdd["wall_thickness_per_angle"]=[]
    for intersection_points in TBRdd["divertor_intersections_per_angle"]
        if length(intersection_points)==0
            divertor_thickness=0.0
        elseif length(intersection_points)==2
            divertor_thickness=norm(vcat(intersection_points[1][1]-intersection_points[2][1],intersection_points[1][2]-intersection_points[2][2]))
        elseif length(intersection_points)==4
            divertor_thickness=norm(vcat(intersection_points[1][1]-intersection_points[2][1],intersection_points[1][2]-intersection_points[2][2]))+norm(vcat(intersection_points[3][1]-intersection_points[4][1],intersection_points[3][2]-intersection_points[4][2]))
        end
        append!(TBRdd["divertor_thickness_per_angle"],divertor_thickness)
    end
    for intersection_points in TBRdd["blanket_intersections_per_angle"]
        if length(intersection_points)==0
            blanket_thickness=0.0
        elseif length(intersection_points)==2
            blanket_thickness=norm(vcat(intersection_points[1][1]-intersection_points[2][1],intersection_points[1][2]-intersection_points[2][2]))
        elseif length(intersection_points)==4
            blanket_thickness=norm(vcat(intersection_points[1][1]-intersection_points[2][1],intersection_points[1][2]-intersection_points[2][2]))+norm(vcat(intersection_points[3][1]-intersection_points[4][1],intersection_points[3][2]-intersection_points[4][2]))
        end
        append!(TBRdd["blanket_thickness_per_angle"],blanket_thickness)
    end
    for intersection_points in TBRdd["wall_intersections_per_angle"]
        if length(intersection_points)==0
            wall_thickness=0.0
        elseif length(intersection_points)==2
            wall_thickness=norm(vcat(intersection_points[1][1]-intersection_points[2][1],intersection_points[1][2]-intersection_points[2][2]))
        elseif length(intersection_points)==4
            wall_thickness=norm(vcat(intersection_points[1][1]-intersection_points[2][1],intersection_points[1][2]-intersection_points[2][2]))+norm(vcat(intersection_points[3][1]-intersection_points[4][1],intersection_points[3][2]-intersection_points[4][2]))
        end
        append!(TBRdd["wall_thickness_per_angle"],wall_thickness)
    end

    wall_r=[r*100-plasma_midpoint[1] for r in dd.neutronics.time_slice[1].wall_loading.flux_r]
    wall_z=[z*100-plasma_midpoint[2] for z in dd.neutronics.time_slice[1].wall_loading.flux_z]
    wall_angles=[]
    for (r,z) in zip(wall_r,wall_z)
        if r<0 && z<0
            append!(wall_angles, atan(z/r)+pi)
        elseif r>=0 && z<0
            append!(wall_angles, 2*pi+atan(z/r))
        elseif r<0 && z>=0
            append!(wall_angles, pi+atan(z/r))
        else
            append!(wall_angles, atan(z/r))
        end
    end


    fw_total_power=sum(dd.neutronics.time_slice[1].wall_loading.power)
    fw_fractional_power_list=[0.0 for x in TBRdd["angles"]]
    angle_bounds=TBRdd["angle_boundaries"]
    for (wall_angle_index, wall_angle) in enumerate(wall_angles)
        if wall_angle>=2*pi
            error("Angle exceeds 2*pi, problem somewhere")
        end
        for (angle_bin_index, angle_bin_lower_bound) in enumerate(angle_bounds)
            if angle_bin_lower_bound == last(angle_bounds)
                if wall_angle >= angle_bin_lower_bound && wall_angle < first(angle_bounds)+2*pi
                    fw_fractional_power_list[angle_bin_index]+=fw_power_list[wall_angle_index]/fw_total_power
                end
            elseif angle_bin_lower_bound < 0.0 
                angle_bin_upper_bound=angle_bounds[angle_bin_index+1]  
                if wall_angle >= 2*pi+angle_bin_lower_bound || wall_angle < angle_bin_upper_bound
                    fw_fractional_power_list[angle_bin_index]+=fw_power_list[wall_angle_index]/fw_total_power
                end
            else
                angle_bin_upper_bound=angle_bounds[angle_bin_index+1]
                if wall_angle >= angle_bin_lower_bound && wall_angle < angle_bin_upper_bound
                    fw_fractional_power_list[angle_bin_index]+=fw_power_list[wall_angle_index]/fw_total_power
                end
            end
        end
    end

    TBRdd["power fraction per angle"]=fw_fractional_power_list
    tbr=0.0
    blanket_model_section=NNeutronics.Blanket()
    for (index, power_frac) in enumerate(TBRdd["power fraction per angle"])
        if TBRdd["blanket_thickness_per_angle"][index] == 0
            local_tbr=0.0
        else
            local_tbr=NNeutronics.TBR(blanket_model_section,TBRdd["wall_thickness_per_angle"][index]/100,TBRdd["blanket_thickness_per_angle"][index]/100, TBRdd["divertor_thickness_per_angle"][index]/100, Li6)
        end
        tbr+=local_tbr*power_frac
        # println("wall:",TBRdd["wall_thickness_per_angle"][index]/100,"  blanket:",TBRdd["blanket_thickness_per_angle"][index]/100," divertor:", TBRdd["divertor_thickness_per_angle"][index]/100," TBR:",local_tbr)
    end
    # println("Total TBR: ", tbr)

    if do_plot
        plot=Plots.plot(legend=false, xlabel="R (m)", ylabel="Z (m)",size=(600,600))
        xvals=[]
        yvals=[]
        for intersections in TBRdd["divertor_intersections_per_angle"]
            for intersection in intersections
                append!(xvals,intersection[1]/100)
                append!(yvals,intersection[2]/100)
            end
        end
        for intersections in TBRdd["blanket_intersections_per_angle"]
            for intersection in intersections
                append!(xvals,intersection[1]/100)
                append!(yvals,intersection[2]/100)
            end
        end
        for intersections in TBRdd["wall_intersections_per_angle"]
            for intersection in intersections
                append!(xvals,intersection[1]/100)
                append!(yvals,intersection[2]/100)
            end
        end

        plot!(plot, xvals, yvals, seriestype = :scatter)
        for layer in vcat(wall_layer_names, vcat(blanket_layer_names, divertor_layer_names)) 
            rcoords=TBRdd[layer]["r_coordinates"]./100
            zcoords=TBRdd[layer]["z_coordinates"]./100
            plot!(plot, rcoords,zcoords)
        end
        display(plot)
    end

    return actor
end