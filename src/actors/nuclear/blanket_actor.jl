import NNeutronics

#= ============ =#
#  ActorBlanket  #
#= ============ =#
Base.@kwdef mutable struct FUSEparameters__ActorBlanket{T<:Real} <: ParametersActor{T}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    _time::Float64 = NaN
    minimum_first_wall_thickness::Entry{T} = Entry{T}("m", "Minimum first wall thickness"; default=0.02)
    blanket_multiplier::Entry{T} = Entry{T}("-", "Neutron thermal power multiplier in blanket"; default=1.2)
    thermal_power_extraction_efficiency::Entry{T} = Entry{T}("-",
        "Fraction of thermal power that is carried out by the coolant at the blanket interface, rather than being lost in the surrounding strutures.";
        default=1.0)
    verbose::Entry{Bool} = act_common_parameters(verbose=false)
end

mutable struct ActorBlanket{D,P} <: CompoundAbstractActor{D,P}
    dd::IMAS.dd{D}
    par::FUSEparameters__ActorBlanket{P}
    act::ParametersAllActors
    function ActorBlanket(dd::IMAS.dd{D}, par::FUSEparameters__ActorBlanket{P}, act::ParametersAllActors; kw...) where {D<:Real,P<:Real}
        logging_actor_init(ActorBlanket)
        par = par(kw...)
        return new{D,P}(dd, par, act)
    end
end

"""
    ActorBlanket(dd::IMAS.dd, act::ParametersAllActors; kw...)

Evaluates blankets tritium breeding ratio (TBR), heat deposition, and neutron leakage

!!! note 
    Stores data in `dd.blanket`
"""
function ActorBlanket(dd::IMAS.dd, act::ParametersAllActors; kw...)
    actor = ActorBlanket(dd, act.ActorBlanket, act; kw...)
    step(actor)
    finalize(actor)
    return actor
end

function _step(actor::ActorBlanket)
    dd = actor.dd
    empty!(dd.blanket)
    
    blankets = IMAS.get_build_layers(dd.build.layer, type=_blanket_)
    if isempty(blankets)
        @warn "No blanket present for ActorBlanket to do anything"
        return actor
    end

    dd = actor.dd
    par = actor.par
    act = actor.act

    eqt = dd.equilibrium.time_slice[]
    ntt = dd.neutronics.time_slice[]

    raxis = eqt.global_quantities.magnetic_axis.r
    zaxis = eqt.global_quantities.magnetic_axis.z

    total_power_neutrons = sum(ntt.wall_loading.power)
    total_power_radiated = 0.0 # IMAS.radiative_power(dd.core_profiles.profiles_1d[])

    # note flux is defined at the cells, not the nodes
    wall_r = (dd.neutronics.first_wall.r[1:end-1] .+ dd.neutronics.first_wall.r[2:end]) ./ 2.0
    wall_z = (dd.neutronics.first_wall.z[1:end-1] .+ dd.neutronics.first_wall.z[2:end]) ./ 2.0

    # all blanket modules
    blankets = [structure for structure in dd.build.structure if structure.type == Int(_blanket_)]
    resize!(dd.blanket.module, length(blankets))

    # Optimize midplane layers thicknesses that:
    # - achieve target TBR
    # - minimize leakage
    # - while fitting in current blanket thickness
    modules_relative_thickness13 = Float64[]
    modules_effective_thickness = Matrix{Float64}[]
    modules_wall_loading_power = Vector{Float64}[]
    modules_flux_geometric_scale = Vector{Float64}[]
    blanket_model_1d = NNeutronics.Blanket()
    for (istructure, structure) in enumerate(blankets)
        bm = dd.blanket.module[istructure]
        bm.name = structure.name

        # get the relevant blanket layers (d1, d2, d3)
        if sum(structure.outline.r) / length(structure.outline.r) < eqt.boundary.geometric_axis.r && length(blankets) > 1 # HFS
            fs = _hfs_
        else
            fs = _lfs_
        end
        d1 = IMAS.get_build_layers(dd.build.layer, type=_wall_, fs=fs)
        if !isempty(d1)
            # if there are multiple walls we choose the one closest to the plasma
            if fs == _hfs_
                d1 = d1[end]
            else
                d1 = d1[1]
            end
        end
        d2 = IMAS.get_build_layer(dd.build.layer, type=_blanket_, fs=fs)
        d3 = IMAS.get_build_layers(dd.build.layer, type=_shield_, fs=fs)
        if !isempty(d3)
            # if there are multiple shields we choose the one closest to the plasma
            if fs == _hfs_
                d3 = d3[end]
            else
                d3 = d3[1]
            end
        end
        append!(modules_relative_thickness13, [1.0, 1.0])

        # assign blanket module layers (designed to handle missing wall and/or missing shield)
        resize!(bm.layer, 3)
        for (kl, dl) in enumerate(bm.layer)
            dl.name = "dummy layer $kl"
            dl.midplane_thickness = 0.0
            dl.material = "vacuum"
        end
        for (kl, dl) in enumerate([d1, d2, d3])
            if !isempty(dl)
                bm.layer[kl].name = dl.name
                bm.layer[kl].midplane_thickness = dl.thickness
                bm.layer[kl].material = dl.material
            end
        end

        # identify first wall portion of the blanket module
        tmp = convex_hull(vcat(eqt.boundary.outline.r, structure.outline.r), vcat(eqt.boundary.outline.z, structure.outline.z); closed_polygon=true)
        index = findall(x -> x == 1, [IMAS.PolygonOps.inpolygon((r, z), tmp) for (r, z) in zip(wall_r, wall_z)])
        istart = argmin(abs.(wall_z[index]))
        if IMAS.getindex_circular(wall_z[index], istart + 1) > wall_z[index][istart]
            i = argmin(wall_z[index])
            index = vcat(index[i:end], index[1:i-1])
        else
            i = argmin(wall_z[index])
            index = vcat(index[i+1:end], index[1:i])
        end

        # calculate effective thickness of the layers based on the direction of the average impinging neutron flux
        effective_thickness = zeros(length(index), 3)
        flux_geometric_scale = zeros(length(index))
        r_coords = zeros(4)
        z_coords = zeros(4)
        for (k, (r0, z0, fr, fz)) in enumerate(zip(wall_r[index], wall_z[index], ntt.wall_loading.flux_r[index], ntt.wall_loading.flux_z[index]))
            fn = norm(fr, fz)
            r_coords[1] = r0
            z_coords[1] = z0
            for (ilayer, layer) in enumerate([d1, d2, d3])
                if isempty(layer)
                    hit = []
                else
                    _, hit = IMAS.intersection([r0, fr / fn * 1000 + r0], [z0, fz / fn * 1000 + z0], layer.outline.r, layer.outline.z)
                end
                if isempty(hit)
                    r_coords[ilayer+1] = r_coords[ilayer]
                    z_coords[ilayer+1] = z_coords[ilayer]
                else
                    r_coords[ilayer+1] = hit[1][1]
                    z_coords[ilayer+1] = hit[1][2]
                end
            end
            effective_thickness[k, :] .= sqrt.(diff(r_coords[:]) .^ 2.0 .+ diff(z_coords[:]) .^ 2.0)

            # approximate geometric scale of neutron flux from fist wall to the back of the blanket
            R1 = r_coords[1]
            a1 = sqrt((r_coords[1] - raxis)^2 + (z_coords[1] - zaxis)^2)
            R3 = r_coords[end]
            a3 = sqrt((r_coords[end] - raxis)^2 + (z_coords[end] - zaxis)^2)
            flux_geometric_scale[k] = (R1 * a1) / (R3 * a3) # 4π²Rr
        end

        # assign totals in IDS
        total_neutron_capture_fraction = sum(ntt.wall_loading.power[index]) / total_power_neutrons
        # evaluate radiative_capture_fraction (could be done better, since ratiation source may not be coming from core)
        total_radiative_capture_fraction = total_neutron_capture_fraction
        bmt = resize!(bm.time_slice)
        bmt.power_incident_neutrons = total_power_neutrons .* total_neutron_capture_fraction
        bmt.power_incident_radiated = total_power_radiated .* total_radiative_capture_fraction
        bmt.power_thermal_neutrons = bmt.power_incident_neutrons * par.blanket_multiplier
        bmt.power_thermal_radiated = bmt.power_incident_radiated
        bmt.power_thermal_extracted = par.thermal_power_extraction_efficiency * (bmt.power_thermal_neutrons + bmt.power_thermal_radiated)

        push!(modules_effective_thickness, effective_thickness)
        push!(modules_wall_loading_power, ntt.wall_loading.power[index])
        push!(modules_flux_geometric_scale, flux_geometric_scale)
    end

    # Optimize layers thicknesses and minimize Li6 enrichment needed to match
    # target TBR and minimize neutron leakage (realistic geometry and wall loading)
    function target_TBR2D(blanket_model::NNeutronics.Blanket, modules_relative_thickness13::Vector{<:Real}, Li6::Real, dd::IMAS.dd, modules_effective_thickness::Vector{Matrix{Float64}}, modules_wall_loading_power::Vector{<:Any}, total_power_neutrons::Real, min_d1::Float64=0.02, target::Float64=0.0)
        energy_grid = NNeutronics.energy_grid()
        total_tritium_breeding_ratio = 0.0
        Li6 = min(max(abs(Li6), 0.0), 100.0)
        modules_peak_wall_flux = zeros(length(dd.blanket.module))
        modules_peak_escape_flux = zeros(length(dd.blanket.module))
        extra_cost = 0.0
        for (ibm, bm) in enumerate(dd.blanket.module)
            bmt = bm.time_slice[]
            bm.layer[2].material = @sprintf("lithium-lead: Li6/7=%3.3f", Li6)
            module_tritium_breeding_ratio = 0.0
            module_wall_loading_power = sum(modules_wall_loading_power[ibm])
            d1 = bm.layer[1].midplane_thickness * abs(modules_relative_thickness13[1+(ibm-1)*2])
            if d1 < min_d1
                extra_cost += (min_d1 - d1)
                d1 = max(d1, min_d1)
            end
            d2 = bm.layer[2].midplane_thickness
            d3 = bm.layer[3].midplane_thickness * abs(modules_relative_thickness13[2+(ibm-1)*2])
            dtot = d1 + d2 + d3
            x1 = d1 / dtot
            x2 = d2 / dtot
            x3 = d3 / dtot
            for k in eachindex(modules_wall_loading_power[ibm])
                ed1 = modules_effective_thickness[ibm][k, 1] * x1
                ed2 = modules_effective_thickness[ibm][k, 2] * x2
                ed3 = modules_effective_thickness[ibm][k, 3] * x3
                module_tritium_breeding_ratio += (NNeutronics.TBR(blanket_model, ed1, ed2, ed3, Li6) * modules_wall_loading_power[ibm][k] / module_wall_loading_power)
                #NOTE: leakeage_energy is total number of neutrons in each energy bin, so just a sum is correct
                LE = NNeutronics.leakeage_energy(blanket_model, ed1, ed2, ed3, Li6, energy_grid)::Vector{Float64}
                escape_flux = sum(LE) * modules_wall_loading_power[ibm][k] * modules_flux_geometric_scale[ibm][k]
                if escape_flux > modules_peak_escape_flux[ibm]
                    modules_peak_escape_flux[ibm] = escape_flux
                end
            end
            modules_peak_wall_flux[ibm] = maximum(modules_wall_loading_power[ibm])
            total_tritium_breeding_ratio += module_tritium_breeding_ratio * module_wall_loading_power / total_power_neutrons
            if target === 0.0
                bmt.tritium_breeding_ratio = module_tritium_breeding_ratio
                bmt.peak_wall_flux = modules_peak_wall_flux[ibm]
                bmt.peak_escape_flux = modules_peak_escape_flux[ibm]
                bm.layer[1].midplane_thickness = d1
                bm.layer[2].midplane_thickness = d2
                bm.layer[3].midplane_thickness = d3
            end
        end
        if target === 0.0
            @ddtime(dd.blanket.tritium_breeding_ratio = total_tritium_breeding_ratio)
        else
            cost = norm([(total_tritium_breeding_ratio - target), maximum(modules_peak_escape_flux) / sum(modules_peak_wall_flux), extra_cost])
            return cost
        end
    end

    res = Optim.optimize(x -> target_TBR2D(blanket_model_1d, x[1:end-1], x[end], dd, modules_effective_thickness, modules_wall_loading_power, total_power_neutrons, par.minimum_first_wall_thickness, dd.requirements.tritium_breeding_ratio), vcat(modules_relative_thickness13, 50.0), Optim.NelderMead())#; autodiff=:forward)#, rel_tol=1E-6)
    total_tritium_breeding_ratio = target_TBR2D(blanket_model_1d, res.minimizer[1:end-1], abs(res.minimizer[end]), dd, modules_effective_thickness, modules_wall_loading_power, total_power_neutrons, par.minimum_first_wall_thickness)
    if par.verbose
        println(res)
    end

    # rebuild geometry
    ActorCXbuild(dd, act)

    return actor
end
