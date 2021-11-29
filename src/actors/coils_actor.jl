@Base.kwdef mutable struct PFcoilsOptTrace
    coils::Vector = []
    currents::Vector = []
    λ_regularize::Vector = []
    cost_ψ::Vector = []
    cost_currents::Vector = []
    cost_bound::Vector = []
    cost_spacing::Vector = []
    cost_total::Vector = []
end

mutable struct PFcoilsOptActor <: CoilsActor
    eq_in::IMAS.equilibrium
    eq_out::IMAS.equilibrium
    time::Real
    pf_active::IMAS.pf_active
    radial_build::IMAS.radial_build
    symmetric::Bool
    λ_regularize::Real
    λ_norm::Real
    trace::PFcoilsOptTrace
end

#= == =#
# INIT #
#= == =#
using Equilibrium
using PolygonOps
using StaticArrays
using DSP
using Random
using Interpolations
using Optim
using AD_GS
using LinearAlgebra
using Statistics
using Plots
import Contour

# NOTE on nomenclature:
# - optim: coils that have theri position and current optimized
# - pinned: coisl with fixed position but current is optimized
# - fixed: fixed position and current
  
"""
    initialize_coils(rb::IMAS.radial_build, ncoils_OH::Int, n_pf_coils_per_gap_region::Vector{Int})

Use radial build layers outline to initialize PF coils distribution
"""
function initialize_coils(rb::IMAS.radial_build, ncoils_OH::Int, n_pf_coils_per_gap_region::Vector)

    resolution = 257
    rmask, zmask, mask = IMAS.structures_mask(rb, resolution=resolution)

    pf_active = IMAS.pf_active()

    resize!(rb.pf_coils_rail, length(n_pf_coils_per_gap_region) + 1)

    # OH coils are distributed on a rail within the OH region
    OH_layer = IMAS.get_radial_build(rb, type=1)
    r_ohcoils = ones(ncoils_OH) .* (sum(extrema(OH_layer.outline.r)) / 2.)
    z_ohcoils = collect(range(minimum(OH_layer.outline.z), maximum(OH_layer.outline.z), length=ncoils_OH))    
    oh_coils = [PointCoil(r, z) for (r, z) in zip(r_ohcoils, z_ohcoils)]
    rb.pf_coils_rail[1].name = "OH"
    rb.pf_coils_rail[1].coils_number = ncoils_OH
    rb.pf_coils_rail[1].outline.r = r_ohcoils
    rb.pf_coils_rail[1].outline.z = z_ohcoils
    rb.pf_coils_rail[1].outline.distance = range(-1, 1, length=ncoils_OH)
    for c in oh_coils
        k = length(pf_active.coil) + 1
        resize!(pf_active.coil, k)
        resize!(pf_active.coil[k].element, 1)
        pf_active.coil[k].identifier = "optim"
        pf_active.coil[k].element[1].geometry.rectangle.r = c.R
        pf_active.coil[k].element[1].geometry.rectangle.z = c.Z
        pf_active.coil[k].element[1].geometry.rectangle.width = maximum(c.R) - minimum(c.R)
        pf_active.coil[k].element[1].geometry.rectangle.height = maximum(c.Z) - minimum(c.Z)
        set_field_time_array(pf_active.coil[k].current, :time, 1, 0.0)
        set_field_time_array(pf_active.coil[k].current, :data, 1, 0.0)
    end

    # Now add actual PF coils to regions of vacuum
    krail = 0
    for (k, layer) in enumerate(rb.layer)
        if (layer.hfs == 1 || k == length(rb.layer)) && ! is_missing(layer.outline, :r)
            if ! is_missing(layer, :material) && layer.material == "vacuum"

                krail += 1
                if isa(n_pf_coils_per_gap_region[krail], Int)
                    ncoils = n_pf_coils_per_gap_region[krail]
                else
                    ncoils = length(n_pf_coils_per_gap_region[krail])
                end

                # add rail info to radial_build IDS
                rb.pf_coils_rail[1 + krail].name = replace(replace(layer.name, "hfs " => ""), "lfs " => "")
                rb.pf_coils_rail[1 + krail].coils_number = ncoils

                if ncoils == 0
                    rb.pf_coils_rail[1 + krail].outline.r = Float64[]
                    rb.pf_coils_rail[1 + krail].outline.z = Float64[]
                    rb.pf_coils_rail[1 + krail].outline.distance = Float64[]
                    continue
                end

                # pick layers with outline information
                if layer.hfs == 1
                    outer_layer = IMAS.get_radial_build(rb, identifier=rb.layer[k].identifier, hfs=1)
                    inner_layer = IMAS.get_radial_build(rb, identifier=rb.layer[k + 1].identifier, hfs=[1,0])
                else
                    inner_layer = IMAS.get_radial_build(rb, identifier=rb.layer[k - 1].identifier, hfs=1)
                    outer_layer = IMAS.get_radial_build(rb, identifier=rb.layer[k].identifier, hfs=[1,0])
                end

                # take two outlines and interpolate them on the same θ
                # inner_r, inner_z, outer_r, outer_z, θ = two_curves_same_θ(inner_layer.outline.r, inner_layer.outline.z, outer_layer.outline.r, outer_layer.outline.z)

                clerance = 7 * length(rmask) / 257 # this is reasonable on a 257 mask grid
                buff = (clerance * 2) * (rmask[2] - rmask[1])

                # generate rail between the two layers where coils will be placed and will be able to slide during the `optimization` phase
                poly = LibGEOS.buffer(xy_polygon(inner_layer.outline.r, inner_layer.outline.z), buff)
                mid_r = [v[1] for v in LibGEOS.coordinates(poly)[1]]
                mid_z = [v[2] for v in LibGEOS.coordinates(poly)[1]]

                # mark what regions on that rail do not intersect solid structures and can hold coils
                clerance = Int(ceil(clerance))
                valid_k = []
                for (k, (r, z)) in enumerate(zip(mid_r, mid_z))
                    ir = argmin(abs.(rmask .- r))
                    iz = argmin(abs.(zmask .- z))
                    if (ir - clerance) < 1 || (ir + clerance) > length(rmask) || (iz - clerance) < 1 || (iz + clerance) > length(zmask)
                        continue
                    end
                    if all(mask[(-clerance:clerance) .+ ir,(-clerance:clerance) .+ iz] .== 0)
                        push!(valid_k, k)
                    end
                end
                istart = argmax(diff(valid_k))
                valid_r = fill(NaN, size(mid_r)...)
                valid_z = fill(NaN, size(mid_z)...)
                valid_r[valid_k] = mid_r[valid_k]
                valid_z[valid_k] = mid_z[valid_k]
                valid_r = vcat(valid_r[istart + 1:end], valid_r[1:istart])
                valid_z = vcat(valid_z[istart + 1:end], valid_z[1:istart])

                if isa(n_pf_coils_per_gap_region[krail], Int)
                    coils_distance = range(-(1-1/ncoils),1-1/ncoils,length=ncoils)
                else
                    coils_distance = n_pf_coils_per_gap_region[krail]
                end

                # evaluate distance along rail
                d_distance = sqrt.(diff(vcat(valid_r, valid_r[1])).^2.0 .+ diff(vcat(valid_z, valid_z[1])).^2.0)
                d_distance[isnan.(d_distance)] .= 0.0
                distance = cumsum(d_distance)
                valid_z = valid_z[d_distance .!= 0]
                valid_r = valid_r[d_distance .!= 0]
                distance = distance[d_distance .!= 0]
                distance = (distance .- distance[1])
                distance = (distance ./ distance[end]).*2.0.-1.0

                # add rail info to radial_build IDS
                rb.pf_coils_rail[1 + krail].outline.r = valid_r
                rb.pf_coils_rail[1 + krail].outline.z = valid_z
                rb.pf_coils_rail[1 + krail].outline.distance = distance

                # uniformely distribute coils
                r_coils = IMAS.interp(distance, valid_r)(coils_distance)
                z_coils = IMAS.interp(distance, valid_z)(coils_distance)
                coils = [PointCoil(r, z) for (r, z) in zip(r_coils, z_coils)]

                # populate IMAS data structure
                for c in coils
                    k = length(pf_active.coil) + 1
                    resize!(pf_active.coil, k)
                    resize!(pf_active.coil[k].element, 1)
                    pf_active.coil[k].identifier = "optim"
                    pf_active.coil[k].element[1].geometry.rectangle.r = c.R
                    pf_active.coil[k].element[1].geometry.rectangle.z = c.Z
                    pf_active.coil[k].element[1].geometry.rectangle.width = maximum(c.R) - minimum(c.R)
                    pf_active.coil[k].element[1].geometry.rectangle.height = maximum(c.Z) - minimum(c.Z)
                    set_field_time_array(pf_active.coil[k].current, :time, 1, 0.0)
                    set_field_time_array(pf_active.coil[k].current, :data, 1, 0.0)
                end
            end
        end
    end
    return pf_active
end

function PFcoilsOptActor(eq_in::IMAS.equilibrium, rb::IMAS.radial_build, ncoils_OH::Int, ncoils_per_region::Vector, λ_regularize=1E-13)
    # initialize coils location
    pf_active = initialize_coils(rb, ncoils_OH, ncoils_per_region)

    # basic constructors
    eq_out = deepcopy(eq_in)
    symmetric = false
    λ_norm = 1.0
    time_index = 1
    time = eq_in.time[time_index]

    # constructor
    actor = PFcoilsOptActor(eq_in, eq_out, time, pf_active, rb, symmetric, λ_regularize, λ_norm, PFcoilsOptTrace())

    # calculate initial current distribution (NOTE: this will assign actor.λ_norm)
    step(actor, maxiter=0)

    return actor
end

function mask_interpolant_function(rb::IMAS.radial_build)
    # generate mask
    rmask, zmask, mask = IMAS.structures_mask(rb)

    # Cubic spline interpolation on the log to ensure positivity of the cost
    mask_log_interpolant_raw = Interpolations.CubicSplineInterpolation((rmask, zmask), log10.(1 .+ mask))
    mask_log_interpolant_raw = Interpolations.extrapolate(mask_log_interpolant_raw.itp, Interpolations.Flat());
    function mask_log_interpolant(r, z)
        return (10.0.^(mask_log_interpolant_raw(r, z)) .- 1)
    end
    return mask_log_interpolant
end

#= == =#
# STEP #
#= == =#
# utility functions for packing and unpacking info in/out of optimization function

function pack_mask(optim_coils::Vector, λ_regularize::Float64, symmetric::Bool)::Vector{Float64}
    coilz = vcat([c.R for c in optim_coils if ((! symmetric) || c.Z >= 0)], [c.Z for c in optim_coils if ((! symmetric) || c.Z >= 0)])
    packed = vcat(coilz, log10(λ_regularize))
    return packed
end

function unpack_mask(packed::Vector, symmetric::Bool)
    coilz = packed[1:end - 1]
    λ_regularize = packed[end]
    optim_coils =  [PointCoil(coilz[k], coilz[k + Int(length(coilz) / 2)]) for (k, c) in enumerate(coilz[1:Int(end / 2)])]
    if symmetric
        optim_coils = vcat(optim_coils, [PointCoil(c.R, -c.Z) for c in optim_coils if c.Z != 0])
    end
    return optim_coils, 10^λ_regularize
end

function optimize_coils_mask(EQfixed::Equilibrium.AbstractEquilibrium; pinned_coils::Vector, optim_coils::Vector, fixed_coils::Vector, fixed_currents::Vector, symmetric::Bool, λ_regularize::Real, λ_ψ::Real, λ_currents::Real, rb::IMAS.radial_build, maxiter::Int, verbose::Bool)
    mask_interpolant = mask_interpolant_function(rb)
    fixed_eq = ψp_on_fixed_eq_boundary(EQfixed, fixed_coils, fixed_currents)
    packed = pack_mask(optim_coils, λ_regularize, symmetric)

    trace = PFcoilsOptTrace()
    packed_tmp = []
    function placement_cost(packed; do_trace=false)
        push!(packed_tmp, packed)
        (optim_coils, λ_regularize) = unpack_mask(packed, symmetric)
        coils = vcat(pinned_coils, optim_coils)
        currents, cost_ψ = currents_to_match_ψp(fixed_eq..., coils, λ_regularize=λ_regularize, return_cost=true)
        cost_ψ = cost_ψ / λ_ψ
        cost_currents = norm(currents) / length(currents) / λ_currents
        cost_bound = norm(mask_interpolant.([c.R for c in optim_coils], [c.Z for c in optim_coils]))
        cost_spacing = 0
        for (k1, c1) in enumerate(optim_coils)
            for (k2, c2) in enumerate(optim_coils)
                if k1 == k2
                    continue
                end
                cost_spacing += 1 / sqrt((c1.R - c2.R)^2 + (c1.Z - c2.Z)^2)
            end
        end
        cost_spacing = cost_spacing / sum([rail.coils_number for rail in rb.pf_coils_rail])^2
        cost = sqrt(cost_ψ^2 + cost_currents^2 + cost_bound^2)# + cost_spacing^2)
        if do_trace
            push!(trace.currents, [no_Dual(c) for c in currents])
            push!(trace.coils, vcat(pinned_coils, [PointCoil(no_Dual(c.R), no_Dual(c.Z)) for c in optim_coils]))
            push!(trace.λ_regularize, no_Dual(λ_regularize))
            push!(trace.cost_ψ, no_Dual(cost_ψ))
            push!(trace.cost_currents, no_Dual(cost_currents))
            push!(trace.cost_bound, no_Dual(cost_bound))
            push!(trace.cost_spacing, no_Dual(cost_spacing))
            push!(trace.cost_total, no_Dual(cost))
        end
        return cost
    end

    function clb(x)
        placement_cost(packed_tmp[end]; do_trace=true)
        false
    end
    
    # use NelderMead() ; other optimizer that works is Newton(), others have trouble
    res = Optim.optimize(placement_cost, packed, Optim.Newton(), Optim.Options(time_limit=60 * 2, iterations=maxiter, callback=clb); autodiff=:forward)

    if verbose println(res) end
    packed = Optim.minimizer(res)

    (optim_coils, λ_regularize) = unpack_mask(packed, symmetric)

    pinned_currents = [trace.currents[end][k] for k in 1:length(pinned_coils)]
    optim_currents = [trace.currents[end][k] for k in (length(pinned_coils)+1):(length(pinned_coils)+length(optim_coils))]
    
    return (pinned_coils, pinned_currents), (optim_coils, optim_currents), (fixed_coils, fixed_currents), λ_regularize, trace
end

function pack_rail(rb::IMAS.radial_build, λ_regularize::Float64, symmetric::Bool)::Vector{Float64}
    distances = []
    for rail in rb.pf_coils_rail
        # not symmetric
        if ! symmetric
            coil_distances = collect(range(-1.0, 1.0, length=rail.coils_number + 2))[2:end - 1]
        # even symmetric
        elseif mod(rail.coils_number, 2) == 0
            coil_distances = collect(range(-1.0, 1.0, length=rail.coils_number + 2))[2 + Int(rail.coils_number // 2):end - 1]
        # odd symmetric
        else
            coil_distances = collect(range(-1.0, 1.0, length=rail.coils_number + 2))[2 + Int((rail.coils_number - 1) // 2) + 1:end - 1]
        end
        append!(distances, coil_distances)
    end
    packed = vcat(distances, log10(λ_regularize))
    return packed
end

function unpack_rail(packed::Vector, symmetric::Bool, rb::IMAS.radial_build)
    distances = packed[1:end - 1]
    λ_regularize = packed[end]
    optim_coils = []
    kcoil = 0
    for rail in rb.pf_coils_rail
        r_interp = IMAS.interp(rail.outline.distance, rail.outline.r, extrapolation_bc=:flat)
        z_interp = IMAS.interp(rail.outline.distance, rail.outline.z, extrapolation_bc=:flat)
        # not symmetric
        if ! symmetric
            dkcoil = rail.coils_number
            coil_distances = distances[kcoil + 1:kcoil + dkcoil]
        # even symmetric
        elseif mod(rail.coils_number, 2) == 0
            dkcoil = Int(rail.coils_number // 2)
            coil_distances = distances[kcoil + 1:kcoil + dkcoil]
            coil_distances = vcat(- reverse(coil_distances), coil_distances)
        # odd symmetric
        else
            dkcoil = Int((rail.coils_number - 1) // 2)
            coil_distances = distances[kcoil + 1:kcoil + dkcoil]
            coil_distances = vcat(- reverse(coil_distances), 0.0, coil_distances)
        end
        kcoil += dkcoil

        # mirror coil position when they reach the end of the rail
        while any(coil_distances .< -1) || any(coil_distances .> 1)
            coil_distances[coil_distances .< -1] = -2.0 .- coil_distances[coil_distances .< -1]
            coil_distances[coil_distances .> 1] = 2.0 .- coil_distances[coil_distances .> 1]
        end

        # do not make coils cross, keep their order
        coil_distances = sort(coil_distances)

        r_coils = r_interp.(coil_distances)
        z_coils = z_interp.(coil_distances)
        append!(optim_coils, [PointCoil(r, z) for (r, z) in zip(r_coils, z_coils)])
        
    end

    return optim_coils, 10^λ_regularize
end

function optimize_coils_rail(EQfixed::Equilibrium.AbstractEquilibrium; pinned_coils::Vector, optim_coils::Vector, fixed_coils::Vector, fixed_currents::Vector, symmetric::Bool, λ_regularize::Real, λ_ψ::Real, λ_currents::Real, rb::IMAS.radial_build, maxiter::Int, verbose::Bool)
    fixed_eq = ψp_on_fixed_eq_boundary(EQfixed, fixed_coils, fixed_currents)
    packed = pack_rail(rb, λ_regularize, symmetric)

    trace = PFcoilsOptTrace()
    packed_tmp = []
    function placement_cost(packed; do_trace=false)
        push!(packed_tmp, packed)
        (optim_coils, λ_regularize) = unpack_rail(packed, symmetric, rb)
        coils = vcat(pinned_coils, optim_coils)
        currents, cost_ψ = currents_to_match_ψp(fixed_eq..., coils, λ_regularize=λ_regularize, return_cost=true)
        cost_ψ = cost_ψ / λ_ψ
        cost_currents = norm(currents) / length(currents) / λ_currents
        cost_spacing = 0
        for (k1, c1) in enumerate(optim_coils)
            for (k2, c2) in enumerate(optim_coils)
                if k1 == k2
                    continue
                end
                cost_spacing += 1 / sqrt((c1.R - c2.R)^2 + (c1.Z - c2.Z)^2)
            end
        end
        cost_spacing = cost_spacing / sum([rail.coils_number for rail in rb.pf_coils_rail])^2
        cost = sqrt(cost_ψ^2 + cost_currents^2 + cost_spacing^2)
        if do_trace
            push!(trace.currents, [no_Dual(c) for c in currents])
            push!(trace.coils, vcat(pinned_coils, [PointCoil(no_Dual(c.R), no_Dual(c.Z)) for c in optim_coils]))
            push!(trace.λ_regularize, no_Dual(λ_regularize))
            push!(trace.cost_ψ, no_Dual(cost_ψ))
            push!(trace.cost_currents, no_Dual(cost_currents))
            push!(trace.cost_bound, NaN)
            push!(trace.cost_spacing, no_Dual(cost_spacing))
            push!(trace.cost_total, no_Dual(cost))
        end
        return cost
    end

    function clb(x)
        placement_cost(packed_tmp[end]; do_trace=true)
        false
    end
    
    # use NelderMead() ; other optimizer that works is Newton(), others have trouble
    res = Optim.optimize(placement_cost, packed, Optim.NelderMead(), Optim.Options(time_limit=60 * 2, iterations=maxiter, callback=clb); autodiff=:forward)
    if verbose println(res) end
    packed = Optim.minimizer(res)

    (optim_coils, λ_regularize) = unpack_rail(packed, symmetric, rb)

    pinned_currents = [trace.currents[end][k] for k in 1:length(pinned_coils)]
    optim_currents = [trace.currents[end][k] for k in (length(pinned_coils)+1):(length(pinned_coils)+length(optim_coils))]

    return (pinned_coils, pinned_currents), (optim_coils, optim_currents), (fixed_coils, fixed_currents), λ_regularize, trace
end

function step(actor::PFcoilsOptActor;
              symmetric=actor.symmetric,
              λ_regularize=actor.λ_regularize,
              λ_ψ=actor.λ_norm,
              λ_currents=1E7,
              maxiter=10000,
              optimization_scheme=:mask,
              verbose=false)

    # convert equilibrium to Equilibrium.jl format, since this is what AD_GS uses
    time_index = 1
    EQfixed = IMAS2Equilibrium(actor.eq_in.time_slice[time_index])

    # generate coils structure as accepted by AD_GS
    fixed_coils = []
    fixed_currents = []
    pinned_coils = []
    pinned_currents = []
    optim_coils = []
    optim_currents = []
    for coil in actor.pf_active.coil
        c = PointCoil(coil.element[1].geometry.rectangle.r, coil.element[1].geometry.rectangle.z)
        if coil.identifier == "pinned"
            push!(pinned_coils, c)
            push!(pinned_currents, coil.current.data[time_index])
        elseif coil.identifier == "optim"
            push!(optim_coils, c)
            push!(optim_currents, coil.current.data[time_index])
        elseif coil.identifier == "fixed"
            push!(fixed_coils, c)
            push!(fixed_currents, coil.current.data[time_index])
        else
            error("Accepted type of coil.identifier are only \"optim\", \"pinned\", or \"fixed\"")
        end
    end
    pinned = (pinned_coils, pinned_currents)
    optim = (optim_coils, optim_currents)
    fixed = (fixed_coils, fixed_currents)

    # do nothing, simply evaluate equilibrium given existing coil currents
    if maxiter < 0
        # pass

    # find coil currents for this initial configuration
    elseif maxiter == 0
        currents, λ_norm = AD_GS.fixed_eq_currents(EQfixed, vcat(pinned_coils, optim_coils), fixed_coils, fixed_currents, λ_regularize=λ_regularize, return_cost=true)
        pinned_currents = [currents[k] for k in 1:length(pinned_coils)]
        optim_currents = [currents[k] for k in (length(pinned_coils)+1):(length(pinned_coils)+length(optim_coils))]
        pinned = (pinned_coils, pinned_currents)
        optim = (optim_coils, optim_currents)
        fixed = (fixed_coils, fixed_currents)
        actor.λ_norm = λ_norm
        trace = actor.trace

    # run optimization
    elseif maxiter > 0
        λ_ψ = actor.λ_norm
        rb = actor.radial_build
        # run mask type optimizer
        if optimization_scheme == :mask
            (pinned, optim, fixed, λ_regularize, trace) = optimize_coils_mask(EQfixed; pinned_coils, optim_coils, fixed_coils, fixed_currents, symmetric, λ_regularize, λ_ψ, λ_currents, rb, maxiter, verbose)
        # run rail type optimizer
        elseif optimization_scheme == :rail
            (pinned, optim, fixed, λ_regularize, trace) = optimize_coils_rail(EQfixed; pinned_coils, optim_coils, fixed_coils, fixed_currents, symmetric, λ_regularize, λ_ψ, λ_currents, rb, maxiter, verbose)
        else
            error("Supported PFcoilsOptActor optimization_scheme are `:mask` and `:rail`")
        end
        actor.λ_regularize = λ_regularize
        actor.trace = trace
    end

    # update ψ map
    ψ_f2f = fixed2free(EQfixed, vcat(pinned[1], optim[1], fixed[1]), vcat(pinned[2], optim[2], fixed[2]), EQfixed.r, EQfixed.z)
    actor.eq_out.time_slice[time_index].profiles_2d[1].psi = transpose(ψ_f2f)
    # IMAS.flux_surfaces(actor.eq_out.time_slice[time_index]) #### PROBLEM

    # fill coil information in actor.pf_active IDS
    pf_active = actor.pf_active
    k = 1
    for (coiltype, (coils, currents)) in [("pinned",pinned), ("optim",optim), ("fixed",fixed)]
        for (coil, current) in zip(coils, currents)
            pf_active.coil[k].identifier = coiltype
            pf_active.coil[k].element[1].geometry.rectangle.r = coil.R
            pf_active.coil[k].element[1].geometry.rectangle.z = coil.Z
            pf_active.coil[k].element[1].geometry.rectangle.width = maximum(coil.R) - minimum(coil.R)
            pf_active.coil[k].element[1].geometry.rectangle.height = maximum(coil.Z) - minimum(coil.Z)
            set_field_time_array(pf_active.coil[k].current, :time, 1, 0.0)
            set_field_time_array(pf_active.coil[k].current, :data, 1, current)
            k+=1
        end
    end

    # update PFcoilsOptActor fields
    actor.pf_active = pf_active

    return actor
end

#= ====== =#
# PLOTTING #
#= ====== =#
"""
    function plot_pfcoilsactor_cx(pfactor::PFcoilsOptActor; trace=false, mask=false)

Plot PFcoilsOptActor optimization cross-section
"""
@recipe function plot_pfcoilsactor_cx(pfactor::PFcoilsOptActor; equilibrium=true, trace=false, mask=false, rail=true)

    # setup plotting area
    rb = pfactor.radial_build
    xlim --> [0.0,maximum(rb.layer[end].outline.r) * 2.0]
    ylim --> [minimum(rb.layer[end].outline.z),maximum(rb.layer[end].outline.z)]
    aspect_ratio --> :equal

    # plot optimization mask
    if mask
        rmask, zmask, mask = IMAS.structures_mask(pfactor.radial_build)
        cl = Contour.contour(rmask, zmask, mask, 0.5)
        for line in Contour.lines(cl)
            @series begin
                label --> ""
                seriescolor --> :gray
                linewidth --> 3
                Contour.coordinates(line)
            end
        end
    end

    # plot optimization rails
    if rail
        for (krail, rail) in enumerate(rb.pf_coils_rail)
            if ! is_missing(rail.outline,:r)
                @series begin
                    label --> "coil optimization rail"
                    primary --> krail == 1 ? true : false
                    color --> :gray
                    linestyle --> :dash
                    rail.outline.r, rail.outline.z
                end
            end
        end
    end

    # plot pf_active coils
    @series begin
        pfactor.pf_active
    end

    # plot final equilibrium
    if equilibrium
        @series begin
            label --> "Final"
            seriescolor --> :red
            pfactor.eq_out.time_slice[1]
        end

        # plot target equilibrium
        @series begin
            label --> "Target"
            seriescolor --> :black
            lcfs --> true
            linestyle --> :dash
            pfactor.eq_in.time_slice[1]
        end
    end

    if trace
        if length(pfactor.trace.coils) > 0
            for c in 1:length(pfactor.trace.coils[1])
            @series begin
                    label --> ""
                    linewidth --> 1
                    seriesalpha --> 0.5
                    primary --> false
                    seriescolor --> :magenta
                    [no_Dual(k[c].R) for k in pfactor.trace.coils], [no_Dual(k[c].Z) for k in pfactor.trace.coils]
                end
            end
        end
    end

end

"""
    plot_pfcoilsactor_trace(trace::PFcoilsOptTrace, what::Symbol=:cost; start_at::Int=1)

Plot PFcoilsOptActor optimization trace

Attributes:
- what::Symbol=:cost or :currents or individual fields of the PFcoilsOptTrace structure
- start_at=::Int=1 index of the first element of the trace to start plotting
"""
@recipe function plot_pfcoilsactor_trace(trace::PFcoilsOptTrace, what::Symbol=:cost; start_at=1)
    x = (start_at:length(trace.cost_total))
    legend --> :bottomleft
    if what == :cost
        @series begin
            label --> "ψ"
            yscale --> :log10
            x, trace.cost_ψ[start_at:end]
        end
        @series begin
            label --> "currents"
            yscale --> :log10
            x, trace.cost_currents[start_at:end]
        end
        @series begin
            label --> "bounds"
            yscale --> :log10
            x, trace.cost_bound[start_at:end]
        end
        @series begin
            label --> "spacing"
            yscale --> :log10
            x, trace.cost_spacing[start_at:end]
        end
        @series begin
            label --> "total"
            yscale --> :log10
            linestyle --> :dash
            color --> :black
            # ylim --> [minimum(trace.cost_total[start_at:end]) / 10,maximum(trace.cost_total[start_at:end])]
            x, trace.cost_total[start_at:end]
        end

    elseif what == :starting_currents
        @series begin
        label --> "Starting"
        getfield(trace, :currents)[start_at:end][1,:]
    end

    elseif what == :final_currents
            @series begin
            label --> "Final"
            getfield(trace, :currents)[start_at:end][end,:]
        end

    elseif what == :currents
        @series begin
            label --> "Starting"
            getfield(trace, what)[start_at:end][1,:]
        end
        @series begin
            label --> "Final"
            getfield(trace, what)[start_at:end][end,:]
        end

    else
        @series begin
            if occursin("cost_", String(what))
                yscale --> :log10
            end
            label --> String(what)
            x, getfield(trace, what)[start_at:end]
        end
    end
end
