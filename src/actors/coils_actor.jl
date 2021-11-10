@Base.kwdef mutable struct PFcoilsOptTrace
    coils::Vector = []
    currents::Vector = []
    λ_regularize::Vector = []
    cost_ψ::Vector = []
    cost_currents::Vector = []
    cost_bound::Vector = []
    cost_distance::Vector = []
    cost_total::Vector = []
end

mutable struct PFcoilsOptActor <: CoilsActor
    eq_in::IMAS.equilibrium
    eq_out::IMAS.equilibrium
    time::Real
    pf_active::IMAS.pf_active
    rmask::AbstractRange
    zmask::AbstractRange
    mask_log_interpolant
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

function unwrap(v, inplace=false)
    # currently assuming an array
    unwrapped = inplace ? v : copy(v)
    for i in 2:length(v)
        while (unwrapped[i] - unwrapped[i - 1] >= pi)
            unwrapped[i] -= 2pi
        end
        while (unwrapped[i] - unwrapped[i - 1] <= -pi)
            unwrapped[i] += 2pi
    end
    end
    return unwrapped
end
  
function atan_eq(r, z, r0, z0)
    if r[1] == r[end] && z[1] == z[end]
        r = r[1:end - 1]
        z = z[1:end - 1]
    end
    θ = unwrap(atan.(z .- z0, r .- r0))
    if θ[2] < θ[1]
        r = reverse(r)
        z = reverse(z)
        θ = reverse(θ)
    end
    return r, z, θ
end
  
function two_curves_same_θ(r1, z1, r2, z2, scheme=:cubic)
    r0 = (sum(r1) / length(r1) + sum(r2) / length(r2)) / 2.0
    z0 = (sum(z1) / length(z1) + sum(z2) / length(z2)) / 2.0
    r1, z1, θ1 = atan_eq(r1, z1, r0, z0)
    r2, z2, θ2 = atan_eq(r2, z2, r0, z0)
    if length(θ2) > length(θ1)
        r1 = IMAS.interp(vcat(θ1 .- 2 * π, θ1, θ1 .+ 2 * π), vcat(r1, r1, r1), scheme=scheme).(θ2)
        z1 = IMAS.interp(vcat(θ1 .- 2 * π, θ1, θ1 .+ 2 * π), vcat(z1, z1, z1), scheme=scheme).(θ2)
        θ = θ2
    else
        r2 = IMAS.interp(vcat(θ2 .- 2 * π, θ2, θ2 .+ 2 * π), vcat(r2, r2, r2), scheme=scheme).(θ1)
        z2 = IMAS.interp(vcat(θ2 .- 2 * π, θ2, θ2 .+ 2 * π), vcat(z2, z2, z2), scheme=scheme).(θ1)
        θ = θ1
    end
    return r1, z1, r2, z2, θ
end
  
  
function equispace_coils(rb, rmask, zmask, mask, ncoils_per_region)
    coils_in_regions = []
    region = 0
    for (k, layer) in enumerate(vcat(rb.layer[1:end]))
        if (layer.hfs == -1 || k == length(rb.layer)) && ! is_missing(layer.outline, :r)
            if ! is_missing(layer, :material) && layer.material == "vacuum"
                # pick layers with outline information
                if layer.hfs == -1
                    inner_layer = IMAS.get_radial_build(rb, identifier=rb.layer[k].identifier, hfs=-1)
                    outer_layer = IMAS.get_radial_build(rb, identifier=rb.layer[k + 1].identifier, hfs=[-1,0])
                else
                    inner_layer = IMAS.get_radial_build(rb, identifier=rb.layer[k - 1].identifier, hfs=-1)
                    outer_layer = IMAS.get_radial_build(rb, identifier=rb.layer[k].identifier, hfs=[-1,0])
                end

                # take two outlines and interpolate them on the same θ
                inner_r, inner_z, outer_r, outer_z, θ = two_curves_same_θ(inner_layer.outline.r, inner_layer.outline.z, outer_layer.outline.r, outer_layer.outline.z)

                # generate the average ouline between the two layers
                mid_r = (inner_r .+ outer_r) ./ 2
                mid_z = (inner_z .+ outer_z) ./ 2

                # mark what regions on that mid-line are valid to hold coils
                n = 7 # this was found to be reasonable on a 257 mask grid
                n = Int(ceil(n * length(rmask) / 257))
                valid_k = []
                for (k, (r, z)) in enumerate(zip(mid_r, mid_z))
                    ir = argmin(abs.(rmask .- r))
                    iz = argmin(abs.(zmask .- z))
                    if (ir - n) < 1 || (ir + n) > length(rmask) || (iz - n) < 1 || (iz + n) > length(zmask)
                        continue
                    end
                    if all(mask[(-n:n) .+ ir,(-n:n) .+ iz] .== 1)
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

                region += 1
                ncoils = ncoils_per_region[region]

                # evaluate distance along valid mid-line
                d_distance = sqrt.(diff(vcat(valid_r, valid_r[1])).^2.0 .+ diff(vcat(valid_z, valid_z[1])).^2.0)
                d_distance[isnan.(d_distance)] .= 1E-6
                distance = cumsum(d_distance)
                distance = (distance .- distance[1])
                distance = (distance ./ distance[end]) .* (ncoils)

                # uniformely distribute coils
                r_coils = IMAS.interp(distance, valid_r)((1:ncoils) .- 0.5)
                z_coils = IMAS.interp(distance, valid_z)((1:ncoils) .- 0.5)
                push!(coils_in_regions, collect(zip(r_coils, z_coils)))

            end
        end 
    end
    # p=plot(rb)
    # for coils_in_region in coils_in_regions
    #     for coil in coils_in_region
    #         scatter!([coil[1]],[coil[2]],primary=false)
    #     end
    # end
    # display(p)
    return coils_in_regions
end

function PFcoilsOptActor(eq_in::IMAS.equilibrium, rb::IMAS.radial_build, ncoils_OH::Int, ncoils_per_region::Vector{Int}, λ_regularize=1E-13, )

    # time_index = get_time_index(eq_in.time_slice, time)
    # eqt = eq_in.time_slice[time_index]

    time_index = 1
    time = eq_in.time[time_index]
    eqt = eq_in.time_slice[time_index]

    OH_layer = IMAS.get_radial_build(rb, type=1)
    r_ohcoils = ones(ncoils_OH) .* (sum(extrema(OH_layer.outline.r)) / 2.)
    z_ohcoils = collect(range(minimum(OH_layer.outline.z), maximum(OH_layer.outline.z), length=ncoils_OH))    
    fixed_coils = [PointCoil(r, z) for (r, z) in zip(r_ohcoils, z_ohcoils)]

    rmask, zmask, mask = IMAS.structures_mask(rb)

    coils_in_regions = equispace_coils(rb, rmask, zmask, mask, ncoils_per_region)
    optim_coils = []
    for coils_in_region in coils_in_regions
        for coil in coils_in_region
            push!(optim_coils, PointCoil(coil[1], coil[2]))
        end
    end

    coils = vcat(fixed_coils, optim_coils)

    mask = 1.0 .- mask

    # # apply filtering to smooth transition between allowed and forbidden regions
    # # nx = 1
    # # ny = Int(ceil(nx * size(mask)[1] / size(mask)[2]))
    # # filterx = exp.(-(range(-1, 1, length=2 * nx + 1) / nx).^2)
    # # filterx = filterx ./ sum(filterx)
    # # filtery = exp.(-(range(-1, 1, length=2 * ny + 1) / ny).^2)
    # # filtery = filtery ./ sum(filtery)
    # # mask = DSP.conv(filterx, filtery, mask)[nx + 1:end - nx, ny + 1:end - ny]

    # Cubic spline interpolation on the log to ensure positivity of the cost
    mask_log_interpolant_raw = Interpolations.CubicSplineInterpolation((rmask, zmask), log10.(1E-1 .+ mask))
    mask_log_interpolant_raw = Interpolations.extrapolate(mask_log_interpolant_raw.itp, Interpolations.Flat());
    function mask_log_interpolant(r, z)
        return (10.0.^(mask_log_interpolant_raw(r, z)) .- 1E-1)
    end
    
    # find coil currents for this initial configuration
    EQfixed = IMAS2Equilibrium(eqt)
    currents, λ_norm = AD_GS.fixed_eq_currents(EQfixed, coils, λ_regularize=λ_regularize, λ_minimize=0.0, λ_zerosum=0.0, return_cost=true)

    # update psirz based on coil configuration
    eq_out = deepcopy(eq_in)

    # ψ from fixed-boundary gEQDSK
    # make ψ at boundary zero, and very small value outside for plotting
    ψ0_fix, ψb_fix = psi_limits(EQfixed)
    σ₀ = sign(ψ0_fix - ψb_fix)
    ψ_fix = [EQfixed(r, z) for z in EQfixed.z, r in EQfixed.r] .- ψb_fix
    ψ_fix = ifelse.(σ₀ * ψ_fix .> 0, ψ_fix, 1e-6 * ψ_fix)
    
    # ψ at the boundary is determined by the value of the currents
    # calculated in fixed_eq_currents
    ψ_f2f = fixed2free(EQfixed, coils, currents, EQfixed.r, EQfixed.z)
    eq_out.time_slice[time_index].profiles_2d[1].psi = transpose(ψ_f2f)

    # populate IMAS data structure
    pf_active = IMAS.pf_active()
    resize!(pf_active.coil, length(coils))
    for (k, c) in enumerate(coils)
        resize!(pf_active.coil[k].element, 1)
        if c in fixed_coils
            pf_active.coil[k].identifier = "fixed"
        else
            pf_active.coil[k].identifier = "optim"
        end
        pf_active.coil[k].element[1].geometry.rectangle.r = c.R
        pf_active.coil[k].element[1].geometry.rectangle.z = c.Z
        pf_active.coil[k].element[1].geometry.rectangle.width = maximum(c.R) - minimum(c.R)
        pf_active.coil[k].element[1].geometry.rectangle.height = maximum(c.Z) - minimum(c.Z)
        set_field_time_array(pf_active.coil[k].current, :time, 1, 0.0)
        set_field_time_array(pf_active.coil[k].current, :data, 1, currents[k])
    end

    # detect if initial coil configuration was symmetric
    symmetric = sum([c.Z for c in coils]) / length(coils) < 1E-3

    # constructor
    PFcoilsOptActor(eq_in, eq_out, time, pf_active, rmask, zmask, mask_log_interpolant, symmetric, λ_regularize, λ_norm, PFcoilsOptTrace())
end

#= == =#
# STEP #
#= == =#
function step(actor::PFcoilsOptActor;
              symmetric=actor.symmetric,
              λ_regularize=actor.λ_regularize,
              λ_currents=1E7,
              verbose=false)

    # generate coils structure as accepted by AD_GS
    coils = []
    for coil in actor.pf_active.coil
        push!(coils, PointCoil(coil.element[1].geometry.rectangle.r, coil.element[1].geometry.rectangle.z))
    end

    # utility functions for packing and unpacking info in/out of optimization function
    function pack(coils, λ_regularize)
        coilz = vcat([c.R for c in coils if ((! symmetric) || c.Z >= 0)], [c.Z for c in coils if ((! symmetric) || c.Z >= 0)])
        packed = vcat(coilz, log10(λ_regularize))
        return packed
    end
    function unpack(packed)
        coilz = packed[1:end - 1]
        λ_regularize = packed[end]
        coils =  [PointCoil(coilz[k], coilz[k + Int(length(coilz) / 2)]) for (k, c) in enumerate(coilz[1:Int(end / 2)])]
        if symmetric
            coils = vcat(coils, [PointCoil(c.R, -c.Z) for c in coils if c.Z != 0])
        end
        return coils, 10^λ_regularize
    end

    # apply packing/unpacknig function
    packed = pack(coils, λ_regularize)
    (coils, λ_regularize) = unpack(packed)

    function optimize_coils(S, coils, λ_regularize, ψp_cost_norm, currents_cost_norm, mask_interpolant)
        fixed_eq = ψp_on_fixed_eq_boundary(S)
        packed = pack(coils, λ_regularize)
        
        trace = PFcoilsOptTrace()
        packed_tmp = []
        function placement_cost(packed; do_trace=false)
            push!(packed_tmp, packed)
            (coils, λ_regularize) = unpack(packed)
            currents, cost = currents_to_match_ψp(fixed_eq..., coils, λ_regularize=λ_regularize, λ_minimize=0.0, λ_zerosum=0.0, return_cost=true)
            cost_ψ = cost / ψp_cost_norm
            cost_currents = norm(currents) / length(currents) / currents_cost_norm
            cost_bound = norm(mask_interpolant.([c.R for c in coils], [c.Z for c in coils]))
            cx = [c.R for c in coils]
            cy = [c.Z for c in coils]
            distance_matrix = sqrt.((repeat(cx, 1, length(cx)) .- repeat(transpose(cx), length(cx), 1)).^2.0 .+ (repeat(cy, 1, length(cy)) .- repeat(transpose(cy), length(cy), 1)).^2.0)
            cost_distance = norm(1.0 ./ (distance_matrix + LinearAlgebra.I * 1E3)) / length(coils).^2
            cost = sqrt(cost_ψ^2 + cost_currents^2 + cost_bound^2 + cost_distance^2)
            if do_trace
                push!(trace.currents, [no_Dual(c) for c in currents])
                push!(trace.coils, [PointCoil(no_Dual(c.R), no_Dual(c.Z)) for c in coils])
                push!(trace.λ_regularize, no_Dual(λ_regularize))
                push!(trace.cost_ψ, no_Dual(cost_ψ))
                push!(trace.cost_currents, no_Dual(cost_currents))
                push!(trace.cost_bound, no_Dual(cost_bound))
                push!(trace.cost_distance, no_Dual(cost_distance))
                push!(trace.cost_total, no_Dual(cost))
            end
            return cost
        end

        function clb(x)
            placement_cost(packed_tmp[end]; do_trace=true)
            false
        end
        
        # use NelderMead() ; other optimizer that works is Newton(), others have trouble
        res = Optim.optimize(placement_cost, packed, Optim.NelderMead(), Optim.Options(time_limit=30, iterations=10000, allow_f_increases=true, successive_f_tol=100, callback=clb); autodiff=:forward)

        if verbose println(res) end
        packed = Optim.minimizer(res)
        (coils, λ_regularize) = unpack(packed)

        return coils, λ_regularize, trace
    end

    # run optimization
    EQfixed = IMAS2Equilibrium(actor.eq_in.time_slice[time_index])
    (coils, λ_regularize, trace) = optimize_coils(EQfixed, coils, λ_regularize, actor.λ_norm, λ_currents, actor.mask_log_interpolant)
    currents = [trace.currents[end][k] for (k, c) in enumerate(coils)]

    # ψ from fixed-boundary gEQDSK
    # make ψ at boundary zero, and very small value outside for plotting
    ψ0_fix, ψb_fix = psi_limits(EQfixed)
    σ₀ = sign(ψ0_fix - ψb_fix)
    ψ_fix = [EQfixed(r, z) for z in EQfixed.z, r in EQfixed.r] .- ψb_fix
    ψ_fix = ifelse.(σ₀ * ψ_fix .> 0, ψ_fix, 1e-6 * ψ_fix)
    
    # ψ at the boundary is determined by the value of the currents
    # calculated in fixed_eq_currents
    ψ_f2f = fixed2free(EQfixed, coils, currents, EQfixed.r, EQfixed.z)
    actor.eq_out.time_slice[time_index].profiles_2d[1].psi = transpose(ψ_f2f)

    # populate IMAS data structure
    pf_active = IMAS.pf_active()
    resize!(pf_active.coil, length(coils))
    for (k, c) in enumerate(coils)
        resize!(pf_active.coil[k].element, 1)
        pf_active.coil[k].element[1].geometry.rectangle.r = c.R
        pf_active.coil[k].element[1].geometry.rectangle.z = c.Z
        pf_active.coil[k].element[1].geometry.rectangle.width = maximum(c.R) - minimum(c.R)
        pf_active.coil[k].element[1].geometry.rectangle.height = maximum(c.Z) - minimum(c.Z)
        set_field_time_array(pf_active.coil[k].current, :time, 1, 0.0)
        set_field_time_array(pf_active.coil[k].current, :data, 1, currents[k])
    end

    # update PFcoilsOptActor fields
    actor.pf_active = pf_active
    actor.λ_regularize = λ_regularize
    actor.trace = trace

    return trace
end

#= ====== =#
# PLOTTING #
#= ====== =#
"""
    plot_pfcoilsactor_cx(pfactor::PFcoilsOptActor)

Plot PFcoilsOptActor optimization cross-section
"""
@recipe function plot_pfcoilsactor_cx(pfactor::PFcoilsOptActor, trace=true)
    # plot mask
    rmask = pfactor.rmask
    zmask = pfactor.zmask
    dst = pfactor.mask_log_interpolant(rmask, zmask)

    xlims --> [rmask[1],rmask[end] * 2.0]
    ylims --> [zmask[1],zmask[end]]
    aspect_ratio --> :equal

    cl = Contour.contour(rmask, zmask, dst, 0.5)
    for line in Contour.lines(cl)
        @series begin
            label --> ""
            seriescolor --> :gray
            linewidth --> 3
            Contour.coordinates(line)
        end
    end

    # plot pf_active coils
    @series pfactor.pf_active
    
    # plot target equilibrium
    @series begin
        label --> "Target"
        seriescolor --> :black
        lcfs --> true
        pfactor.eq_in.time_slice[1]
    end

    # plot final equilibrium
    @series begin
        label --> "Final"
        seriescolor --> :red
        pfactor.eq_out.time_slice[1]
    end

    if trace
        if length(pfactor.trace.coils) > 0
            for c in 1:length(pfactor.trace.coils[1])
                @series begin
                    label --> ""
                    linewidth --> 1
                    seriesalpha --> 0.5
                    primary --> false
                    seriescolor --> :blue
                    [no_Dual(k[c].R) for k in pfactor.trace.coils], [no_Dual(k[c].Z) for k in pfactor.trace.coils]
                end
            end
        end
    end

end

"""
    function plot_pfcoilsactor_trace(trace::PFcoilsOptTrace, what::Symbol=:cost; start_at::Int=1)

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
            label --> "distance"
            yscale --> :log10
            x, trace.cost_distance[start_at:end]
        end
        @series begin
            label --> "total"
            yscale --> :log10
            x, trace.cost_total[start_at:end]
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
            label --> String(what)
            x, getfield(trace, what)[start_at:end]
        end
    end
end
