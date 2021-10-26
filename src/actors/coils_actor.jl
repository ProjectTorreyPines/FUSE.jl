@Base.kwdef mutable struct PFcoilsOptTrace
    coils::Vector = []
    currents::Vector = []
    λ_regularize::Vector = []
    cost_ψ::Vector = []
    cost_currents::Vector = []
    cost_bound::Vector = []
    cost::Vector = []
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


# The PFcoilsOptActor should eventually also take IMAS.wall, IMAS.tf, IMAS.cryostat IDSs as an inputs
function PFcoilsOptActor(eq_in::IMAS.equilibrium,
                         time::Real,
                         ncoils::Int,
                         coils_outsideof=1.3,
                         coils_insideof=nothing,
                         λ_regularize=1E-13,)
    time_index = get_time_index(eq_in.time_slice, time)
    eqt = eq_in.time_slice[time_index]
    
    if coils_insideof === nothing
        # define extent of computation domain where the placing of pf coils could take place
        xlim = extrema(eqt.profiles_2d[1].grid.dim1)
        ylim = extrema(eqt.profiles_2d[1].grid.dim2)
        coils_insideof = zip([xlim[1],xlim[2],xlim[2],xlim[1],xlim[1]], [ylim[2],ylim[2],ylim[1],ylim[1],ylim[2]])
    else
        xlim = extrema([p[1] for p in coils_insideof])
        ylim = extrema([p[2] for p in coils_insideof])
    end

    # grid the coputation domain 
    resolution = 257
    rmask = range(xlim[1], xlim[2], length=resolution)
    zmask = range(ylim[1], ylim[2], length=resolution * Int(round((ylim[2] - ylim[1] / (xlim[2] - xlim[1])))))
    mask = ones(length(rmask), length(zmask))

    # outer domain (this will be either the cryostat or the TF coils depending on whether the PFs are inside or outside the TFs)
    coils_insideof_array = StaticArrays.SVector.([p[1] for p in coils_insideof], [p[2] for p in coils_insideof])
    for (kr, rr) in enumerate(rmask)
        for (kz, zz) in enumerate(zmask)
            if PolygonOps.inpolygon((rr, zz), coils_insideof_array) == 1
                mask[kr,kz] = 0.0
            end
        end
    end
    
    # geometric center
    R0 = eqt.boundary.geometric_axis.r[end]

    if typeof(coils_outsideof) <: Number
        coils_outsideof = (coils_outsideof, coils_outsideof)
    end
    if typeof(coils_outsideof) <: Tuple
        # plasma boundary
        ψb = eqt.global_quantities.psi_boundary
        ψ0 = eqt.global_quantities.psi_axis
        pr, pz = IMAS.flux_surface(eqt, (ψb - ψ0) * 0.5 + ψ0)

        mr = Statistics.mean(pr)
        dx = maximum(pr) - minimum(pr)
        dz = maximum(pz) - minimum(pz)

        # coils will not be closer to the boundary than this
        coils_outsideof = zip((pr .- mr) .* (1 + coils_outsideof[1] / dx) .+ R0, pz .* (1 + coils_outsideof[2] / dx))
    end
    
    # forbidden plasma area
    coils_outsideof_array = StaticArrays.SVector.([p[1] for p in coils_outsideof], [p[2] for p in coils_outsideof])
    for (kr, rr) in enumerate(rmask)
        for (kz, zz) in enumerate(zmask)
            if PolygonOps.inpolygon((rr, zz), coils_outsideof_array) == 1
                mask[kr,kz] = 1.0
            end
        end
    end

    # apply filtering to smooth transition between allowed and forbidden regions
    nx = 7
    ny = Int(ceil(nx * size(mask)[1] / size(mask)[2]))
    filterx = exp.(-(range(-1, 1, length=2 * nx + 1) / nx).^2)
    filterx = filterx ./ sum(filterx)
    filtery = exp.(-(range(-1, 1, length=2 * ny + 1) / ny).^2)
    filtery = filtery ./ sum(filtery)
    mask = DSP.conv(filterx, filtery, mask)[nx + 1:end - nx, ny + 1:end - ny]

    # never allow the coils to leave the computation domain
    mask[1,2:end] .= 1.0
    mask[end,2:end] .= 1.0
    mask[2:end,1] .= 1.0
    mask[2:end,end] .= 1.0

    # Cubic spline interpolation on the log to ensure positivity of the cost
    mask_log_interpolant_raw = Interpolations.CubicSplineInterpolation((rmask, zmask), log10.(1.0 .+ mask))
    mask_log_interpolant_raw = Interpolations.extrapolate(mask_log_interpolant_raw.itp, Interpolations.Flat());
    function mask_log_interpolant(r, z)
        return 10.0.^(mask_log_interpolant_raw(r, z)) .- 1
    end
    
    # initial position of coils
    pr_coil = ([p[1] for p in coils_outsideof] .- R0) .* 1.1 .+ R0
    pz_coil = [p[2] for p in coils_outsideof] .* 1.2

    # initial coils placement is uniformely spaced along a surface conforming to plasma
    coils = []
    θ_coil = atan.(pz_coil, pr_coil .- R0)
    θ_coil = DSP.unwrap(θ_coil)
    if θ_coil[2] < θ_coil[1]
        θ_coil = -θ_coil
    end
    pr_coilθ = Interpolations.extrapolate(Interpolations.interpolate((θ_coil,), pr_coil, Interpolations.Gridded(Interpolations.Linear())), Interpolations.Periodic())
    pz_coilθ = Interpolations.extrapolate(Interpolations.interpolate((θ_coil,), pz_coil, Interpolations.Gridded(Interpolations.Linear())), Interpolations.Periodic())
    θ = range(π / 4, 2π - π / 4, length=1001)[1:end - 1]
    l2p = vcat(0, cumsum(sqrt.(diff(pr_coilθ(θ)).^2.0 .+ diff(pz_coilθ(θ)).^2.0)))
    lcoils = range(0, l2p[end], length=ncoils + 1)[1:end - 1]
    pr_coilL = Interpolations.extrapolate(Interpolations.interpolate((l2p,), pr_coilθ(θ), Interpolations.Gridded(Interpolations.Linear())), Interpolations.Periodic())
    pz_coilL = Interpolations.extrapolate(Interpolations.interpolate((l2p,), pz_coilθ(θ), Interpolations.Gridded(Interpolations.Linear())), Interpolations.Periodic())
    if mod(ncoils, 2) == 0
        start_coil_L = x -> abs.(pz_coilL.(lcoils[1] .+ x) .+ pz_coilL.(lcoils[end] .+ x))[1]
    else
        start_coil_L = x -> abs.(pz_coilL.(lcoils[1] .+ x))[1]
    end
    L0 = Optim.minimizer(Optim.optimize(start_coil_L, 0.0, l2p[end]))
    coils = [(pr_coilL(L), pz_coilL(L)) for L in lcoils .+ L0]
    if mod(ncoils, 2) == 1
        coils[1] = (coils[1][1], 0.0)
    end
    
    # find coil currents for this initial configuration
    EQfixed = IMAS2Equilibrium(eq_in, time)
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
    eq_out.time_slice[1].profiles_2d[1].psi = transpose(ψ_f2f)

    # populate IMAS data structure
    pf_active = IMAS.pf_active()
    resize!(pf_active.coil, length(coils))
    for (k, c) in enumerate(coils)
        resize!(pf_active.coil[k].element, 1)
        pf_active.coil[k].element[1].geometry.rectangle.r = c[1]
        pf_active.coil[k].element[1].geometry.rectangle.z = c[2]
        pf_active.coil[k].element[1].geometry.rectangle.height = 0.0
        pf_active.coil[k].element[1].geometry.rectangle.width = 0.0
        set_field_time_array(pf_active.coil[k].current, :time, 1, 0.0)
        set_field_time_array(pf_active.coil[k].current, :data, 1, currents[k])
    end

    # detect if initial coil configuration was symmetric
    symmetric = sum([c[2] for c in coils]) / length(coils) < 1E-3

    # constructor
    PFcoilsOptActor(eq_in, eq_out, time, pf_active, rmask, zmask, mask_log_interpolant, symmetric, λ_regularize, λ_norm, PFcoilsOptTrace())
end

#= == =#
# STEP #
#= == =#
function Base.step(actor::PFcoilsOptActor;
                   symmetric=actor.symmetric,
                   λ_regularize=actor.λ_regularize,
                   λ_currents=1E7,
                   verbose=false)

    # generate coils structure as accepted by AD_GS
    coils = []
    for coil in actor.pf_active.coil
        push!(coils, (coil.element[1].geometry.rectangle.r, coil.element[1].geometry.rectangle.z))
    end

    # utility functions for packing and unpacking info in/out of optimization function
    function pack(coils, λ_regularize)
        coilz = vcat([c[1] for c in coils if ((! symmetric) || c[2] >= 0)], [c[2] for c in coils if ((! symmetric) || c[2] >= 0)])
        packed = vcat(coilz, log10(λ_regularize))
        return packed
    end
    function unpack(packed)
        coilz = packed[1:end - 1]
        λ_regularize = packed[end]
        coils =  [(coilz[k], coilz[k + Int(length(coilz) / 2)]) for (k, c) in enumerate(coilz[1:Int(end / 2)])]
        if symmetric
            coils = vcat(coils, [(c[1], -c[2]) for c in coils if c[2] != 0])
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
            cost_bound = norm(mask_interpolant.([c[1] for c in coils], [c[2] for c in coils])) * 10.0
            cost = sqrt.(cost_ψ.^2 + cost_currents.^2 + cost_bound.^2)
            cx = [c[1] for c in coils]
            cy = [c[2] for c in coils]
            coil_distance_cost = 0.1 / minimum(sqrt.((repeat(cx, 1, length(cx)) .- repeat(transpose(cx), length(cx), 1)).^2.0 .+ (repeat(cy, 1, length(cx)) .- repeat(transpose(cy), length(cx), 1)).^2) + LinearAlgebra.I * 1E3)
            cost += coil_distance_cost
            if do_trace
                push!(trace.currents, currents)
                push!(trace.coils, coils)
                push!(trace.λ_regularize, λ_regularize)
                push!(trace.cost_ψ, cost_ψ)
                push!(trace.cost_currents, cost_currents)
                push!(trace.cost_bound, cost_bound)
                push!(trace.cost, cost)
            end
            return cost
        end

        function clb(x)
            placement_cost(packed_tmp[end];do_trace=true)
            false
        end
        
        # use NelderMead() ; other optimizer that works is Newton(), others have trouble
        res = Optim.optimize(placement_cost, packed, Optim.NelderMead(), Optim.Options(time_limit=30, iterations=10000, allow_f_increases=true, successive_f_tol=100, callback=clb); autodiff=:forward)
        packed = Optim.minimizer(res)
        (coils, λ_regularize) = unpack(packed)
        if verbose println(res) end

        return coils, λ_regularize, trace
    end

    # run optimization
    EQfixed = IMAS2Equilibrium(actor.eq_in, actor.time)
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
    actor.eq_out.time_slice[1].profiles_2d[1].psi = transpose(ψ_f2f)

    # populate IMAS data structure
    pf_active = IMAS.pf_active()
    resize!(pf_active.coil, length(coils))
    for (k, c) in enumerate(coils)
        resize!(pf_active.coil[k].element, 1)
        pf_active.coil[k].element[1].geometry.rectangle.r = c[1]
        pf_active.coil[k].element[1].geometry.rectangle.z = c[2]
        pf_active.coil[k].element[1].geometry.rectangle.height = 0.0
        pf_active.coil[k].element[1].geometry.rectangle.width = 0.0
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

    xlims --> [rmask[1],rmask[end]]
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
                    [FUSE.no_Dual(k[c][1]) for k in pfactor.trace.coils], [FUSE.no_Dual(k[c][2]) for k in pfactor.trace.coils]
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
@recipe function plot_pfcoilsactor_trace(trace::PFcoilsOptTrace, what::Symbol=:cost; start_at::Int=1)
    x = (start_at:length(trace.cost))
    if what == :cost
        @series begin
            label --> "ψ"
            yscale --> :log10
            x, [FUSE.no_Dual(y) for y in trace.cost_ψ[start_at:end]]
        end
        @series begin
            label --> "currents"
            yscale --> :log10
            x, [FUSE.no_Dual(y) for y in trace.cost_currents[start_at:end]]
        end
        @series begin
            label --> "bounds"
            yscale --> :log10
            x, [FUSE.no_Dual(y) for y in trace.cost_bound[start_at:end]]
        end
        @series begin
            label --> "total"
            yscale --> :log10
            x, [FUSE.no_Dual(y) for y in trace.cost[start_at:end]]
        end
            
    elseif what == :currents
        @series begin
            label --> "Starting"
        [FUSE.no_Dual(y) for y in getfield(trace, what)[start_at:end]][1,:]
        end
        @series begin
            label --> "Final"
            [FUSE.no_Dual(y) for y in getfield(trace, what)[start_at:end]][end,:]
        end

    else
        @series begin
            label --> String(what)
            x, [FUSE.no_Dual(y) for y in getfield(trace, what)[start_at:end]]
        end
    end
end
