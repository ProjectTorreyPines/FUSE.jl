mutable struct PFcoilsActor <: CoilsActor
    equilibrium::IMAS.equilibrium
    pf_active::IMAS.pf_active
    mask_log_interpolant
end

#= == =#
# INIT #
#= == =#
using PolygonOps
using StaticArrays
using DSP
using Random
using Interpolations
using Optim

# The PFcoilsActor should eventually take IMAS.wall, IMAS.tf, IMAS.cryostat IDSs as an inputs
function PFcoilsActor(equilibrium::IMAS.equilibrium, ncoils::Int)
    eqt = equilibrium.time_slice[1]
    
    # define extent of computation domain where the placing of pf coils could take place
    xlim = eqt.profiles_2d[1].grid.dim1
    ylim = eqt.profiles_2d[1].grid.dim2
    dx = 0.3
    xlim = (maximum([0,minimum(xlim) - (xlim[2] - xlim[1]) * dx]), maximum(xlim) + (xlim[2] - xlim[1]) * dx)
    ylim = (minimum(ylim) - (ylim[2] - ylim[1]) * dx, maximum(ylim) + (ylim[2] - ylim[1]) * dx);

    # grid the coputation domain 
    resolution = 257
    r = range(xlim[1] - dx, xlim[2] + dx, length=resolution)
    z = range(ylim[1] - dx, ylim[2] + dx, length=resolution * Int(round((ylim[2] - ylim[1] / (xlim[2] - xlim[1])))))
    pts = [((kr, kz), (rr, zz)) for (kz, zz) in enumerate(z), (kr, rr) in enumerate(r)]
    mask = zeros(size(pts)...) .+ 1
    
    # outer domain (this will be either the cryostat or the TF coils depending on whether the PFs are inside or outside the TFs)
    pr_outer = [xlim[1],xlim[2],xlim[2],xlim[1],xlim[1]]
    pz_outer = [ylim[2],ylim[2],ylim[1],ylim[1],ylim[2]]
    for ((kr, kz), (rr, zz)) in hcat(pts...)
        point = SA[rr,zz]
        if PolygonOps.inpolygon((rr, zz), StaticArrays.SVector.(pr_outer, pz_outer)) == 1
            mask[kz,kr] = 0
        end
    end
    
    # plasma boundary
    ψb = eqt.global_quantities.psi_boundary
    ψ0 = eqt.global_quantities.psi_axis
    pr, pz = IMAS.flux_surface(eqt, (ψb - ψ0) * 0.999 + ψ0)
    
    R0 = eqt.boundary.geometric_axis.r
    
    # coils will not be closer to the boundary than this
    ml = 1.2
    pr_inner = (pr .- R0) .* ml .+ R0
    pz_inner = pz .* ml * 1.1

    # initial position of coils
    ml = ml + 0.1
    pr_coil = (pr .- R0) .* ml .+ R0
    pz_coil = pz .* ml * 1.1

    # forbidden plasma area
    for ((kr, kz), (rr, zz)) in hcat(pts...)
        point = SA[rr,zz]
        if PolygonOps.inpolygon((rr, zz), StaticArrays.SVector.(pr_inner, pz_inner)) == 1
            mask[kz,kr] = 1
        end
    end
    
    # apply filtering to smooth transition between allowed and forbidden regions
    n = 7
    filter = exp.(-(range(-1, 1, length=2 * n + 1) / n).^2)
    filter = filter ./ sum(filter)
    mask = DSP.conv(filter, filter, mask)[n + 1:end - n,n + 1:end - n]

    # Cubic spline interpolation on the log to ensure positivity of the cost
    mask_log_interpolant_raw = Interpolations.CubicSplineInterpolation((z, r), log10.(1.0 .+ mask))
    mask_log_interpolant_raw = Interpolations.extrapolate(mask_log_interpolant_raw.itp, Interpolations.Flat());
    function mask_log_interpolant(z, r)
        return 10.0.^(mask_log_interpolant_raw(z, r)) .- 1
    end
    
    # initial coils placement is uniformely spaced along a surface conforming to plasma
    coils = []
    θ_coil = atan.(pz_coil, pr_coil .- R0)
    θ_coil = DSP.unwrap(θ_coil)
    if θ_coil[2] < θ_coil[1]
        θ_coil = -θ_coil
    end
    pr_coilθ = Interpolations.extrapolate(Interpolations.interpolate((θ_coil,), pr_coil, Interpolations.Gridded(Interpolations.Linear())), Interpolations.Periodic())
    pz_coilθ = Interpolations.extrapolate(Interpolations.interpolate((θ_coil,), pz_coil, Interpolations.Gridded(Interpolations.Linear())), Interpolations.Periodic())
    θ = range(0, 2π, length=1001)[1:end - 1]
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
    coils = [[pr_coilL(L),pz_coilL(L)] for L in lcoils .+ L0]
    if mod(ncoils, 2) == 1
        coils[1][2] = 0.0
    end
    
    # populate IMAS data structure
    pf_active = IMAS.pf_active()
    resize!(pf_active.coil, length(coils))
    for (k, c) in enumerate(coils)
        resize!(pf_active.coil[k].element, 1)
        pf_active.coil[k].element[1].geometry.rectangle.r = c[1]
        pf_active.coil[k].element[1].geometry.rectangle.z = c[2]
        pf_active.coil[k].element[1].geometry.rectangle.height = 0.0
        pf_active.coil[k].element[1].geometry.rectangle.width = 0.0
    end
    
    PFcoilsActor(equilibrium, pf_active, mask_log_interpolant)
end
