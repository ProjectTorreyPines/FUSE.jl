
using Equilibrium
using Contour
using PolygonOps
using StaticArrays
using Optim
using Interpolations

mutable struct SolovevEquilibriumActor <: EquilibriumActor
    eq_in::IMAS.equilibrium
    time::Real
    S::SolovevEquilibrium
    eq_out::IMAS.equilibrium
end

#= == =#
# INIT #
#= == =#
function SolovevEquilibriumActor(equilibrium::IMAS.equilibrium, time::Real)
    time_index = get_time_index(equilibrium.time_slice, time)
    eqt = equilibrium.time_slice[time_index]
    a = eqt.profiles_1d.r_outboard[end] - eqt.profiles_1d.r_inboard[end]
    R0 = eqt.profiles_1d.r_outboard[end] + eqt.profiles_1d.r_inboard[end]
    κ = eqt.profiles_1d.elongation[end]
    δ = (eqt.profiles_1d.triangularity_lower[end] + eqt.profiles_1d.triangularity_upper[end]) / 2.0
    ϵ = a / R0
    B0 = abs(equilibrium.vacuum_toroidal_field.b0[time_index])
    B0_dir = Int(sign(B0))
    R0 = equilibrium.vacuum_toroidal_field.r0
    qstar = eqt.profiles_1d.q[end]
    Ip_dir = Int(sign(qstar) * B0_dir)
    alpha = 0.0
    S0 = solovev(B0, R0, ϵ, δ, κ, alpha, qstar, B0_dir=B0_dir, Ip_dir=Ip_dir)
    SolovevEquilibriumActor(equilibrium, time, S0, IMAS.equilibrium())
end

#= == =#
# STEP #
#= == =#
function Base.step(actor::SolovevEquilibriumActor; abs_error=1E-3, max_iter=100, verbose=false)
    # non-linear optimization to obtain a target `beta_t`
    S0 = actor.S
    time_index = get_time_index(actor.eq_in.time_slice, actor.time)
    target_beta = actor.eq_in.time_slice[time_index].global_quantities.beta_tor

    function opti(x)
        S = solovev(S0.B0, S0.R0, S0.epsilon, S0.delta, S0.kappa, x[1], S0.qstar, B0_dir=S0.sigma_B0, Ip_dir=S0.sigma_Ip)
        v[1] = S
        precision = abs(S.beta_t - target_beta)
        cost = precision.^2
        return cost
    end
    
    x = [S0.alpha]
    v = Any[S0]
    for k in 1:max_iter
        g = ForwardDiff.gradient(opti, x, )
        x[1] -= (g[1] / abs(S0.beta_t))
        precision = abs(v[1].beta_t.value - target_beta)
        if verbose
            @printf("α=%3.3f β_t=%3.3e precision=%3.3e\n", x[1], v[1].beta_t.value, precision)
        end
        if precision < abs_error
            break
        end
        if k == max_iter
            error("Current β_t=$(v[1].beta_t.value) is not β_t,target $(target_beta)")
        end
    end
    
    actor.S = solovev(S0.B0, S0.R0, S0.epsilon, S0.delta, S0.kappa, x[1], S0.qstar, B0_dir=S0.sigma_B0, Ip_dir=S0.sigma_Ip)
end

#= ====== =#
# FINALIZE #
#= ====== =#
function finalize(actor::SolovevEquilibriumActor, n=129)
    equilibrium = actor.eq_out
    time_index = get_time_index(equilibrium.time_slice, actor.time)
    eqt = equilibrium.time_slice[time_index]

    eqt.profiles_1d.psi = collect(range(Equilibrium.psi_limits(actor.S)..., length=n))

    set_field_time_array(equilibrium.vacuum_toroidal_field, :b0, time_index, actor.S.B0)
    equilibrium.vacuum_toroidal_field.r0 = actor.S.R0

    eqt.profiles_1d.pressure = Equilibrium.pressure(actor.S, eqt.profiles_1d.psi)
    eqt.profiles_1d.dpressure_dpsi = Equilibrium.pressure_gradient(actor.S, eqt.profiles_1d.psi)

    eqt.profiles_1d.f = Equilibrium.poloidal_current(actor.S, eqt.profiles_1d.psi)
    eqt.profiles_1d.f_df_dpsi = eqt.profiles_1d.f .* Equilibrium.poloidal_current_gradient(actor.S, eqt.profiles_1d.psi)

    resize!(eqt.profiles_2d, 1)
    eqt.profiles_2d[1].grid_type.index = 1
    rlims, zlims = Equilibrium.limits(actor.S)
    eqt.profiles_2d[1].grid.dim1 = range(rlims..., length=n)
    eqt.profiles_2d[1].grid.dim2 = range(zlims..., length=n)
    eqt.profiles_2d[1].psi = [actor.S(rr, zz) for rr in eqt.profiles_2d[1].grid.dim1, zz in eqt.profiles_2d[1].grid.dim2]

    eqt.global_quantities.magnetic_axis.r, eqt.global_quantities.magnetic_axis.z = Equilibrium.magnetic_axis(actor.S)

    eqt.profiles_1d.elongation = zero(eqt.profiles_1d.psi)
    eqt.profiles_1d.triangularity_lower = zero(eqt.profiles_1d.psi)
    eqt.profiles_1d.triangularity_upper = zero(eqt.profiles_1d.psi)
    eqt.profiles_1d.r_inboard = zero(eqt.profiles_1d.psi)
    eqt.profiles_1d.r_outboard = zero(eqt.profiles_1d.psi)
    for (k, psi_level) in enumerate(eqt.profiles_1d.psi)
        pr, pz = flux_surface(eqt, psi_level)
        tmp = fluxGeo(pr, pz)
        eqt.profiles_1d.elongation[k] = tmp["kappa"]
        eqt.profiles_1d.triangularity_upper[k] = tmp["delu"]
        eqt.profiles_1d.triangularity_lower[k] = tmp["dell"]
        eqt.profiles_1d.r_outboard[k] = tmp["max_r"]
        eqt.profiles_1d.r_inboard[k] = tmp["min_r"]
        eqt.boundary.elongation_upper = tmp["kapu"]
        eqt.boundary.elongation_lower = tmp["kapl"]
    end
    eqt.profiles_1d.triangularity_upper[1] = 0.0
    eqt.profiles_1d.triangularity_lower[1] = 0.0
    eqt.profiles_1d.elongation[1]= eqt.profiles_1d.elongation[2]

    # eqt.profiles_1d.r_outboard =
    # eqt.profiles_1d.r_inboard = 
    # eqt.profiles_1d.elongation = 
    # eqt.profiles_1d.triangularity_upper = 
    # eqt.profiles_1d.triangularity_lower = 
    # equilibrium.vacuum_toroidal_field.b0 = 
    # equilibrium.vacuum_toroidal_field.r0 = 
    # ...
    return actor.eq_out
end


#= ====================== =#
# IMAS PROESSING FUNCTIONS #
#= ====================== =#
"""
    flux_surface(eqt, psi_level)

returns r,z coordiates of flux surface at given psi_level
"""
function flux_surface(eqt::IMAS.equilibrium__time_slice, psi_level::Real)
    if psi_level == eqt.profiles_1d.psi[1]
        psi_level = eqt.profiles_1d.psi[1] * 0.9 + eqt.profiles_1d.psi[2] * 0.1
    end
    cl = Contour.contour(eqt.profiles_2d[1].grid.dim1,
                         eqt.profiles_2d[1].grid.dim2,
                         eqt.profiles_2d[1].psi,
                         psi_level)
    for line in Contour.lines(cl)
        pr, pz = Contour.coordinates(line)
        if !((pr[1] == pr[end]) & (pz[1] == pz[end]))
            continue
        end
        polygon = SVector.(pr, pz)
        if PolygonOps.inpolygon((eqt.global_quantities.magnetic_axis.r, eqt.global_quantities.magnetic_axis.z), polygon) == 1
            return pr, pz
        end
    end
end

function flux_surfaces_geo(eqt::IMAS.equilibrium__time_slice)
end

"""
    fluxGeo(inputR::Vector{Real}, inputZ::Vector{Real})::Dict

Recturns dictionary with geometric properties of a given flux surface
"""
function fluxGeo(inputR::Vector{T} where T <: Real, inputZ::Vector{T} where T <: Real)::Dict

    # concatenate inputs 3 times to avoid bound errors in minimization
    if inputR[1] == inputR[2]
        inputR = inputR[1:end - 1]
        inputR = inputZ[1:end - 1]
    end
    inputRclose = vcat(inputR, inputR, inputR)
    inputZclose = vcat(inputZ, inputZ, inputZ)

    # These are the extrema indices
    _, imaxr = findmax(inputR)
    _, iminr = findmin(inputR)
    _, imaxz = findmax(inputZ)
    _, iminz = findmin(inputZ)

    t = 1:length(inputRclose)
    
    interpR = Interpolations.CubicSplineInterpolation(t, inputRclose)
    optR = Optim.optimize(x -> -interpR(x), length(inputR) + imaxr - 1, length(inputR) + imaxr + 1)
    tmaxr = Optim.minimizer(optR)
    optR = Optim.optimize(x -> interpR(x), length(inputR) + iminr - 1, length(inputR) + iminr + 1)
    tminr = Optim.minimizer(optR)

    interpZ = Interpolations.CubicSplineInterpolation(t, inputZclose)
    optZ = Optim.optimize(x -> -interpZ(x), length(inputZ) + imaxz - 1, length(inputZ) + imaxz + 1)
    tmaxz = Optim.minimizer(optZ)
    optZ = Optim.optimize(x -> interpZ(x), length(inputZ) + iminz - 1, length(inputZ) + iminz + 1)
    tminz = Optim.minimizer(optZ)

    r_at_max_z, max_z = interpR(tmaxz), interpZ(tmaxz)
    r_at_min_z, min_z = interpR(tminz), interpZ(tminz)
    z_at_max_r, max_r = interpZ(tmaxr), interpR(tmaxr)
    z_at_min_r, min_r = interpZ(tminr), interpR(tminr)

    dl = vcat([0], sqrt.(diff(inputR).^2 + diff(inputZ).^2))

    geo = Dict()
    geo["r_at_max_z"] = r_at_max_z
    geo["r_at_min_z"] = r_at_min_z
    geo["z_at_max_r"] = z_at_max_r
    geo["z_at_min_r"] = z_at_min_r
    geo["max_z"] = max_z
    geo["min_z"] = min_z
    geo["max_r"] = max_r
    geo["min_r"] = min_r
    geo["R"] = 0.5 * (max_r + min_r)
    geo["Z"] = 0.5 * (max_z + min_z)
    geo["a"] = 0.5 * (max_r - min_r)
    geo["eps"] = geo["a"] / geo["R"]
    geo["per"] = sum(dl)
    geo["surfArea"] = 2 * pi * sum(inputR .* dl)
    geo["kappa"] = 0.5 * ((max_z - min_z) / geo["a"])
    geo["kapu"] = (max_z - z_at_max_r) / geo["a"]
    geo["kapl"] = (z_at_max_r - min_z) / geo["a"]
    geo["delu"] = (geo["R"] - r_at_max_z) / geo["a"]
    geo["dell"] = (geo["R"] - r_at_min_z) / geo["a"]
    geo["delta"] = 0.5 * (geo["dell"] + geo["delu"])
    geo["zoffset"] = z_at_max_r

    # NOTE: lonull, upnull, squareness, centroid, zeta have not been translated from OMFIT fluxGeo

    return geo
end