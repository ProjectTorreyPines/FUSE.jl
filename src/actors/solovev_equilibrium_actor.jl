
using Equilibrium
using Printf
import ForwardDiff
import Optim

mutable struct SolovevEquilibriumActor <: EquilibriumActor
    eq_in::IMAS.equilibrium
    time::Real
    S::SolovevEquilibrium
    eq_out::IMAS.equilibrium
end

function IMAS2Equilibrium(equilibrium::IMAS.equilibrium,
                          time::Real)
    time_index = get_time_index(equilibrium.time_slice, time)
    eqt = equilibrium.time_slice[time_index]

    dim1=range(eqt.profiles_2d[1].grid.dim1[1],eqt.profiles_2d[1].grid.dim1[end],length=length(eqt.profiles_2d[1].grid.dim1))
    @assert collect(dim1) ≈ eqt.profiles_2d[1].grid.dim1
    dim2=range(eqt.profiles_2d[1].grid.dim2[1],eqt.profiles_2d[1].grid.dim2[end],length=length(eqt.profiles_2d[1].grid.dim2))
    @assert collect(dim2) ≈ eqt.profiles_2d[1].grid.dim2
    psi=range(eqt.profiles_1d.psi[1],eqt.profiles_1d.psi[end],length=length(eqt.profiles_1d.psi))
    @assert collect(psi) ≈ eqt.profiles_1d.psi

    Equilibrium.efit(Equilibrium.cocos(11), # COCOS
                     dim1, # Radius/R range
                     dim2, # Elevation/Z range
                     psi, # Polodial Flux range (polodial flux from magnetic axis)
                     eqt.profiles_2d[1].psi, # Polodial Flux on RZ grid (polodial flux from magnetic axis)
                     eqt.profiles_1d.f, # Polodial Current
                     eqt.profiles_1d.pressure, # Plasma pressure
                     eqt.profiles_1d.q, # Q profile
                     eqt.profiles_1d.psi.*0, # Electric Potential
                     (eqt.global_quantities.magnetic_axis.r, eqt.global_quantities.magnetic_axis.z), # Magnetic Axis (raxis,zaxis)
                     Int(sign(equilibrium.vacuum_toroidal_field.b0[time_index])*sign(eqt.global_quantities.ip)) # sign(dot(J,B))
                     )
end

#= == =#
# INIT #
#= == =#
"""
    function SolovevEquilibriumActor(equilibrium::IMAS.equilibrium, time::Real)

Constructor for the SolovevEquilibriumActor structure
“One size fits all” analytic solutions to the Grad–Shafranov equation
Phys. Plasmas 17, 032502 (2010); https://doi.org/10.1063/1.3328818

- qstar: Kink safety factor

- alpha: Constant affecting the pressure
"""
function SolovevEquilibriumActor(equilibrium::IMAS.equilibrium,
                                 time::Real;
                                 qstar = 1.5,
                                 alpha = 0.0)
    time_index = get_time_index(equilibrium.time_slice, time)
    eqt = equilibrium.time_slice[time_index]

    a = eqt.boundary.minor_radius
    R0 = eqt.boundary.geometric_axis.r
    κ = eqt.boundary.elongation
    δ = eqt.boundary.triangularity
    ϵ = a / R0
    B0 = abs(equilibrium.vacuum_toroidal_field.b0[time_index])
    B0_dir = Int(sign(B0))
    R0 = equilibrium.vacuum_toroidal_field.r0
    Ip_dir = Int(sign(qstar) * B0_dir)

    if length(eqt.boundary.x_point)>0
        xpoint = (eqt.boundary.x_point[1].r, eqt.boundary.x_point[1].z)
    else
        xpoint = nothing
    end

    S0 = solovev(B0, R0, ϵ, δ, κ, alpha, qstar, B0_dir=B0_dir, Ip_dir=Ip_dir, xpoint=xpoint)

    SolovevEquilibriumActor(equilibrium, time, S0, IMAS.equilibrium())
end

#= == =#
# STEP #
#= == =#
"""
    Base.step(actor::SolovevEquilibriumActor; verbose=false)

Non-linear optimization to obtain a target `ip` and `beta_normal`
"""
function Base.step(actor::SolovevEquilibriumActor; verbose=false)
    S0 = actor.S
    time_index = get_time_index(actor.eq_in.time_slice, actor.time)

    target_ip = actor.eq_in.time_slice[time_index].global_quantities.ip
    target_beta = actor.eq_in.time_slice[time_index].global_quantities.beta_normal

    B0, R0, epsilon, delta, kappa, alpha, qstar, target_ip, target_beta = promote(S0.B0, S0.R0, S0.epsilon, S0.delta, S0.kappa, S0.alpha, S0.qstar, target_ip, target_beta)

    # NOTE: some problems when running with xpoint, that I suspect are due to issues with flux surface tracing in Equilibrium.jl
    function opti(x)
        S = solovev(B0, R0, epsilon, delta, kappa, x[1], x[2], B0_dir=S0.sigma_B0, Ip_dir=S0.sigma_Ip, xpoint=S0.xpoint)
        beta_cost = abs((Equilibrium.beta_n(S) - target_beta)^2/target_beta)
        ip_cost = abs((Equilibrium.plasma_current(S) - target_ip)^2/target_ip)
        return (beta_cost + ip_cost)^2
    end

    if isa(B0, ForwardDiff.Dual)
        optim_method = Optim.NelderMead()
    else
        optim_method = Optim.Newton()
    end
    res = Optim.optimize(opti, [alpha, qstar], optim_method, Optim.Options(g_tol=1E-1); autodiff=:forward)
    
    if verbose
        println(res)
    end

    actor.S = solovev(B0, R0, epsilon, delta, kappa, res.minimizer[1], res.minimizer[2], B0_dir=S0.sigma_B0, Ip_dir=S0.sigma_Ip, xpoint=S0.xpoint)
    return res
end

#= ====== =#
# FINALIZE #
#= ====== =#
"""
    finalize(actor::SolovevEquilibriumActor, n::Integer=129)::IMAS.equilibrium

Store SolovevEquilibriumActor data in IMAS.equilibrium
"""
function finalize(actor::SolovevEquilibriumActor,
                  resolution::Integer=129,
                  rlims::NTuple{2,<:Real} = Equilibrium.limits(actor.S)[1],
                  zlims::NTuple{2,<:Real} = Equilibrium.limits(actor.S)[2])::IMAS.equilibrium

    equilibrium = actor.eq_out
    time_index = get_time_index(equilibrium.time_slice, actor.time)
    eqt = equilibrium.time_slice[time_index]

    eqt.profiles_1d.psi = collect(range(Equilibrium.psi_limits(actor.S)..., length=resolution))

    set_field_time_array(equilibrium.vacuum_toroidal_field, :b0, time_index, actor.S.B0)
    equilibrium.vacuum_toroidal_field.r0 = actor.S.R0

    eqt.profiles_1d.pressure = Equilibrium.pressure(actor.S, eqt.profiles_1d.psi)
    eqt.profiles_1d.dpressure_dpsi = Equilibrium.pressure_gradient(actor.S, eqt.profiles_1d.psi)

    eqt.profiles_1d.f = Equilibrium.poloidal_current(actor.S, eqt.profiles_1d.psi)
    eqt.profiles_1d.f_df_dpsi = eqt.profiles_1d.f .* Equilibrium.poloidal_current_gradient(actor.S, eqt.profiles_1d.psi)

    # magnetic axis
    eqt.global_quantities.magnetic_axis.r, eqt.global_quantities.magnetic_axis.z = Equilibrium.magnetic_axis(actor.S)

    # generate grid with vertex on magnetic axis
    resize!(eqt.profiles_2d, 1)
    eqt.profiles_2d[1].grid_type.index = 1
    if true
        eqt.profiles_2d[1].grid.dim1 = range(rlims[1], rlims[2], length=resolution)
        eqt.profiles_2d[1].grid.dim2 = range(zlims[1], zlims[2], length=resolution)#Int(ceil(resolution*actor.S.kappa)))
    else
        dr = (eqt.global_quantities.magnetic_axis.r - (rlims[2] + rlims[1]) / 2.0)
        dr0 = (rlims[2] - rlims[1]) / resolution
        ddr = mod(dr, dr0)
        eqt.profiles_2d[1].grid.dim1 = range(rlims[1] - ddr, rlims[2] + ddr, length=resolution) .+ ddr
        _, i = findmin(abs.(eqt.profiles_2d[1].grid.dim1 .- eqt.global_quantities.magnetic_axis.r))
        eqt.profiles_2d[1].grid.dim1 = eqt.profiles_2d[1].grid.dim1 .- eqt.profiles_2d[1].grid.dim1[i] .+ eqt.global_quantities.magnetic_axis.r
        dz = (eqt.global_quantities.magnetic_axis.z - (zlims[2] + zlims[1]) / 2.0)
        dz0 = (zlims[2] - zlims[1]) / resolution
        ddz = mod(dz, dz0)
        eqt.profiles_2d[1].grid.dim2 = range(zlims[1] - ddz, zlims[2] + ddz, length=resolution) .+ ddz
        _, i = findmin(abs.(eqt.profiles_2d[1].grid.dim2 .- eqt.global_quantities.magnetic_axis.z))
        eqt.profiles_2d[1].grid.dim2 = eqt.profiles_2d[1].grid.dim2 .- eqt.profiles_2d[1].grid.dim2[i] .+ eqt.global_quantities.magnetic_axis.z
    end

    eqt.profiles_2d[1].psi = [actor.S(rr, zz) for rr in eqt.profiles_2d[1].grid.dim1, zz in eqt.profiles_2d[1].grid.dim2]

    eqt.profiles_2d[1].b_field_r = zeros(size(eqt.profiles_2d[1].psi)...)
    eqt.profiles_2d[1].b_field_tor = zeros(size(eqt.profiles_2d[1].psi)...)
    eqt.profiles_2d[1].b_field_z = zeros(size(eqt.profiles_2d[1].psi)...)
    for (kr, rr) in enumerate(eqt.profiles_2d[1].grid.dim1), (kz, zz) in enumerate(eqt.profiles_2d[1].grid.dim2)
        (eqt.profiles_2d[1].b_field_r[kr,kz], eqt.profiles_2d[1].b_field_tor[kr,kz], eqt.profiles_2d[1].b_field_z[kr,kz]) = Bfield(actor.S, rr, zz)
    end

    IMAS.flux_surfaces(eqt, actor.S.B0, actor.S.R0)

    return actor.eq_out
end