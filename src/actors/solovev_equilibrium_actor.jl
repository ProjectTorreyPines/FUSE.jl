
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
    qstar = 1.5
    Ip_dir = Int(sign(qstar) * B0_dir)
    alpha = 0.0
    S0 = solovev(B0, R0, ϵ, δ, κ, alpha, qstar, B0_dir=B0_dir, Ip_dir=Ip_dir)
    SolovevEquilibriumActor(equilibrium, time, S0, IMAS.equilibrium())
end

#= == =#
# STEP #
#= == =#
function Base.step(actor::SolovevEquilibriumActor; verbose=false)
    # non-linear optimization to obtain a target `beta_normal`
    S0 = actor.S
    time_index = get_time_index(actor.eq_in.time_slice, actor.time)

    function opti(x)
        S = solovev(S0.B0, S0.R0, S0.epsilon, S0.delta, S0.kappa, x[1], x[2], B0_dir=S0.sigma_B0, Ip_dir=S0.sigma_Ip)
        beta_cost=(Equilibrium.beta_n(S) - actor.eq_in.time_slice[time_index].global_quantities.beta_normal).^2
        ip_cost=((Equilibrium.plasma_current(S)- actor.eq_in.time_slice[time_index].global_quantities.ip)/1E6).^2
        return beta_cost+ip_cost
    end

    res = Optim.optimize(opti, [S0.alpha,S0.qstar], Optim.Newton(); autodiff=:forward)
    
    if verbose
        println(res)
    end

    return actor.S = solovev(S0.B0, S0.R0, S0.epsilon, S0.delta, S0.kappa, res.minimizer[1], res.minimizer[2], B0_dir=S0.sigma_B0, Ip_dir=S0.sigma_Ip)
end

#= ====== =#
# FINALIZE #
#= ====== =#
function finalize(actor::SolovevEquilibriumActor, n::Integer=129)
    @assert mod(n, 2) == 1 "`n` in finalize SolovevEquilibriumActor must be a odd number"
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

    # magnetic axis
    eqt.global_quantities.magnetic_axis.r, eqt.global_quantities.magnetic_axis.z = Equilibrium.magnetic_axis(actor.S)

    # generate grid with vertex on magnetic axis
    resize!(eqt.profiles_2d, 1)
    eqt.profiles_2d[1].grid_type.index = 1
    rlims, zlims = Equilibrium.limits(actor.S)
    dr = (eqt.global_quantities.magnetic_axis.r - (rlims[2] + rlims[1]) / 2.0)
    dr0 = (rlims[2] - rlims[1]) / n
    ddr = mod(dr, dr0)
    eqt.profiles_2d[1].grid.dim1 = range(rlims[1] - ddr, rlims[2] + ddr, length=n) .+ ddr
    _, i = findmin(abs.(eqt.profiles_2d[1].grid.dim1 .- eqt.global_quantities.magnetic_axis.r))
    eqt.profiles_2d[1].grid.dim1 = eqt.profiles_2d[1].grid.dim1 .- eqt.profiles_2d[1].grid.dim1[i] .+ eqt.global_quantities.magnetic_axis.r
    dz = (eqt.global_quantities.magnetic_axis.z - (zlims[2] + zlims[1]) / 2.0)
    dz0 = (zlims[2] - zlims[1]) / n
    ddz = mod(dz, dz0)
    eqt.profiles_2d[1].grid.dim2 = range(zlims[1] - ddz, zlims[2] + ddz, length=n) .+ ddz
    _, i = findmin(abs.(eqt.profiles_2d[1].grid.dim2 .- eqt.global_quantities.magnetic_axis.z))
    eqt.profiles_2d[1].grid.dim2 = eqt.profiles_2d[1].grid.dim2 .- eqt.profiles_2d[1].grid.dim2[i] .+ eqt.global_quantities.magnetic_axis.z

    eqt.profiles_2d[1].psi = [actor.S(rr, zz) for rr in eqt.profiles_2d[1].grid.dim1, zz in eqt.profiles_2d[1].grid.dim2]

    eqt.profiles_2d[1].b_field_r = zeros(size(eqt.profiles_2d[1].psi)...)
    eqt.profiles_2d[1].b_field_tor = zeros(size(eqt.profiles_2d[1].psi)...)
    eqt.profiles_2d[1].b_field_z = zeros(size(eqt.profiles_2d[1].psi)...)
    for (kr, rr) in enumerate(eqt.profiles_2d[1].grid.dim1), (kz, zz) in enumerate(eqt.profiles_2d[1].grid.dim2)
        (eqt.profiles_2d[1].b_field_r[kr,kz], eqt.profiles_2d[1].b_field_tor[kr,kz], eqt.profiles_2d[1].b_field_z[kr,kz]) = Bfield(actor.S, rr, zz)
    end



    IMAS.flux_surfaces(eqt, actor.S.B0, actor.S.R0)

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