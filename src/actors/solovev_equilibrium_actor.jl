
using Equilibrium
using Printf
import ForwardDiff

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

    eqt.profiles_2d[1].b_field_r = zeros(size(eqt.profiles_2d[1].psi)...)
    eqt.profiles_2d[1].b_field_tor = zeros(size(eqt.profiles_2d[1].psi)...)
    eqt.profiles_2d[1].b_field_z = zeros(size(eqt.profiles_2d[1].psi)...)
    for (kr, rr) in enumerate(eqt.profiles_2d[1].grid.dim1), (kz, zz) in enumerate(eqt.profiles_2d[1].grid.dim2)
        (eqt.profiles_2d[1].b_field_r[kr,kz], eqt.profiles_2d[1].b_field_tor[kr,kz], eqt.profiles_2d[1].b_field_z[kr,kz]) = Bfield(actor.S, rr, zz)
    end

    # magnetic axis
    eqt.global_quantities.magnetic_axis.r, eqt.global_quantities.magnetic_axis.z = Equilibrium.magnetic_axis(actor.S)

    IMAS.flux_surfaces(eqt)

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