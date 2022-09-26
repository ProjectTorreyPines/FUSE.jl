import TGLFNN: flux_solution
#= ===================== =#
#  ActorNeoclassical      #
#= ===================== =#
mutable struct ActorNeoclassical <: PlasmaAbstractActor
    dd::IMAS.dd
    par::ParametersActor
    flux_solutions::AbstractVector{<:flux_solution}
end

function ParametersActor(::Type{Val{:ActorNeoclassical}})
    par = ParametersActor(nothing)
    par.neoclassical_model = Switch([:changhinton], "", "Neoclassical model to run"; default=:changhinton)
    par.rho_transport = Entry(AbstractVector{<:Real}, "", "rho_tor_norm values to compute neoclassical fluxes on"; default=0.2:0.1:0.8)
    return par
end

"""
    ActorNeoclassical(dd::IMAS.dd, act::ParametersAllActors; kw...)

The ActorNeoclassical evaluates the neoclassical predicted turbulence at a set of rho_tor_norm grid points
"""
function ActorNeoclassical(dd::IMAS.dd, act::ParametersAllActors; kw...)
    par = act.ActorNeoclassical(kw...)
    actor = ActorNeoclassical(dd, par)
    step(actor)
    finalize(actor)
    return actor
end

function ActorNeoclassical(dd::IMAS.dd, par::ParametersActor; kw...)
    par = par(kw...)
    return ActorNeoclassical(dd, par, flux_solution[])
end

"""
    step(actor::ActorNeoclassical)

Runs Neoclassical actor to evaluate the turbulence flux on a Vector of gridpoints
"""
function step(actor::ActorNeoclassical)
    par = actor.par
    dd = actor.dd
    model = resize!(dd.core_transport.model, "identifier.index" => 5)
    model.identifier.name = string(par.neoclassical_model)
    m1d = resize!(model.profiles_1d)
    m1d.grid_flux.rho_tor_norm = par.rho_transport

    if par.neoclassical_model == :changhinton
        actor.flux_solutions = [neoclassical_changhinton(dd, rho, 1) for rho in par.rho_transport]
    end

    return actor
end

"""
    function finalize(actor::ActorNeoclassical)

Writes results to dd.core_transport
"""
function finalize(actor::ActorNeoclassical)
    dd = actor.dd
    cp1d = dd.core_profiles.profiles_1d[]
    eqt = dd.equilibrium.time_slice[]

    model = dd.core_transport.model[[idx for idx in keys(actor.dd.core_transport.model) if actor.dd.core_transport.model[idx].identifier.name == string(actor.par.neoclassical_model)][1]]
    model.profiles_1d[].total_ion_energy.flux = zeros(length(actor.par.rho_transport))
    for (neoclassical_idx, rho) in enumerate(actor.par.rho_transport)
        rho_transp_idx = findfirst(i -> i == rho, model.profiles_1d[].grid_flux.rho_tor_norm)
        rho_cp_idx = argmin(abs.(cp1d.grid.rho_tor_norm .- rho))
        model.profiles_1d[].total_ion_energy.flux[rho_transp_idx] = actor.flux_solutions[neoclassical_idx].ENERGY_FLUX_i * IMAS.gyrobohm_energy_flux(cp1d, eqt)[rho_cp_idx] # W / m^2
    end
    return actor
end

"""
    calc_neo(dd::IMAS.dd, rho_fluxmatch::Float64, iion::Integer, AS_ion::Float64)

This function calculates the neoclassical flux using Chang-Hinton model which has has been modified assuming Zi = 1, and ni=ne 
"""
function neoclassical_changhinton(dd::IMAS.dd, rho_fluxmatch::Real, iion::Integer)
    eq = dd.equilibrium
    eqt = eq.time_slice[]
    eq1d = eqt.profiles_1d
    prof1d = dd.core_profiles.profiles_1d[]
    rmin = 0.5 * (eq1d.r_outboard - eq1d.r_inboard)

    m_to_cm = 1e2
    Rmaj0 = eq1d.geometric_axis.r[1] * m_to_cm

    a = dd.equilibrium.time_slice[].boundary.minor_radius * m_to_cm
    rho_cp = prof1d.grid.rho_tor_norm
    rho_eq = eq1d.rho_tor_norm
    gridpoint_eq = argmin(abs.(rho_eq .- rho_fluxmatch))
    gridpoint_cp = argmin(abs.(rho_cp .- rho_fluxmatch))
    eps = rmin[gridpoint_eq] * m_to_cm / Rmaj0

    q = eq1d.q[gridpoint_eq]

    ne = prof1d.electrons.density_thermal * 1e-6
    Te = prof1d.electrons.temperature[gridpoint_cp]
    Ti = prof1d.ion[iion].temperature

    rmin = IMAS.interp1d(rho_eq, rmin).(rho_cp)
    Rmaj = IMAS.interp1d(eq1d.rho_tor_norm, eq1d.gm8).(prof1d.grid.rho_tor_norm)
    drmaj = IMAS.gradient(rmin, Rmaj)
    shift = -drmaj[gridpoint_cp]

    dlntdr = -IMAS.gradient(rmin, Ti)[gridpoint_cp] / Ti[gridpoint_cp]
    dlntdr = dlntdr * a / m_to_cm
    dlnndr = -IMAS.gradient(rmin, ne)[gridpoint_cp] / ne[gridpoint_cp]
    dlnndr = dlnndr * a / m_to_cm
    Ti = Ti[gridpoint_cp]
    ne = ne[gridpoint_cp]
    k = 1.6e-12
    e = 4.8e-10

    mi = 1.6726e-24 * prof1d.ion[iion].element[1].a

    c_s = sqrt(k * Te / mi)
    loglam = 24.0 - log(sqrt(ne / Te))

    k0 = 0.66
    a0 = 1.03
    b0 = 0.31
    c0 = 0.74

    Zi = prof1d.ion[iion].element[1].z_n
    alpha = Zi - 1
    nui = sqrt(2.0) * pi * ne * (Zi * e)^4 * loglam / sqrt(mi) / (k * Ti)^1.5
    nu = nui * a / c_s / sqrt(Ti / Ti) * (Ti / Te)^1.5
    nui_HH = nu * (4.0 / 3.0) / sqrt(2.0 * pi)
    nui_star_HH = nui_HH * (Rmaj0 / a) * abs(q / (sqrt(eps) * sqrt(eps) * sqrt(eps)))
    mu_star = (1.0 + 1.54 * alpha) * nui_star_HH

    CH_Bmag2inv_avg = ((1.0 + 1.5 * (eps * eps + eps * shift)
                        + 0.375 * eps * eps * eps * shift)
                       /
                       (1.0 + 0.5 * eps * shift))
    CH_Bmag2avg_inv = ((sqrt(1.0 - eps * eps) * (1.0 + 0.5 * eps * shift))
                       /
                       (1.0 + (shift / eps) * (sqrt(1.0 - eps * eps) - 1.0)))
    CH_I_div_psip = q / eps

    F2 = (0.5 / sqrt(eps)) * (CH_Bmag2inv_avg - CH_Bmag2avg_inv)

    K1 = (-alpha * (0.83 + 0.42 * alpha) / (0.58 + alpha)
          *
          (c0 * mu_star * sqrt(eps * eps * eps) * F2)
          /
          (1.0 + c0 * mu_star * sqrt(eps * eps * eps)))

    K2 = ((k0 * (1.0 + 1.54 * alpha)
           +
           (1.88 * sqrt(eps) - 1.54 * eps) * (1.0 + 3.75 * alpha)) * CH_Bmag2inv_avg
          /
          (1 + a0 * sqrt(mu_star) + b0 * mu_star)
          +
          k0 * sqrt(eps * eps * eps) * (c0 * c0 / b0) * mu_star * F2
          * (1.0 + 1.33 * alpha * (1.0 + 0.6 * alpha) / (1.0 + 1.79 * alpha))
          /
          (1.0 + c0 * sqrt(eps * eps * eps) * mu_star))

    neo_rho_star_in = 0.001
    #    @show Ti/Te
    efluxi = (CH_I_div_psip^2
              *
              Ti / Te
              * (prof1d.ion[iion].element[1].a / Zi^2
                 * neo_rho_star_in^2 * sqrt(eps) * nui_HH)
              * ((K2 + K1) * dlntdr + K1 * dlnndr))

    qneo_gb = neo_rho_star_in^2
    sol = flux_solution(0.0, 0.0, 0.0, efluxi / qneo_gb)
    return sol
end