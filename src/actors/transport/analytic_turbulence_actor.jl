import GACODE
import IMAS

#= ======================= =#
#  ActorAnalyticTurbulence  #
#= ======================= =#
@actor_parameters_struct ActorAnalyticTurbulence{T} begin
    model::Switch{Symbol} = Switch{Symbol}([:GyroBohm, :BgB], "-", "Analytic transport model"; default=:GyroBohm)
    αBgB::Entry{T} = Entry{T}("-", "Scale factor for BgB transport"; default=0.01)
    χeB_coefficient::Entry{T} = Entry{T}("-", "Coefficient of Bohm component in χe. χe=αBgB*(χeB_coefficient*χeB+χeGB_coefficient*χeGB)"; default=0.01)
    χeGB_coefficient::Entry{T} = Entry{T}("-", "Coefficient of gyro-Bohm component in χe. χe=αBgB*(χeB_coefficient*χeB+χeGB_coefficient*χeGB)"; default=50.0)
    χiB_coefficient::Entry{T} = Entry{T}("-", "Coefficient of Bohm component in χi. χi=αBgB*(χiB_coefficient*χiB+χiGB_coefficient*χiGB)"; default=0.001)
    χiGB_coefficient::Entry{T} = Entry{T}("-", "Coefficient of gyro-Bohm component in χi. χi=αBgB*(χiB_coefficient*χiB+χiGB_coefficient*χiGB)"; default=1.0)
    rho_transport::Entry{AbstractVector{T}} = Entry{AbstractVector{T}}("-", "rho_tor_norm values to compute fluxes on"; default=0.25:0.1:0.85)
    do_plot::Entry{Bool} = act_common_parameters(; do_plot=false)
end

mutable struct ActorAnalyticTurbulence{D,P} <: SingleAbstractActor{D,P}
    dd::IMAS.dd{D}
    par::OverrideParameters{P,FUSEparameters__ActorAnalyticTurbulence{P}}
    flux_solutions::Vector{GACODE.FluxSolution{D}}
end

"""
    ActorAnalyticTurbulence(dd::IMAS.dd, act::ParametersAllActors; kw...)

Evaluates analytic turbulent transport models including GyroBohm and Bohm+gyro-Bohm (BgB) models.

The actor supports two transport models:

  - `:GyroBohm`: Simple gyro-Bohm scaling returning unit fluxes
  - `:BgB`: Detailed Bohm + gyro-Bohm model calculating electron and ion energy diffusivities (χe, χi)
    and particle flux (Γe) based on local plasma parameters, pressure gradients, and magnetic geometry
    Tholerus, Emmi, et al. Nuclear Fusion 64.10 (2024): 106030

The BgB model computes transport coefficients using local temperature and density gradients,
safety factor profiles, and magnetic field geometry. Results are normalized to gyro-Bohm
and stored in `dd.core_transport` as anomalous transport fluxes.
"""
function ActorAnalyticTurbulence(dd::IMAS.dd, act::ParametersAllActors; kw...)
    actor = ActorAnalyticTurbulence(dd, act.ActorAnalyticTurbulence; kw...)
    step(actor)
    finalize(actor)
    return actor
end

function ActorAnalyticTurbulence(dd::IMAS.dd{D}, par::FUSEparameters__ActorAnalyticTurbulence{P}; kw...) where {D<:Real,P<:Real}
    logging_actor_init(ActorAnalyticTurbulence)
    par = OverrideParameters(par; kw...)
    return ActorAnalyticTurbulence(dd, par, GACODE.FluxSolution{D}[])
end

"""
    _step(actor::ActorAnalyticTurbulence)

Runs the selected analytic turbulent transport model on specified radial grid points.

For the BgB model, calculates:

  - χeB, χeGB: Bohm and gyro-Bohm electron heat diffusivities
  - χiB, χiGB: Bohm and gyro-Bohm ion heat diffusivities
  - Γe: Electron particle flux based on density gradients

All transport coefficients are scaled by user-specified coefficients and stored as
FluxSolution objects normalized to gyro-Bohm units.
"""
function _step(actor::ActorAnalyticTurbulence{D,P}) where {D<:Real,P<:Real}
    dd = actor.dd
    par = actor.par

    cp1d = dd.core_profiles.profiles_1d[]
    eqt = dd.equilibrium.time_slice[]
    eqt1d = eqt.profiles_1d

    if actor.par.model == :GyroBohm
        flux_solution = GACODE.FluxSolution{D}(1.0, 1.0, 1.0, [1.0 for ion in cp1d.ion], 1.0)
        actor.flux_solutions = GACODE.FluxSolution{D}[flux_solution for irho in par.rho_transport]

    elseif actor.par.model == :BgB
        A1 = 1.0
        A2 = 0.3
        αgB = 5e-6
        αB = 2e-4

        rho_cp = cp1d.grid.rho_tor_norm
        rho_eq = eqt1d.rho_tor_norm

        bax = eqt.global_quantities.magnetic_axis.b_field_tor
        q_profile = IMAS.interp1d(rho_eq, eqt1d.q).(rho_cp)
        volume = IMAS.interp1d(rho_eq, eqt1d.volume).(rho_cp)
        dvoldrho = IMAS.interp1d(rho_eq, eqt1d.dvolume_drho_tor).(rho_cp)*eqt1d.rho_tor[end]
        vprime_miller = IMAS.interp1d(rho_eq, GACODE.volume_prime_miller_correction(eqt)).(rho_cp)
        surf = IMAS.interp1d(rho_eq, eqt1d.surface).(rho_cp)

        rmin = GACODE.r_min_core_profiles(eqt1d, rho_cp)
        Te = cp1d.electrons.temperature
        dlntedr = .-IMAS.calc_z(rho_cp, Te, :backward)
        ne = cp1d.electrons.density_thermal
        dlnnedr = .-IMAS.calc_z(rho_cp, ne, :backward)

        Ti = cp1d.ion[1].temperature
        zeff = cp1d.zeff
        dlntidr = .-IMAS.calc_z(rho_cp, Ti, :backward)

        Q_GB = GACODE.gyrobohm_energy_flux(cp1d, eqt)
        Γ_GB = GACODE.gyrobohm_particle_flux(cp1d, eqt)

        dpe = @. ne * IMAS.mks.e * Te * (dlntedr .+ dlnnedr)
        dpi = @. ne * IMAS.mks.e * Ti * (dlntidr .+ dlnnedr) / zeff
        χeB = @. αB * rmin[end] * q_profile^2 / abs(bax) * Te * (dlntedr .+ dlnnedr)
        χeGB = @. αgB * sqrt(Te / bax^2) * Te * dlntedr

        χiB = 2.0 * χeB
        χiGB = 0.5 * χeGB

        χe = @. actor.par.αBgB * (actor.par.χeB_coefficient * χeB + actor.par.χeGB_coefficient * χeGB)
        χi = @. actor.par.αBgB * (actor.par.χiB_coefficient * χiB + actor.par.χiGB_coefficient * χiGB)

        Dp = @. (A1 + (A2 - A1) * rho_cp) * χe * χi / (χe + χi)
        vin = @. 0.5 * Dp * (surf ^ 2 / volume) / dvoldrho
        vin[1] = 0.0
        Γe = @. (Dp * dlnnedr - vin) * ne

        gridpoint_cp = [argmin_abs(cp1d.grid.rho_tor_norm, ρ) for ρ in par.rho_transport]

        actor.flux_solutions = [
            GACODE.FluxSolution{D}(
                χe[irho] * dpe[irho] / Q_GB[irho] / vprime_miller[irho],
                χi[irho] * dpi[irho] / Q_GB[irho] / vprime_miller[irho],
                Γe[irho] / Γ_GB[irho] / vprime_miller[irho],
                [Γe[irho] / Γ_GB[irho] / vprime_miller[irho] / ion.element[1].z_n for ion in cp1d.ion],
                1.0
            )
            for irho in gridpoint_cp
        ]

        if par.do_plot
            plot(rho_cp, χeB; ylabel="χ", xlabel="ρ", label="χeB")
            display(plot!(rho_cp, 5e3 * χeGB; ylabel="χ", xlabel="ρ", label="5000 χeGB"))
        end
    end

    return actor
end

"""
    _finalize(actor::ActorAnalyticTurbulence)

Writes calculated transport fluxes to `dd.core_transport.model[:anomalous]`.

The model identifier is set to the selected transport model name (:GyroBohm or :BgB)
and flux results are converted from normalized GACODE format to IMAS format using
`GACODE.flux_gacode_to_imas` for electron/ion energy, particle, and momentum fluxes.
"""
function _finalize(actor::ActorAnalyticTurbulence)
    dd = actor.dd
    par = actor.par
    cp1d = dd.core_profiles.profiles_1d[]
    eqt = dd.equilibrium.time_slice[]

    model = resize!(dd.core_transport.model, :anomalous; wipe=false)
    model.identifier.name = string(par.model)
    m1d = resize!(model.profiles_1d)
    m1d.grid_flux.rho_tor_norm = par.rho_transport

    GACODE.flux_gacode_to_imas((:electron_energy_flux, :ion_energy_flux, :electron_particle_flux, :ion_particle_flux, :momentum_flux), actor.flux_solutions, m1d, eqt, cp1d)

    return actor
end
