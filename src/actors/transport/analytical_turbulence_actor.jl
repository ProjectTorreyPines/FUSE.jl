import GACODE
import IMAS

#= ======================= =#
#  ActorAnalyticTurbulence  #
#= ======================= =#
Base.@kwdef mutable struct FUSEparameters__ActorAnalyticTurbulence{T<:Real} <: ParametersActor{T}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    _time::Float64 = NaN
    model::Switch{Symbol} = Switch{Symbol}([:GyroBohm,:BgB], "-", "Analytic transport model"; default=:GyroBohm)
    αBgB::Entry{T} = Entry{T}("-", "Scale factor for BgB transport"; default=0.01)
    rho_transport::Entry{AbstractVector{T}} = Entry{AbstractVector{T}}("-", "rho_tor_norm values to compute fluxes on"; default=0.25:0.1:0.85)
    do_plot::Entry{Bool} = act_common_parameters(; do_plot=false)
end

mutable struct ActorAnalyticTurbulence{D,P} <: SingleAbstractActor{D,P}
    dd::IMAS.dd{D}
    par::OverrideParameters{P,FUSEparameters__ActorAnalyticTurbulence{P}}
    flux_solutions::Vector{<:GACODE.FluxSolution}
end

"""
    ActorAnalyticTurbulence(dd::IMAS.dd, act::ParametersAllActors; kw...)

Evaluates analytic turbulence models
"""
function ActorAnalyticTurbulence(dd::IMAS.dd, act::ParametersAllActors; kw...)
    actor = ActorAnalyticTurbulence(dd, act.ActorAnalyticTurbulence; kw...)
    step(actor)
    finalize(actor)
    return actor
end

function ActorAnalyticTurbulence(dd::IMAS.dd, par::FUSEparameters__ActorAnalyticTurbulence; kw...)
    logging_actor_init(ActorAnalyticTurbulence)
    par = OverrideParameters(par; kw...)
    return ActorAnalyticTurbulence(dd, par, GACODE.FluxSolution[])
end

"""
    _step(actor::ActorAnalyticTurbulence)

Runs analytic turbulent transport model on a vector of gridpoints
"""
function _step(actor::ActorAnalyticTurbulence)
    dd = actor.dd
    par = actor.par

    cp1d = dd.core_profiles.profiles_1d[]
    eqt = dd.equilibrium.time_slice[]
    eqt1d = eqt.profiles_1d
    if actor.par.model == :GyroBohm
        flux_solution = GACODE.FluxSolution(1.0, 1.0, 1.0, [1.0 for ion in cp1d.ion], 1.0)
        actor.flux_solutions = [flux_solution for irho in par.rho_transport]
    elseif actor.par.model == :BgB
        A1  = 1.0
        A2 = 0.3
        αgB = 5e-6
        αB = 2e-4
    
        rho_cp = cp1d.grid.rho_tor_norm
        rho_eq = eqt1d.rho_tor_norm
        
        bax = eqt.global_quantities.magnetic_axis.b_field_tor
        q_profile = IMAS.interp1d(rho_eq, eqt1d.q).(rho_cp)
        vprime_miller = IMAS.interp1d(rho_eq, GACODE.volume_prime_miller_correction(eqt)).(rho_cp)

        rmin = GACODE.r_min_core_profiles(eqt1d, rho_cp)
        Te = cp1d.electrons.temperature
        dlntedr = .-IMAS.calc_z(rho_cp, Te, :backward)
        ne = cp1d.electrons.density_thermal
        dlnnedr = .-IMAS.calc_z(rho_cp, ne, :backward)

        Ti = cp1d.ion[1].temperature
        zeff = cp1d.zeff
        dlntidr = .-IMAS.calc_z(rho_cp, Ti, :backward)

        Q_GB = GACODE.gyrobohm_energy_flux(cp1d,eqt)
        Γ_GB = GACODE.gyrobohm_particle_flux(cp1d,eqt)

        dpe = @. ne*IMAS.mks.e*Te*(dlntedr .+ dlnnedr) 
        dpi = @. ne*IMAS.mks.e*Ti*(dlntidr .+ dlnnedr) / zeff
        χeB = @.  αB * rmin[end] * q_profile^2 / abs(bax) * Te * (dlntedr .+ dlnnedr)
        χeGB = @. αgB * sqrt(Te/bax^2) * Te * dlntedr 

        χiB = 2.0*χeB
        χiGB = 0.5*χeGB

        χe = @. actor.par.αBgB * (0.01*χeB+50.0*χeGB) 
        χi = @. actor.par.αBgB * (0.001*χiB+χiGB)
        Γe = @. (A1+(A2-A1)*rho_cp)*χe*χi/(χe+χi) * dlnnedr  * ne
        gridpoint_cp = [argmin_abs(cp1d.grid.rho_tor_norm, ρ) for ρ in par.rho_transport]

        actor.flux_solutions = [
            GACODE.FluxSolution(
                χe[irho] * dpe[irho] / Q_GB[irho] / vprime_miller[irho] ,
                χi[irho] * dpi[irho] / Q_GB[irho] / vprime_miller[irho],
                Γe[irho]/Γ_GB[irho] / vprime_miller[irho],
                [Γe[irho]/Γ_GB[irho] / vprime_miller[irho]  / ion.element[1].z_n for ion in cp1d.ion],
                1.0
            )
            for irho in gridpoint_cp
        ]

        if par.do_plot
            plot(rho_cp, χeB,ylabel="χ", xlabel="ρ", label="χeB")
            display(plot!(rho_cp, 5e3*χeGB,ylabel="χ", xlabel="ρ", label="5000 χeGB"))
        end
    end

    return actor
end

"""
    _finalize(actor::ActorAnalyticTurbulence)

Writes results to dd.core_transport
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

    GACODE.flux_gacode_to_imas((:electron_energy_flux, :ion_energy_flux,  :electron_particle_flux, :ion_particle_flux, :momentum_flux), actor.flux_solutions, m1d, eqt, cp1d)

    return actor
end
