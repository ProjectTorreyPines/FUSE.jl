#= ============== =#
#  ActorDivertors  #
#= ============== =#
import BoundaryPlasmaModels

Base.@kwdef mutable struct FUSEparameters__ActorDivertors{T<:Real} <: ParametersActor{T}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    _time::Float64 = NaN
    heat_flux_model::Switch{Symbol} = Switch{Symbol}([:lengyel, :stangeby], "-", "Divertor heat flux model"; default=:lengyel)
    impurities::Entry{Vector{Symbol}} = Entry{Vector{Symbol}}("-", "Vector of impurity species"; default=Symbol[])
    impurities_fraction::Entry{Vector{T}} = Entry{Vector{T}}("-", "Vector of impurity fractions"; default=T[])
    heat_spread_factor::Entry{T} = Entry{T}("-", "Heat flux expansion factor in the private flux region (eg. due to transport) should be >= 1.0"; default=1.0)
    thermal_power_extraction_efficiency::Entry{T} = Entry{T}("-", "Fraction of thermal power that is carried out by the coolant at the divertor interface, rather than being lost in the surrounding strutures."; default=1.0)
    verbose::Entry{Bool} = act_common_parameters(verbose=false)
end

mutable struct ActorDivertors{D,P} <: SingleAbstractActor{D,P}
    dd::IMAS.dd{D}
    par::OverrideParameters{P,FUSEparameters__ActorDivertors{P}}
    boundary_plasma_models::Vector{BoundaryPlasmaModels.DivertorHeatFluxModel}
end

"""
    ActorDivertors(dd::IMAS.dd, act::ParametersAllActors; kw...)

Evaluates divertor loading and deposited power

!!! note

    Stores data in `dd.divertors`
"""
function ActorDivertors(dd::IMAS.dd, act::ParametersAllActors; kw...)
    actor = ActorDivertors(dd, act.ActorDivertors; kw...)
    step(actor)
    finalize(actor)
    return actor
end

function ActorDivertors(dd::IMAS.dd, par::FUSEparameters__ActorDivertors; kw...)
    logging_actor_init(ActorDivertors)
    par = OverrideParameters(par; kw...)
    return ActorDivertors(dd, par, Vector{BoundaryPlasmaModels.DivertorHeatFluxModel}())
end

function _step(actor::ActorDivertors)
    dd = actor.dd
    par = actor.par

    eqt = dd.equilibrium.time_slice[]
    cp1d = dd.core_profiles.profiles_1d[]

    sol1 = IMAS.sol(eqt, dd.wall; levels=1)[:lfs][1] # first SOL open field line on the lfs
    Psol = IMAS.power_sol(dd.core_sources, cp1d)
    λ_omp = IMAS.widthSOL_eich(eqt, Psol)

    identifiers = IMAS.identify_strike_surface(sol1, dd.divertors)
    identifiers = [(k_divertor, k_target) for (k_divertor, k_target, k_tile) in identifiers]

    for (k_divertor, divertor) in enumerate(dd.divertors.divertor)
        for (k_target, target) in enumerate(divertor.target)

            # identify which end of the field line strikes the divertor/target that we are considering
            id = (k_divertor, k_target)
            if id == identifiers[1]
                strike_index = 1
            elseif id == identifiers[2]
                strike_index = length(sol1.r)
            else
                @ddtime(target.power_conducted.data = 0.0)
                @ddtime(target.power_convected.data = 0.0)
                @ddtime(target.power_incident.data = 0.0)
                continue
            end

            # Power conducted/convected by the plasma on this divertor target  [W]
            @ddtime(target.power_conducted.data = Psol / 2.0)
            @ddtime(target.power_convected.data = 0.0)

            # λq at the outer midplane [m]
            resize!(target.two_point_model)
            target.two_point_model[].sol_heat_decay_length = λ_omp

            # target flux expansion [-] and wetted area [m²]
            flxexp = @ddtime(target.flux_expansion.data = sol1.total_flux_expansion[strike_index])
            λ_target = flxexp * λ_omp
            wetted_area = @ddtime(target.wetted_area.data = λ_target * 2π * sol1.r[strike_index])

            # strike angles [rad]
            @ddtime(target.tilt_angle_tor.data = atan(sol1.Bp[strike_index] / sol1.Bt[strike_index]))
            @ddtime(target.tilt_angle_pol.data = sol1.strike_angles[strike_index == 1 ? 1 : 2])

            # Setup model based on dd and run
            boundary_plasma_model = BoundaryPlasmaModels.DivertorHeatFluxModel(par.heat_flux_model)
            push!(actor.boundary_plasma_models, boundary_plasma_model)
            BoundaryPlasmaModels.setup_model(boundary_plasma_model, target, eqt, cp1d, sol1; par.impurities, par.impurities_fraction, par.heat_spread_factor)
            boundary_plasma_model() # run
            if par.verbose
                println("      == $(uppercase(target.name)) ==")
                BoundaryPlasmaModels.summary(boundary_plasma_model)
                println()
            end

            # get the peak power heatflux [W.m²]
            power_flux_peak = @ddtime(target.power_flux_peak.data = boundary_plasma_model.results.q_perp_target_spread)

            # get the total power on the divertor [W]
            power_incident = power_flux_peak * wetted_area
            @ddtime(target.power_incident.data = power_incident)
        end

        # update totals
        IMAS.divertor_totals_from_targets!(divertor)

        @ddtime(divertor.power_thermal_extracted.data = par.thermal_power_extraction_efficiency * @ddtime(divertor.power_incident.data))
    end

    return actor
end
