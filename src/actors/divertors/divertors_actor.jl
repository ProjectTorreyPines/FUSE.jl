#= ============== =#
#  ActorDivertors  #
#= ============== =#
import BoundaryPlasmaModels

Base.@kwdef mutable struct FUSEparameters__ActorDivertors{T<:Real} <: ParametersActor{T}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    _time::Float64 = NaN
    boundary_model::Switch{Symbol} = Switch{Symbol}([:slcoupled, :lengyel, :stangeby], "-", "SOL boundary model"; default=:slcoupled)
    impurities::Entry{Vector{Symbol}} = Entry{Vector{Symbol}}("-", "Vector of impurity species"; default=Symbol[])
    impurities_fraction::Entry{Vector{T}} = Entry{Vector{T}}("-", "Vector of impurity fractions"; default=T[])
    heat_spread_factor::Entry{T} = Entry{T}("-", "Heat flux expansion factor in the private flux region (eg. due to transport) should be >= 1.0"; default=one(T))
    thermal_power_extraction_efficiency::Entry{T} = Entry{T}("-", "Fraction of thermal power that is carried out by the coolant at the divertor interface, rather than being lost in the surrounding strutures."; default=one(T))
    verbose::Entry{Bool} = act_common_parameters(verbose=false)
end

mutable struct ActorDivertors{D,P} <: SingleAbstractActor{D,P}
    dd::IMAS.dd{D}
    par::OverrideParameters{P,FUSEparameters__ActorDivertors{P}}
    boundary_plasma_models::Vector{BoundaryPlasmaModels.SOLBoundaryModel}
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
    return ActorDivertors(dd, par, Vector{BoundaryPlasmaModels.SOLBoundaryModel}())
end

function _step(actor::ActorDivertors)
    dd = actor.dd
    par = actor.par

    eqt = dd.equilibrium.time_slice[]
    cp1d = dd.core_profiles.profiles_1d[]
    cs = dd.core_sources

    sol1 = IMAS.sol(eqt, dd.wall; levels=1)[:lfs][1] # first SOL open field line on the lfs
    Psol = IMAS.power_sol(dd.core_sources, cp1d)
    λ_omp = IMAS.widthSOL_eich(eqt, Psol)
    index_inner = 1
    index_outer = length(sol1.r)
    if sol1.r[1] < sol1.r[end]
        index_inner = 1
        index_outer = length(sol1.r)
    else
        index_inner = length(sol1.r)
        index_outer = 1
    end

    strike_indices = (index_outer, index_inner)
    empty!(actor.boundary_plasma_models)
    for (i, strike_index) in enumerate(strike_indices)
        boundary_plasma_model = BoundaryPlasmaModels.SOLBoundaryModel(par.boundary_model)
        push!(actor.boundary_plasma_models, boundary_plasma_model)
        BoundaryPlasmaModels.setup_model(boundary_plasma_model, λ_omp, eqt, cp1d, cs, sol1; par.impurities, par.impurities_fraction, par.heat_spread_factor,strike_index)
        boundary_plasma_model() # run
    end
    # --- Decide how to write back results --- #
    if !_has_any_divertor(dd)
        write_divertor_totals_only_from_sol!(dd, sol1, λ_omp, actor.boundary_plasma_models, par)
    return actor
    end

    if !_has_any_target(dd)
        write_divertor_totals_only_from_sol!(dd, sol1, λ_omp, actor.boundary_plasma_models, par)
        return actor
    end

    populate_divertor_targets_from_sol!(dd, sol1, λ_omp, actor.boundary_plasma_models, Psol, par, strike_indices)
    return actor
end

function populate_divertor_targets_from_sol!(
    dd,
    sol1,
    λ_omp,
    boundary_models::Vector,
    Psol,
    par,
    strike_indices::Tuple{Int,Int};
)
    identifiers = try
        IMAS.identify_strike_surface(sol1, dd.divertors)
    catch err
        @warn "identify_strike_surface failed: $err"
        Tuple{Int,Int,Int}[(0,0,0),(0,0,0)]
    end

    if length(identifiers) < 2 || any(==( (0,0,0) ), identifiers)
        @warn "No valid strike surfaces; skip per-target write-back."
        write_divertor_totals_only_from_sol!(dd, sol1, λ_omp, boundary_models, par)
        return nothing
    end

    strike_angle_slot = (1, 2)

    for i in 1:2
        (k_div, k_tgt) = identifiers[i]
        target = dd.divertors.divertor[k_div].target[k_tgt]

        strike_index = strike_indices[i]
        bpm          = boundary_models[i]

        # λq @ OMP
        resize!(target.two_point_model)
        target.two_point_model[].sol_heat_decay_length = λ_omp

        # flux expansion & wetted area
        flxexp = @ddtime(target.flux_expansion.data = sol1.total_flux_expansion[strike_index])
        λ_target = flxexp * λ_omp
        wetted_area = @ddtime(target.wetted_area.data = λ_target * 2π * sol1.r[strike_index])

        # strike angles
        @ddtime(target.tilt_angle_tor.data = atan(sol1.Bp[strike_index] / sol1.Bt[strike_index]))
        @ddtime(target.tilt_angle_pol.data = sol1.strike_angles[strike_angle_slot[i]])

        @ddtime(target.power_flux_peak.data = bpm.results.q_perp_target_spread)

        power_incident = @ddtime(target.power_flux_peak.data) * wetted_area
        @ddtime(target.power_incident.data = power_incident)

        @ddtime(target.power_conducted.data = Psol / 2.0)
        @ddtime(target.power_convected.data = 0.0)
    end

    for divertor in dd.divertors.divertor
        IMAS.divertor_totals_from_targets!(divertor)
        @ddtime(divertor.power_thermal_extracted.data =
            par.thermal_power_extraction_efficiency * @ddtime(divertor.power_incident.data))
    end

    return nothing
end
_has_any_divertor(dd) = !isempty(dd.divertors.divertor)

function _has_any_target(dd)
    for d in dd.divertors.divertor
        if !isempty(d.target); return true; end
    end
    return false
end
"""
Write total divertor quantities using only SOL data
(used when no target geometry is available)
"""
function write_divertor_totals_only_from_sol!(
    dd,
    sol1,
    λ_omp,
    boundary_models::Vector,
    par::OverrideParameters
)
    function compute_incident_power(strike_index::Int, boundary_model)::Float64
        flux_expansion = sol1.total_flux_expansion[strike_index]
        λ_target = flux_expansion * λ_omp
        wetted_area = λ_target * 2π * sol1.r[strike_index]
        q_peak = boundary_model.results.q_perp_target_spread
        return q_peak * wetted_area
    end

    index_inner = 1
    index_outer = length(sol1.r)
    if !(sol1.r[1] < sol1.r[end])
        index_inner, index_outer = length(sol1.r), 1
    end
    strike_indices = (index_outer, index_inner)

    power_outer = compute_incident_power(strike_indices[1], boundary_models[1])
    power_inner = compute_incident_power(strike_indices[2], boundary_models[2])
    total_power_incident = power_outer + power_inner

    for divertor in dd.divertors.divertor
        @ddtime(divertor.power_incident.data = total_power_incident)
        @ddtime(divertor.power_thermal_extracted.data =
            par.thermal_power_extraction_efficiency * @ddtime(divertor.power_incident.data))
    end
    return nothing
end
