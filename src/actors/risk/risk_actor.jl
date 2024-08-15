using FusionMaterials 

#= ============ =#
#  ActorRisk  #
#= ============ =#
Base.@kwdef mutable struct FUSEparameters__ActorRisk{T<:Real} <: ParametersActor{T}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    _time::Float64 = NaN
    models::Switch{Symbol} = Switch{Symbol}([:engineering, :plasma, :both], "-", "Risk model to run"; default = :both)
    trl_to_risk::Switch{Symbol} = Switch{Symbol}([:exponential, :linear], "-", "Relationship between increasing TRL and decreasing risk"; default = :exponential)
    fraction_disruptions_mitigated::Entry{Real} = Entry{Real}("-", "Fraction of disruptions that are successfully mitigated"; default = 0.95)
    dwell_time::Entry{Real} = Entry{Real}("-", "Time between pulses as a fraction of flattop duration"; default = 0.2)
end

mutable struct ActorRisk{D,P} <: CompoundAbstractActor{D,P}
    dd::IMAS.dd{D}
    par::FUSEparameters__ActorRisk{P}
    act::ParametersAllActors
end

"""
	ActorRisk(dd::IMAS.dd, act::ParametersAllActors; kw...)

Quantifies risk associated to engineering and/or plasma design choices

!!! note 
	Stores data in `dd.risk`
"""
function ActorRisk(dd::IMAS.dd, act::ParametersAllActors; kw...)
    actor = ActorRisk(dd, act.ActorRisk, act; kw...)
    step(actor)
    finalize(actor)
    return actor
end

function _step(actor::ActorRisk)
    dd = actor.dd
    par = actor.par
    rsk = dd.risk 
    cst = dd.costing 
    bd = dd.build 

    empty!(rsk)

    ###############
    # Engineering #
    ###############
    n_engineering = 0
    if !isempty(dd.blanket)
        n_engineering += 1
    end

    if cst.model == "ARIES"
        n_engineering += 3 # each magnet system costed individually in ARIES
    elseif cst.model == "Sheffield"
        n_engineering += 1 # magnets costed as a single system in Sheffield
    end
    
    loss = resize!(rsk.engineering.loss, n_engineering)

    # Coils

    coil_types = Dict("TF" => bd.tf.technology, "OH" => bd.oh.technology, "PF" => bd.pf_active.technology)
    if cst.model == "ARIES"
        
        for (i,coil_type) in enumerate(keys(coil_types))
            loss[i].description = "$coil_type failure"
            loss[i].probability = magnet_risk(coil_types[coil_type], :exponential)
            loss[i].severity = IMAS.select_direct_captial_cost(dd, coil_type)
        end
    elseif cst.model == "Sheffield"
        loss[1].description = "Magnet failure"
        loss[1].probability = 0.0
        for coil_type in keys(coil_types)
            loss[1].probability += 1/3 * magnet_risk(coil_types[coil_type], :exponential)
        end
        tokamak_idx = findfirst(x->x.name == "tokamak", cst.cost_direct_capital.system)
        pc_idx = findfirst(x->x.name == "primary coils", cst.cost_direct_capital.system[tokamak_idx].subsystem)
        pc_cost = cst.cost_direct_capital.system[tokamak_idx].subsystem[pc_idx].cost 
        loss[1].severity = pc_cost
    end

    # Blanket 
    if !isempty(dd.blanket)
        next_idx = findfirst(isempty(loss) for loss in dd.risk.engineering.loss)
        loss[next_idx].description = "Blanket failure"
        loss[next_idx].probability = blanket_risk(dd, :exponential)
        if cst.model == "ARIES"
            loss[next_idx].severity = IMAS.select_direct_captial_cost(dd, "blanket")
        elseif cst.model == "Sheffield"
            da = DollarAdjust(dd)
            blanket_cost = cost_direct_capital_Sheffield(:blanket, true, dd, da)
            loss[next_idx].severity = blanket_cost
        end
    end

    ##########
    # Plasma #
    ##########
    if isempty(dd.stability.model)
        @warn "Plasma risk won't be calculated; must run ActorStabilityLimits before ActorRisk to compute plasma risk"
    end

    plasma_failure_probability(dd)
    plasma_failure_severity(dd)

    recovery_time = disruption_recovery_time(dd.risk.plasma.loss[1].severity) # units are years

    # here, disruption probability is expected to be in units of disruptions per year 

    if !isempty(rsk.plasma.loss)
        for loss in rsk.plasma.loss
            disrupt_per_year = disruptions_per_year(loss.probability, par.fraction_disruptions_mitigated, par.dwell_time, dd)
            adjusted_LCoE = cst.levelized_CoE / (1 - fraction_recovery(disrupt_per_year, recovery_time))
            loss.risk = adjusted_LCoE - cst.levelized_CoE
        end
    end
    
    return actor
end

function _finalize(actor::ActorRisk)
    sort!(actor.dd.risk.engineering.loss; by=x -> x.risk, rev=true)
    sort!(actor.dd.risk.plasma.loss; by=x -> x.risk, rev=true)
    return actor 
end 

function magnet_risk(coil_tech::Union{Missing,FusionMaterials.IMAS_build_coil_techs}, function_type::Symbol)
    conductor_trl = Dict("copper" => 9, "nbti" => 8, "rebco" => 3, "nb3sn" => 7, "nb3sn_iter" => 7, "nb3sn_kdemo" => 3)
    return trl_to_risk(function_type, conductor_trl[coil_tech.material]) 
end

function blanket_risk(dd::IMAS.dd, function_type::Symbol)
    bd = dd.build
    blanket_trl = Dict("lithium_lead" => 4, "flibe" => 3)

    if !isempty(dd.blanket)
        blanket_layer = IMAS.get_build_layer(bd.layer; type=IMAS._blanket_, fs=IMAS._hfs_)
    end

    return trl_to_risk(function_type, blanket_trl[blanket_layer.material])
end

##########
# Plasma #
##########

function disruption_probability(x::Real)
# define a tanh function fitted such that:
#   - at 0% of the stability limit, the probability of disruption is 0% 
#   - at 95% of the stability limit, the probability of disruption is 50% 
#   - at 100% of the stability limit, the probability of disruption is 100% 
    p = [52.43288920416563, 6.29025922319696, -8.619606592129893, 52.435212811655994]
    if 0.0 ≤ x ≤ 1.0
        probability = p[1] * tanh.(p[2] * x .+ p[3]) .+ p[4]
    elseif x > 1.0
        probability = 1.0
    end
    return probability
end

function disruptions_per_year(probability_per_pulse::Real, fraction_mitigated::Real, dwell_time::Real, dd::IMAS.dd)
    pulse_length_hours = dd.build.oh.flattop_duration / 3600
    pulses_per_year = (24 / ((1 + dwell_time) * pulse_length_hours)) * 365 * dd.costing.availability # 1.2 accounts for dwell time between pulses
    disruptions_per_year = pulses_per_year * probability_per_pulse * (1 - fraction_mitigated)
    return disruptions_per_year 
end

function plasma_failure_probability(dd::IMAS.dd)
    rsk = dd.risk 
    if isempty(dd.mhd_linear)
        resize!(rsk.plasma.loss, length(dd.stability.model))
    else
        resize!(rsk.plasma.loss, length(dd.stability.model) + 1)
        # Include vertical stability
        mode = dd.mhd_linear.time_slice[].toroidal_mode[1]
        rsk.plasma.loss[end].description = mode.perturbation_type.description

        if mode.growthrate ≥ 0.15 || mode.growthrate < 0.0
            rsk.plasma.loss[end].probability = 0.0
        elseif mode.growthrate < 0.15 
            rsk.plasma.loss[end].probability = disruption_probability(mode.growthrate / 0.15)
        end
    end

    for (i,model) in enumerate(dd.stability.model)
        rsk.plasma.loss[i].description = model.identifier.description
        rsk.plasma.loss[i].probability = disruption_probability(model.fraction[1])
    end    
end

function plasma_failure_severity(dd::IMAS.dd)
    for loss in dd.risk.plasma.loss
        loss.severity = disruption_severity(dd)
    end    
end

function disruption_severity(dd::IMAS.dd)
    eq = dd.equilibrium
    eqt = eq.time_slice[]
   
    ip = eqt.global_quantities.ip
    b0 = eqt.global_quantities.vacuum_toroidal_field.b0
    beta_tor = eqt.global_quantities.beta_tor
    a = eqt.boundary.minor_radius
    ϵ = a / eqt.boundary.geometric_axis.r

    iter_severity_vertical_force = 117
    severity_vertical_force = ((ip / 1e6) * abs(b0) / (1 - ϵ)) / iter_severity_vertical_force

    iter_severity_current_quench = 56
    severity_current_quench = ((ip / 1e6)^2 / a^2 ) / iter_severity_current_quench

    iter_severity_thermal_quench = 1.19
    severity_thermal_quench = sqrt(a) * b0^2 * beta_tor / iter_severity_thermal_quench

    return (severity_vertical_force + severity_current_quench + severity_thermal_quench) / 3

end

function disruption_recovery_time(severity::Real)
    # normalize relationship between disruption severity and recovery time such that for a disruption 20% more severe than an ITER full power shot, 
    # recovery time is on the order of a year 
    recovery_time = (exp.(severity) .- 1).* (1 ./ (exp.(1.2) .- 1))
    return recovery_time
end

function fraction_recovery(disruption_probability::Real, recovery_time::Real)
    # from Maris et al, the amount of time spent recovering from disruptions as a function of disruption probability and severity 
    # disruption probability expected in units of disruptions per year 
    fraction_recovery = disruption_probability * recovery_time / (1 + (disruption_probability * recovery_time))
    return fraction_recovery
end

function fraction_damage(disruption_probability::Real, average_damage::Real, plant_target_lifetime::Real)
    # also from Maris et al, the probability that a plant will have to close early due to damage from disruptions
    proportion_damage_budget = disruption_probability * average_damage
    probability_early_shutdown = 1 - exp(-plant_target_lifetime * proportion_damage_budget)
    return (1 - probability_early_shutdown)
end

###############
# Engineering # 
###############

function trl_to_risk(function_type::Symbol, trl::Real)
    if function_type == :exponential 
        k = 1/(exp(-1) + exp(-9))
        risk = k * exp(-(trl)) + (k*(exp(-9)))
    elseif function_type == :linear 
        risk = -(0.1*trl) + 1.1
    else
        error("trl_to_risk type not defined")
    end

    return risk 
end


