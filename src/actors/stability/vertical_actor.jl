import VacuumFields

#= ====================== =#
#  ActorVerticalStability  #
#= ====================== =#
@actor_parameters_struct ActorVerticalStability{T} begin
    #== actor parameters ==#
    model::Entry{Bool} = Entry{Bool}("-", "Tunr on/off model of vertical stability"; default=true)
    #== display and debugging parameters ==#
    do_plot::Entry{Bool} = act_common_parameters(; do_plot=false)
end

mutable struct ActorVerticalStability{D,P} <: CompoundAbstractActor{D,P}
    dd::IMAS.dd{D}
    par::OverrideParameters{P,FUSEparameters__ActorVerticalStability{P}}
    act::ParametersAllActors{P}
    stability_margin::D
    normalized_growth_rate::D
    passive_coils::Vector{VacuumFields.MultiCoil}
end

"""
    ActorVerticalStability(dd::IMAS.dd, act::ParametersAllActors; kw...)

Computes vertical stability metrics for tokamak plasmas using the equilibrium and coil configuration.

This actor calculates two key stability metrics:
- Stability margin: measures the plasma's resilience to vertical displacement events (VDEs)
- Normalized growth rate (γτ): quantifies how rapidly vertical instabilities develop

The analysis considers both active PF coils and passive conducting structures. Results are stored in
`dd.mhd_linear.time_slice[].toroidal_mode` where stability margin > 0.15 and γτ < 10 indicate stability.
"""
function ActorVerticalStability(dd::IMAS.dd, act::ParametersAllActors; kw...)
    actor = ActorVerticalStability(dd, act.ActorVerticalStability, act; kw...)
    step(actor)
    finalize(actor)
    return actor
end

function ActorVerticalStability(dd::IMAS.dd{D}, par::FUSEparameters__ActorVerticalStability{P}, act::ParametersAllActors{P}; kw...) where {D<:Real,P<:Real}
    logging_actor_init(ActorVerticalStability)
    par = OverrideParameters(par; kw...)
    return ActorVerticalStability(dd, par, act, D(NaN), D(NaN), VacuumFields.MultiCoil[])
end

"""
    _step(actor::ActorVerticalStability)

Calculates vertical stability metrics using the VacuumFields package.

The calculation involves:
1. Loading active coils from `dd.pf_active` with their current distribution
2. Loading passive conducting structures from `dd.pf_passive`
3. Creating a plasma image from the equilibrium
4. Computing stability margin using the method of images
5. Computing normalized growth rate considering circuit time constants

The method validates that all coils have positive resistance before calculating the normalized growth rate.
"""
function _step(actor::ActorVerticalStability)
    dd = actor.dd
    par = actor.par

    # Defaults
    actor.stability_margin = NaN
    actor.normalized_growth_rate = NaN

    if !par.model
        return actor
    end

    active_coils = VacuumFields.IMAS_pf_active__coils(dd; actor.act.ActorPFactive.green_model)
    # BCL 4/24/25: Could use MultiCoils, but would need to account for default resistances like IMAS_pf_active__coils
    # active_coils = VacuumFields.MultiCoils(dd.pf_active; active_only=true)
    if all(VacuumFields.current_per_turn(coil) == 0.0 for coil in active_coils)
        @warn "Active coils have no current. Can't compute vertical stability metrics"
        return actor
    end

    # load passive structures from pf_passive
    actor.passive_coils = VacuumFields.MultiCoils(dd.pf_passive)

    eqt = dd.equilibrium.time_slice[]
    Ip = eqt.global_quantities.ip
    image = VacuumFields.Image(eqt)
    coils = vcat(active_coils, actor.passive_coils)
    actor.stability_margin = VacuumFields.stability_margin(image, coils, Ip)

    for (k, coil) in enumerate(active_coils)
        if VacuumFields.resistance(coil) <= 0.0
            @warn "Active coil #$(k) has invalid resistance: $(VacuumFields.resistance(coil)). Can't compute normalized growth rate.\nOffending coil: $(repr(coil))"
            return actor
        end
    end

    for (k, coil) in enumerate(actor.passive_coils)
        if VacuumFields.resistance(coil) <= 0.0
            @warn "Passive coil #$(k) has invalid resistance: $(VacuumFields.resistance(coil)). Can't compute normalized growth rate.\nOffending coil: $(repr(coil))"
            return actor
        end
    end

    _, _, actor.normalized_growth_rate = VacuumFields.normalized_growth_rate(image, coils, Ip)

    return actor
end

"""
    _finalize(actor::ActorVerticalStability)

Stores computed stability metrics in the IMAS data structure.

Stability metrics are stored in `dd.mhd_linear.time_slice[].toroidal_mode`:
- Stability margin: stored with perturbation_type "m_s" and n_tor=0
- Normalized growth rate: stored with perturbation_type "γτ" and n_tor=0

If plotting is enabled, generates a visualization showing passive structures considered in the analysis.
"""
function _finalize(actor::ActorVerticalStability)
    dd = actor.dd
    par = actor.par

    mhd = resize!(dd.mhd_linear.time_slice; wipe=false)

    # Stability margin
    if !isnan(actor.stability_margin)
        mode = resize!(mhd.toroidal_mode, "perturbation_type.name" => "m_s", "n_tor" => 0)
        mode.perturbation_type.description = "Vertical stability margin > 0.15 for stability"
        mode.stability_metric = actor.stability_margin
    else
        deleteat!(mhd.toroidal_mode, "perturbation_type.name" => "m_s", "n_tor" => 0)
    end

    # Normalized growth rate
    if !isnan(actor.normalized_growth_rate)
        mode = resize!(mhd.toroidal_mode, "perturbation_type.name" => "γτ", "n_tor" => 0)
        mode.perturbation_type.description = "Normalized vertical growth rate < 10 for stability"
        mode.stability_metric = actor.normalized_growth_rate
    else
        deleteat!(mhd.toroidal_mode, "perturbation_type.name" => "γτ", "n_tor" => 0)
    end

    # plot
    if par.do_plot
        plot(dd.build; wireframe=true, title="Passive structures considered for vertical stability")
        colors = distinguishable_colors(length(actor.passive_coils))
        for (k, coil) in enumerate(actor.passive_coils)
            plot!(coil; color=colors[k], alpha=0.5)
        end
        display(plot!(; aspect_ratio=:equal))
    end

    return actor
end
