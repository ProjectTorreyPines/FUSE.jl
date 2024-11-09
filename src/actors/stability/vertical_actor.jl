import VacuumFields

#= ====================== =#
#  ActorVerticalStability  #
#= ====================== =#
Base.@kwdef mutable struct FUSEparameters__ActorVerticalStability{T<:Real} <: ParametersActor{T}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    _time::Float64 = NaN
    #== actor parameters ==#
    model::Entry{Bool} = Entry{Bool}("-", "Tunr on/off model of vertical stability"; default=true)
    #== display and debugging parameters ==#
    do_plot::Entry{Bool} = act_common_parameters(; do_plot=false)
end

mutable struct ActorVerticalStability{D,P} <: CompoundAbstractActor{D,P}
    dd::IMAS.dd{D}
    par::FUSEparameters__ActorVerticalStability{P}
    act::ParametersAllActors{P}
    stability_margin::D
    normalized_growth_rate::D
    passive_coils::Vector{VacuumFields.MultiCoil}
end

"""
    ActorVerticalStability(dd::IMAS.dd, act::ParametersAllActors; kw...)

Compute vertical stability metrics
"""
function ActorVerticalStability(dd::IMAS.dd, act::ParametersAllActors; kw...)
    par = act.ActorVerticalStability(kw...)
    actor = ActorVerticalStability(dd, par, act)
    step(actor)
    finalize(actor)
    return actor
end

function ActorVerticalStability(dd::IMAS.dd{D}, par::FUSEparameters__ActorVerticalStability{P}, act::ParametersAllActors{P}; kw...) where {D<:Real,P<:Real}
    logging_actor_init(ActorVerticalStability)
    par = par(kw...)
    return ActorVerticalStability(dd, par, act, D(NaN), D(NaN), VacuumFields.MultiCoil[])
end

"""
    _step(actor::ActorVerticalStability)

Compute vertical stability metrics
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
    if all(VacuumFields.current(coil) == 0.0 for coil in active_coils)
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

Store vertical stability metrics
"""
function _finalize(actor::ActorVerticalStability)
    dd = actor.dd
    par = actor.par

    mhd = resize!(dd.mhd_linear.time_slice)
    resize!(mhd.toroidal_mode, 2)

    # Stability margin
    mode = mhd.toroidal_mode[1]
    mode.perturbation_type.description = "Vertical stability margin, > 0.15 for stability (N.B., not in Hz)"
    mode.perturbation_type.name = "m_s"
    mode.n_tor = 0
    mode.growthrate = actor.stability_margin # not in Hz

    # Normalized growth rate
    mode = mhd.toroidal_mode[2]
    mode.perturbation_type.description = "Normalized vertical growth rate, < 10 for stability (N.B., not in Hz)"
    mode.perturbation_type.name = "γτ"
    mode.n_tor = 0
    mode.growthrate = actor.normalized_growth_rate # not in Hz

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
