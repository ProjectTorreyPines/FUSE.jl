#= ========================= =#
#  ActorEquilibriumTransport  #
#= ========================= =#
Base.@kwdef mutable struct FUSEparameters__ActorEquilibriumTransport{T} <: ParametersActor where {T<:Real}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    do_plot::Entry{Bool} = Entry(Bool, "-", "plot"; default=false)
    max_iter::Entry{Int} = Entry(Int, "-", "max number of transport-equilibrium iterations"; default=5)
    convergence_error::Entry{T} = Entry(T, "-", "Convergence error threshold"; default=1E-2)
end

mutable struct ActorEquilibriumTransport <: PlasmaAbstractActor
    dd::IMAS.dd
    par::FUSEparameters__ActorEquilibriumTransport
    act::ParametersAllActors
    actor_jt::ActorSteadyStateCurrent
    actor_eq::ActorEquilibrium
    actor_tr::ActorTauenn
end

"""
    ActorEquilibriumTransport(dd::IMAS.dd, act::ParametersAllActors; kw...)

Compound actor that runs the following actors in succesion:
* ActorSteadyStateCurrent
* ActorTauenn
* ActorEquilibrium

!!! note 
    Stores data in `dd.equilibrium, dd.core_profiles, dd.core_sources`
"""
function ActorEquilibriumTransport(dd::IMAS.dd, act::ParametersAllActors; kw...)
    par = act.ActorEquilibriumTransport(kw...)
    actor = ActorEquilibriumTransport(dd, par, act)
    step(actor)
    finalize(actor)
    return actor
end

function ActorEquilibriumTransport(dd::IMAS.dd, par::FUSEparameters__ActorEquilibriumTransport, act::ParametersAllActors; kw...)
    logging_actor_init(ActorEquilibriumTransport)
    par = par(kw...)
    actor_jt = ActorSteadyStateCurrent(dd, act.ActorSteadyStateCurrent)
    actor_eq = ActorEquilibrium(dd, act.ActorEquilibrium, act)
    actor_tr = ActorTauenn(dd, act.ActorTauenn)
    return ActorEquilibriumTransport(dd, par, act, actor_jt, actor_eq, actor_tr)
end

function _step(actor::ActorEquilibriumTransport)
    dd = actor.dd
    par = actor.par
    act = actor.act

    if par.do_plot
        pe = plot(dd.equilibrium; color=:gray, label=" (before)", coordinate=:rho_tor_norm)
        pp = plot(dd.core_profiles; color=:gray, label=" (before)")
        ps = plot(dd.core_sources; color=:gray, label=" (before)")
    end

    # Set j_ohmic to steady state
    finalize(step(actor.actor_jt))

    # set CHEASE switches specific to this workflow
    if actor.actor_eq.par.model == :CHEASE
        chease_par = actor.actor_eq.eq_actor.par
        orig_par_chease = deepcopy(chease_par)
        chease_par.free_boundary = false
        chease_par.rescale_eq_to_ip = true
    end

    try
        iter = 1
        total_error = 0.0
        while (iter == 1) || (total_error > par.convergence_error)
            # get current and pressure profiles before updating them
            j_tor_before = dd.core_profiles.profiles_1d[].j_tor
            pressure_before = dd.core_profiles.profiles_1d[].pressure

            # run transport actor
            finalize(step(actor.actor_tr))

            # Set j_ohmic to steady state
            finalize(step(actor.actor_jt))

            # prepare equilibrium input based on transport core_profiles output
            prepare(dd, :ActorEquilibrium, act)

            # run equilibrium actor with the updated beta
            finalize(step(actor.actor_eq))

            # evaluate change in current and pressure profiles after the update
            j_tor_after = dd.core_profiles.profiles_1d[].j_tor
            pressure_after = dd.core_profiles.profiles_1d[].pressure
            error_jtor = sum((j_tor_after .- j_tor_before) .^ 2) / sum(j_tor_before .^ 2)
            error_pressure = sum((pressure_after .- pressure_before) .^ 2) / sum(pressure_before .^ 2)
            total_error = sqrt(error_jtor + error_pressure) / 2.0

            if act.ActorEquilibrium.model == :Solovev
                total_error = 0 # temporary fix to force Solovev to run exactly once
            end

            if iter == par.max_iter
                @warn "Max number of iterations ($(par.max_iter)) has been reached with convergence error of $(round(total_error,digits = 3)) compared to threshold of $(par.convergence_error)"
                break
            end

            if par.do_plot
                println("Iteration = $iter , convergence error = $(round(total_error,digits = 3)), threshold = $(par.convergence_error)")
            end
            iter += 1
        end

    finally
        if actor.actor_eq.par.model == :CHEASE
            actor.actor_eq.eq_actor.par = orig_par_chease
            if orig_par_chease.free_boundary
                finalize(step(actor.actor_eq))
            end
        end
    end

    if par.do_plot
        display(plot!(pe, dd.equilibrium, coordinate=:rho_tor_norm, label=" (after)"))
        display(plot!(pp, dd.core_profiles, label=" (after)"))
        display(plot!(ps, dd.core_sources, label=" (after)"))
    end

    return actor
end


#= ================== =#
#  ActorWholeFacility  #
#= ================== =#
Base.@kwdef mutable struct FUSEparameters__ActorWholeFacility{T} <: ParametersActor where {T<:Real}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
end

mutable struct ActorWholeFacility <: FacilityAbstractActor
    dd::IMAS.dd
    par::FUSEparameters__ActorWholeFacility
    act::ParametersAllActors
    EquilibriumTransport::Union{Nothing,ActorEquilibriumTransport}
    PlasmaLimits::Union{Nothing,ActorPlasmaLimits}
    HFSsizing::Union{Nothing,ActorHFSsizing}
    LFSsizing::Union{Nothing,ActorLFSsizing}
    CXbuild::Union{Nothing,ActorCXbuild}
    PFcoilsOpt::Union{Nothing,ActorPFcoilsOpt}
    PassiveStructures::Union{Nothing,ActorPassiveStructures}
    Neutronics::Union{Nothing,ActorNeutronics}
    Blanket::Union{Nothing,ActorBlanket}
    Divertors::Union{Nothing,ActorDivertors}
    BalanceOfPlant::Union{Nothing,ActorBalanceOfPlant}
    Costing::Union{Nothing,ActorCosting}
end

"""
    ActorWholeFacility(dd::IMAS.dd, act::ParametersAllActors; kw...)

Compound actor that runs all the physics, engineering and costing actors needed to model the whole plant:
* ActorEquilibriumTransport
* ActorPlasmaLimits
* ActorHFSsizing
* ActorLFSsizing
* ActorCXbuild
* ActorPFcoilsOpt
* ActorPassiveStructures
* ActorNeutronics
* ActorBlanket
* ActorDivertors
* ActorBalanceOfPlant
* ActorCosting

!!! note 
    Stores data in `dd`
"""
function ActorWholeFacility(dd::IMAS.dd, act::ParametersAllActors; kw...)
    par = act.ActorWholeFacility(kw...)
    actor = ActorWholeFacility(dd, par, act)
    step(actor)
    finalize(actor)
    return actor
end

function ActorWholeFacility(dd::IMAS.dd, par::FUSEparameters__ActorWholeFacility, act::ParametersAllActors; kw...)
    logging_actor_init(ActorWholeFacility)
    par = par(kw...)

    ActorWholeFacility(dd, par, act,
        nothing,
        nothing,
        nothing,
        nothing,
        nothing,
        nothing,
        nothing,
        nothing,
        nothing,
        nothing,
        nothing,
        nothing)
end

function _step(actor::ActorWholeFacility)
    dd = actor.dd
    act = actor.act
    actor.EquilibriumTransport = ActorEquilibriumTransport(dd, act)
    actor.PlasmaLimits == ActorPlasmaLimits(dd, act)
    actor.HFSsizing = ActorHFSsizing(dd, act)
    actor.LFSsizing = ActorLFSsizing(dd, act)
    actor.CXbuild = ActorCXbuild(dd, act)
    actor.PFcoilsOpt = ActorPFcoilsOpt(dd, act)
    if act.ActorPFcoilsOpt.update_equilibrium && act.ActorCXbuild.rebuild_wall
        actor.CXbuild = ActorCXbuild(dd, act)
        actor.PFcoilsOpt = ActorPFcoilsOpt(dd, act; update_equilibrium=false)
    end
    actor.PassiveStructures = ActorPassiveStructures(dd, act)
    actor.Neutronics = ActorNeutronics(dd, act)
    actor.Blanket = ActorBlanket(dd, act)
    actor.Divertors = ActorDivertors(dd, act)
    actor.BalanceOfPlant = ActorBalanceOfPlant(dd, act)
    actor.Costing = ActorCosting(dd, act)
    return actor
end