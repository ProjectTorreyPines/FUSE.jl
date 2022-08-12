#= ========================= =#
#  ActorEquilibriumTransport  #
#= ========================= =#
mutable struct ActorEquilibriumTransport <: PlasmaAbstractActor
    dd::IMAS.dd
    par::ParametersActor
    act::ParametersAllActors
    actor_jt::ActorSteadyStateCurrent
    actor_eq::ActorEquilibrium
    actor_tr::ActorTauenn
end

function ParametersActor(::Type{Val{:ActorEquilibriumTransport}})
    par = ParametersActor(nothing)
    par.do_plot = Entry(Bool, "", "plot"; default=false)
    par.iterations = Entry(Int, "", "transport-equilibrium iterations"; default=1)
    return par
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

function ActorEquilibriumTransport(dd::IMAS.dd, par::ParametersActor, act::ParametersAllActors; kw...)
    par = par(kw...)
    actor_jt = ActorSteadyStateCurrent(dd, act.ActorSteadyStateCurrent)
    actor_eq = ActorEquilibrium(dd, act.ActorEquilibrium, act)
    actor_tr = ActorTauenn(dd, act.ActorTauenn)
    return ActorEquilibriumTransport(dd, par, act, actor_jt, actor_eq, actor_tr)
end

function step(actor::ActorEquilibriumTransport)
    dd = actor.dd
    par = actor.par
    act = actor.act

    if par.do_plot
        pe = plot(dd.equilibrium; color=:gray, label="old", coordinate=:rho_tor_norm)
        pp = plot(dd.core_profiles; color=:gray, label=" (old)")
        ps = plot(dd.core_sources; color=:gray, label=" (old)")
    end

    # Set j_ohmic to steady state
    finalize(step(actor.actor_jt))

    for iteration in 1:par.iterations
        # run transport actor
        finalize(step(actor.actor_tr))

        # prepare equilibrium input based on transport core_profiles output
        prepare(dd, :ActorEquilibrium, act)

        # run equilibrium actor with the updated beta
        finalize(step(actor.actor_eq))

        # Set j_ohmic to steady state
        finalize(step(actor.actor_jt))
    end

    if par.do_plot
        display(plot!(pe, dd.equilibrium, coordinate=:rho_tor_norm))
        display(plot!(pp, dd.core_profiles))
        display(plot!(ps, dd.core_sources))
    end

    return actor
end


#= ================== =#
#  ActorWholeFacility  #
#= ================== =#
mutable struct ActorWholeFacility <: FacilityAbstractActor
    dd::IMAS.dd
end

function ParametersActor(::Type{Val{:ActorWholeFacility}})
    par = ParametersActor(nothing)
    return par
end

"""
    ActorWholeFacility(dd::IMAS.dd, act::ParametersAllActors; kw...)

Compound actor that runs all the actors needed to model the whole plant

!!! note 
    Stores data in `dd`
"""
function ActorWholeFacility(dd::IMAS.dd, act::ParametersAllActors; kw...)
    par = act.ActorWholeFacility(kw...)
    actor = ActorWholeFacility(dd)
    step(actor; act)
    finalize(actor)
    return actor
end

function step(actor::ActorWholeFacility; act::Union{Missing,ParametersAllActors}=missing, iterations::Int=1, do_plot::Bool=false)
    dd = actor.dd
    ActorEquilibriumTransport(dd, act)
    ActorHFSsizing(dd, act)
    ActorLFSsizing(dd, act)
    ActorCXbuild(dd, act)
    ActorPFcoilsOpt(dd, act)
    ActorNeutronics(dd, act)
    ActorBlanket(dd, act)
    ActorDivertors(dd, act)
    ActorBalanceOfPlant(dd, act)
    ActorCosting(dd, act)
    return actor
end