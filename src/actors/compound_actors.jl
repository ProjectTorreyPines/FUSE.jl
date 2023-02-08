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

Base.@kwdef struct FUSEparameters__ActorEquilibriumTransport{T} <: ParametersActor where {T<:Real}
    do_plot = Entry(Bool, "", "plot"; default=false)
    iterations = Entry(Int, "", "transport-equilibrium iterations"; default=1)
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
        pe = plot(dd.equilibrium; color=:gray, label="old", coordinate=:rho_tor_norm)
        pp = plot(dd.core_profiles; color=:gray, label=" (old)")
        ps = plot(dd.core_sources; color=:gray, label=" (old)")
    end

    # Set j_ohmic to steady state
    finalize(step(actor.actor_jt))

    @assert act.ActorCHEASE.rescale_eq_to_ip "Running CHEASE with ActorEquilibriumTransport requires ActorCHEASE.rescale_eq_to_ip = true"
    act_chease = deepcopy(act.ActorCHEASE)

    actor.actor_eq.par = act_chease
    act_chease.free_boundary = false

    max_iter = 5
    iter = 0
    ip = dd.equilibrium.time_slice[].global_quantities.ip
    conv_criteria = ip / 1e2
    avg_diff = conv_criteria + 1

    while avg_diff > conv_criteria
        # run transport actor
        finalize(step(actor.actor_tr))

        # Set j_ohmic to steady state
        finalize(step(actor.actor_jt))

        j_tor_before = dd.core_profiles.profiles_1d[].j_tor

        # prepare equilibrium input based on transport core_profiles output
        prepare(dd, :ActorEquilibrium, act)

        # run equilibrium actor with the updated beta
        finalize(step(actor.actor_eq))

        j_tor_after = dd.core_profiles.profiles_1d[].j_tor

        avg_diff = sum(abs.(j_tor_after .- j_tor_before)) / length(j_tor_after)

        if act.ActorEquilibrium.model == :Solovev
            avg_diff = 0 # temporary fix to force Solovev to run exactly once
        end

        if iter == max_iter
            @warn "Max number of iterations has been reached ($max_iter), current difference is $(round(avg_diff,digits = 3)), convergence criteria is $conv_criteria"
            break
        end
        iter += 1

    end

    if act.ActorCHEASE.free_boundary
        act_chease.free_boundary = true
        finalize(step(actor.actor_eq))
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
    par::ParametersActor
    act::ParametersAllActors
    EquilibriumTransport::Union{Nothing,ActorEquilibriumTransport}
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

Base.@kwdef struct FUSEparameters__ActorWholeFacility{T} <: ParametersActor where {T<:Real}
end

"""
    ActorWholeFacility(dd::IMAS.dd, act::ParametersAllActors; kw...)

Compound actor that runs all the physics, engineering and costing actors needed to model the whole plant:
* ActorEquilibriumTransport
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

function ActorWholeFacility(dd::IMAS.dd, par::ParametersActor, act::ParametersAllActors; kw...)
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
        nothing)
end

function _step(actor::ActorWholeFacility)
    dd = actor.dd
    act = actor.act
    actor.EquilibriumTransport = ActorEquilibriumTransport(dd, act)
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