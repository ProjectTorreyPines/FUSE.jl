#= ========================= =#
#  ActorEquilibriumTransport  #
#= ========================= =#
mutable struct ActorEquilibriumTransport <: AbstractActor
    dd::IMAS.dd
end

function ParametersActor(::Type{Val{:ActorEquilibriumTransport}})
    par = ParametersActor(nothing)
    par.do_plot = Entry(Bool, "", "plot"; default=false)
    par.iterations = Entry(Int, "", "transport-equilibrium iterations"; default=1)
    return par
end

"""
    ActorEquilibriumTransport(dd::IMAS.dd, act::ParametersActor; kw...)

Compound actor that runs the following actors in succesion:
```julia
ActorSteadyStateCurrent(dd, act)    # Current evolution to steady-state 
ActorTauenn(dd, act)                # For transport
ActorSolovev(dd, act)               # Equilibrium
ActorSteadyStateCurrent(dd, act)    # Consistent current
```

!!! note 
    Stores data in `dd.equilibrium, dd.core_profiles, dd.core_sources`
"""
function ActorEquilibriumTransport(dd::IMAS.dd, act::ParametersActor; kw...)
    par = act.ActorEquilibriumTransport(kw...)
    actor = ActorEquilibriumTransport(dd)
    step(actor; act, iterations=par.iterations, do_plot=par.do_plot)
    finalize(actor)
    return actor
end

function step(actor::ActorEquilibriumTransport; act::Union{Missing,ParametersActor}=missing, iterations::Int=1, do_plot::Bool=false)
    dd = actor.dd
    if act === missing
        act = ParametersActor()
    end

    if do_plot
        pe = plot(dd.equilibrium; color=:gray, label="")
        pp = plot(dd.core_profiles; color=:gray, label="")
        ps = plot(dd.core_sources; color=:gray)
    end

    # Set j_ohmic to steady state
    ActorSteadyStateCurrent(dd, act)

    for iteration in 1:iterations
        # run transport actor
        ActorTauenn(dd, act)

        # Set beta_normal from equilbrium to the kinetic beta_n
        if !isempty(dd.core_profiles.profiles_1d)
            dd.equilibrium.time_slice[].global_quantities.beta_normal = @ddtime(dd.summary.global_quantities.beta_tor_thermal_norm.value)
        end

        # run equilibrium actor with the updated beta
        ActorSolovev(dd, act)

        # Set j_ohmic to steady state
        ActorSteadyStateCurrent(dd, act)
    end

    if do_plot
        display(plot!(pe, dd.equilibrium))
        display(plot!(pp, dd.core_profiles))
        display(plot!(ps, dd.core_sources))
    end

    return dd
end

#= ================ =#
#  ActorWholeDevice  #
#= ================ =#
mutable struct ActorWholeDevice <: AbstractActor
    dd::IMAS.dd
end

function ParametersActor(::Type{Val{:ActorWholeDevice}})
    par = ParametersActor(nothing)
    return par
end

"""
    ActorWholeDevice(dd::IMAS.dd, act::ParametersActor; kw...)

Compound actor that runs the following actors in succesion:
```julia
ActorEquilibriumTransport(dd, act)
ActorHFSsizing(dd, act)
ActorLFSsizing(dd, act)
ActorCXbuild(dd, act)
ActorPFcoilsOpt(dd, act)
ActorCosting(dd, act)
ActorBalanceOfPlant(dd,act)    
```

!!! note 
    Stores data in `dd`
"""
function ActorWholeDevice(dd::IMAS.dd, act::ParametersActor; kw...)
    par = act.ActorWholeDevice(kw...)
    actor = ActorWholeDevice(dd)
    step(actor; act)
    finalize(actor)
    return actor
end

function step(actor::ActorWholeDevice; act::Union{Missing,ParametersActor}=missing, iterations::Int=1, do_plot::Bool=false)
    dd = actor.dd
    ActorEquilibriumTransport(dd, act)
    ActorHFSsizing(dd, act)
    ActorLFSsizing(dd, act)
    ActorCXbuild(dd, act)
    #ActorNeutronics(dd, act) # not really connected to anything at this point
    ActorPFcoilsOpt(dd, act)
    ActorCosting(dd, act)
    ActorBalanceOfPlant(dd,act)
end