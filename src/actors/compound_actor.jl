
#= ========================= =#
#  ActorEquilibriumTransport  #
#= ========================= =#
mutable struct ActorEquilibriumTransport <: ActorAbstract
    dd::IMAS.dd
end

function ParametersActor(::Type{Val{:ActorEquilibriumTransport}})
    par = ParametersActor(nothing)
    par.do_plot = Entry(Bool, "", "plot"; default=false)
    par.iterations = Entry(Int, "", "transport-equilibrium iterations"; default=1)
    return par
end

function ActorEquilibriumTransport(dd::IMAS.dd, act::ParametersActor; kw...)
    par = act.ActorEquilibriumTransport(kw...)
    actor = ActorEquilibriumTransport(dd)
    step(actor; act, iterations=par.iterations, do_plot=par.do_plot)
    finalize(actor)
    return actor
end

function step(actor::ActorEquilibriumTransport;  act::Union{Missing,ParametersActor}=missing, iterations::Int=1, do_plot::Bool=false)
    dd = actor.dd
    if act === missing
        act = ParametersActor()
    end

    # Set j_ohmic to steady state
    IMAS.j_ohmic_steady_state!(dd.equilibrium.time_slice[], dd.core_profiles.profiles_1d[])

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
        IMAS.j_ohmic_steady_state!(dd.equilibrium.time_slice[], dd.core_profiles.profiles_1d[])
    end

    if do_plot
        display(plot(dd.equilibrium, psi_levels_out=[], label=ini.general.casename))
        display(plot(dd.core_profiles))
        display(plot(dd.core_sources))
    end

    return dd
end