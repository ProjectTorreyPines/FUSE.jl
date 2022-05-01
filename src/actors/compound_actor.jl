
#= ========================= =#
#  EquilibriumTransportActor  #
#= ========================= =#
mutable struct EquilibriumTransportActor <: AbstractActor
    dd::IMAS.dd
end

function ActorParameters(::Type{Val{:EquilibriumTransportActor}})
    par = ActorParameters(nothing)
    par.do_plot = Entry(Bool, "", "plot"; default=false)
    par.iterations = Entry(Int, "", "transport-equilibrium iterations"; default=1)
    return par
end

function EquilibriumTransportActor(dd::IMAS.dd, act::ActorParameters; kw...)
    par = act.EquilibriumTransportActor(kw...)
    actor = EquilibriumTransportActor(dd)
    step(actor; act, iterations=par.iterations, do_plot=par.do_plot)
    finalize(actor)
    return actor
end

function step(actor::EquilibriumTransportActor;  act::Union{Missing,ActorParameters}=missing, iterations::Int=1, do_plot::Bool=false)
    dd = actor.dd
    if act === missing
        act = ActorParameters()
    end

    # Set j_ohmic to steady state
    IMAS.j_ohmic_steady_state!(dd.equilibrium.time_slice[], dd.core_profiles.profiles_1d[])

    for iteration in 1:iterations
        # run transport actor
        TauennActor(dd, act)

        # Set beta_normal from equilbrium to the kinetic beta_n
        if !isempty(dd.core_profiles.profiles_1d)
            dd.equilibrium.time_slice[].global_quantities.beta_normal = @ddtime(dd.summary.global_quantities.beta_tor_thermal_norm.value)
        end

        # run equilibrium actor with the updated beta
        SolovevActor(dd, act)

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