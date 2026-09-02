#= ================== =#
#  ActorBetaMatch  #
#= ================== =#
@actor_parameters_struct ActorBetaMatch{T} begin
    βn_target::Entry{T} = Entry{T}("-", "Target value for normalized toroidal beta")
    T_shaping::Entry{T} = Entry{T}("-", "Shaping coefficient for the temperature profile")
    thermal::Entry{Bool} = Entry{Bool}("-", "Match thermal beta or total beta"; default=true)
    norm::Entry{Bool} = Entry{Bool}("-", "Match normalized beta or toroidal beta"; default=true)
    fix_pedestal::Entry{Bool} = Entry{Bool}("-", "Only scale core temperature using Hmode_Profiles"; default=true)
end

mutable struct ActorBetaMatch{D,P} <: CompoundAbstractActor{D,P}
    dd::IMAS.DD{D}
    par::OverrideParameters{P,FUSEparameters__ActorBetaMatch{P}}
    act::ParametersAllActors{P}
    function ActorBetaMatch(dd::IMAS.DD{D}, par::FUSEparameters__ActorBetaMatch{P}, act::ParametersAllActors{P}; kw...) where {D<:Real,P<:Real}
        logging_actor_init(ActorBetaMatch)
        par = OverrideParameters(par; kw...)
        return new{D,P}(dd, par, act)
    end
end

"""
    ActorBetaMatch(dd::IMAS.dd, act::ParametersAllActors; kw...)

Updates core temperature profiles in order to match target betaN.
Uses Hmode_profiles() function to scale core temperature profiles while keeping pedestal unchanged.
Does not change density profiles.
Does not change pedestal parameters.

"""
function ActorBetaMatch(dd::IMAS.DD, act::ParametersAllActors; kw...)
    actor = ActorBetaMatch(dd, act.ActorBetaMatch, act; kw...)
    step(actor)
    finalize(actor)
    return actor
end

"""
    _step(actor::ActorBetaMatch)
"""
function _step(actor::ActorBetaMatch)
    dd = actor.dd
    par = actor.par
    cp1d = dd.core_profiles.profiles_1d[]

    # get profile parameters
    wped = 1 - dd.summary.local.pedestal.position.rho_tor_norm[1]
    te = dd.core_profiles.profiles_1d[].electrons.temperature
    ti = dd.core_profiles.profiles_1d[].ion[1].temperature
    te_sep = te[end]
    ti_sep = ti[end]
    te_ped = dd.summary.local.pedestal.t_e.value[1]
    ti_ped = dd.summary.local.pedestal.t_i_average.value[1]
    t_ratio = ti[1]/te[1]
    
    function cost_function(x)
        if par.fix_pedestal
            te_new = IMAS.Hmode_profiles(te_sep, te_ped, x[1], length(cp1d.grid.rho_tor_norm), par.T_shaping, par.T_shaping, wped)
            ti_new = IMAS.Hmode_profiles(ti_sep, ti_ped, x[1]*t_ratio, length(cp1d.grid.rho_tor_norm), par.T_shaping, par.T_shaping, wped)
        else
            te_new  = x[1] / te[1] * te
            ti_new = te_new * t_ratio
        end
        cp1d.electrons.temperature = te_new
        for ion in cp1d.ion
            ion.temperature = ti_new
        end
        βn_new = IMAS.beta_tor(dd.equilibrium.time_slice[], cp1d, norm=par.norm, thermal=par.thermal)
        return abs(par.βn_target - βn_new)
    end
    
    # run optimizer to find Te0 that gives target betaN
    res = Optim.optimize(x -> cost_function(x), 0, 1e6)
    #display(res)
    Te0 = res.minimizer[1]
    
    # if optimized Te0 is less than te_ped, increase to be equal to te_ped
    if par.fix_pedestal
        if Te0 < te_ped
            Te0 = te_ped
            te_new = IMAS.Hmode_profiles(te_sep, te_ped, Te0, length(cp1d.grid.rho_tor_norm), par.T_shaping, par.T_shaping, wped)
            ti_new = IMAS.Hmode_profiles(ti_sep, ti_ped, Te0*t_ratio, length(cp1d.grid.rho_tor_norm), par.T_shaping, par.T_shaping, wped)
            cp1d.electrons.temperature = te_new
            for ion in cp1d.ion
                ion.temperature = ti_new
            end
        end
    end

    return actor
end

"""
    _finalize(actor::ActorBetaMatch)

Updates IMAS.core_sources
"""
function _finalize(actor::ActorBetaMatch)
    dd = actor.dd

    # update sources
    IMAS.intrinsic_sources!(dd)

    return actor
end