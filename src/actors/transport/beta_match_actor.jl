#= ================== =#
#  ActorBetaMatch  #
#= ================== =#
Base.@kwdef mutable struct FUSEparameters__ActorBetaMatch{T<:Real} <: ParametersActor{T}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    _time::Float64 = NaN
    βn_target::Entry{T} = Entry{T}("-", "Target value for normalized toroidal beta")
    T_shaping::Entry{T} = Entry{T}("-", "Shaping coefficient for the temperature profile")
end

mutable struct ActorBetaMatch{D,P} <: CompoundAbstractActor{D,P}
    dd::IMAS.dd{D}
    par::FUSEparameters__ActorBetaMatch{P}
end

"""
    ActorBetaMatch(dd::IMAS.dd, act::ParametersAllActors; kw...)

Updates core temperature profiles in order to match target betaN.
Uses Hmode_profiles() function to scale core temperature profiles while keeping pedestal unchanged.
Does not change density profiles.
Does not change pedestal parameters.

"""
function ActorBetaMatch(dd::IMAS.dd, act::ParametersAllActors; kw...)
    actor = ActorBetaMatch(dd, act.ActorBetaMatch, act; kw...)
    step(actor)
    finalize(actor)
    return actor
end

function ActorBetaMatch(dd::IMAS.dd, par::FUSEparameters__ActorBetaMatch, act::ParametersAllActors; kw...)
    par = par(kw...)
    return ActorBetaMatch(dd, par)
end

"""
    _step(actor::ActorBetaMatch)
"""
function _step(actor::ActorBetaMatch)
    par = actor.par
    dd = actor.dd
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
        te_new = IMAS.Hmode_profiles(te_sep, te_ped, x[1], length(cp1d.grid.rho_tor_norm), par.T_shaping, par.T_shaping, wped)
        ti_new = IMAS.Hmode_profiles(ti_sep, ti_ped, x[1]*t_ratio, length(cp1d.grid.rho_tor_norm), par.T_shaping, par.T_shaping, wped)
        cp1d.electrons.temperature = te_new
        for ion in cp1d.ion
            ion.temperature = ti_new
        end
        βn_new = IMAS.beta_tor_norm(dd.equilibrium, cp1d)
        return abs(par.βn_target - βn_new)
    end
    
    # run optimizer to find Te0 that gives target betaN
    res = Optim.optimize(x -> cost_function(x), 0, 1e6)
    #display(res)
    Te0 = res.minimizer[1]
    
    # if optimized Te0 is less than te_ped, increase to be equal to te_ped
    if Te0 < te_ped
        Te0 = te_ped
        te_new = IMAS.Hmode_profiles(te_sep, te_ped, Te0, length(cp1d.grid.rho_tor_norm), par.T_shaping, par.T_shaping, wped)
        ti_new = IMAS.Hmode_profiles(ti_sep, ti_ped, Te0*t_ratio, length(cp1d.grid.rho_tor_norm), par.T_shaping, par.T_shaping, wped)
        cp1d.electrons.temperature = te_new
        for ion in cp1d.ion
            ion.temperature = ti_new
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
    IMAS.sources!(dd)

    return actor
end
