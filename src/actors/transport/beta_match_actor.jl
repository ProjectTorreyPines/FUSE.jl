#= ================== =#
#  ActorBetaMatch  #
#= ================== =#
Base.@kwdef mutable struct FUSEparameters__ActorBetaMatch{T<:Real} <: ParametersActor{T}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    _time::Float64 = NaN
    βn_target::Entry{T} = Entry{T}("-", "Target value for normalized toroidal beta")
    Te_shaping::Entry{T} = Entry{T}("-", "Shaping coefficient for the temperature profile"; default=1.8)
    ne_shaping::Entry{T} = Entry{T}("-", "Shaping coefficient for the density profile"; default=1.8)
    T_ratio_pedestal::Entry{T} = Entry{T}("-", "Ion to electron temperature ratio in the pedestal"; default=1.0)
    T_ratio_core::Entry{T} = Entry{T}("-", "Ion to electron temperature ratio in the core"; default=1.0)
    do_plot::Entry{Bool} = Entry{Bool}("-", "Display optimization details and plot equilibrium, current profiles"; default=false)
end

mutable struct ActorBetaMatch{D,P} <: SingleAbstractActor{D,P}
    dd::IMAS.dd{D}
    par::FUSEparameters__ActorBetaMatch{P}
end

"""
    ActorBetaMatch(dd::IMAS.dd, act::ParametersAllActors; kw...)

Updates core temperature profiles in order to match target betaN.
Uses Hmode_profiles() function to scale core temperature profiles.

"""
function ActorBetaMatch(dd::IMAS.dd, act::ParametersAllActors; kw...)
    actor = ActorBetaMatch(dd, act.ActorBetaMatch, act; kw...)
    step(actor)
    finalize(actor)
    return actor
end

function ActorBetaMatch(dd::IMAS.dd, par::FUSEparameters__ActorBetaMatch, act::ParametersAllActors; kw...)
    logging_actor_init(ActorBetaMatch)
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

    # run ActorPedestal to update pedestal (without blend_core_edge)
    #finalize(step(actor.actor_ped))

    # get profile parameters
    ne = cp1d.electrons.density_thermal
    te = cp1d.electrons.temperature
    ti = cp1d.ion[1].temperature
    wped = 1 - dd.summary.local.pedestal.position.rho_tor_norm[1]
    ne_ped = dd.summary.local.pedestal.n_e.value[1]
    te_ped = dd.summary.local.pedestal.t_e.value[1]

    # update electron density profile with new pedestal width
    ne_new = IMAS.Hmode_profiles(ne[end], ne_ped, ne[1], length(cp1d.grid.rho_tor_norm), par.ne_shaping, par.ne_shaping, wped)
    cp1d.electrons.density_thermal = ne_new
    if any(ne_new .< 0)
        error("ne profile is negative for n0=$(ne[1]) /m^3 and Tped=$(ne_ped) /m^3")
    end
    # update ion density profiles using
    # * existing ratios of electron to ion densities
    ion_fractions = zeros(Float64, length(cp1d.ion), length(ne))
    for (ii, ion) in enumerate(cp1d.ion)
        ion_fractions[ii, :] = ion.density_thermal ./ ne
        ion.density_thermal = ion_fractions[ii, :] .* cp1d.electrons.density_thermal
    end

    function cost_function(x)
        te_new = IMAS.Hmode_profiles(te[end], te_ped, x[1], length(cp1d.grid.rho_tor_norm), par.Te_shaping, par.Te_shaping, wped)
        ti_new = IMAS.Hmode_profiles(ti[end], te_ped*par.T_ratio_pedestal, x[1]*par.T_ratio_core, length(cp1d.grid.rho_tor_norm), par.Te_shaping, par.Te_shaping, wped)
        cp1d.electrons.temperature = te_new
        for ion in cp1d.ion
            ion.temperature = ti_new
        end
        # update sources
        #IMAS.sources!(dd)
        # calc new βn
        βn_new = IMAS.beta_tor_norm(dd.equilibrium, cp1d)
        return abs(par.βn_target - βn_new) / par.βn_target
    end
    
    # run optimizer to find Te0 that gives target betaN
    res = Optim.optimize(x -> cost_function(x), 0, 1e6)
    Te0 = res.minimizer[1]
    
    # if optimized Te0 is less than te_ped, increase to be equal to te_ped
    if Te0 < te_ped
        display("Te0 ($(Te0) eV) is less than te_ped ($(te_ped) eV)")
        Te0 = te_ped
        te_new = IMAS.Hmode_profiles(te[end], te_ped, Te0, length(cp1d.grid.rho_tor_norm), par.Te_shaping, par.Te_shaping, wped)
        ti_new = IMAS.Hmode_profiles(ti[end], te_ped*par.T_ratio_pedestal, Te0*par.T_ratio_core, length(cp1d.grid.rho_tor_norm), par.Te_shaping, par.Te_shaping, wped)
        cp1d.electrons.temperature = te_new
        for ion in cp1d.ion
            ion.temperature = ti_new
        end
    end

    # plot results
    if par.do_plot
        pressure = cp1d.pressure
        volume = cp1d.grid.volume
        pressure_avg = IMAS.trapz(volume, pressure) / volume[end]
        display("<P> = $(pressure_avg)")
        display("Ip = $(IMAS.Ip(cp1d))")
        display("βn = $(IMAS.beta_tor_norm(dd.equilibrium, cp1d))")
        display(res)
        display(plot(dd.core_profiles))
        display(plot(dd.equilibrium))
        display(plot(dd.core_sources; only=5))
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
