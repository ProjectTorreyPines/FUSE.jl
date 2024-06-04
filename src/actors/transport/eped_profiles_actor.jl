#= ================== =#
#  ActorEPEDProfiles  #
#= ================== =#
Base.@kwdef mutable struct FUSEparameters__ActorEPEDProfiles{T<:Real} <: ParametersActorPlasma{T}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    _time::Float64 = NaN
    T_shaping::Entry{T} = Entry{T}("-", "Shaping coefficient for the temperature profile"; default=1.8)
    n_shaping::Entry{T} = Entry{T}("-", "Shaping coefficient for the density profile"; default=1.8)
    T_ratio_pedestal::Entry{T} = Entry{T}("-", "Ion to electron temperature ratio in the pedestal"; default=1.0)
    T_ratio_core::Entry{T} = Entry{T}("-", "Ion to electron temperature ratio in the core"; default=1.0)
end

mutable struct ActorEPEDProfiles{D,P} <: CompoundAbstractActor{D,P}
    dd::IMAS.dd{D}
    par::FUSEparameters__ActorEPEDProfiles{P}
end

"""
    ActorEPEDProfiles(dd::IMAS.dd, act::ParametersAllActors; kw...)

Updates pedestal height and width and blends with core profiles that are defined by shaping factors.

Does not change on-axis values of plasma profiles.
"""
function ActorEPEDProfiles(dd::IMAS.dd, act::ParametersAllActors; kw...)
    actor = ActorEPEDProfiles(dd, act.ActorEPEDProfiles, act; kw...)
    step(actor)
    finalize(actor)
    return actor
end

function ActorEPEDProfiles(dd::IMAS.dd, par::FUSEparameters__ActorEPEDProfiles, act::ParametersAllActors; kw...)
    par = par(kw...)
    return ActorEPEDProfiles(dd, par)
end

"""
    _step(actor::ActorEPEDProfiles)
"""
function _step(actor::ActorEPEDProfiles)
    par = actor.par
    dd = actor.dd

    eqt = dd.equidistant.time_slice[]
    cp1d = dd.core_profiles.profiles_1d[]

    sol = run_EPED(eqt, cp1d; ne_ped_from=:pulse_schedule, zeff_ped_from=:pulse_schedule, βn_from=:equilibrium,  ip_from=:pulse_schedule, only_powerlaw=false, warn_nn_train_bounds=false)

    ne_ped = IMAS.get_from(dd, Val{:ne_ped}, :pulse_schedule)
    zeff_ped =  IMAS.get_from(dd, Val{:zeff_ped}, :pulse_schedule)

    pped = sol.pressure.GH.H
    w_ped = IMAS.interp1d(cp1d.grid.psi_norm, cp1d.grid.rho_tor_norm).(1 - sol.width.GH.H)

    impurity = [ion.element[1].z_n for ion in cp1d.ion if Int(floor(ion.element[1].z_n)) != 1][1]
    zi = sum(impurity) / length(impurity)
    nival =ne_ped  * (zeff_ped - 1) / (zi^2 - zi)
    nval = ne_ped  - zi * nival
    nsum = ne_ped  + nval + nival
    tped = (pped * 1e6) / nsum / constants.e

    # update electron temperature profile using
    # * new pedestal height & width
    # * existing Te0 & T_shaping 
    Te = cp1d.electrons.temperature
    Te_ped = 2.0 * tped / (1.0 + par.T_ratio_pedestal)
    w_ped = @ddtime(dd.summary.local.pedestal.position.rho_tor_norm)
    tval = IMAS.Hmode_profiles(Te[end], Te_ped, Te[1], length(cp1d.grid.rho_tor_norm), par.T_shaping, par.T_shaping, 1.0 - w_ped)
    cp1d.electrons.temperature = tval
    if any(tval .< 0)
        throw("Te profile is negative for T0=$(Te[1]) eV and Tped=$(Te_ped) eV")
    end

    # update ion temperature profiles using:
    # * new pedestal height & width
    # * existing Te0 & T_shaping, and T_ratio_pedestal & T_ratio_core
    tval_ions = IMAS.Hmode_profiles(
        Te[end] * par.T_ratio_pedestal,
        Te_ped * par.T_ratio_pedestal,
        Te[1] * par.T_ratio_core,
        length(cp1d.grid.rho_tor_norm),
        par.T_shaping,
        par.T_shaping,
        1.0 - w_ped
    )
    if any(tval_ions .< 0)
        throw("Ti profile is negative for T0=$(Te[1]*par.T_ratio_core) eV and Tped=$(Te_ped*par.T_ratio_pedestal) eV")
    end
    for ion in cp1d.ion
        ion.temperature = tval_ions
    end

    # update electron density profile using
    # * new pedestal height & width
    # * existing ne0 & n_shaping 
    ne = cp1d.electrons.density_thermal
    # first store ratios of electron density to ion densities
    ion_fractions = zeros(Float64, length(cp1d.ion), length(ne))
    for (ii, ion) in enumerate(cp1d.ion)
        ion_fractions[ii, :] = ion.density_thermal ./ ne
    end
    nval = IMAS.Hmode_profiles(ne[end], ne_ped, ne[1], length(cp1d.grid.rho_tor_norm), par.n_shaping, par.n_shaping, 1.0 - w_ped)
    cp1d.electrons.density_thermal = nval
    if any(nval .< 0)
        throw("ne profile is negative for n0=$(ne[1]) /m^3 and Tped=$(ne_ped) /m^3")
    end

    # update ion density profiles using
    # * existing ratios of electron to ion densities
    for (ii, ion) in enumerate(cp1d.ion)
        ion.density_thermal = ion_fractions[ii, :] .* cp1d.electrons.density_thermal
    end

    return actor
end

"""
    _finalize(actor::ActorEPEDProfiles)

Updates IMAS.core_sources
"""
function _finalize(actor::ActorEPEDProfiles)
    dd = actor.dd

    # update sources
    IMAS.sources!(dd)

    return actor
end
