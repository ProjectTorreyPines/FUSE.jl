#= ================== =#
#  ActorEPEDprofiles  #
#= ================== =#
Base.@kwdef mutable struct FUSEparameters__ActorEPEDprofiles{T<:Real} <: ParametersActor{T}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    _time::Float64 = NaN
    Te_shaping::Entry{T} = Entry{T}("-", "Shaping coefficient for the temperature profile")
    ne_shaping::Entry{T} = Entry{T}("-", "Shaping coefficient for the density profile")
    T_ratio_pedestal::Entry{T} = Entry{T}("-", "Ion to electron temperature ratio in the pedestal")
    T_ratio_core::Entry{T} = Entry{T}("-", "Ion to electron temperature ratio in the core")
end

mutable struct ActorEPEDprofiles{D,P} <: SingleAbstractActor{D,P}
    dd::IMAS.dd{D}
    par::FUSEparameters__ActorEPEDprofiles{P}
    function ActorEPEDprofiles(dd::IMAS.dd{D}, par::FUSEparameters__ActorEPEDprofiles{P}; kw...) where {D<:Real,P<:Real}
        par = par(kw...)
        return new{D,P}(dd, par)
    end
end

"""
    ActorEPEDprofiles(dd::IMAS.dd, act::ParametersAllActors; kw...)

Updates pedestal height and width and blends with core profiles that are defined by shaping factors.

Does not change on-axis values of plasma profiles.
"""
function ActorEPEDprofiles(dd::IMAS.dd, act::ParametersAllActors; kw...)
    actor = ActorEPEDprofiles(dd, act.ActorEPEDprofiles, act; kw...)
    step(actor)
    finalize(actor)
    return actor
end

"""
    _step(actor::ActorEPEDprofiles)
"""
function _step(actor::ActorEPEDprofiles)
    par = actor.par
    dd = actor.dd
    cp1d = dd.core_profiles.profiles_1d[]

    sol = run_EPED(dd; ne_from=:pulse_schedule, zeff_from=:pulse_schedule, βn_from=:equilibrium, ip_from=:pulse_schedule, only_powerlaw=false, warn_nn_train_bounds=false)
    pped = sol.pressure.GH.H
    wped = sol.width.GH.H

    rho_ped = IMAS.interp1d(cp1d.grid.psi_norm, cp1d.grid.rho_tor_norm).(1 - sol.width.GH.H)

    ne_ped = IMAS.get_from(dd, Val{:ne_ped}, :pulse_schedule, rho_ped)
    zeff_ped = IMAS.get_from(dd, Val{:zeff_ped}, :pulse_schedule, rho_ped)

    impurity = [ion.element[1].z_n for ion in cp1d.ion if Int(floor(ion.element[1].z_n)) != 1][1]
    zi = sum(impurity) / length(impurity)
    nival = ne_ped * (zeff_ped - 1) / (zi^2 - zi)
    nval = ne_ped - zi * nival
    nsum = ne_ped + nval + nival
    tped = (pped * 1e6) / nsum / IMAS.mks.e

    # update electron temperature profile using
    # * new pedestal height & width
    # * existing Te0 & Te_shaping 
    Te = cp1d.electrons.temperature
    Te_ped = 2.0 * tped / (1.0 + par.T_ratio_pedestal)
    tval = IMAS.Hmode_profiles(Te[end], Te_ped, Te[1], length(cp1d.grid.rho_tor_norm), par.Te_shaping, par.Te_shaping, wped)
    cp1d.electrons.temperature = tval
    if any(tval .< 0)
        error("Te profile is negative for T0=$(Te[1]) eV and Tped=$(Te_ped) eV")
    end

    # update ion temperature profiles using:
    # * new pedestal height & width
    # * existing Te0 & Te_shaping, and T_ratio_pedestal & T_ratio_core
    tval_ions = IMAS.Hmode_profiles(
        Te[end] * par.T_ratio_pedestal,
        Te_ped * par.T_ratio_pedestal,
        Te[1] * par.T_ratio_core,
        length(cp1d.grid.rho_tor_norm),
        par.Te_shaping,
        par.Te_shaping,
        wped)
    if any(tval_ions .< 0)
        error("Ti profile is negative for T0=$(Te[1]*par.T_ratio_core) eV and Tped=$(Te_ped*par.T_ratio_pedestal) eV")
    end
    for ion in cp1d.ion
        ion.temperature = tval_ions
    end

    # update electron density profile using
    # * new pedestal height & width
    # * existing ne0 & ne_shaping 
    ne = cp1d.electrons.density_thermal
    # first store ratios of electron density to ion densities
    ion_fractions = zeros(Float64, length(cp1d.ion), length(ne))
    for (ii, ion) in enumerate(cp1d.ion)
        ion_fractions[ii, :] = ion.density_thermal ./ ne
    end
    nval = IMAS.Hmode_profiles(ne[end], ne_ped, ne[1], length(cp1d.grid.rho_tor_norm), par.ne_shaping, par.ne_shaping, wped)
    cp1d.electrons.density_thermal = nval
    if any(nval .< 0)
        error("ne profile is negative for n0=$(ne[1]) /m^3 and Tped=$(ne_ped) /m^3")
    end

    # update ion density profiles using
    # * existing ratios of electron to ion densities
    for (ii, ion) in enumerate(cp1d.ion)
        ion.density_thermal = ion_fractions[ii, :] .* cp1d.electrons.density_thermal
    end

    return actor
end

"""
    _finalize(actor::ActorEPEDprofiles)

Updates IMAS.core_sources
"""
function _finalize(actor::ActorEPEDprofiles)
    dd = actor.dd

    # update sources
    IMAS.sources!(dd)

    return actor
end
