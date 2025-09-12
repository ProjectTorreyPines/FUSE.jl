#= ================== =#
#  ActorEPEDprofiles  #
#= ================== =#
@actor_parameters_struct ActorEPEDprofiles{T} begin
    Te_shaping::Entry{T} = Entry{T}("-", "Shaping coefficient for the temperature profile")
    ne_shaping::Entry{T} = Entry{T}("-", "Shaping coefficient for the density profile")
    T_ratio_pedestal::Entry{T} = Entry{T}("-", "Ion to electron temperature ratio in the pedestal")
    T_ratio_core::Entry{T} = Entry{T}("-", "Ion to electron temperature ratio in the core")
end

mutable struct ActorEPEDprofiles{D,P} <: CompoundAbstractActor{D,P}
    dd::IMAS.dd{D}
    par::OverrideParameters{P,FUSEparameters__ActorEPEDprofiles{P}}
    act::ParametersAllActors{P}
    function ActorEPEDprofiles(dd::IMAS.dd{D}, par::FUSEparameters__ActorEPEDprofiles{P}, act::ParametersAllActors{P}; kw...) where {D<:Real,P<:Real}
        logging_actor_init(ActorEPEDprofiles)
        par = OverrideParameters(par; kw...)
        return new{D,P}(dd, par, act)
    end
end

"""
    ActorEPEDprofiles(dd::IMAS.dd, act::ParametersAllActors; kw...)

Constructs complete plasma profiles by blending EPED pedestal predictions with shaped core profiles.

This actor combines EPED-predicted pedestal conditions with parametrically shaped core 
profiles to create complete, consistent temperature and density profiles across the entire 
plasma. It maintains fixed on-axis values while ensuring smooth transitions between core 
and pedestal regions.

Key operations:
- Runs EPED model to predict pedestal pressure and width
- Applies H-mode profile functions with configurable shaping parameters
- Maintains ion-electron temperature ratio consistency (core vs. pedestal)
- Preserves particle density ratios between species
- Ensures proper pedestal positioning at ρ=0.9 reference location

Profile construction workflow:
1. **EPED Prediction**: Calculates pedestal pressure and width from current plasma state
2. **Temperature Profiles**: Constructs Te and Ti using H-mode functions with shaping
3. **Density Profiles**: Updates ne and ni maintaining species fraction consistency  
4. **Source Updates**: Recalculates all power and particle sources for updated profiles

Shaping parameters:
- `Te_shaping`: Controls electron temperature profile curvature
- `ne_shaping`: Controls electron density profile curvature  
- `T_ratio_pedestal`: Ti/Te ratio in pedestal region
- `T_ratio_core`: Ti/Te ratio in core region

The actor preserves the on-axis plasma values while optimizing the profile shapes 
to match both EPED pedestal predictions and prescribed core-pedestal transitions.

!!! note

    Does not modify on-axis values of plasma profiles, only shapes and pedestal conditions
"""
function ActorEPEDprofiles(dd::IMAS.dd, act::ParametersAllActors; kw...)
    actor = ActorEPEDprofiles(dd, act.ActorEPEDprofiles, act; kw...)
    step(actor)
    finalize(actor)
    return actor
end

"""
    _step(actor::ActorEPEDprofiles)

Executes EPED model and constructs shaped plasma profiles.

The step function runs the EPED prediction, extracts pedestal boundary conditions,
and applies H-mode profile functions to create consistent temperature and density 
profiles that smoothly transition from shaped cores to EPED-predicted pedestals.
"""
function _step(actor::ActorEPEDprofiles)
    dd = actor.dd
    par = actor.par
    act = actor.act

    cp1d = dd.core_profiles.profiles_1d[]

    sol = run_EPED(dd; ne_from=:pulse_schedule, zeff_from=:pulse_schedule, βn_from=:equilibrium, ip_from=:pulse_schedule, act.ActorEPED.only_powerlaw, act.ActorEPED.warn_nn_train_bounds)
    pped = sol.pressure.GH.H
    wped = sol.width.GH.H

    rho_ped = IMAS.interp1d(cp1d.grid.psi_norm, cp1d.grid.rho_tor_norm).(1 - sol.width.GH.H)

    ne_ped = IMAS.get_from(dd, Val(:ne_ped), :pulse_schedule, rho_ped)
    zeff_ped = IMAS.get_from(dd, Val(:zeff_ped), :pulse_schedule, rho_ped)

    impurity = [ion.element[1].z_n for ion in cp1d.ion if !IMAS.is_hydrogenic(ion)][1]
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
        error("Te profile is negative for T0=$(Te[1]) [eV] and Tped=$(Te_ped) [eV]")
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
        error("Ti profile is negative for T0=$(Te[1]*par.T_ratio_core) [eV] and Tped=$(Te_ped*par.T_ratio_pedestal) [eV]")
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
    nval = IMAS.Hmode_profiles(ne[end], ne_ped, ne[1], length(cp1d.grid.rho_tor_norm), par.ne_shaping, par.ne_shaping, 0.05)
    cp1d.electrons.density_thermal = IMAS.ped_height_at_09(cp1d.grid.rho_tor_norm, nval, ne_ped)
    if any(nval .< 0)
        error("ne profile is negative for n0=$(ne[1]) [m⁻³] and ne_ped=$(ne_ped) [m⁻³]")
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

Updates plasma sources to be consistent with the new profile shapes.

Recalculates all power, particle, and momentum sources in dd.core_sources 
to ensure consistency with the updated temperature and density profiles.
"""
function _finalize(actor::ActorEPEDprofiles)
    dd = actor.dd

    # update sources
    IMAS.sources!(dd)

    return actor
end
