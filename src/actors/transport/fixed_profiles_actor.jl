#= ================== =#
#  ActorFixedProfiles  #
#= ================== =#
Base.@kwdef mutable struct FUSEparameters__ActorFixedProfiles{T} <: ParametersActor where {T<:Real}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    T_shaping::Entry{T} = Entry{T}("-", "Shaping coefficient for the temperature profile"; default=1.8)
    n_shaping::Entry{T} = Entry{T}("-", "Shaping coefficient for the density profile"; default=1.8)
    T_ratio_pedestal::Entry{T} = Entry{T}("-", "Ion to electron temperature ratio in the pedestal"; default=1.0)
    T_ratio_core::Entry{T} = Entry{T}("-", "Ion to electron temperature ratio in the core"; default=1.0)
    n_ITB_height::Entry{T} = Entry{T}("-", "Density ITB height as a ratio of core-ped different"; default=0.0)
    T_ITB_height::Entry{T} = Entry{T}("-", "Temperature ITB height as a ratio of core-ped different"; default=0.0)
    ITB_width::Entry{T} = Entry{T}("rho", "ITB width"; default=0.1)
    ITB_radius::Entry{T} = Entry{T}("rho", "ITB radial location"; default=0.6)
    update_pedestal::Entry{Bool} = Entry{Bool}("-", "Update pedestal with EPED-NN and modify profiles accordingly"; default=true)
end

mutable struct ActorFixedProfiles{D,P} <: PlasmaAbstractActor{D,P}
    dd::IMAS.dd{D}
    par::FUSEparameters__ActorFixedProfiles{P}
    ped_actor::ActorPedestal{D,P}
end

"""
    ActorFixedProfiles(dd::IMAS.dd, act::ParametersAllActors; kw...)

Updates pedestal height and width and blends with core profiles that are defined by shaping factors.

Does not change on-axis values of plasma profiles.
"""
function ActorFixedProfiles(dd::IMAS.dd, act::ParametersAllActors; kw...)
    actor = ActorFixedProfiles(dd, act.ActorFixedProfiles, act; kw...)
    step(actor)
    finalize(actor)
    return actor
end

function ActorFixedProfiles(dd::IMAS.dd, par::FUSEparameters__ActorFixedProfiles, act::ParametersAllActors; kw...)
    par = par(kw...)
    ped_actor = ActorPedestal(dd, act.ActorPedestal; update_core_profiles=false, par.T_ratio_pedestal, ip_from=:core_profiles, βn_from=:equilibrium)
    return ActorFixedProfiles(dd, par, ped_actor)
end

"""
    _step(actor::ActorFixedProfiles)
"""
function _step(actor::ActorFixedProfiles)
    par = actor.par
    dd = actor.dd
    cp1d = dd.core_profiles.profiles_1d[]

    # update pedestal
    if par.update_pedestal
        finalize(step(actor.ped_actor))
    end

    # update electron temperature profile using
    # * new pedestal height & width
    # * existing Te0 & T_shaping 
    # * ITB parameters 

    Te = cp1d.electrons.temperature
    Te_ped = @ddtime(dd.summary.local.pedestal.t_e.value)
    w_ped = @ddtime(dd.summary.local.pedestal.position.rho_tor_norm)

    Te_ITB_height = (Te[1] - Te_ped) * par.T_ITB_height

    tval = IMAS.ITB_profiles(
        Te[end], 
        Te_ped, 
        Te[1], 
        length(cp1d.grid.rho_tor_norm), 
        par.T_shaping, 
        par.T_shaping, 
        1.0 - w_ped,
        par.ITB_radius,
        par.ITB_width,
        Te_ITB_height,
    )
    cp1d.electrons.temperature = tval
    if any(tval .< 0)
        throw("Te profile is negative for T0=$(Te[1]) eV and Tped=$(Te_ped) eV")
    end

    # update ion temperature profiles using:
    # * new pedestal height & width
    # * existing Te0 & T_shaping, and T_ratio_pedestal & T_ratio_core
    # * ITB parameters 

    Ti_ITB_height = (Te[1] * par.T_ratio_core - Te_ped * par.T_ratio_pedestal) * par.T_ITB_height

    tval_ions = IMAS.ITB_profiles(
        Te[end] * par.T_ratio_pedestal,
        Te_ped * par.T_ratio_pedestal,
        Te[1] * par.T_ratio_core,
        length(cp1d.grid.rho_tor_norm),
        par.T_shaping,
        par.T_shaping,
        1.0 - w_ped,
        par.ITB_radius,
        par.ITB_width,
        Ti_ITB_height,
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
    ne_ped = @ddtime(dd.summary.local.pedestal.n_e.value)
    # first store ratios of electron density to ion densities
    ion_fractions = zeros(Float64, length(cp1d.ion), length(ne))
    for (ii, ion) in enumerate(cp1d.ion)
        ion_fractions[ii, :] = ion.density_thermal ./ ne
    end

    ne_ITB_height = (ne[1] - ne_ped) * par.n_ITB_height

    nval = IMAS.ITB_profiles(
        ne[end], 
        ne_ped, 
        ne[1], 
        length(cp1d.grid.rho_tor_norm), 
        par.n_shaping, 
        par.n_shaping, 
        1.0 - w_ped,
        par.ITB_radius,
        par.ITB_width,
        ne_ITB_height,
        )
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
    _finalize(actor::ActorFixedProfiles)

Updates IMAS.core_sources
"""
function _finalize(actor::ActorFixedProfiles)
    dd = actor.dd

    # update sources
    IMAS.sources!(dd)

    return actor
end
