import TAUENN

#= ===================== =#
#  ActorFixedProfiles     #
#= ===================== =#
Base.@kwdef mutable struct FUSEparameters__ActorFixedProfiles{T} <: ParametersActor where {T<:Real}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    T_shaping::Entry{T} = Entry{T}("-", "Shaping coefficient for the temperature profile"; default=1.8)
    n_shaping::Entry{T} = Entry{T}("-", "Shaping coefficient for the density profile"; default=1.8)
    T_ratio_pedestal::Entry{T} = Entry{T}("-", "Ion to electron temperature ratio in the pedestal"; default=1.0)
    T_ratio_core::Entry{T} = Entry{T}("-", "Ion to electron temperature ratio in the core"; default=1.0)
    update_pedestal::Entry{Bool} = Entry{Bool}("-","update pedestal with eped_nn" ;default=true)
    verbose::Entry{Bool} = Entry{Bool}("-", "Verbose"; default=false)
end

mutable struct ActorFixedProfiles{D,P} <: PlasmaAbstractActor
    dd::IMAS.dd{D}
    par::FUSEparameters__ActorFixedProfiles{P}
    ped_actor::ActorPedestal{D,P}
end

"""
    ActorFixedProfiles(dd::IMAS.dd, act::ParametersAllActors; kw...)

Stand-in actor for Fixed Profile method in EquilibriumTransport workflow.
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
    ped_actor = ActorPedestal(dd, act.ActorPedestal)
    return ActorFixedProfiles(dd, par, ped_actor)
end

"""
    step(actor::ActorFixedProfiles)

"""
function _step(actor::ActorFixedProfiles)
    par = actor.par
    dd = actor.dd

    # update pedestal with EPED-NN
    if par.update_pedestal

        # run EPED-NN
        finalize(step(actor.ped_actor))

        # update electron temperature profile
        # uses new pedestal height & width, existing Te0 & T_shaping 
        cp1d = dd.core_profiles.profiles_1d[]
        Te = cp1d.electrons.temperature
        tval = IMAS.Hmode_profiles(Te[end], dd.summary.local.pedestal.t_e.value[], Te[1], length(cp1d.grid.rho_tor_norm), par.T_shaping, par.T_shaping, 1 - dd.summary.local.pedestal.position.rho_tor_norm[])
    
        cp1d.electrons.temperature = tval
        if any(tval .< 0)
            throw("Te profile is negative for T0=$(Te[1]) eV and Tped=$(dd.summary.local.pedestal.t_e.value[]) eV")
        end
    
        # update ion temperature profiles
        # uses new pedestal height & width, existing Te0 & T_shaping, and T_ratio_pedestal & T_ratio_core
        tval_ions = IMAS.Hmode_profiles(Te[end]*par.T_ratio_pedestal, dd.summary.local.pedestal.t_e.value[]*par.T_ratio_pedestal, Te[1]*par.T_ratio_core, length(cp1d.grid.rho_tor_norm), par.T_shaping, par.T_shaping, 1 - dd.summary.local.pedestal.position.rho_tor_norm[])
        if any(tval_ions .< 0)
            throw("Ti profile is negative for T0=$(Te[1]*par.T_ratio_core) eV and Tped=$(dd.summary.local.pedestal.t_e.value[]*par.T_ratio_pedestal) eV")
        end
        for ion in cp1d.ion
            ion.temperature = tval_ions
        end

        # update electron density profile
        # uses new pedestal height & width, existing ne0 & n_shaping 

        # first store ratios of electron density to ion densities
        ion_fractions = zeros(Float64, length(cp1d.ion), length(cp1d.electrons.density))
        for (ii, ion) in enumerate(cp1d.ion)
            ion_fractions[ii, :] = ion.density_thermal ./ cp1d.electrons.density
        end

        ne = cp1d.electrons.density
        nval = IMAS.Hmode_profiles(ne[end], dd.summary.local.pedestal.n_e.value[], ne[1], length(cp1d.grid.rho_tor_norm), par.n_shaping, par.n_shaping, 1 - dd.summary.local.pedestal.position.rho_tor_norm[])

        cp1d.electrons.density = nval
        if any(nval .< 0)
            throw("ne profile is negative for n0=$(ne[1]) /m^3 and Tped=$(dd.summary.local.pedestal.n_e.value[]) /m^3")
        end

        # update ion density profiles
        # uses stored ratios of electron density to ion densities
        for (ii, ion) in enumerate(cp1d.ion)
            ion.density_thermal = ion_fractions[ii, :] .* cp1d.electrons.density
        end

    end


    return actor
end

"""
    finalize(actor::ActorFixedProfiles)

Updates IMAS.sources
"""
function _finalize(actor::ActorFixedProfiles)
    dd = actor.dd

    # update sources

    IMAS.sources!(dd)

    return actor
end
