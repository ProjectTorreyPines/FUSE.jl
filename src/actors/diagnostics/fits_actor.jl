#= ================ =#
#  ActorFitProfiles  #
#= ================ =#
Base.@kwdef mutable struct FUSEparameters__ActorFitProfiles{T<:Real} <: ParametersActor{T}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    _time::Float64 = NaN
    #== actor parameters ==#
    time_averaging::Entry{Float64} = Entry{Float64}("s", "Time averaging window")
    rho_averaging::Entry{Float64} = Entry{Float64}("", "rho averaging window")
    rho_grid::Entry{Int} = Entry{Int}("-", "Number of points in rho"; default=101)
    time_basis_ids::Switch{Symbol} = Switch{Symbol}([:equilibrium, :core_profiles], "-", "Time basis to use"; default=:core_profiles)
    #== display and debugging parameters ==#
end

mutable struct ActorFitProfiles{D,P} <: CompoundAbstractActor{D,P}
    dd::IMAS.dd{D}
    par::OverrideParameters{P,FUSEparameters__ActorFitProfiles{P}}
    act::ParametersAllActors{P}
    function ActorFitProfiles(dd::IMAS.dd{D}, par::FUSEparameters__ActorFitProfiles{P}, act::ParametersAllActors{P}; kw...) where {D<:Real,P<:Real}
        logging_actor_init(ActorFitProfiles)
        par = OverrideParameters(par; kw...)
        return new{D,P}(dd, par, act)
    end
end

"""
    ActorFitProfiles(dd::IMAS.dd, act::ParametersAllActors; kw...)

FitProfiles experimental data
"""
function ActorFitProfiles(dd::IMAS.dd, act::ParametersAllActors; kw...)
    actor = ActorFitProfiles(dd, act.ActorFitProfiles, act; kw...)
    step(actor)
    finalize(actor)
    return actor
end

function _step(actor::ActorFitProfiles{D,P}) where {D<:Real,P<:Real}
    dd = actor.dd
    par = actor.par

    # parameters and switches
    time_basis = getproperty(dd, par.time_basis_ids).time
    rho_tor_norm = range(0.0, 1.0, par.rho_grid)
    rho_tor_norm12 = 0.05:(rho_tor_norm[2]-rho_tor_norm[1]):1.2
    smooth1 = 1.0
    smooth2 = par.rho_averaging

    # set aside original raw data, since we'll be overwriting it internally 
    dd1 = IMAS.dd{D}()
    for field in (:thomson_scattering, :charge_exchange)
        setproperty!(dd1, field, deepcopy(getproperty(dd, field)))
    end

    # identify outliers on raw data
    for experimental_ids in (dd.thomson_scattering, dd.charge_exchange)
        tg = IMAS.time_groups(experimental_ids; min_channels=5)
        IMAS.adaptive_outlier_removal!(tg)
    end

    # use the same time-basis
    for experimental_ids in (dd.thomson_scattering, dd.charge_exchange)
        tg = IMAS.time_dependent_leaves(experimental_ids)
        times_coords = Dict{IMAS.IDS,IMAS.Coordinate}()
        for group in values(tg)
            for leaf in group
                data = IMAS.smooth_by_convolution(leaf.ids, leaf.field, time_basis; window_size=par.time_averaging)
                setproperty!(leaf.ids, leaf.field, data)
                time_coord = IMAS.time_coordinate(leaf.ids, leaf.field)
                times_coords[time_coord.ids] = time_coord
            end
        end
        for time_coord in values(times_coords)
            setproperty!(time_coord.ids, time_coord.field, time_basis)
        end
    end

    # empty core_profiles
    empty!(dd.core_profiles)
    for time0 in time_basis
        cp1d = resize!(dd.core_profiles.profiles_1d, time0)
        cp1d.grid.rho_tor_norm = rho_tor_norm
        bulk_ion, imp_ion = resize!(cp1d.ion, 2)
        IMAS.ion_element!(bulk_ion, :D)
        IMAS.ion_element!(imp_ion, :C)
    end
    dd.core_profiles.time = time_basis

    # get space-time dependent data and return interpolators
    itp_te = IMAS.fit2d(Val(:t_e), dd; transform=sqrt)
    itp_ne = IMAS.fit2d(Val(:n_e), dd; transform=sqrt)
    itp_zeff = IMAS.fit2d(Val(:zeff), dd; transform=x -> sqrt(max(x, 1.0) - 1.0))
    itp_nimp = IMAS.fit2d(Val(:n_i_over_n_e), dd; transform=sqrt)
    itp_ti = IMAS.fit2d(Val(:t_i), dd; transform=sqrt)

    # fit Te
    for (k, time0) in enumerate(time_basis)
        cp1d = dd.core_profiles.profiles_1d[k]
        data = itp_te(rho_tor_norm12, range(time0, time0, length(rho_tor_norm12))) .^ 2
        cp1d.electrons.temperature = IMAS.fit1d(rho_tor_norm12, data, rho_tor_norm; smooth1, smooth2).fit
    end

    # fit ne
    for (k, time0) in enumerate(time_basis)
        cp1d = dd.core_profiles.profiles_1d[k]
        data = itp_ne(rho_tor_norm12, range(time0, time0, length(rho_tor_norm12))) .^ 2
        cp1d.electrons.density_thermal = IMAS.fit1d(rho_tor_norm12, data, rho_tor_norm; smooth1, smooth2).fit
    end

    # fit Zeff
    for (k, time0) in enumerate(time_basis)
        cp1d = dd.core_profiles.profiles_1d[k]
        data = itp_zeff(rho_tor_norm12, range(time0, time0, length(rho_tor_norm12))) .^ 2.0 .+ 1.0
        cp1d.zeff = IMAS.fit1d(rho_tor_norm12, data, rho_tor_norm; smooth1, smooth2).fit
    end

    # fit ni
    for (k, time0) in enumerate(time_basis)
        cp1d = dd.core_profiles.profiles_1d[k]
        bulk_ion = cp1d.ion[1]
        imp_ion = cp1d.ion[2]

        data = itp_nimp(rho_tor_norm12, range(time0, time0, length(rho_tor_norm12))) .^ 2
        data .*= itp_ne(rho_tor_norm12, range(time0, time0, length(rho_tor_norm12))) .^ 2
        n_i = IMAS.fit1d(rho_tor_norm12, data, rho_tor_norm; smooth1, smooth2).fit

        bulk_ion.density_thermal = zero(rho_tor_norm)
        imp_ion.density_thermal = n_i
    end

    # fit Ti
    for (k, time0) in enumerate(time_basis)
        cp1d = dd.core_profiles.profiles_1d[k]
        data = itp_ti(rho_tor_norm12, range(time0, time0, length(rho_tor_norm12))) .^ 2
        ti = IMAS.fit1d(rho_tor_norm12, data, rho_tor_norm; smooth1, smooth2).fit
        for ion in cp1d.ion
            ion.temperature = ti
        end
    end

    # quasi neutrality
    for (k, time0) in enumerate(time_basis)
        cp1d = dd.core_profiles.profiles_1d[k]
        IMAS.enforce_quasi_neutrality!(cp1d, :D)
    end

    # # rotation
    # index = [IMAS.hasdata(cp1d, :rotation_frequency_tor_sonic) for cp1d in dd1.core_profiles.profiles_1d]
    # min_k_orig = findfirst(index)
    # if min_k_orig !== nothing
    #     for (k, time0) in enumerate(time_basis)
    #         k_orig = max(min_k_orig, IMAS.nearest_causal_time(dd1.core_profiles.time, time0; bounds_error=false).index)
    #         dd.core_profiles.profiles_1d[k].rotation_frequency_tor_sonic = dd1.core_profiles.profiles_1d[k_orig].rotation_frequency_tor_sonic
    #     end
    # end

    # restore the original raw data
    for field in (:thomson_scattering, :charge_exchange)
        setproperty!(dd, field, getproperty(dd1, field))
    end

    return actor
end
