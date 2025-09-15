#= ================ =#
#  ActorFitProfiles  #
#= ================ =#
Base.@kwdef mutable struct FUSEparameters__ActorFitProfiles{T<:Real} <: ParametersActor{T}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    _time::Float64 = NaN
    #== actor parameters ==#
    time_averaging::Entry{Float64} = Entry{Float64}("s", "Time averaging window")
    rho_averaging::Entry{Float64} = Entry{Float64}("-", "rho averaging window")
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
    rho_tor_norm12 = 0.1:(rho_tor_norm[2]-rho_tor_norm[1]):1.2
    smooth1 = 0.1
    smooth2 = par.rho_averaging

    # # set aside original raw data, since we'll be overwriting it internally 
    # dd1 = IMAS.dd{D}()
    # for field in (:thomson_scattering, :charge_exchange)
    #     setproperty!(dd1, field, deepcopy(getproperty(dd, field)))
    # end

    # identify outliers on raw data
    for experimental_ids in (dd.thomson_scattering, dd.charge_exchange)
        tg = IMAS.time_groups(experimental_ids; min_channels=0)
        IMAS.adaptive_outlier_removal!(tg; min_channels=5)
    end

    # disregard divertor thomson (D3D specific! need to generalize)
    for (k,ch) in reverse!(collect(enumerate(dd.thomson_scattering.channel)))
        if contains(lowercase(ch.name), "divertor")
            popat!(dd.thomson_scattering.channel, k)
        end
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

    # fit Te
    itp_te = IMAS.fit2d(Val(:t_e), dd; transform=abs)
    for (kt, time0) in enumerate(time_basis)
        cp1d = dd.core_profiles.profiles_1d[kt]
        data = itp_te(rho_tor_norm12, range(time0, time0, length(rho_tor_norm12)))
        cp1d.electrons.temperature = IMAS.fit1d(rho_tor_norm12, data, rho_tor_norm; smooth1, smooth2).fit
    end

    # scale thomson scattering density based on interferometer measurements
    if !isempty(dd.interferometer.channel)
        n_points = 101
        interferometer_calibration_times = time_basis[1:2:end]

        # identify thomson scattering subsystems
        ts_subsystems_mapper = Dict{String,Vector{Int}}()
        for (kch, ch) in enumerate(dd.thomson_scattering.channel)
            subsystem_name = IMAS.extract_subsystem_from_channel_name(ch.name)
            if subsystem_name ∉ keys(ts_subsystems_mapper)
                ts_subsystems_mapper[subsystem_name] = Int[]
            end
            push!(ts_subsystems_mapper[subsystem_name], kch)
        end

        # interpolate experimental measurements at calibration times
        experiments =
            [IMAS.linear_interp1d(ch_data.n_e_line_average.time, ch_data.n_e_line_average.data).(interferometer_calibration_times) for ch_data in dd.interferometer.channel]

        # evaluate rho,ne at the calibration times for each thomson scattering channel
        itp_ne = IMAS.fit2d(Val(:n_e), dd; transform=abs)
        nes = Vector{D}[]
        chρs = Vector{D}[]
        chr = D[ch.position.r for ch in dd.thomson_scattering.channel]
        chz = D[ch.position.z for ch in dd.thomson_scattering.channel]
        for time0 in interferometer_calibration_times
            i = IMAS.nearest_causal_time(dd.equilibrium.time, time0; bounds_error=false).index
            eqt = dd.equilibrium.time_slice[i]
            r, z, RHO_interpolant = IMAS.ρ_interpolant(eqt)
            ρ = RHO_interpolant.(chr, chz)
            data = itp_ne(ρ, range(time0, time0, length(ρ)))
            push!(chρs, ρ)
            push!(nes, data)
        end

        # optimization scales density for groups of thomson scattering channels
        function cost(scales)
            scales = abs.(scales)
            data = deepcopy(nes)
            c = Float64[]
            for (kt, time0) in enumerate(interferometer_calibration_times)
                for (kch, subsystem) in enumerate(keys(ts_subsystems_mapper))
                    data[kt][ts_subsystems_mapper[subsystem]] .*= scales[kch]
                end
                eqt = dd.equilibrium.time_slice[time0]
                index = sortperm(chρs[kt])
                for (kch, ch) in enumerate(dd.interferometer.channel)
                    density_thermal = IMAS.line_average(eqt, data[kt][index], chρs[kt][index], ch.line_of_sight; n_points)
                    simulation = density_thermal.line_average
                    push!(c, norm((experiments[kch][kt] .- simulation) ./ experiments[kch][kt]))
                end
            end
            return norm(c)
        end
        res = Optim.optimize(cost, fill(1.0, length(ts_subsystems_mapper)), Optim.NelderMead())

        # scale raw data in thomson_scattering IDS
        scales = res.minimizer
        scales_string = join(["$subsystem_name=$(@sprintf("%3.3f",scale))" for (subsystem_name,scale) in zip(keys(ts_subsystems_mapper),scales)], ", ")
        @info "Thomson subsystems scaled to match interferometer measurements: $scales_string"
        for (kch, subsystem) in enumerate(keys(ts_subsystems_mapper))
            for chnum in ts_subsystems_mapper[subsystem]
                dd.thomson_scattering.channel[chnum].n_e.data .*= scales[kch]
            end
        end
    end

    # fit ne
    itp_ne = IMAS.fit2d(Val(:n_e), dd; transform=abs)
    for (k, time0) in enumerate(time_basis)
        cp1d = dd.core_profiles.profiles_1d[k]
        data = itp_ne(rho_tor_norm12, range(time0, time0, length(rho_tor_norm12)))
        cp1d.electrons.density_thermal = IMAS.fit1d(rho_tor_norm12, data, rho_tor_norm; smooth1, smooth2).fit
    end

    # # fit Zeff
    # itp_zeff = IMAS.fit2d(Val(:zeff), dd; transform=x -> abs(max(x, 1.0) - 1.0))
    # for (k, time0) in enumerate(time_basis)
    #     cp1d = dd.core_profiles.profiles_1d[k]
    #     data = itp_zeff(rho_tor_norm12, range(time0, time0, length(rho_tor_norm12))) .+ 1.0
    #     cp1d.zeff = IMAS.fit1d(rho_tor_norm12, data, rho_tor_norm; smooth1, smooth2).fit
    # end

    # fit ni
    itp_nimp = IMAS.fit2d(Val(:n_i_over_n_e), dd; transform=abs)
    for (k, time0) in enumerate(time_basis)
        cp1d = dd.core_profiles.profiles_1d[k]
        bulk_ion = cp1d.ion[1]
        imp_ion = cp1d.ion[2]

        data = itp_nimp(rho_tor_norm12, range(time0, time0, length(rho_tor_norm12)))
        data .*= itp_ne(rho_tor_norm12, range(time0, time0, length(rho_tor_norm12)))
        n_i = IMAS.fit1d(rho_tor_norm12, data, rho_tor_norm; smooth1, smooth2).fit

        bulk_ion.density_thermal = zero(rho_tor_norm)
        imp_ion.density_thermal = n_i
    end

    # fit Ti
    itp_ti = IMAS.fit2d(Val(:t_i), dd; transform=abs)
    for (k, time0) in enumerate(time_basis)
        cp1d = dd.core_profiles.profiles_1d[k]
        data = itp_ti(rho_tor_norm12, range(time0, time0, length(rho_tor_norm12)))
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

    # fit ωtor
    itp_ωtor = IMAS.fit2d(Val(:ω_tor), dd)
    for (k, time0) in enumerate(time_basis)
        cp1d = dd.core_profiles.profiles_1d[k]
        data = itp_ωtor(rho_tor_norm12, range(time0, time0, length(rho_tor_norm12)))
        ωtor = IMAS.fit1d(rho_tor_norm12, data, rho_tor_norm; smooth1, smooth2).fit
        for ion in cp1d.ion
            ion.rotation_frequency_tor = ωtor
        end
    end

    # restore the original raw data
    # for field in (:thomson_scattering, :charge_exchange)
    #     setproperty!(dd, field, getproperty(dd1, field))
    # end

    return actor
end
