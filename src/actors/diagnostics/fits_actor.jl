#= ================ =#
#  ActorFitProfiles  #
#= ================ =#
Base.@kwdef mutable struct FUSEparameters__ActorFitProfiles{T<:Real} <: ParametersActor{T}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    _time::Float64 = NaN
    #== actor parameters ==#

    #== display and debugging parameters ==#

end

mutable struct ActorFitProfiles{D,P} <: CompoundAbstractActor{D,P}
    dd::IMAS.dd{D}
    par::OverrideParameters{P,FUSEparameters__ActorFitProfiles{P}}
    act::ParametersAllActors{P}
    stage::Dict{Symbol,IMAS.dd{D}}
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

function ActorFitProfiles(dd::IMAS.dd{D}, par::FUSEparameters__ActorFitProfiles{P}, act::ParametersAllActors{P}; kw...) where {D<:Real,P<:Real}
    logging_actor_init(ActorFitProfiles)
    par = OverrideParameters(par; kw...)
    return ActorFitProfiles(dd, par, act, Dict{Symbol,IMAS.dd{D}}())
end

function _step(actor::ActorFitProfiles{D,P}) where {D<:Real,P<:Real}
    dd = actor.dd
    stage = actor.stage

    par = (time_average_window=0.1, spatial_average_smoothing=0.1, rho_grid=101, time_basis=dd.equilibrium.time)
    rho_tor_norm = range(0.0, 1.0, par.rho_grid)

    dd0 = IMAS.dd{D}()
    dd0.global_time = sum(extrema(par.time_basis)) / 2.0
    dd0.equilibrium = dd.equilibrium
    dd0.thomson_scattering = deepcopy(dd.thomson_scattering)
    dd0.charge_exchange = deepcopy(dd.charge_exchange)

    # identify outliers on raw data
    stage[:outliers] = dd0
    tg = IMAS.time_groups(dd0; min_channels=5)
    IMAS.adaptive_outlier_removal!(tg)

    # use the same time-basis
    stage[:retimed] = dd1 = dd0
    dd1.equilibrium = IMAS.equilibrium{D}()
    tg = IMAS.time_dependent_data(dd1)
    times_coords = Dict{IMAS.IDS,IMAS.IMASdd.Coordinate}()
    for group in values(tg)
        for leaf in group
            data = IMAS.smooth_by_convolution(leaf.ids, leaf.field, par.time_basis; window_size=0.05, no_nan=true)
            setproperty!(leaf.ids, leaf.field, data)
            time_coord = IMAS.time_coordinate(leaf.ids, leaf.field)
            times_coords[time_coord.ids] = time_coord
        end
    end
    for time_coord in values(times_coords)
        setproperty!(time_coord.ids, time_coord.field, par.time_basis)
    end
    dd1.equilibrium = dd.equilibrium

    # # identify outliers on raw data
    # stage[:no_nans] = dd2 = dd1
    # tg = IMAS.time_groups(dd2; min_channels=5)
    # IMAS.adaptive_outlier_removal!(tg)

    # fitting
    for time0 in par.time_basis
        cp1d = resize!(dd1.core_profiles.profiles_1d, time0)
        cp1d.grid.rho_tor_norm = rho_tor_norm

        t_e = fit1d(Val(:t_e), dd1, time0, rho_tor_norm; smooth1=0.5, smooth2=par.spatial_average_smoothing)
        cp1d.electrons.temperature = t_e.fit

        n_e = fit1d(Val(:n_e), dd1, time0, rho_tor_norm; smooth1=0.5, smooth2=par.spatial_average_smoothing)
        cp1d.electrons.density_thermal = n_e.fit

        zeff = fit1d(Val(:zeff), dd1, time0, rho_tor_norm; smooth1=0.5, smooth2=par.spatial_average_smoothing)
        cp1d.zeff = zeff.fit

        t_i = fit1d(Val(:t_i), dd1, time0, rho_tor_norm; smooth1=0.5, smooth2=par.spatial_average_smoothing)
        cp1d.t_i_average = t_i.fit

        n_imp = fit1d(Val(:n_imp), dd1, time0, rho_tor_norm; smooth1=0.5, smooth2=par.spatial_average_smoothing)
        bulk_ion, imp_ion = resize!(cp1d.ion, 2)
        IMAS.ion_element!(bulk_ion, :D)
        IMAS.ion_element!(imp_ion, :C)
        bulk_ion.density_thermal = zero(rho_tor_norm)
        bulk_ion.temperature = t_i.fit
        imp_ion.density_thermal = n_imp.fit
        imp_ion.temperature = t_i.fit
        IMAS.enforce_quasi_neutrality!(cp1d, :D)
    end

    #dd.core_profiles = dd1.core_profiles

    return actor
end

function fit1d(what::Any, dd::IMAS.dd{T}, time0::Float64, rho_tor_norm::AbstractVector{T}; smooth1::Float64, smooth2::Float64) where {T<:Real}
    # get data
    rho, data = getdata(what, dd, time0)

    # remove any NaN
    index = .!isnan.(rho) .&& .!isnan.(data)

    if sum(.!index) == length(data)
        return (rho=rho, data=data, fit=rho_tor_norm .* NaN, rho_tor_norm=rho_tor_norm, rho_linearized=rho)
    else
        rho = rho[index]
        data = data[index]
    end

    # sort by rho
    index = sortperm(rho)
    data .= data[index]
    rho .= rho[index]

    # linearize space
    result = IMAS.smooth_by_convolution(IMAS.Measurements.value.(data); xi=rho, xo=rho_tor_norm, window_size=smooth1, no_nan=true)
    g = abs.(IMAS.gradient(rho_tor_norm, result))
    g .= cumsum(g)
    g .= cumsum(g)
    rho_inverse = @. (g - g[1]) / (g[end] - g[1]) * (rho_tor_norm[end] - rho_tor_norm[1]) + rho_tor_norm[1]
    rho_linearized = IMAS.interp1d(rho_tor_norm, rho_inverse).(rho)

    #rho_linearized, data = symmetrize(rho_linearized, data, :even)
    result = IMAS.smooth_by_convolution(data; xi=rho_linearized, xo=rho_inverse, window_size=smooth2, interpolate=10, no_nan=true)

    return (rho=rho, data=data, fit=result, rho_tor_norm=rho_tor_norm, rho_linearized=rho_linearized)
end

function fit1d(what::Val{:n_imp}, dd::IMAS.dd{T}, time0::Float64, rho_tor_norm::AbstractVector{T}; smooth1::Float64, smooth2::Float64) where {T<:Real}
    n_i_over_n_e = fit1d(Val(:n_i_over_n_e), dd, time0, rho_tor_norm; smooth1, smooth2)
    ne = dd.core_profiles.profiles_1d[time0].electrons.density_thermal
    result = ne .* n_i_over_n_e.fit
    data = n_i_over_n_e.data .* IMAS.interp1d(rho_tor_norm, ne).(n_i_over_n_e.rho) / 6.0
    return (rho=n_i_over_n_e.rho, data=data, fit=result, rho_tor_norm=rho_tor_norm, rho_linearized=n_i_over_n_e.rho_linearized)
end

function getdata(what_val::Union{Val{:t_e},Val{:n_e}}, dd::IMAS.dd{T}, time0::Float64) where {T<:Real}
    what = typeof(what_val).parameters[1]
    ts = dd.thomson_scattering

    eqt = dd.equilibrium.time_slice[time0]
    r, z, RHO_interpolant = IMAS.ρ_interpolant(eqt)
    first = :nan
    last = :nan
    data = [IMAS.get_time_array(getproperty(ch, what), :data, time0, :linear; first, last) for ch in ts.channel]
    data_σ = [IMAS.get_time_array(getproperty(ch, what), :data_σ, time0, :linear; first, last) for ch in ts.channel if IMAS.hasdata(getproperty(ch, what), :data_σ)]
    rho = [RHO_interpolant(ch.position.r, ch.position.z) for ch in ts.channel]

    return rho, IMAS.Measurements.measurement.(data, data_σ)
end

function getdata(what_val::Union{Val{:t_i},Val{:n_i_over_n_e},Val{:zeff}}, dd::IMAS.dd{T}, time0::Float64) where {T<:Real}
    what = typeof(what_val).parameters[1]
    cer = dd.charge_exchange

    eqt = dd.equilibrium.time_slice[time0]
    r, z, RHO_interpolant = IMAS.ρ_interpolant(eqt)
    first = :nan
    last = :nan

    if what == :zeff
        data = T[IMAS.get_time_array(getproperty(ch, what), :data, time0, :linear; first, last) for ch in cer.channel]
        data_σ = T[IMAS.get_time_array(getproperty(ch, what), :data_σ, time0, :linear; first, last) for ch in cer.channel if IMAS.hasdata(getproperty(ch, what), :data_σ)]
    else
        data = T[IMAS.get_time_array(getproperty(ch.ion[1], what), :data, time0, :linear; first, last) for ch in cer.channel]
        data_σ =
            T[IMAS.get_time_array(getproperty(ch.ion[1], what), :data_σ, time0, :linear; first, last) for ch in cer.channel if IMAS.hasdata(getproperty(ch.ion[1], what), :data_σ)]
    end
    rho = [
        RHO_interpolant(IMAS.get_time_array(ch.position.r, :data, time0, :linear; first, last), IMAS.get_time_array(ch.position.z, :data, time0, :linear; first, last)) for
        ch in cer.channel
    ]

    if isempty(data_σ)
        return rho, data
    else
        return rho, IMAS.Measurements.measurement.(data, data_σ)
    end
end

function symmetrize(rho, data, symmetry)
    @assert symmetry in (:none, :even, :odd)
    if symmetry == :even
        rho = [-reverse(rho); rho]
        data = [reverse(data); data]
    elseif symmetry == :odd
        rho = [-reverse(rho); rho]
        data = [-reverse(data); data]
    end
    return rho, data
end
