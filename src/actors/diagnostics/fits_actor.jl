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
            data = IMAS.smooth_by_convolution(leaf.ids, leaf.field, par.time_basis; window_size=par.time_average_window)
            setproperty!(leaf.ids, leaf.field, data)
            time_coord = IMAS.time_coordinate(leaf.ids, leaf.field)
            times_coords[time_coord.ids] = time_coord
        end
    end
    for time_coord in values(times_coords)
        setproperty!(time_coord.ids, time_coord.field, par.time_basis)
    end
    dd1.equilibrium = dd.equilibrium

    # fitting
    for time0 in par.time_basis
        cp1d = resize!(dd1.core_profiles.profiles_1d, time0)
        cp1d.grid.rho_tor_norm = rho_tor_norm
    end
    dd1.core_profiles.time = par.time_basis

    rho_tor_norm12 = range(0.0, 1.2, Int(ceil(par.rho_grid * 1.2)))

    # Te
    # task1 = Threads.@spawn begin
    if true
        index = [IMAS.hasdata(cp1d.electrons, :temperature) for cp1d in dd.core_profiles.profiles_1d]
        min_k_orig = findfirst(index)
        for (k, time0) in enumerate(par.time_basis)
            k_orig = max(min_k_orig, IMAS.nearest_causal_time(dd.core_profiles.time, time0; bounds_error=false).index)
            dd1.core_profiles.profiles_1d[k].electrons.temperature = dd.core_profiles.profiles_1d[k_orig].electrons.temperature
        end
    else
        itp = fit2d(Val(:t_e), dd1; transform=sqrt)
        for (k, time0) in enumerate(par.time_basis)
            cp1d = dd1.core_profiles.profiles_1d[k]
            data = itp.(rho_tor_norm12, range(time0, time0, length(rho_tor_norm12))) .^ 2
            cp1d.electrons.temperature = smooth_with_edge(rho_tor_norm12, data, rho_tor_norm; smooth1=0.5, smooth2=par.spatial_average_smoothing).fit
        end
    end
    # end

    # ne
    # task2 = Threads.@spawn begin
    if true
        index = [IMAS.hasdata(cp1d.electrons, :density_thermal) for cp1d in dd.core_profiles.profiles_1d]
        min_k_orig = findfirst(index)
        for (k, time0) in enumerate(par.time_basis)
            k_orig = max(min_k_orig, IMAS.nearest_causal_time(dd.core_profiles.time, time0; bounds_error=false).index)
            dd1.core_profiles.profiles_1d[k].electrons.density_thermal = dd.core_profiles.profiles_1d[k_orig].electrons.density_thermal
        end
    else
        itp = itp_ne = fit2d(Val(:n_e), dd1; transform=sqrt)
        for (k, time0) in enumerate(par.time_basis)
            cp1d = dd1.core_profiles.profiles_1d[k]
            data = itp.(rho_tor_norm12, range(time0, time0, length(rho_tor_norm12))) .^ 2
            cp1d.electrons.density_thermal = smooth_with_edge(rho_tor_norm12, data, rho_tor_norm; smooth1=0.5, smooth2=par.spatial_average_smoothing).fit
        end
    end
    # end

    # # zeff
    # # task3 = Threads.@spawn begin
    # if true
    #     index = [!IMAS.ismissing(cp1d, :zeff) for cp1d in dd.core_profiles.profiles_1d]
    #     min_k_orig = findfirst(index)
    #     for (k, time0) in enumerate(par.time_basis)
    #         k_orig = max(min_k_orig, IMAS.nearest_causal_time(dd.core_profiles.time, time0; bounds_error=false).index)
    #         dd1.core_profiles.profiles_1d[k].zeff = dd.core_profiles.profiles_1d[k_orig].zeff
    #     end
    # else
    #     itp = fit2d(Val(:zeff), dd1; transform=x -> sqrt(max(x, 1.0) - 1.0))
    #     for (k, time0) in enumerate(par.time_basis)
    #         cp1d = dd1.core_profiles.profiles_1d[k]
    #         data = itp.(rho_tor_norm12, range(time0, time0, length(rho_tor_norm12))) .^ 2.0 .+ 1.0
    #         cp1d.zeff = smooth_with_edge(rho_tor_norm12, data, rho_tor_norm; smooth1=0.5, smooth2=par.spatial_average_smoothing).fit
    #     end
    # end
    # # end

    # ti
    # task4 = Threads.@spawn begin
    if true
        index = [!IMAS.ismissing(cp1d, :t_i_average) for cp1d in dd.core_profiles.profiles_1d]
        min_k_orig = findfirst(index)
        for (k, time0) in enumerate(par.time_basis)
            k_orig = max(min_k_orig, IMAS.nearest_causal_time(dd.core_profiles.time, time0; bounds_error=false).index)
            dd1.core_profiles.profiles_1d[k].t_i_average = dd.core_profiles.profiles_1d[k_orig].t_i_average
        end
    else
        itp = fit2d(Val(:t_i), dd1; transform=sqrt)
        for (k, time0) in enumerate(par.time_basis)
            cp1d = dd1.core_profiles.profiles_1d[k]
            data = itp.(rho_tor_norm12, range(time0, time0, length(rho_tor_norm12))) .^ 2
            cp1d.t_i_average = smooth_with_edge(rho_tor_norm12, data, rho_tor_norm; smooth1=0.5, smooth2=par.spatial_average_smoothing).fit
        end
    end
    # end

    # wait(task1)
    # wait(task2)
    # wait(task3)
    # wait(task4)

    # rotation
    # task5 = Threads.@spawn begin
    if true
        index = [IMAS.hasdata(cp1d.ion[1], :density_thermal) for cp1d in dd.core_profiles.profiles_1d]
        min_k_orig = findfirst(index)
        for (k, time0) in enumerate(par.time_basis)
            cp1d = dd1.core_profiles.profiles_1d[k]
            bulk_ion, imp_ion = resize!(cp1d.ion, 2)
            IMAS.ion_element!(bulk_ion, :D)
            IMAS.ion_element!(imp_ion, :C)
            k_orig = max(min_k_orig, IMAS.nearest_causal_time(dd.core_profiles.time, time0; bounds_error=false).index)
            bulk_ion.density_thermal = dd.core_profiles.profiles_1d[k_orig].ion[1].density_thermal
            imp_ion.density_thermal = dd.core_profiles.profiles_1d[k_orig].ion[2].density_thermal
            bulk_ion.temperature = dd.core_profiles.profiles_1d[k_orig].t_i_average
            imp_ion.temperature = dd.core_profiles.profiles_1d[k_orig].t_i_average
        end
    else
        itp = fit2d(Val(:n_i_over_n_e), dd1; transform=sqrt)
        for (k, time0) in enumerate(par.time_basis)
            cp1d = dd1.core_profiles.profiles_1d[k]
            data = itp.(rho_tor_norm12, range(time0, time0, length(rho_tor_norm12))) .^ 2
            data .*= itp_ne.(rho_tor_norm12, range(time0, time0, length(rho_tor_norm12))) .^ 2

            n_i = smooth_with_edge(rho_tor_norm12, data, rho_tor_norm; smooth1=0.5, smooth2=par.spatial_average_smoothing).fit

            bulk_ion, imp_ion = resize!(cp1d.ion, 2)
            IMAS.ion_element!(bulk_ion, :D)
            IMAS.ion_element!(imp_ion, :C)
            bulk_ion.density_thermal = zero(rho_tor_norm)
            bulk_ion.temperature = cp1d.t_i_average
            imp_ion.density_thermal = n_i
            imp_ion.temperature = cp1d.t_i_average
            IMAS.enforce_quasi_neutrality!(cp1d, :D)
        end
    end
    # end

    index = [IMAS.hasdata(cp1d, :rotation_frequency_tor_sonic) for cp1d in dd.core_profiles.profiles_1d]
    min_k_orig = findfirst(index)
    if min_k_orig !== nothing
        for (k, time0) in enumerate(par.time_basis)
            k_orig = max(min_k_orig, IMAS.nearest_causal_time(dd.core_profiles.time, time0; bounds_error=false).index)
            dd1.core_profiles.profiles_1d[k].rotation_frequency_tor_sonic = dd.core_profiles.profiles_1d[k_orig].rotation_frequency_tor_sonic
        end
    end

    dd.core_profiles = dd1.core_profiles

    return actor
end

function fit2d(what::Any, dd::IMAS.dd{T}; transform::F=x -> x) where {T<:Real,F}
    # get data
    time, rho, data = FUSE.getdata(what, dd)

    # remove any NaN
    index = .!isnan.(data)
    if sum(.!index) != length(data)
        time = @views time[index]
        rho = @views rho[index]
        data = @views data[index]
    end

    return IMAS.NaturalNeighbours.interpolate(rho, time, transform.(data))
end

function smooth_with_edge(rho, data, rho_tor_norm; smooth1::Float64, smooth2::Float64) where {T<:Real}
    # linearize space
    result = IMAS.smooth_by_convolution(IMAS.Measurements.value.(data); xi=rho, xo=rho_tor_norm, window_size=smooth1)
    g = abs.(IMAS.gradient(rho_tor_norm, result))
    g .= cumsum(g)
    g .= cumsum(g)
    rho_inverse = @. (g - g[1]) / (g[end] - g[1]) * (rho_tor_norm[end] - rho_tor_norm[1]) + rho_tor_norm[1]
    rho_linearized = IMAS.interp1d(rho_tor_norm, rho_inverse).(rho)

    #rho_linearized, data = symmetrize(rho_linearized, data, :even)
    result = IMAS.smooth_by_convolution(data; xi=rho_linearized, xo=rho_inverse, window_size=smooth2)

    return (rho=rho, data=data, fit=result, rho_tor_norm=rho_tor_norm, rho_linearized=rho_linearized)
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

function getdata(what_val::Union{Val{:t_i},Val{:n_i_over_n_e},Val{:zeff}}, dd::IMAS.dd{T}) where {T<:Real}
    what = typeof(what_val).parameters[1]
    cer = dd.charge_exchange

    data = T[]
    data_σ = T[]
    time = T[]
    chr = T[]
    chz = T[]
    for ch in cer.channel
        if what == :zeff
            ch_data = getproperty(ch, what)
        else
            ch_data = getproperty(ch.ion[1], what)
        end
        append!(data, getproperty(ch_data, :data))
        if IMAS.hasdata(ch_data, :data_σ)
            append!(data_σ, getproperty(ch_data, :data_σ))
        end
        append!(time, getproperty(ch_data, :time))
        append!(chr, ch.position.r.data)
        append!(chz, ch.position.z.data)
    end

    time0 = unique(time)
    rho = chz .* 0.0
    for time0 in unique(time)
        i = IMAS.nearest_causal_time(dd.equilibrium.time, time0; bounds_error=false).index
        eqt = dd.equilibrium.time_slice[i]
        r, z, RHO_interpolant = IMAS.ρ_interpolant(eqt)
        index = time .== time0
        rho[index] = RHO_interpolant.(chr[index], chz[index])
    end

    if isempty(data_σ)
        return (time=time, rho=rho, data=data)
    else
        return (time=time, rho=rho, data=IMAS.Measurements.measurement.(data, data_σ))
    end
end


function getdata(what_val::Union{Val{:t_e},Val{:n_e}}, dd::IMAS.dd{T}) where {T<:Real}
    what = typeof(what_val).parameters[1]
    ts = dd.thomson_scattering

    data = T[]
    data_σ = T[]
    time = T[]
    chr = T[]
    chz = T[]
    for ch in ts.channel
        ch_data = getproperty(ch, what)
        _data = getproperty(ch_data, :data)
        append!(data, _data)
        if IMAS.hasdata(ch_data, :data_σ)
            append!(data_σ, getproperty(ch_data, :data_σ))
        end
        append!(time, getproperty(ch_data, :time))
        append!(chr, fill(ch.position.r, size(_data)))
        append!(chz, fill(ch.position.z, size(_data)))
    end

    time0 = unique(time)
    rho = chz .* 0.0
    for time0 in unique(time)
        i = IMAS.nearest_causal_time(dd.equilibrium.time, time0; bounds_error=false).index
        eqt = dd.equilibrium.time_slice[i]
        r, z, RHO_interpolant = IMAS.ρ_interpolant(eqt)
        index = time .== time0
        rho[index] = RHO_interpolant.(chr[index], chz[index])
    end

    if isempty(data_σ)
        return (time=time, rho=rho, data=data)
    else
        return (time=time, rho=rho, data=IMAS.Measurements.measurement.(data, data_σ))
    end
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

    if isempty(data_σ)
        return rho, data
    else
        return rho, IMAS.Measurements.measurement.(data, data_σ)
    end
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
