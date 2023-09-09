import NumericalIntegration: cumul_integrate

"""
    init_nbi(dd::IMAS.dd, ini::ParametersAllInits, act::ParametersAllActors)

Initialize `dd.nbi` starting from `ini` and `act` parameters
"""
function init_nbi(dd::IMAS.dd, ini::ParametersAllInits, act::ParametersAllActors)
    init_nbi(dd,
        ini.nbi.power_launched,
        ini.nbi.beam_energy,
        ini.nbi.beam_mass,
        ini.nbi.toroidal_angle,
        getproperty(ini.nbi, :efficiency_conversion, missing),
        getproperty(ini.nbi, :efficiency_transmission, missing))
end

function init_nbi(
    dd::IMAS.dd,
    power_launched::Union{Real,AbstractVector{<:Real}},
    beam_energy::Union{Real,AbstractVector{<:Real}},
    beam_mass::Union{Real,AbstractVector{<:Real}},
    toroidal_angle::Union{Real,AbstractVector{<:Real}},
    efficiency_conversion::Union{Missing,Real,AbstractVector{<:Real}},
    efficiency_transmission::Union{Missing,Real,AbstractVector{<:Real}})

    power_launched, beam_energy, beam_mass, toroidal_angle, efficiency_conversion, efficiency_transmission = same_length_vectors(power_launched, beam_energy, beam_mass, toroidal_angle, efficiency_conversion, efficiency_transmission)
    empty!(dd.nbi.unit)
    resize!(dd.nbi.unit, length(power_launched))
    for (idx, nbu) in enumerate(dd.nbi.unit)
        nbu.name = length(power_launched) > 1 ? "nbi_$idx" : "nbi"
        @ddtime(nbu.energy.data = beam_energy[idx])
        @ddtime(nbu.power_launched.data = power_launched[idx])
        nbu.available_launch_power = power_launched[idx]
        nbu.species.a = beam_mass[idx]
        # 1 beamlet
        beamlet = resize!(nbu.beamlets_group, 1)[1]
        beamlet.angle = toroidal_angle[idx] / 360 * 2pi
        # Efficiencies
        if efficiency_conversion[idx] !== missing
            nbu.efficiency.conversion = efficiency_conversion[idx]
        end
        if efficiency_transmission[idx] !== missing
            nbu.efficiency.transmission = efficiency_transmission[idx]
        end
    end

    return dd
end

"""
    init_ec_launchers(dd::IMAS.dd, ini::ParametersAllInits, act::ParametersAllActors)

Initialize `dd.ec_launchers` starting from `ini` and `act` parameters
"""
function init_ec_launchers(dd::IMAS.dd, ini::ParametersAllInits, act::ParametersAllActors)
    init_ec_launchers(dd,
        ini.ec_launchers.power_launched,
        getproperty(ini.ec_launchers, :efficiency_conversion, missing),
        getproperty(ini.ec_launchers, :efficiency_transmission, missing))
end

function init_ec_launchers(
    dd::IMAS.dd,
    power_launched::Union{Real,AbstractVector{<:Real}},
    efficiency_conversion::Union{Real,AbstractVector{<:Real},Missing},
    efficiency_transmission::Union{Real,AbstractVector{<:Real},Missing})

    (power_launched, efficiency_conversion, efficiency_transmission) = same_length_vectors(power_launched, efficiency_conversion, efficiency_transmission)
    empty!(dd.ec_launchers.beam)
    resize!(dd.ec_launchers.beam, length(power_launched))
    for (idx, ecb) in enumerate(dd.ec_launchers.beam)
        ecb.name = length(power_launched) > 1 ? "ec_$idx" : "ec"
        @ddtime(ecb.power_launched.data = power_launched[idx])
        ecb.available_launch_power = power_launched[idx]
        if efficiency_conversion[idx] !== missing
            ecb.efficiency.conversion = efficiency_conversion[idx]
        end
        if efficiency_transmission[idx] !== missing
            ecb.efficiency.transmission = efficiency_transmission[idx]
        end
    end
end

"""
    init_ic_antennas(dd::IMAS.dd, ini::ParametersAllInits, act::ParametersAllActors)

Initialize `dd.ic_antennas` starting from `ini` and `act` parameters
"""
function init_ic_antennas(dd::IMAS.dd, ini::ParametersAllInits, act::ParametersAllActors)
    init_ic_antennas(dd,
        ini.ic_antennas.power_launched,
        getproperty(ini.ic_antennas, :efficiency_conversion, missing),
        getproperty(ini.ic_antennas, :efficiency_transmission, missing),
        getproperty(ini.ic_antennas, :efficiency_coupling, missing))
end

function init_ic_antennas(
    dd::IMAS.dd,
    power_launched::Union{Real,AbstractVector{<:Real}},
    efficiency_conversion::Union{Real,AbstractVector{<:Real},Missing},
    efficiency_transmission::Union{Real,AbstractVector{<:Real},Missing},
    efficiency_coupling::Union{Real,AbstractVector{<:Real},Missing})

    (power_launched, efficiency_conversion, efficiency_transmission, efficiency_coupling) = same_length_vectors(power_launched, efficiency_conversion, efficiency_transmission, efficiency_coupling)
    empty!(dd.ic_antennas.antenna)
    resize!(dd.ic_antennas.antenna, length(power_launched))
    for (idx, ica) in enumerate(dd.ic_antennas.antenna)
        ica.name = length(power_launched) > 1 ? "ic_$idx" : "ic"
        @ddtime(ica.power_launched.data = power_launched[idx])
        ica.available_launch_power = power_launched[idx]
        if efficiency_coupling[idx] !== missing
            ica.efficiency.coupling = efficiency_coupling[idx]
        end
        if efficiency_conversion[idx] !== missing
            ica.efficiency.conversion = efficiency_conversion[idx]
        end
        if efficiency_transmission[idx] !== missing
            ica.efficiency.transmission = efficiency_transmission[idx]
        end
    end
end

"""
    init_lh_antennas(dd::IMAS.dd, ini::ParametersAllInits, act::ParametersAllActors)

Initialize `dd.lh_antennas` starting from `ini` and `act` parameters
"""
function init_lh_antennas(dd::IMAS.dd, ini::ParametersAllInits, act::ParametersAllActors)
    init_lh_antennas(dd,
        ini.lh_antennas.power_launched,
        getproperty(ini.lh_antennas, :efficiency_conversion, missing),
        getproperty(ini.lh_antennas, :efficiency_transmission, missing),
        getproperty(ini.lh_antennas, :efficiency_coupling, missing))
end

function init_lh_antennas(
    dd::IMAS.dd,
    power_launched::Union{Real,AbstractVector{<:Real}},
    efficiency_conversion::Union{Real,AbstractVector{<:Real},Missing},
    efficiency_transmission::Union{Real,AbstractVector{<:Real},Missing},
    efficiency_coupling::Union{Real,AbstractVector{<:Real},Missing})
    (power_launched, efficiency_conversion, efficiency_transmission, efficiency_coupling) = same_length_vectors(power_launched, efficiency_conversion, efficiency_transmission, efficiency_coupling)
    empty!(dd.lh_antennas.antenna)
    resize!(dd.lh_antennas.antenna, length(power_launched))
    for (idx, lha) in enumerate(dd.lh_antennas.antenna)
        lha = resize!(dd.lh_antennas.antenna, idx)[idx]
        lha.name = length(power_launched) > 1 ? "lh_$idx" : "lh"
        @ddtime(lha.power_launched.data = power_launched[idx])
        lha.available_launch_power = power_launched[idx]
        if efficiency_coupling[idx] !== missing
            lha.efficiency.coupling = efficiency_coupling[idx]
        end
        if efficiency_conversion[idx] !== missing
            lha.efficiency.conversion = efficiency_conversion[idx]
        end
        if efficiency_transmission[idx] !== missing
            lha.efficiency.transmission = efficiency_transmission[idx]
        end
    end
end

"""
    unique_core_sources_names!(cs::IMAS.core_sources)

add trailing number to sources that have the same name
"""
function unique_core_sources_names!(cs::IMAS.core_sources)
    tmp = Dict{String,Int}()
    for source in cs.source
        if source.identifier.name ∉ keys(tmp)
            tmp[source.identifier.name] = 1
        else
            tmp[source.identifier.name] += 1
        end
    end
    for source in cs.source
        if tmp[source.identifier.name] == 1
            tmp[source.identifier.name] = 0
        end
    end
    for source in reverse(cs.source)
        if source.identifier.name ∈ keys(tmp) && tmp[source.identifier.name] > 0
            tmp[source.identifier.name] -= 1
            source.identifier.name = "$(source.identifier.name)_$(tmp[source.identifier.name]+1)"
        end
    end
end

"""
    init_core_sources!(dd::IMAS.dd, ini::ParametersAllInits, act::ParametersAllActors, dd1::IMAS.dd)

Initialize `dd.nbi`, `dd.ec_launchers`, `dd.ic_antennas`, `dd.lh_antennas` starting from `ini` and `act` parameters
"""
function init_core_sources!(dd::IMAS.dd, ini::ParametersAllInits, act::ParametersAllActors, dd1::IMAS.dd)
    TimerOutputs.reset_timer!("init_core_sources")
    TimerOutputs.@timeit timer "init_core_sources" begin
        init_from = ini.general.init_from

        if init_from == :ods
            if !ismissing(dd1.core_sources, :time) && length(dd1.core_sources.time) > 0
                dd.core_sources = dd1.core_sources
                unique_core_sources_names!(dd.core_sources)
                if isempty(dd1.ec_launchers.beam) && findfirst(:ec, dd.core_sources.source) !== missing
                    @assert act.ActorHCD.ec_model == :none "Init using sources from ODS requires `act.ActorHCD.ec_model = :none` or data in `dd.ec_launchers.beam`"
                end
                if isempty(dd1.ic_antennas.antenna) && findfirst(:ic, dd.core_sources.source) !== missing
                    @assert act.ActorHCD.ic_model == :none "Init using sources from ODS requires `act.ActorHCD.ic_model = :none` or data in `dd.ic_launchers.antenna`"
                end
                if isempty(dd1.lh_antennas.antenna) && findfirst(:lh, dd.core_sources.source) !== missing
                    @assert act.ActorHCD.lh_model == :none "Init using sources from ODS requires `act.ActorHCD.lc_model = :none` or data in `dd.lh_antennas.antenna`"
                end
                if isempty(dd1.nbi.unit) && findfirst(:nbi, dd.core_sources.source) !== missing
                    @assert act.ActorHCD.nb_model == :none "Init using sources from ODS requires `act.ActorHCD.nb_model = :none` or data in `dd.nbi.unit`"
                end
            else
                init_from = :scalars
            end
        end

        if init_from == :scalars
            if !ismissing(ini.ec_launchers, :power_launched) && any(ini.ec_launchers.power_launched .> 0)
                init_ec_launchers(dd, ini, act)
            end
            if !ismissing(ini.ic_antennas, :power_launched) && any(ini.ic_antennas.power_launched .> 0)
                init_ic_antennas(dd, ini, act)
            end
            if !ismissing(ini.lh_antennas, :power_launched) && any(ini.lh_antennas.power_launched .> 0)
                init_lh_antennas(dd, ini, act)
            end
            if !ismissing(ini.nbi, :power_launched) && any(ini.nbi.power_launched .> 0)
                init_nbi(dd, ini, act)
            end
        end

        ActorHCD(dd, act)

        IMAS.sources!(dd)

        return dd
    end
end

