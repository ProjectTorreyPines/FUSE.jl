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

    for idx in eachindex(power_launched)
        nbu = resize!(dd.nbi.unit, idx)[idx]
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
    for idx in eachindex(power_launched)
        ecl = resize!(dd.ec_launchers.beam, idx)[idx]
        ecl.name = length(power_launched) > 1 ? "ec_$idx" : "ec"
        @ddtime(ecl.power_launched.data = power_launched[idx])
        ecl.available_launch_power = power_launched[idx]
        if efficiency_conversion[idx] !== missing
            ecl.efficiency.conversion = efficiency_conversion[idx]
        end
        if efficiency_transmission[idx] !== missing
            ecl.efficiency.transmission = efficiency_transmission[idx]
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
    for idx in eachindex(power_launched)
        ica = resize!(dd.ic_antennas.antenna, idx)[idx]
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
    for idx in eachindex(power_launched)
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
    init_core_sources(dd::IMAS.dd, ini::ParametersAllInits, act::ParametersAllActors)

Initialize `dd.nbi`, `dd.ec_launchers`, `dd.ic_antennas`, `dd.lh_antennas` starting from `ini` and `act` parameters
"""
function init_core_sources(dd::IMAS.dd, ini::ParametersAllInits, act::ParametersAllActors)
    TimerOutputs.reset_timer!("init_core_sources")
    TimerOutputs.@timeit timer "init_core_sources" begin
        init_from = ini.general.init_from

        if init_from == :ods
            dd1 = IMAS.json2imas(ini.ods.filename)
            if !ismissing(dd1.core_sources, :time) && length(keys(dd1.core_sources.time)) > 0
                dd.global_time = max(dd.global_time, maximum(dd1.core_sources.time))
                dd.core_sources = dd1.core_sources
            else
                init_from = :scalars
            end
        end

        if init_from == :scalars || (init_from == :ods && isempty(dd1.nbi))
            if !ismissing(ini.nbi, :power_launched) && any(ini.nbi.power_launched .> 0)
                init_nbi(dd, ini, act)
            end
        end
        if init_from == :scalars || (init_from == :ods && isempty(dd1.ec_launchers))
            if !ismissing(ini.ec_launchers, :power_launched) && any(ini.ec_launchers.power_launched .> 0)
                init_ec_launchers(dd, ini, act)
            end
        end
        if init_from == :scalars || (init_from == :ods && isempty(dd1.ic_antennas))
            if !ismissing(ini.ic_antennas, :power_launched) && any(ini.ic_antennas.power_launched .> 0)
                init_ic_antennas(dd, ini, act)
            end
        end
        if init_from == :scalars || (init_from == :ods && isempty(dd1.lh_antennas))
            if !ismissing(ini.lh_antennas, :power_launched) && any(ini.lh_antennas.power_launched .> 0)
                init_lh_antennas(dd, ini, act)
            end
        end

        ActorHCD(dd, act)
        
        IMAS.sources!(dd)

        return dd
    end
end

