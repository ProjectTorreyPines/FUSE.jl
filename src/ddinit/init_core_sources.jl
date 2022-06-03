import NumericalIntegration: cumul_integrate

"""
    init_nbi(dd::IMAS.dd, ini::ParametersInit, act::ParametersActor)

Initialize `dd.nbi` starting from 0D `ini` parameters and `act` actor parameters.
"""
function init_nbi(dd::IMAS.dd, ini::ParametersInit, act::ParametersActor)
    init_nbi(dd,
        ini.nbi.power_launched,
        ini.nbi.beam_energy,
        ini.nbi.beam_mass,
        ini.nbi.toroidal_angle,
        evalmissing(ini.nbi, :efficiency_conversion),
        evalmissing(ini.nbi, :efficiency_transmission))
    ActorNBIsimple(dd, act)
end

function init_nbi(
    dd::IMAS.dd,
    power_launched::Union{Real,Vector},
    beam_energy::Union{Real,Vector},
    beam_mass::Union{Real,Vector},
    toroidal_angle::Union{Real,Vector},
    efficiency_conversion::Union{Missing,Real,Vector},
    efficiency_transmission::Union{Missing,Real,Vector})

    power_launched, beam_energy, beam_mass, toroidal_angle, efficiency_conversion, efficiency_transmission = same_length_vectors(power_launched, beam_energy, beam_mass, toroidal_angle, efficiency_conversion, efficiency_transmission)

    for idx in 1:length(power_launched)
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
        nbu.efficiency.conversion = efficiency_conversion[idx]
        nbu.efficiency.transmission = efficiency_transmission[idx]
    end

    return dd
end

"""
    init_ec_launchers(dd::IMAS.dd, ini::ParametersInit, act::ParametersActor)

Initialize `dd.ec_launchers` starting from 0D `ini` parameters and `act` actor parameters.
"""
function init_ec_launchers(dd::IMAS.dd, ini::ParametersInit, act::ParametersActor)
    init_ec_launchers(dd,
        ini.ec_launchers.power_launched,
        evalmissing(ini.ec_launchers, :efficiency_conversion),
        evalmissing(ini.ec_launchers, :efficiency_transmission))
    ActorECsimple(dd, act)
end

function init_ec_launchers(
    dd::IMAS.dd,
    power_launched::Union{Real,Vector},
    efficiency_conversion::Union{Real,Vector,Missing},
    efficiency_transmission::Union{Real,Vector,Missing})

    (power_launched, efficiency_conversion, efficiency_transmission) = same_length_vectors(power_launched, efficiency_conversion, efficiency_transmission)
    for idx in 1:length(power_launched)
        ecl = resize!(dd.ec_launchers.launcher, idx)[idx]
        ecl.name = length(power_launched) > 1 ? "ec_$idx" : "ec"
        @ddtime(ecl.power_launched.data = power_launched[idx])
        ecl.available_launch_power = power_launched[idx]
        ecl.efficiency.conversion = efficiency_conversion[idx]
        ecl.efficiency.transmission = efficiency_transmission[idx]
    end
end

"""
    init_ic_antennas(dd::IMAS.dd, ini::ParametersInit, act::ParametersActor)

Initialize `dd.ic_antennas` starting from 0D `ini` parameters and `act` actor parameters.
"""
function init_ic_antennas(dd::IMAS.dd, ini::ParametersInit, act::ParametersActor)
    init_ic_antennas(dd,
        ini.ic_antennas.power_launched,
        evalmissing(ini.ic_antennas, :efficiency_conversion),
        evalmissing(ini.ic_antennas, :efficiency_transmission),
        evalmissing(ini.ic_antennas, :efficiency_coupling))
    ActorICsimple(dd, act)
end

function init_ic_antennas(
    dd::IMAS.dd,
    power_launched::Union{Real,Vector},
    efficiency_conversion::Union{Real,Vector,Missing},
    efficiency_transmission::Union{Real,Vector,Missing},
    efficiency_coupling::Union{Real,Vector,Missing})

    (power_launched, efficiency_conversion, efficiency_transmission, efficiency_coupling) = same_length_vectors(power_launched, efficiency_conversion, efficiency_transmission, efficiency_coupling)
    for idx in 1:length(power_launched)
        ica = resize!(dd.ic_antennas.antenna, idx)[idx]
        ica.name = length(power_launched) > 1 ? "ic_$idx" : "ic"
        @ddtime(ica.power_launched.data = power_launched[idx])
        ica.available_launch_power = power_launched[idx]
        ica.efficiency.conversion = efficiency_conversion[idx]
        ica.efficiency.transmission = efficiency_transmission[idx]
        ica.efficiency.coupling = efficiency_coupling[idx]
    end
end

"""
    init_lh_antennas(dd::IMAS.dd, ini::ParametersInit, act::ParametersActor)

Initialize `dd.lh_antennas` starting from 0D `ini` parameters and `act` actor parameters.
"""
function init_lh_antennas(dd::IMAS.dd, ini::ParametersInit, act::ParametersActor)
    init_lh_antennas(dd,
        ini.lh_antennas.power_launched,
        evalmissing(ini.lh_antennas, :efficiency_conversion),
        evalmissing(ini.lh_antennas, :efficiency_transmission),
        evalmissing(ini.lh_antennas, :efficiency_coupling))
    ActorLHsimple(dd, act)
end

function init_lh_antennas(
    dd::IMAS.dd,
    power_launched::Union{Real,Vector},
    efficiency_conversion::Union{Real,Vector,Missing},
    efficiency_transmission::Union{Real,Vector,Missing},
    efficiency_coupling::Union{Real,Vector,Missing})
    (power_launched, efficiency_conversion, efficiency_transmission, efficiency_coupling) = same_length_vectors(power_launched, efficiency_conversion, efficiency_transmission, efficiency_coupling)
    for idx in 1:length(power_launched)
        lha = resize!(dd.lh_antennas.antenna, idx)[idx]
        lha.name = length(power_launched) > 1 ? "lh_$idx" : "lh"
        @ddtime(lha.power_launched.data = power_launched[idx])
        lha.available_launch_power = power_launched[idx]
        lha.efficiency.conversion = efficiency_conversion[idx]
        lha.efficiency.transmission = efficiency_transmission[idx]
        lha.efficiency.coupling = efficiency_coupling[idx]
    end
end

"""
    init_core_sources(dd::IMAS.dd, ini::ParametersInit, act::ParametersActor)

Initialize `dd.nbi`, `dd.ec_launchers`, `dd.ic_antennas`, `dd.lh_antennas` starting from 0D `ini` parameters and `act` actor parameters.
"""
function init_core_sources(dd::IMAS.dd, ini::ParametersInit, act::ParametersActor)
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

    if init_from == :scalars
        if !ismissing(ini.nbi, :power_launched) && any(ini.nbi.power_launched .> 0)
            init_nbi(dd, ini, act)
        end
        if !ismissing(ini.ec_launchers, :power_launched) && any(ini.ec_launchers.power_launched .> 0)
            init_ec_launchers(dd, ini, act)
        end
        if !ismissing(ini.ic_antennas, :power_launched) && any(ini.ic_antennas.power_launched .> 0)
            init_ic_antennas(dd, ini, act)
        end
        if !ismissing(ini.lh_antennas, :power_launched) && any(ini.lh_antennas.power_launched .> 0)
            init_lh_antennas(dd, ini, act)
        end
    end

    return dd
end

