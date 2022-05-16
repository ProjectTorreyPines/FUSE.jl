import NumericalIntegration: cumul_integrate

function init_nbi(dd::IMAS.dd, ini::ParametersInit)
    return init_nbi(dd,
        ini.nbi.power_launched,
        ini.nbi.beam_energy,
        ini.nbi.beam_mass,
        ini.nbi.toroidal_angle,
        evalmissing(ini.nbi, :efficiency_conversion),
        evalmissing(ini.nbi, :efficiency_transmission))
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
        nbu = resize!(dd.nbi.unit, idx)
        nbu.name = length(power_launched) > 1 ? "nbi_$idx" : "nbi"
        @ddtime(nbu.energy.data = beam_energy[idx])
        @ddtime(nbu.power_launched.data = power_launched[idx])
        nbu.available_launch_power = power_launched[idx]
        nbu.species.a = beam_mass[idx]
        # 1 beamlet
        beamlet = resize!(nbu.beamlets_group, 1)
        beamlet.angle = toroidal_angle[idx] / 360 * 2pi
        # Efficiencies
        nbu.efficiency.conversion = efficiency_conversion[idx]
        nbu.efficiency.transmission = efficiency_transmission[idx]
    end

    return dd
end

function init_ec_launchers(dd::IMAS.dd, ini::ParametersInit)
    return init_ec_launchers(dd, ini.ec.power_launched, evalmissing(ini.ec.efficiency_conversion), evalmissing(ini.ec.efficiency_transmission))
end

function init_ec_launchers(
    dd::IMAS.dd,
    power_launched::Union{Real,Vector},
    efficiency_conversion::Union{Real,Vector},
    efficiency_transmission::Union{Real,Vector})

    (power_launched, efficiency_conversion, efficiency_transmission) = same_length_vectors(power_launched, efficiency_conversion, efficiency_transmission)
    for idx in 1:length(power_launched)
        ecl = resize!(dd.ec_launchers.launcher, idx)
        ecl.name = length(power_launched) > 1 ? "ec_$idx" : "ec"
        @ddtime(ecl.power_launched.data = power_launched[idx])
        ecl.available_launch_power = power_launched[idx]
        ecl.efficiency.conversion = efficiency_conversion[idx]
        ecl.efficiency.transmission = efficiency_transmission[idx]
    end
end

function init_ic_antennas(dd::IMAS.dd, ini::ParametersInit)
    return init_ic_antennas(dd, ini.ic.power_launched, evalmissing(ini.ic.efficiency_conversion), evalmissing(ini.ic.efficiency_transmission), evalmissing(ini.ic.efficiency_coupling))
end

function init_ic_antennas(
    dd::IMAS.dd,
    power_launched::Union{Real,Vector},
    efficiency_conversion::Union{Real,Vector},
    efficiency_transmission::Union{Real,Vector},
    efficiency_coupling::Union{Real,Vector})

    (power_launched, efficiency_conversion, efficiency_transmission, efficiency_coupling) = same_length_vectors(power_launched, efficiency_conversion, efficiency_transmission, efficiency_coupling)
    for idx in 1:length(power_launched)
        ica = resize!(dd.ic_antennas.antenna, idx)
        ica.name = length(power_launched) > 1 ? "ic_$idx" : "ic"
        @ddtime(ica.power_launched.data = power_launched[idx])
        ica.available_launch_power = power_launched[idx]
        ica.efficiency.conversion = efficiency_conversion[idx]
        ica.efficiency.transmission = efficiency_transmission[idx]
        ica.efficiency.coupling = efficiency_coupling[idx]
    end
end

function init_lh_antennas(dd::IMAS.dd, ini::ParametersInit)
    return init_lh_antennas(dd, ini.lh.power_launched, evalmissing(ini.lh.efficiency_conversion), evalmissing(ini.lh.efficiency_transmission), evalmissing(ini.lh.efficiency_coupling))
end

function init_lh_antennas(
    dd::IMAS.dd,
    power_launched::Union{Real,Vector},
    efficiency_conversion::Union{Real,Vector},
    efficiency_transmission::Union{Real,Vector},
    efficiency_coupling::Union{Real,Vector})
    (power_launched, efficiency_conversion, efficiency_transmission, efficiency_coupling) = same_length_vectors(power_launched, efficiency_conversion, efficiency_transmission, efficiency_coupling)
    for idx in 1:length(power_launched)
        lha = resize!(dd.lh_antennas.antenna, idx)
        lha.name = length(power_launched) > 1 ? "lh_$idx" : "lh"
        @ddtime(lha.power_launched.data = power_launched[idx])
        lha.available_launch_power = power_launched[idx]
        lha.efficiency.conversion = efficiency_conversion[idx]
        lha.efficiency.transmission = efficiency_transmission[idx]
        lha.efficiency.coupling = efficiency_coupling[idx]
    end
end

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
            init_nbi(dd, ini)
            SimpleNBIactor(dd, act)
        end
        if !ismissing(ini.ec, :power_launched) && any(ini.ec.power_launched .> 0)
            init_ec_launchers(dd, ini)
            SimpleECactor(dd, act)
        end
        if !ismissing(ini.ic, :power_launched) && any(ini.ic.power_launched .> 0)
            init_ic_antennas(dd, ini)
            SimpleICactor(dd, act)
        end
        if !ismissing(ini.lh, :power_launched) && any(ini.lh.power_launched .> 0)
            init_lh_antennas(dd, ini)
            SimpleLHactor(dd, act)
        end
    end

    return dd
end

