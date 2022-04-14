import NumericalIntegration: cumul_integrate

function init_nbi(dd::IMAS.dd, ini::InitParameters)
    return init_nbi(dd, ini.nbi.power_launched, ini.nbi.beam_energy, ini.nbi.beam_mass, ini.nbi.toroidal_angle)
end

function init_nbi(
    dd::IMAS.dd,
    power_launched::Union{Real,Vector},
    beam_energy::Union{Real,Vector},
    beam_mass::Union{Real,Vector},
    toroidal_angle::Union{Real,Vector})

    power_launched, beam_energy, beam_mass, toroidal_angle = same_length_vectors(power_launched, beam_energy, beam_mass, toroidal_angle)

    for idx in 1:length(power_launched)
        resize!(dd.nbi.unit, idx)
        nbi_u = dd.nbi.unit[idx]
        @ddtime(nbi_u.energy.data = beam_energy[idx])
        @ddtime(nbi_u.power_launched.data = power_launched[idx])
        nbi_u.species.a = beam_mass[idx]
        # 1 beamlet
        resize!(nbi_u.beamlets_group, 1)
        nbi_u.beamlets_group[1].angle = toroidal_angle[idx] / 360 * 2pi
    end

    return dd
end

function init_ec_launchers(dd::IMAS.dd, ini::InitParameters)
    return init_ec_launchers(dd, ini.ec.power_launched)
end

function init_ec_launchers(dd::IMAS.dd, power_launched::Union{Real,Vector})
    (power_launched,) = same_length_vectors(power_launched)
    for idx in 1:length(power_launched)
        resize!(dd.ec_launchers.launcher, idx)
        @ddtime dd.ec_launchers.launcher[idx].power_launched.data = power_launched[idx]
    end
end

function init_ic_antennas(dd::IMAS.dd, ini::InitParameters)
    return init_ic_antennas(dd, ini.ic.power_launched)
end

function init_ic_antennas(dd::IMAS.dd, power_launched::Union{Real,Vector})
    (power_launched,) = same_length_vectors(power_launched)
    for idx in 1:length(power_launched)
        resize!(dd.ic_antennas.antenna, idx)
        @ddtime dd.ic_antennas.antenna[idx].power_launched.data = power_launched[idx]
    end
end

function init_lh_antennas(dd::IMAS.dd, ini::InitParameters)
    return init_lh_antennas(dd, ini.lh.power_launched)
end

function init_lh_antennas(dd::IMAS.dd, power_launched::Union{Real,Vector})
    (power_launched,) = same_length_vectors(power_launched)
    for idx in 1:length(power_launched)
        resize!(dd.lh_antennas.antenna, idx)
        @ddtime dd.lh_antennas.antenna[idx].power_launched.data = power_launched[idx]
    end
end

function init_core_sources(dd::IMAS.dd, ini::InitParameters, act::ActorParameters)
    init_from = ini.general.init_from

    if init_from == :gasc
        init_from = :scalars
    end

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
            finalize(step(SimpleNBIactor(dd)))
        end
        if !ismissing(ini.ec, :power_launched) && any(ini.ec.power_launched .> 0)
            init_ec_launchers(dd, ini)
            finalize(step(SimpleECactor(dd)))
        end
        if !ismissing(ini.ic, :power_launched) && any(ini.ic.power_launched .> 0)
            init_ic_antennas(dd, ini)
            finalize(step(SimpleICactor(dd)))
        end
        if !ismissing(ini.lh, :power_launched) && any(ini.lh.power_launched .> 0)
            init_lh_antennas(dd, ini)
            finalize(step(SimpleLHactor(dd)))
        end
    end

    return dd
end

