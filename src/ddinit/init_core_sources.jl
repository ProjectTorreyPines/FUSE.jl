import NumericalIntegration: cumul_integrate

function init_nbi(dd::IMAS.dd, par::Parameters)
    return init_nbi(dd, par.nbi.beam_power, par.nbi.beam_energy, par.nbi.beam_mass, par.nbi.toroidal_angle)
end

function init_nbi(
    dd::IMAS.dd,
    beam_power::Union{Real,Vector},
    beam_energy::Union{Real,Vector},
    beam_mass::Union{Real,Vector},
    toroidal_angle::Union{Real,Vector})

    beam_power, beam_energy, beam_mass, toroidal_angle = IMAS.same_length_vectors(beam_power, beam_energy, beam_mass, toroidal_angle)

    for idx in 1:length(beam_power)
        resize!(dd.nbi.unit, idx)
        nbi_u = dd.nbi.unit[idx]
        @ddtime(nbi_u.energy.data = beam_energy[idx])
        @ddtime(nbi_u.power_launched.data = beam_power[idx])
        nbi_u.species.a = beam_mass[idx]
        # 1 beamlet
        resize!(nbi_u.beamlets_group, 1)
        nbi_u.beamlets_group[1].angle = toroidal_angle[idx] / 360 * 2pi
    end

    return dd
end

function init_ec_launchers(dd::IMAS.dd, par::Parameters)
    return init_ec_launchers(dd, par.ec.power_launched)
end

function init_ec_launchers(dd::IMAS.dd, power_launched::Union{Real,Vector})
    power_launched = [p for p in power_launched]
    for idx in 1:length(power_launched)
        resize!(dd.ec_launchers.launcher, idx)
        @ddtime dd.ec_launchers.launcher[idx].power_launched.data = power_launched[idx]
    end
end

function init_ic_antennas(dd::IMAS.dd, par::Parameters)
    return init_ic_antennas(dd, par.ic.power_launched)
end

function init_ic_antennas(dd::IMAS.dd, power_launched::Union{Real,Vector})
    power_launched = [p for p in power_launched]
    for idx in 1:length(power_launched)
        resize!(dd.ic_antennas.antenna, idx)
        @ddtime dd.ic_antennas.antenna[idx].power_launched.data = power_launched[idx]
    end
end

function init_lh_antennas(dd::IMAS.dd, par::Parameters)
    return init_lh_antennas(dd, par.lh.power_launched)
end

function init_lh_antennas(dd::IMAS.dd, power_launched::Union{Real,Vector})
    power_launched = [p for p in power_launched]
    for idx in 1:length(power_launched)
        resize!(dd.lh_antennas.antenna, idx)
        @ddtime dd.lh_antennas.antenna[idx].power_launched.data = power_launched[idx]
    end
end


function init_core_sources(dd::IMAS.dd, par::Parameters)
    init_from = par.general.init_from

    if init_from == :gasc
        gasc = GASC(par.gasc.filename, par.gasc.case)
        init_core_sources(dd, gasc, par)

    elseif init_from == :ods
        dd1 = IMAS.json2imas(par.ods.filename)
        if length(keys(dd1.core_sources)) > 0
            dd.global_time = max(dd.global_time, maximum(dd1.core_sources.time))
            dd.core_sources = dd1.core_sources
        else
            init_from = :scalars
        end
    end

    if init_from == :scalars
        if !ismissing(par.nbi, :beam_power) && any(par.nbi.beam_power .> 0)
            init_nbi(dd, par)
            finalize(step(simpleNBIactor(dd)))
        end
        if !ismissing(par.ec, :power_launched) && any(par.ec.power_launched .> 0)
            init_ec_launchers(dd, par)
            finalize(step(simpleECactor(dd)))
        end
        if !ismissing(par.ic, :power_launched) && any(par.ic.power_launched .> 0)
            init_ic_antennas(dd, par)
            finalize(step(simpleICactor(dd)))
        end
        if !ismissing(par.lh, :power_launched) && any(par.lh.power_launched .> 0)
            init_lh_antennas(dd, par)
            finalize(step(simpleLHactor(dd)))
        end
    end

    return dd
end

function init_core_sources(dd::IMAS.dd, gasc::GASC, par::Parameters)
    gasc = gasc.solution

    injected_power = gasc["OUTPUTS"]["current drive"]["powerAux"] * 1E6
    cd_power = injected_power / gasc["INPUTS"]["current drive"]["auxCDPowerFactor"]
    heating_power = injected_power - cd_power
    plug_power = injected_power / gasc["INPUTS"]["power efficiency"]["efficiencyAux"]

    cd_powers = Dict()
    for system in ["NNB", "NB", "LH", "FW", "EC", "HI"]
        cd_powers[system] = cd_power * gasc["INPUTS"]["current drive"]["$(system)CDFraction"]
    end

    par = deepcopy(par)
    par.general.init_from = :scalars
    par.nbi.beam_power = Float64[]
    par.nbi.beam_energy = Float64[]
    if heating_power >0
        push!(par.nbi.beam_power, heating_power)
        push!(par.nbi.beam_energy, 200e3)
    end
    if cd_powers["NB"] >0
        push!(par.nbi.beam_power, cd_powers["NB"])
        push!(par.nbi.beam_energy, 200e3)
    end
    if cd_powers["NNB"] >0
        push!(par.nbi.beam_power, cd_powers["NNB"])
        push!(par.nbi.beam_energy, 1000e3)
    end
    par.lh.power_launched = Float64[]
    if cd_powers["LH"] > 0
        push!(par.lh.power_launched, cd_powers["LH"])
    end
    if cd_powers["HI"] > 0
        push!(par.lh.power_launched, cd_powers["HI"])
    end
    par.ic.power_launched = cd_powers["FW"]
    par.ec.power_launched = cd_powers["EC"]

    init_core_sources(dd, par)
    return dd
end
