import NumericalIntegration: cumul_integrate

function init_nbi(dd::IMAS.dd, par::Parameters)
    init_nbi(dd, par.nbi.beam_energy, par.nbi.beam_mass, par.nbi.beam_power, par.nbi.toroidal_angle)
end

function init_nbi(
    dd::IMAS.dd,
    beam_energy::Union{Real,Vector},
    beam_mass::Union{Real,Vector},
    beam_power::Union{Real,Vector},
    toroidal_angle::Union{Real,Vector})

    beam_energy, beam_mass, beam_power, toroidal_angle = IMAS.same_length_vectors(beam_energy, beam_mass, beam_power, toroidal_angle)

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
    dd
end

function init_core_sources(dd::IMAS.dd, par::Parameters)
    init_from = par.general.init_from

    if init_from == :gasc
        gasc = GASC(par.gasc.filename, par.gasc.case)
        init_core_sources(dd, gasc, par)

    elseif init_from == :ods
        dd1 = IMAS.json2imas(par.ods.filename)
        if length(keys(dd1.core_sources)) > 0
            dd.core_sources = dd1.core_sources
        else
            init_from = :scalars
        end
    end

    if init_from == :scalars
        if par.nbi.beam_power > 0
            init_nbi(dd, par)

            nbiactor = simpleNBIactor(dd)
            FUSE.step(nbiactor)
            FUSE.finalize(nbiactor)
        end
    end

    return dd
end

function init_core_sources(dd::IMAS.dd, gasc::GASC, par::Parameters)
    gasc = gasc.solution

    heating_power = gasc["OUTPUTS"]["current drive"]["powerAux"]
    cd_power = heating_power / gasc["INPUTS"]["current drive"]["auxCDPowerFactor"]
    plug = heating_power / gasc["INPUTS"]["power efficiency"]["efficiencyAux"]

    cd_powers = Dict()
    for system in ["NNB", "NB", "LH", "FW", "EC", "HI"]
        cd_powers[system] = cd_power * gasc["INPUTS"]["current drive"]["$(system)CDFraction"]
    end

    par = deepcopy(par)
    par.general.init_from = :scalars
    par.nbi.beam_power = heating_power # this will need to be expanded
    par.nbi.beam_energy = 200e3
    par.nbi.beam_mass = 2
    par.nbi.toroidal_angle = 0.0

    init_core_sources(dd, par)
    return dd
end
