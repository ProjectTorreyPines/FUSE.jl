"""
    init_ec!(dd::IMAS.dd, ini::ParametersAllInits, act::ParametersAllActors, dd1::IMAS.dd=IMAS.dd())

Initialize `dd.ec_launchers` starting from `ini` and `act` parameters
"""
function init_ec!(dd::IMAS.dd, ini::ParametersAllInits, act::ParametersAllActors, dd1::IMAS.dd=IMAS.dd())
    empty!(dd.ec_launchers.beam)
    resize!(dd.ec_launchers.beam, length(ini.ec_launcher))
    for (idx, (ecb, ini_ecb, ps_ecb)) in enumerate(zip(dd.ec_launchers.beam, ini.ec_launcher, dd.pulse_schedule.ec.beam))
        ecb.name = length(ini.ec_launcher) > 1 ? "ec_$idx" : "ec"
        ecb.available_launch_power = maximum(ps_ecb.power_launched.reference)
        # Efficiencies
        ecb.efficiency.conversion = ini_ecb.efficiency_conversion
        ecb.efficiency.transmission = ini_ecb.efficiency_transmission
    end
    return dd
end

"""
    init_ic!(dd::IMAS.dd, ini::ParametersAllInits, act::ParametersAllActors, dd1::IMAS.dd=IMAS.dd())

Initialize `dd.ic_antennas` starting from `ini` and `act` parameters
"""
function init_ic!(dd::IMAS.dd, ini::ParametersAllInits, act::ParametersAllActors, dd1::IMAS.dd=IMAS.dd())
    empty!(dd.ic_antennas.antenna)
    resize!(dd.ic_antennas.antenna, length(ini.ic_antenna))
    for (idx, (ica, ini_ica, ps_ica)) in enumerate(zip(dd.ic_antennas.antenna, ini.ic_antenna, dd.pulse_schedule.ic.antenna))
        ica.name = length(ini.ic_antenna) > 1 ? "ic_$idx" : "ic"
        ica.available_launch_power = maximum(ps_ica.power.reference)
        # Efficiencies
        ica.efficiency.coupling = ini_ica.efficiency_coupling
        ica.efficiency.conversion = ini_ica.efficiency_conversion
        ica.efficiency.transmission = ini_ica.efficiency_transmission
    end
    return dd
end

"""
    init_lh!(dd::IMAS.dd, ini::ParametersAllInits, act::ParametersAllActors, dd1::IMAS.dd=IMAS.dd())

Initialize `dd.lh_antennas` starting from `ini` and `act` parameters
"""
function init_lh!(dd::IMAS.dd, ini::ParametersAllInits, act::ParametersAllActors, dd1::IMAS.dd=IMAS.dd())
    empty!(dd.lh_antennas.antenna)
    resize!(dd.lh_antennas.antenna, length(ini.lh_antenna))
    for (idx, (lha, ini_lha, ps_lha)) in enumerate(zip(dd.lh_antennas.antenna, ini.lh_antenna, dd.pulse_schedule.lh.antenna))
        lha.name = length(ini.lh_antenna) > 1 ? "lh_$idx" : "lh"
        lha.available_launch_power = maximum(ps_lha.power.reference)
        # Efficiencies
        lha.efficiency.coupling = ini_lha.efficiency_coupling
        lha.efficiency.conversion = ini_lha.efficiency_conversion
        lha.efficiency.transmission = ini_lha.efficiency_transmission
    end
    return dd
end

"""
    init_nb!(dd::IMAS.dd, ini::ParametersAllInits, act::ParametersAllActors, dd1::IMAS.dd=IMAS.dd())

Initialize `dd.nbi` starting from `ini` and `act` parameters
"""
function init_nb!(dd::IMAS.dd, ini::ParametersAllInits, act::ParametersAllActors, dd1::IMAS.dd=IMAS.dd())
    empty!(dd.nbi.unit)
    resize!(dd.nbi.unit, length(ini.nb_unit))
    for (idx, (nbu, ini_nbu, ps_nbu)) in enumerate(zip(dd.nbi.unit, ini.nb_unit, dd.pulse_schedule.nbi.unit))
        nbu.name = length(ini.nb_unit) > 1 ? "nbi_$idx" : "nbi"
        @ddtime(nbu.energy.data = ini_nbu.beam_energy)
        nbu.available_launch_power = maximum(ps_nbu.power.reference)
        nbu.species.a = ini_nbu.beam_mass
        # 1 beamlet
        beamlet = resize!(nbu.beamlets_group, 1)[1]
        beamlet.angle = ini_nbu.toroidal_angle / 360 * 2pi
        # Efficiencies
        nbu.efficiency.conversion = ini_nbu.efficiency_conversion
        nbu.efficiency.transmission = ini_nbu.efficiency_transmission
    end
    return dd
end

"""
    init_pl!(dd::IMAS.dd, ini::ParametersAllInits, act::ParametersAllActors, dd1::IMAS.dd=IMAS.dd())

Initialize `dd.pellet_launcher` starting from `ini` and `act` parameters
"""
function init_pl!(dd::IMAS.dd, ini::ParametersAllInits, act::ParametersAllActors, dd1::IMAS.dd=IMAS.dd())
    empty!(dd.pellets)
    resize!(dd.pellets.launcher, length(ini.pellet_launcher))
    for (idx, (pll, ini_pll)) in enumerate(zip(dd.pellets.launcher, ini.pellet_launcher))
        pll.name = length(ini.pellet_launcher) > 1 ? "pellet_$idx" : "pellet"
        pll.shape.type.name = string(ini_pll.shape)
        pll.shape.type.index = IMAS.name_2_index(pll.shape)[ini_pll.shape]
        pll.shape.size = ini_pll.size
        resize!(pll.species, 1)
        pll.species[1].label = string(ini_pll.species)
    end
    return dd
end

"""
    init_hcd!(dd::IMAS.dd, ini::ParametersAllInits, act::ParametersAllActors, dd1::IMAS.dd=IMAS.dd())

Initialize `dd.nbi`, `dd.ec_launchers`, `dd.ic_antennas`, `dd.lh_antennas` starting from `ini` and `act` parameters
"""
function init_hcd!(dd::IMAS.dd, ini::ParametersAllInits, act::ParametersAllActors, dd1::IMAS.dd=IMAS.dd())
    TimerOutputs.reset_timer!("init_hcd")
    TimerOutputs.@timeit timer "init_hcd" begin
        init_from = ini.general.init_from

        # EC
        if init_from == :ods && IMAS.hasdata(dd1.ec_launchers)
            dd.ec_launchers = deepcopy(dd1.ec_launchers)
        else
            init_ec!(dd, ini, act, dd1)
        end

        # IC
        if init_from == :ods && IMAS.hasdata(dd1.ic_antennas)
            dd.ic_antennas = deepcopy(dd1.ic_antennas)
        else
            init_ic!(dd, ini, act, dd1)
        end

        # LH
        if init_from == :ods && IMAS.hasdata(dd1.lh_antennas)
            dd.lh_antennas = deepcopy(dd1.lh_antennas)
        else
            init_lh!(dd, ini, act, dd1)
        end

        # NB
        if init_from == :ods && IMAS.hasdata(dd1.nbi)
            dd.nbi = deepcopy(dd1.nbi)
        else
            init_nb!(dd, ini, act, dd1)
        end

        # PL
        if init_from == :ods && IMAS.hasdata(dd1.pellets)
            dd.pellets = deepcopy(dd1.pellets)
        else
            init_pl!(dd, ini, act, dd1)
        end

        ActorHCD(dd, act)

        return dd
    end
end

