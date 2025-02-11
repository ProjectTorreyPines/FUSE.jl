"""
    init_ec!(dd::IMAS.dd, ini::ParametersAllInits, act::ParametersAllActors, dd1::IMAS.dd=IMAS.dd())

Initialize `dd.ec_launchers` starting from `ini` and `act` parameters
"""
function init_ec!(dd::IMAS.dd, ini::ParametersAllInits, act::ParametersAllActors, dd1::IMAS.dd=IMAS.dd())
    resize!(dd.ec_launchers.beam, length(ini.ec_launcher); wipe=false)
    @assert length(dd.ec_launchers.beam) == length(ini.ec_launcher) == length(dd.pulse_schedule.ec.beam)
    for (idx, (ecb, ini_ecb, ps_ecb)) in enumerate(zip(dd.ec_launchers.beam, ini.ec_launcher, dd.pulse_schedule.ec.beam))
        if ismissing(ecb, :name)
            ecb.name = length(ini.ec_launcher) > 1 ? "ec_$idx" : "ec"
        end
        ps_ecb.name = ecb.name
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
    resize!(dd.ic_antennas.antenna, length(ini.ic_antenna); wipe=false)
    @assert length(dd.ic_antennas.antenna) == length(ini.ic_antenna) == length(dd.pulse_schedule.ic.antenna)
    for (idx, (ica, ini_ica, ps_ica)) in enumerate(zip(dd.ic_antennas.antenna, ini.ic_antenna, dd.pulse_schedule.ic.antenna))
        if ismissing(ica, :name)
            ica.name = length(ini.ic_antenna) > 1 ? "ic_$idx" : "ic"
        end
        ps_ica.name = ica.name
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
    resize!(dd.lh_antennas.antenna, length(ini.lh_antenna); wipe=false)
    @assert length(dd.lh_antennas.antenna) == length(ini.lh_antenna) == length(dd.pulse_schedule.lh.antenna)
    for (idx, (lha, ini_lha, ps_lha)) in enumerate(zip(dd.lh_antennas.antenna, ini.lh_antenna, dd.pulse_schedule.lh.antenna))
        if ismissing(lha, :name)
            lha.name = length(ini.lh_antenna) > 1 ? "lh_$idx" : "lh"
        end
        ps_lha.name = lha.name
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
    resize!(dd.nbi.unit, length(ini.nb_unit); wipe=false)
    @assert length(dd.nbi.unit) == length(ini.nb_unit) == length(dd.pulse_schedule.nbi.unit)
    for (idx, (nbu, ini_nbu, ps_nbu)) in enumerate(zip(dd.nbi.unit, ini.nb_unit, dd.pulse_schedule.nbi.unit))
        if ismissing(nbu, :name)
            nbu.name = length(ini.nb_unit) > 1 ? "nbi_$idx" : "nbi"
        end
        ps_nbu.name = nbu.name
        @ddtime(nbu.energy.data = ini_nbu.beam_energy)
        nbu.available_launch_power = maximum(ps_nbu.power.reference)
        nbu.species.a = ini_nbu.beam_mass
        nbu.species.z_n = 1.0
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
    resize!(dd.pellets.launcher, length(ini.pellet_launcher); wipe=false)
    @assert length(dd.nbi.unit) == length(ini.nb_unit) == length(dd.pulse_schedule.nbi.unit)
    for (idx, (pll, ini_pll)) in enumerate(zip(dd.pellets.launcher, ini.pellet_launcher))
        if ismissing(pll, :name)
            pll.name = length(ini.pellet_launcher) > 1 ? "pellet_$idx" : "pellet"
        end
        #ps_ppl.name = pll.name
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
        end
        init_ec!(dd, ini, act, dd1)

        # IC
        if init_from == :ods && IMAS.hasdata(dd1.ic_antennas)
            dd.ic_antennas = deepcopy(dd1.ic_antennas)
        end
        init_ic!(dd, ini, act, dd1)

        # LH
        if init_from == :ods && IMAS.hasdata(dd1.lh_antennas)
            dd.lh_antennas = deepcopy(dd1.lh_antennas)
        end
        init_lh!(dd, ini, act, dd1)

        # NB
        if init_from == :ods && IMAS.hasdata(dd1.nbi)
            dd.nbi = deepcopy(dd1.nbi)
        end
        init_nb!(dd, ini, act, dd1)

        # PL
        if init_from == :ods && IMAS.hasdata(dd1.pellets)
            dd.pellets = deepcopy(dd1.pellets)
        end
        init_pl!(dd, ini, act, dd1)

        ActorHCD(dd, act)

        return dd
    end
end

