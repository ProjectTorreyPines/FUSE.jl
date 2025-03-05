"""
    case_parameters(::Type{Val{:STEP}}; init_from::Symbol=:scalars, pf_from::Symbol=init_from)

UKAEA STEP design
"""
function case_parameters(::Type{Val{:STEP}}; init_from::Symbol=:scalars, pf_from::Symbol=init_from)
    @assert init_from in (:scalars, :ods)
    @assert pf_from in (:scalars, :ods)

    ini = ParametersInits()
    act = ParametersActors()
    #### INI ####

    ini.general.casename = "STEP"
    ini.general.init_from = :scalars
    ini.ods.filename = joinpath("__FUSE__", "sample", "STEP_starting_point.json")

    ini.build.layers = OrderedCollections.OrderedDict(
        :gap_OH => 0.233,
        :OH => 0.133,
        :gap_TF_OH => 0.016847172081829065,
        :hfs_TF => 0.4043321299638989,
        :hfs_gap_coils => 0.0,
        :hfs_thermal_shield => 0.03369434416365813,
        :hfs_vacuum_vessel => 0.5559566787003614,
        :hfs_blanket => 0.030541516245487195,
        :hfs_first_wall => 0.02,
        :plasma => 4.380264741275571,
        :lfs_first_wall => 0.02,
        :lfs_blanket => 0.6538868832731644,
        :lfs_vacuum_vessel => 0.6064981949458499,
        :lfs_thermal_shield => 0.13477737665463252 + 0.06738868832731626,
        :lfs_gap_coils => 0.25,
        :lfs_TF => 0.4043321299638989,
        :gap_cryostat => 1.5,
        :cryostat => 0.2
    )
    ini.build.plasma_gap = 0.125
    ini.build.symmetric = true
    ini.build.divertors = :double
    ini.build.n_first_wall_conformal_layers = 1

    ini.equilibrium.B0 = 3.2
    ini.equilibrium.R0 = 3.6
    ini.equilibrium.ϵ = 2.0 / ini.equilibrium.R0
    ini.equilibrium.κ = 2.93
    ini.equilibrium.δ = 0.59
    ini.equilibrium.ip = 21.1e6 # from PyTok
    ini.equilibrium.xpoints = :double
    ini.equilibrium.boundary_from = :scalars

    ini.core_profiles.ne_setting = :greenwald_fraction_ped
    ini.core_profiles.ne_value = 0.85
    ini.core_profiles.ne_shaping = 1.2
    ini.core_profiles.Te_core = 20E3
    ini.core_profiles.Te_shaping = 2.5
    ini.core_profiles.Ti_Te_ratio = 1.0
    ini.core_profiles.bulk = :DT
    ini.core_profiles.zeff = 2.5 # from PyTok
    ini.core_profiles.helium_fraction = 0.01  # No helium fraction in PyTok
    ini.core_profiles.impurity = :Ne #Barium :Ba
    ini.core_profiles.rot_core = 0.0
    ini.core_profiles.ejima = 0.1

    ini.build.layers[:OH].coils_inside = 5
    ini.build.layers[:lfs_vacuum_vessel].coils_inside = [9,12]
    ini.build.layers[:lfs_gap_coils].coils_inside = [6,7, 8, 10, 11, 13, 14, 15]

    ini.oh.technology = :rebco
    ini.pf_active.technology = :nb3sn_iter
    ini.tf.technology = :rebco

    ini.tf.shape = :rectangle
    ini.tf.n_coils = 12
    ini.tf.ripple = 0.005 # this is to avoid the TF coming in too close

    resize!(ini.ec_launcher, 1)
    ini.ec_launcher[1].power_launched = 150.e6
    resize!(act.ActorSimpleEC.actuator, 1)
    act.ActorSimpleEC.actuator[1].rho_0 = 0.0
    act.ActorSimpleEC.actuator[1].width = 0.25
    act.ActorSimpleEC.actuator[1].ηcd_scale = 0.5

    ini.requirements.flattop_duration = 1000.0
    ini.requirements.tritium_breeding_ratio = 1.1
    ini.requirements.power_electric_net = 250e6

    dd = IMAS.dd()
    # if init_from==:ods we need to sanitize the ODS that was given to us
    if init_from == :ods
        # Fix the core profiles
        dd = load_ods(ini.ods.filename)
        cp1d = dd.core_profiles.profiles_1d[]

        rho = cp1d.grid.rho_tor_norm
        cp1d.electrons.density_thermal = ((1.0 .- rho .^ 4) .* 1.6 .+ 0.4) .* 1E20

        cp1d.zeff = fill(ini.core_profiles.zeff, length(rho))
        cp1d.rotation_frequency_tor_sonic = IMAS.Hmode_profiles(0.0, ini.core_profiles.rot_core / 8, ini.core_profiles.rot_core, length(cp1d.grid.rho_tor_norm), 1.4, 1.4, 0.05)

        # Set ions:
        bulk_ion, imp_ion, he_ion = resize!(cp1d.ion, 3)
        # 1. DT
        IMAS.ion_element!(bulk_ion, ini.core_profiles.bulk)
        @assert bulk_ion.element[1].z_n == 1.0 "Bulk ion `$(ini.core_profiles.bulk)` must be a Hydrogenic isotope [:H, :D, :DT, :T]"
        # 2. Impurity
        IMAS.ion_element!(imp_ion, ini.core_profiles.impurity)
        # 3. He
        IMAS.ion_element!(he_ion, :He4)

        # pedestal
        summary = dd.summary
        @ddtime summary.local.pedestal.n_e.value = cp1d.electrons.density_thermal[argmin(abs.(rho .- (1 - ini.core_profiles.w_ped)))]
        @ddtime summary.local.pedestal.position.rho_tor_norm = 1 - ini.core_profiles.w_ped
        @ddtime summary.local.pedestal.zeff.value = ini.core_profiles.zeff

        # Zeff and quasi neutrality for a helium constant fraction with one impurity specie
        niFraction = zeros(3)
        # DT == 1
        # Imp == 2
        # He == 3
        zimp = imp_ion.element[1].z_n
        niFraction[3] = ini.core_profiles.helium_fraction
        niFraction[1] = (zimp - ini.core_profiles.zeff + 4 * niFraction[3] - 2 * zimp * niFraction[3]) / (zimp - 1)
        niFraction[2] = (ini.core_profiles.zeff - niFraction[1] - 4 * niFraction[3]) / zimp^2
        @assert !any(niFraction .< 0.0) "zeff impossible to match for given helium fraction [$(ini.core_profiles.helium_fraction))] and zeff [$(ini.core_profiles.zeff)]"
        ni_core = 0.0
        for i in 1:length(cp1d.ion)
            cp1d.ion[i].density_thermal = cp1d.electrons.density_thermal .* niFraction[i]
            cp1d.ion[i].temperature = cp1d.ion[1].temperature
            ni_core += cp1d.electrons.density_thermal[1] * niFraction[i]
        end

        ini.equilibrium.pressure_core = 1.175e6

        # -----------------
        ini.general.dd = dd
        ini.general.init_from = :ods
    end

    # pf_active
    # if pf_from=:ods we take the pf_coils that were given to us
    if pf_from == :ods
        coils = dd.pf_active.coil
        pf_rz = [
            (2.0429184549356223, 8.6986301369863),
            (4.017167381974248, 9.623287671232877),
            (6.815450643776823, 9.623287671232877),
            (6.832618025751072, 6.386986301369863),
            (8.309012875536478, 2.1061643835616444),
            (8.309012875536478, -2.1061643835616444),
            (6.832618025751072, -6.386986301369863),
            (6.815450643776823, -9.623287671232877),
            (4.017167381974248, -9.623287671232877),
            (2.0429184549356223, -8.6986301369863)]

        oh_zh = [
            (-6.471803481967896, 0.9543053940181108),
            (-5.000022627813203, 0.9677463150606198),
            (0.0, 4.9731407857281855),
            (5.000022627813203, 0.9677463150606198),
            (6.471803481967896, 0.9543053940181108)]

        r_oh = ini.build.layers[1].thickness + ini.build.layers[2].thickness / 2.0
        b = ini.equilibrium.ϵ * ini.equilibrium.R0 * ini.equilibrium.κ
        z_oh = (ini.equilibrium.Z0 - b, ini.equilibrium.Z0 + b)
        z_ohcoils, h_oh = size_oh_coils(z_oh[1], z_oh[2], 0.1, length(oh_zh), 1.0, 0.0)
        oh_zh = [(z, h_oh) for z in z_ohcoils]

        empty!(coils)
        resize!(coils, length(oh_zh) .+ length(pf_rz))

        for (idx, (z, h)) in enumerate(oh_zh)
            resize!(coils[idx].element, 1)
            pf_geo = coils[idx].element[1].geometry
            pf_geo.geometry_type = 2
            pf_geo.rectangle.r = r_oh
            pf_geo.rectangle.z = z
            pf_geo.rectangle.height = h
            pf_geo.rectangle.width = ini.build.layers[2].thickness
        end

        for (idx, (r, z)) in enumerate(pf_rz)
            idx += length(oh_zh)
            resize!(coils[idx].element, 1)
            pf_geo = coils[idx].element[1].geometry
            pf_geo.geometry_type = 2
            pf_geo.rectangle.r = r
            pf_geo.rectangle.z = z
            pf_geo.rectangle.height = 0.61
            pf_geo.rectangle.width = 0.53
        end

        IMAS.set_coils_function(coils, ini.equilibrium.R0)

        # -----------------
        ini.general.dd = dd
        ini.general.init_from = :ods
    end

    #### ACT ####

    act.ActorPFdesign.symmetric = true

    act.ActorCoreTransport.model = :FluxMatcher

    act.ActorFluxMatcher.evolve_densities = :fixed
    act.ActorFluxMatcher.rho_transport = 0.3:0.05:0.8

    act.ActorTGLF.tglfnn_model = "sat0_em_d3d"

    act.ActorPlasmaLimits.models = Symbol[]

    # High-beta ST equilibrium is tricky to converge
    act.ActorTEQUILA.relax = 0.01

    return ini, act
end
