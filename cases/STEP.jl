
function case_parameters(::Type{Val{:STEP}}; init_from::Symbol=:scalars)::Tuple{ParametersAllInits,ParametersAllActors}
    ini, act = case_parameters(:STEP_scalars)

    if init_from == :ods
        # Fix the core profiles
        dd = IMAS.json2imas(ini.ods.filename)
        cp1d = dd.core_profiles.profiles_1d[]

        rho = cp1d.grid.rho_tor_norm

        cp1d.zeff = ones(length(rho)) .* ini.core_profiles.zeff
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

        ini.general.dd = dd
        ini.general.init_from = :ods
    end
    return ini, act
end


"""
    case_parameters(:STEP)

STEP

Arguments:
"""
function case_parameters(::Type{Val{:STEP_scalars}})::Tuple{ParametersAllInits,ParametersAllActors}
    ini = ParametersInits(; n_ec=1)
    act = ParametersActors()
    #### INI ####

    ini.general.casename = "STEP"
    ini.general.init_from = :scalars
    ini.ods.filename = joinpath(@__DIR__, "..", "sample", "STEP_starting_point.json")

    ini.build.layers = OrderedCollections.OrderedDict(
        :gap_OH => 0.2984356197352587,
        :OH => 0.13477737665463296,
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
        :lfs_gap_coils => 3.0,
        :lfs_TF => 0.4043321299638989,
        :gap_cryostat => 1.5,
        :cryostat => 0.2
    )
    ini.build.layers[:cryostat].shape = :rectangle
    ini.build.plasma_gap = 0.125
    ini.build.symmetric = true
    ini.build.divertors = :double
    ini.build.n_first_wall_conformal_layers = 5

    ini.equilibrium.B0 = 3.2
    ini.equilibrium.R0 = 3.6
    ini.equilibrium.ϵ = 1 / 1.8
    ini.equilibrium.κ = 2.93
    ini.equilibrium.δ = 0.59
    ini.equilibrium.ip = 20.9e6
    ini.equilibrium.xpoints = :double
    ini.equilibrium.boundary_from = :scalars

    ini.core_profiles.zeff = 2.0 # unkown
    ini.core_profiles.rot_core = 0.0
    ini.core_profiles.bulk = :DT
    ini.core_profiles.impurity = :Ne
    ini.core_profiles.helium_fraction = 0.01

    ini.oh.n_coils = 8
    ini.oh.technology = :HTS

    ini.pf_active.n_coils_inside = 10
    ini.pf_active.n_coils_outside = 0
    ini.pf_active.technology = :HTS

    ini.tf.n_coils = 12
    ini.tf.technology = :HTS
    ini.tf.shape = :rectangle

    act.ActorPFcoilsOpt.symmetric = true
    act.ActorEquilibrium.symmetrize = true

    ini.ec_launcher[1].power_launched = 150.e6 #  some at rho = 0.7 with a 0.2 width some in core 

    ini.requirements.flattop_duration = 1800.0
    ini.requirements.tritium_breeding_ratio = 1.1
    ini.requirements.power_electric_net = 100e6

    act.ActorFluxMatcher.evolve_densities = :flux_match
    act.ActorTGLF.user_specified_model = "sat1_em_iter"

    act.ActorStabilityLimits.models = Symbol[]

    set_new_base!(ini)
    set_new_base!(act)

    return ini, act
end
