"""
    case_parameters(
    	:ITER,
    	init_from::Symbol,
    	boundary_from::Symbol=:auto,
    	ne_setting::Symbol=:ne_ped,
    	time_dependent::Bool=false)

ITER
"""
function case_parameters(
    ::Type{Val{:ITER}};
    init_from::Symbol,
    boundary_from::Symbol=:auto,
    ne_setting::Symbol=:ne_ped,
    time_dependent::Bool=false)::Tuple{ParametersAllInits,ParametersAllActors}

    ini = ParametersInits(; n_nb=1, n_ec=1, n_ic=1, n_lh=1, n_pl=1)
    act = ParametersActors()

    ini.general.casename = "ITER_$(init_from)"
    ini.general.init_from = init_from

    if init_from == :ods
        wall_ods = joinpath("__FUSE__", "sample", "ITER_wall_ods.json")
        pf_active_ods = joinpath("__FUSE__", "sample", "ITER_pf_active_ods.json")
        equilibrium_ods = joinpath("__FUSE__", "sample", "ITER_equilibrium_ods.json")
        ini.ods.filename = "$(wall_ods),$(pf_active_ods),$(equilibrium_ods)"
        act.ActorCXbuild.rebuild_wall = false
        # act.ActorStabilityLimits.raise_on_breach = false
        if boundary_from == :auto
            boundary_from = :ods
        end
    else
        ini.equilibrium.B0 = -5.3
        act.ActorCXbuild.rebuild_wall = true
        if boundary_from == :auto
            boundary_from = :MXH_params
        end
    end
    act.ActorEquilibrium.model = :TEQUILA

    ini.equilibrium.xpoints = :lower
    ini.equilibrium.boundary_from = boundary_from

    if boundary_from == :scalars
        ini.equilibrium.R0 = 6.2
        ini.equilibrium.Z0 = 0.4
        ini.equilibrium.ϵ = 0.32
        ini.equilibrium.κ = 1.85
        ini.equilibrium.δ = 0.485
        ini.equilibrium.ζ = -0.09583
        ini.equilibrium.tilt = 0.01
        ini.equilibrium.𝚶 = 0.15912
    elseif boundary_from == :MXH_params
        R0 = 6.192066877538616
        Z0 = 0.35168671862415857
        ϵ = 0.32360681350046777
        κ = 1.8457782310407964
        c0 = 0.013312894232172886
        c = [0.18097024640622442, -0.06587365660361746]
        s = [0.4532039691273859, 0.11378936281355961]
        ini.equilibrium.MXH_params = [R0, Z0, ϵ, κ, c0, c..., s...]
    elseif boundary_from == init_from == :ods
        # pass
    else
        error("invalid boundary_from=:$boundary_from")
    end

    ini.equilibrium.pressure_core = 0.643e6
    ini.equilibrium.ip = 15E6

    # explicitly set thickness of radial build layers
    ini.build.layers = layers = OrderedCollections.OrderedDict{Symbol,Float64}()
    layers[:gap_OH] = 1.356
    layers[:OH] = 0.8
    layers[:gap_OH_TF] = 0.1
    layers[:hfs_TF] = 0.882
    layers[:hfs_gap_TF_shield] = 0.016 + 0.1 + 0.032
    layers[:hfs_vacuum_vessel] = 0.338
    layers[:hfs_gap_vacuum_vessel_blanket] = 0.01
    layers[:hfs_shield] = 0.465
    layers[:hfs_wall] = 0.06
    layers[:plasma] = 0.095 + 4.0 + 0.194
    layers[:lfs_wall] = 0.06
    layers[:lfs_shield] = 0.491
    layers[:lfs_gap_blanket_vacuum_vessel] = 0.01
    layers[:lfs_vacuum_vessel] = 0.758
    layers[:lfs_gap_shield_TF] = 0.25 + 0.2 + 0.037
    layers[:lfs_TF] = 0.882
    layers[:gap_cryostat] = 3.343
    layers[:cryostat] = 0.05
    ini.build.layers = layers
    ini.build.n_first_wall_conformal_layers = 4

    ini.build.layers[:OH].coils_inside = 6
    ini.build.layers[:gap_cryostat].coils_inside = 6
    act.ActorPFdesign.symmetric = false

    ini.oh.technology = :nb3sn_iter
    ini.pf_active.technology = :nbti
    ini.tf.technology = :nb3sn_iter

    ini.tf.shape = :circle_ellipse
    ini.tf.n_coils = 18

    ini.requirements.flattop_duration = 500.0 # 500 s for Q=10 scenario

    ini.core_profiles.ne_setting = ne_setting
    if ne_setting == :ne_ped
        ini.core_profiles.ne_value = 0.9 * 0.75 * IMAS.greenwald_density(ini.equilibrium.ip, layers[:plasma] / 2.0)
        act.ActorPedestal.density_match = :ne_ped
    elseif ne_setting == :greenwald_fraction_ped
        ini.core_profiles.ne_value = 0.9 * 0.75
        act.ActorPedestal.density_match = :ne_ped
    elseif ne_setting == :greenwald_fraction
        ini.core_profiles.ne_value = 0.9
        act.ActorPedestal.density_match = :ne_line
    end
    ini.core_profiles.ne_shaping = 1.0
    ini.core_profiles.helium_fraction = 0.01
    ini.core_profiles.T_ratio = 1.0
    ini.core_profiles.T_shaping = 1.8
    ini.core_profiles.zeff = 2.0
    ini.core_profiles.rot_core = 1e4
    ini.core_profiles.bulk = :DT
    ini.core_profiles.impurity = :Ne

    ini.nb_unit[1].power_launched = 33.4e6
    ini.nb_unit[1].beam_energy = 1e6
    ini.nb_unit[1].toroidal_angle = 20.0 * deg

    ini.ec_launcher[1].power_launched = 20E6
    ini.ec_launcher[1].rho_0 = 0.0

    ini.ic_antenna[1].power_launched = 24E6

    ini.lh_antenna[1].power_launched = 10E6

    ini.pellet_launcher[1].shape = :cylindrical
    ini.pellet_launcher[1].species = :T
    ini.pellet_launcher[1].size = Float64[0.003, 0.004] / 2.0
    ini.pellet_launcher[1].frequency = 0.02 # Hz

    if time_dependent
        ini.time.pulse_shedule_time_basis = range(0, 300; step=1.0)
        ini.time.simulation_start = 300.0
        rampup_ends = 10.0

        ini.rampup.side = :lfs
        ini.rampup.ends_at = rampup_ends
        ini.rampup.diverted_at = rampup_ends * 0.8

        ini.equilibrium.pressure_core = t -> ramp(t / rampup_ends) .^ 2 * 0.643e6
        ini.equilibrium.ip = t -> ramp(t / rampup_ends) * 14E6 + ramp((t - 100) / 100) * 1E6

        ini.nb_unit[1].power_launched = t -> (1 .+ ramp((t - 100) / 100)) * 16.7e6
        ini.ec_launcher[1].power_launched = t -> (1 .+ ramp((t - 100) / 100)) * 10E6
        ini.ic_antenna[1].power_launched = t -> (1 .+ ramp((t - 100) / 100)) * 12E6
        ini.lh_antenna[1].power_launched = t -> (1 .+ ramp((t - 100) / 100)) * 5E6
        ini.pellet_launcher[1].frequency = t -> (1 .+ ramp((t - 100) / 100)) * 0.01 # Hz
    end

    #### ACT ####

    act.ActorTGLF.tglfnn_model = "sat1_em_iter"

    act.ActorFluxMatcher.step_size = 0.5

    act.ActorWholeFacility.update_build = false

    return ini, act
end

function TraceCAD(::Type{Val{:ITER}})
    x_length = 16.9
    x_offset = -1.08
    y_offset = 0.1
    return TraceCAD(:ITER, x_length, x_offset, y_offset)
end
