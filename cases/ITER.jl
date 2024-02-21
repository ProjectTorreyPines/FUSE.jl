"""
    case_parameters(:ITER; init_from::Symbol)

ITER

Arguments:

  - `init_from`: `:scalars` or `:ods` (ODS contains equilibrium and wall information)
"""
function case_parameters(::Type{Val{:ITER}}; init_from::Symbol, boundary_from=:MXH_params)::Tuple{ParametersAllInits,ParametersAllActors}
    ini = ParametersInits(; n_nb=1, n_ec=1, n_ic=1, n_lh=1, n_pl=1)
    act = ParametersActors()

    # checking init_from and boundary_from
    @assert boundary_from != :rz_points "boundary from :rz_points isn't supported in the ITER case"

    ini.general.casename = "ITER_$(init_from)"
    ini.general.init_from = init_from

    if init_from == :ods
        ini.ods.filename = joinpath("__FUSE__", "sample", "ITER_eq_ods.json")
        act.ActorCXbuild.rebuild_wall = false
        ini.equilibrium.boundary_from = :ods
        ini.equilibrium.xpoints = :lower
    else
        ini.equilibrium.B0 = -5.3
        ini.equilibrium.ip = 15e6
        ini.equilibrium.pressure_core = 0.643e6

        ini.equilibrium.xpoints = :lower
        ini.equilibrium.boundary_from = boundary_from

        R0 = 6.2
        Z0 = 0.4
        ϵ = 0.32
        κ = 1.85
        δ = 0.485
        ζ = -0.09583
        if boundary_from == :scalars
            ini.equilibrium.R0 = R0
            ini.equilibrium.Z0 = Z0
            ini.equilibrium.ϵ = ϵ
            ini.equilibrium.κ = κ
            ini.equilibrium.δ = δ
            ini.equilibrium.ζ = ζ
        else
            ini.equilibrium.MXH_params = [
                R0, Z0, ϵ, κ, 0.00337,
                0.15912, -0.05842, -0.04573, 0.00694, 0.00614, 0.00183,
                asin(δ), -ζ, -0.05597, -0.01655, 0.00204, 0.00306]
        end

        act.ActorCXbuild.rebuild_wall = true
    end
    act.ActorEquilibrium.model = :TEQUILA

    ini.equilibrium.ip = (t -> ramp(t / 100.0) * 10E6 + ramp((t - 100) / 100.0) * 5E6) ↔ (2, (100.0, 200.0), (10E6, 15E6))

    ini.time.pulse_shedule_time_basis = range(0, 300, 1000)
    ini.time.simulation_start = 300.0

    # explicitly set thickness of radial build layers
    ini.build.layers = layers = OrderedCollections.OrderedDict{Symbol,Float64}()
    layers[:gap_OH] = 1.329
    layers[:OH] = 0.734
    layers[:gap_OH_TF] = 0.112
    layers[:hfs_TF] = 0.909
    layers[:hfs_gap_TF_shield] = 0.016
    layers[:hfs_gap_shield] = 0.1
    layers[:hfs_gap_shield_vacuum_vessel] = 0.032
    layers[:hfs_vacuum_vessel] = 0.338
    layers[:hfs_gap_vacuum_vessel_blanket] = 0.01
    layers[:hfs_shield] = 0.465
    layers[:hfs_wall] = 0.06
    layers[:plasma] = 0.095+4.0+0.194
    layers[:lfs_wall] = 0.06
    layers[:lfs_blanket] = 0.491
    layers[:lfs_gap_blanket_vacuum_vessel] = 0.01
    layers[:lfs_vacuum_vessel] = 0.758
    layers[:lfs_gap_vacuum_vessel_shield] = 0.037
    layers[:lfs_gap_shield] = 0.20
    layers[:lfs_gap_shield_TF] = 0.05+0.2
    layers[:lfs_TF] = 0.882
    layers[:gap_tf_cryostat] = 3.343
    layers[:cryostat] = 0.05
    ini.build.layers = layers
    ini.build.n_first_wall_conformal_layers = 7

    ini.oh.n_coils = 6
    ini.pf_active.n_coils_inside = 0
    ini.pf_active.n_coils_outside = 6
    ini.pf_active.technology = :nbti
    act.ActorPFdesign.symmetric = false

    ini.tf.shape = :double_ellipse
    ini.tf.n_coils = 18
    ini.tf.technology = :nb3sn_iter

    ini.oh.technology = :nb3sn_iter
    ini.requirements.flattop_duration = 500.0 # 500 s for Q=10 scenario

    ini.core_profiles.greenwald_fraction = 0.9
    ini.core_profiles.greenwald_fraction_ped = ini.core_profiles.greenwald_fraction * 0.75
    ini.core_profiles.helium_fraction = 0.01
    ini.core_profiles.T_ratio = 1.0
    ini.core_profiles.T_shaping = 1.8
    ini.core_profiles.n_shaping = 0.9
    ini.core_profiles.zeff = 2.0
    ini.core_profiles.rot_core = 0.0
    ini.core_profiles.bulk = :DT
    ini.core_profiles.impurity = :Ne

    ini.nb_unit[1].power_launched = (t -> 16.7e6 + ramp((t - 100) / 100.0) * 16.7e6) ↔ (2, (100.0, 300.0), (0.0, 33.4E6), (:match, :float))
    ini.nb_unit[1].beam_energy = 1e6

    ini.ec_launcher[1].power_launched = (t -> 10e6 + ramp((t - 100) / 100.0) * 10e6) ↔ (2, (100.0, 300.0), (0.0, 20E6), (:match, :float))
    ini.ec_launcher[1].rho_0 = (t->0.5) ↔ (2, (100.0, 300.0), (0.0, 0.8), (:match, :float))

    ini.ic_antenna[1].power_launched = (t -> 12e6 + ramp((t - 100) / 100.0) * 12e6) ↔ (2, (100.0, 300.0), (0.0, 24E6), (:match, :float))

    ini.lh_antenna[1].power_launched = (t -> 5e6 + ramp((t - 100) / 100.0) * 5e6) ↔ (2, (100.0, 300.0), (0.0, 10E6), (:match, :float))

    ini.pellet_launcher[1].shape = :cylindrical
    ini.pellet_launcher[1].species = :T
    ini.pellet_launcher[1].size = Float64[0.003, 0.004] / 2.0
    ini.pellet_launcher[1].frequency = 0.01 # Hz

    act.ActorFluxMatcher.evolve_densities = :flux_match
    act.ActorTGLF.user_specified_model = "sat1_em_iter"

    act.ActorWholeFacility.update_build = false

    set_new_base!(ini)
    set_new_base!(act)

    return ini, act
end

function TraceCAD(::Type{Val{:ITER}})
    x_length = 16.9
    x_offset = -1.08
    y_offset = 0.1
    TraceCAD(:ITER, x_length, x_offset, y_offset)
end