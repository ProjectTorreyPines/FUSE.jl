"""
    case_parameters(:ITER; init_from::Symbol)

ITER

Arguments:
* `init_from`: `:scalars` or `:ods` (ODS contains equilibrium and wall information)
"""
function case_parameters(::Type{Val{:ITER}}; init_from::Symbol)::Tuple{ParametersAllInits,ParametersAllActors}
    ini = ParametersInits()
    act = ParametersActors()
    ini.general.casename = "ITER_$(init_from)"
    ini.general.init_from = init_from

    ini.build.blanket = 0.0
    ini.build.shield = 0.5
    ini.build.vessel = 0.125
    ini.material.shield = "Tungsten"

    if init_from == :ods
        ini.ods.filename = joinpath(@__DIR__, "..", "sample", "ITER_eq_ods.json")
        act.ActorCXbuild.rebuild_wall = false
        act.ActorHFSsizing.fixed_aspect_ratio = true

        ini.equilibrium.boundary_from = :MXH_params
        ini.equilibrium.xpoints_number = 1
    else
        ini.equilibrium.B0 = -5.3
        ini.equilibrium.ip = 15e6
        ini.equilibrium.pressure_core = 0.643e6

        ini.equilibrium.boundary_from = :MXH_params
        ini.equilibrium.xpoints_number = 1

        R0 = 6.2
        Z0 = 0.0
        ϵ = 0.32
        κ = 1.85
        δ = 0.485
        ζ = -0.09583
        if ini.equilibrium.boundary_from == :scalars
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
        act.ActorHFSsizing.fixed_aspect_ratio = true
        # act.ActorEquilibrium.model = :CHEASE
    end

    # explicitly set thickness of radial build layers
    ini.build.layers = layers = DataStructures.OrderedDict()
    layers[:gap_OH] = 0.80
    layers[:OH] = 1.30
    layers[:hfs_TF] = 1.10
    layers[:hfs_vacuum_vessel] = 0.30
    layers[:hfs_shield] = 0.40
    layers[:hfs_wall] = 0.1
    layers[:plasma] = 4.40
    layers[:lfs_wall] = 0.1
    layers[:lfs_shield] = 0.40
    layers[:lfs_vacuum_vessel] = 1.05
    layers[:lfs_TF] = 1.10
    layers[:gap_cryostat] = 2.34
    layers[:cryostat] = 0.30
    ini.build.n_first_wall_conformal_layers = 3

    ini.pf_active.n_oh_coils = 6
    ini.pf_active.n_pf_coils_inside = 0
    ini.pf_active.n_pf_coils_outside = 6
    ini.pf_active.technology = coil_technology(:ITER, :PF)
    act.ActorPFcoilsOpt.symmetric = false

    ini.tf.shape = :princeton_D_scaled
    ini.tf.n_coils = 18
    ini.tf.technology = coil_technology(:ITER, :TF)

    ini.oh.technology = coil_technology(:ITER, :OH)
    act.ActorFluxSwing.operate_at_j_crit = false

    ini.requirements.flattop_duration = 1800.0

    ini.core_profiles.greenwald_fraction = 0.9
    ini.core_profiles.helium_fraction = 0.01
    ini.core_profiles.T_shaping = 1.8
    ini.core_profiles.zeff = 2.0
    ini.core_profiles.rot_core = 0.0
    ini.core_profiles.bulk = :DT
    ini.core_profiles.impurity = :Ne

    ini.nbi.power_launched = 2 * 16.7e6
    ini.nbi.beam_energy = 1e6
    ini.ec_launchers.power_launched = 2 * 10e6
    ini.ic_antennas.power_launched = 24 * 1e6

    act.ActorTransportSolver.evolve_densities = Dict(
        :Ne        => :match_ne_scale,
        :DT        => :quasi_neutrality,
        :He        => :match_ne_scale,
        :electrons => :flux_match)

    set_new_base!(ini)
    set_new_base!(act)

    return ini, act
end
