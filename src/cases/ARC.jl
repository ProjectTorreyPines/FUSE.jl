"""
    case_parameters(:ARC; init_from::Symbol=:scalars, flux_matcher::Bool=false)

CFS/MIT ARC V3A design

Zero-dimensional parameters follow table 2 and table 4 of Hillesheim et al.,
"Overview of the physics basis for the ARC fusion power plant", J. Plasma Phys. (2026),
doi:10.1017/S0022377826101706.

"""
function case_parameters(::Val{:ARC}; init_from::Symbol=:scalars, flux_matcher::Bool=false)
    ini = ParametersInits()
    act = ParametersActors()
    ini.general.casename = "ARC"
    ini.general.init_from = init_from

    # Keep any wall we are given: the released one is far more detailed than anything FUSE
    # would generate. Safe unconditionally -- `ActorCXbuild` still builds a wall from the
    # equilibrium when `dd.wall` is empty, regardless of this flag.
    act.ActorCXbuild.rebuild_wall = false

    if init_from == :ods
        ini.ods.filename = joinpath("__FUSE__", "sample", "ARCv3a_starting_dd.json")
        act.ActorWholeFacility.update_build = false
    end
    ini.equilibrium.boundary_from = :scalars

    # table 2, with κ/δ/ζ from an order-4 MXH fit to the geqdsk-ARCv3a separatrix
    # (paper quotes κ_sep = 1.8 and δ_sep = 0.65, which use different definitions than MXH)
    ini.equilibrium.R0 = 4.6229
    ini.equilibrium.ϵ = 0.2556
    ini.equilibrium.κ = 1.8915
    ini.equilibrium.δ = 0.558
    ini.equilibrium.ζ = -0.2297
    ini.equilibrium.B0 = -11.4
    ini.equilibrium.Z0 = 0.0
    ini.equilibrium.ip = 12.0e6
    ini.equilibrium.xpoints = :double

    # Radial build of ARC V3A, measured at the machine midplane.
    ini.build.n_first_wall_conformal_layers = 4
    ini.build.scale_layers_to_R0 = false # thicknesses below are absolute, not relative
    ini.build.plasma_gap = 0.017         # separatrix sits only ~2 cm off the first wall
    layers = OrderedCollections.OrderedDict{Symbol,Float64}()
    #                                       thickness      R        source
    layers[:gap_OH] = 0.7269               #            0.7269     [dev] CS inner edge
    layers[:OH] = 0.7562                   #            1.4831     [dev] CS coils
    layers[:hfs_TF] = 0.9069               #            2.3900     [fig2] TF case, bucked against the CS
    layers[:hfs_gap_TF_blanket_tank] = 0.07#            2.4600     [fig2]
    layers[:hfs_tank_wall] = 0.3          #            2.4800     [dev] blanket-tank shell
    layers[:hfs_blanket] = 0.55            #            3.3300     [dev] FLiBe (WC shield sits inside)
    layers[:hfs_vacuum_vessel_outer] = 0.03#            3.3600     [dev]
    layers[:hfs_gap_vacuum_vessel] = 0.045 #            3.4050     [dev] FLiBe cooling channels
    layers[:hfs_vacuum_vessel_inner] = 0.01#            3.4150     [dev]
    layers[:hfs_first_wall] = 0.005        #            3.4200     [dev] tungsten
    layers[:plasma] = 2.40                 #            5.8200     [dev] limiter contour
    layers[:lfs_first_wall] = 0.005        #            5.8250     [dev] tungsten
    layers[:lfs_vacuum_vessel_inner] = 0.01#            5.8350     [dev]
    layers[:lfs_gap_vacuum_vessel] = 0.045 #            5.8800     [dev] FLiBe cooling channels
    layers[:lfs_vacuum_vessel_outer] = 0.03#            5.9100     [dev]
    layers[:lfs_blanket] = 1.29            #            7.2000     [dev] FLiBe
    layers[:lfs_tank_wall] = 0.02          #            7.2200     [dev] blanket-tank shell
    layers[:lfs_gap_blanket_tank_TF] = 0.49#            7.7100     [fig2]
    layers[:lfs_TF] = 1.19                 #            8.9000     [fig2] TF case
    layers[:gap_cryostat] = 1.20           #           10.1000     [assm]
    layers[:cryostat] = 0.20               #           10.3000     [assm]

    ini.build.layers = layers

    ini.build.layers[:hfs_blanket].material = :flibe
    ini.build.layers[:lfs_blanket].material = :flibe
    ini.build.layers[:hfs_first_wall].material = :tungsten
    ini.build.layers[:lfs_first_wall].material = :tungsten

    # 6 CS coils (cs1u/l, cs2u/l, cs3u/l) then 10 PF coils (pf1..pf5, u/l), in the order
    # they appear in the device description
    ini.build.layers[:OH].coils_inside = 6
    ini.build.layers[:gap_cryostat].coils_inside = collect(7:16)

    ini.oh.technology = :rebco
    ini.pf_active.technology = :rebco
    ini.tf.technology = :rebco


    ini.tf.shape = :princeton_D
    ini.tf.n_coils = 18

    ini.center_stack.bucked = true
    ini.center_stack.plug = true

    ini.requirements.power_electric_net = 400e6      # table 2
    ini.requirements.flattop_duration = 900.0        # table 2
    ini.requirements.tritium_breeding_ratio = 1.1
    ini.requirements.coil_j_margin = 0.1
    ini.requirements.coil_stress_margin = 0.1

    # table 3 and table 4, with core values from input.gacode_V3A
    # f_G in the paper is a line-averaged Greenwald fraction, hence :ne_line below
    ini.core_profiles.ne_setting = :greenwald_fraction
    ini.core_profiles.ne_value = 0.9                 # f_G, table 4
    act.ActorPedestal.density_match = :ne_line
    ini.core_profiles.ne_shaping = 0.9
    ini.core_profiles.Te_core = 24.45e3
    ini.core_profiles.Te_shaping = 1.8
    ini.core_profiles.Ti_Te_ratio = 0.9              # table 3
    ini.core_profiles.zeff = 1.52                    # table 3
    ini.core_profiles.bulk = :DT
    ini.core_profiles.impurity = :B          # input.gacode_V3A seeds boron (~2% of ne) plus trace W
    ini.core_profiles.helium_fraction = 0.10
    ini.core_profiles.rot_core = 0.0

    # ICRF is the sole auxiliary heating system (50 MW max coupled; 21.5 MW absorbed
    # in the nominal flattop scenario of table 4)
    resize!(ini.ic_antenna, 1)
    ini.ic_antenna[1].power_launched = 21.5e6

    #### ACT ####

    act.ActorPFdesign.symmetric = true

    if !flux_matcher
        act.ActorCoreTransport.model = :none
    end

    act.ActorFluxMatcher.max_iterations = 50
    act.ActorFluxMatcher.verbose = true

    act.ActorTGLF.electromagnetic = false
    act.ActorTGLF.sat_rule = :sat0
    act.ActorTGLF.model = :TJLF

    return ini, act
end

# NOTE: ARC.jpg is the V2 CAD trace; the offsets below are not calibrated to V3A
function TraceCAD(::Val{:ARC})
    x_length = 7.23
    x_offset = 0.57
    y_offset = 0.05
    return TraceCAD(:ARC, x_length, x_offset, y_offset)
end
