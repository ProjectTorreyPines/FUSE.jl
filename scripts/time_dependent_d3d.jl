using FUSE

shot =  parse(Int64,ARGS[1])
workers = parse(Int64, ARGS[2])
FUSE.parallel_environment("localhost", workers)
ini, act = FUSE.case_parameters(:D3D, shot; fit_profiles=true) #NBI with balanced torque

ini.time.simulation_start = ini.general.dd.equilibrium.time_slice[2].time
dd = IMAS.dd()
FUSE.init!(dd, ini, act)

experiment_LH = FUSE.LH_analysis(dd; do_plot=false)

act.ActorPedestal.model = :dynamic
act.ActorPedestal.tau_n = experiment_LH.tau_n
act.ActorPedestal.tau_t = experiment_LH.tau_t
act.ActorWPED.ped_to_core_fraction = experiment_LH.W_ped_to_core_fraction
act.ActorEPED.ped_factor = 1.0
act.ActorPedestal.T_ratio_pedestal = 1.0 # Ti/Te in the pedestal

act.ActorPedestal.density_ratio_L_over_H = 1.0
act.ActorPedestal.zeff_ratio_L_over_H = 1.0

# LH-transition at user-defined times
act.ActorPedestal.mode_transitions = experiment_LH.mode_transitions
act.ActorPedestal.mode_transitions[5.2] = :L_mode

act.ActorEquilibrium.model = :EGGO # EGGO

act.ActorFRESCO.nR = 65
act.ActorFRESCO.nZ = 65

act.ActorNeutralFueling.τp_over_τe = 0.5

act.ActorFluxMatcher.evolve_plasma_sources = false
act.ActorFluxMatcher.algorithm = :simple
act.ActorFluxMatcher.max_iterations = 10
act.ActorFluxMatcher.verbose = false
act.ActorFluxMatcher.evolve_pedestal = false

act.ActorTGLF.tglfnn_model = "sat1_em_d3d"

dd.global_time = ini.general.dd.equilibrium.time_slice[2].time # start_time
δt = 0.025
final_time = ini.general.dd.equilibrium.time[end]
act.ActorDynamicPlasma.Nt = Int(ceil((final_time - dd.global_time) / δt))
act.ActorDynamicPlasma.Δt = final_time - dd.global_time

act.ActorDynamicPlasma.evolve_current = true
act.ActorDynamicPlasma.evolve_equilibrium = true
act.ActorDynamicPlasma.evolve_transport = true
act.ActorDynamicPlasma.evolve_hcd = true
act.ActorDynamicPlasma.evolve_pf_active = false
act.ActorDynamicPlasma.evolve_pedestal = true

act.ActorDynamicPlasma.ip_controller = false
act.ActorDynamicPlasma.time_derivatives_sources = true

actor = FUSE.ActorDynamicPlasma(dd, act; verbose=false)
result_path = get(ENV, "FUSE_RESULT_ARCHIVE", "")
if length(result_path) > 0
    Base.run(`mkdir -p $result_path`)
end
IMAS.imas2h5i(dd, joinpath(result_path,"fuse_time_dependent_$shot.h5"))
