using FUSE
FUSE.logging(Logging.Info; actors=Logging.Error);

qce_regime = true


### NEW DEFINITIONS ###


# function alpha_edge_crit(dd) # from https://doi.org/10.1088/1741-4326/ad89da

#     kappa = dd.equilibrium.time_slice[1].boundary.elongation
#     delta = dd.equilibrium.time_slice[1].boundary.triangularity

#     α_edge_crit = 0.64*kappa^2.2 * (1+delta)^0.9

#     return α_edge_crit
# end

# function alpha_edge_faitsh(dd) # from ???

#     return α_edge_faitsh
# end




ini,act = FUSE.case_parameters(:KDEMO);

### Act settings
act.ActorEquilibrium.model = :TEQUILA
act.ActorCoreTransport.model = :FluxMatcher
act.ActorFluxMatcher.algorithm = :simple

act.ActorPFdesign.model=:optimal
act.ActorFluxSwing.operate_oh_at_j_crit = true # this maximizes flattop inside the fluxswing actor

act.ActorWholeFacility.update_plasma = true

# This is handeled by the constraint funcitons
act.ActorHFSsizing.error_on_performance = false
act.ActorHFSsizing.error_on_technology = false

act.ActorFluxMatcher.evolve_densities = :flux_match

### Ini settings act.ActorEPED.scaling = 0.8 [0.4, 0.6, 0.8]

resize!(ini.ec_launcher,1)
ini.ec_launcher[1].power_launched = 1e7 ↔ [1e4,5e7]
ini.ec_launcher[1].efficiency_conversion = 0.45
ini.ec_launcher[1].efficiency_transmission = 0.8
act.ActorSimpleEC.actuator[1].rho_0 = 0.2 #↔ [0.1,0.9]

ini.ic_antenna[1].power_launched = 5.0e7  ↔ [1e4,5e7]


ini.core_profiles.impurity = :Kr  #↔ (:Kr, :Ar, :Ne, :Xe)
act.ActorFluxMatcher.max_iterations = 500


#ini.equilibrium.ζ = 0.1 ↔ [0,0.2]
ini.equilibrium.B0 = 10.0 ↔ [8.0, 15.0]
ini.equilibrium.ip = 11.0e6 ↔ [10.0e6, 22e6]
ini.equilibrium.R0 = 5.5 ↔ [4.5, 8.0]
ini.equilibrium.pressure_core = missing

ini.bop.cycle_type = :brayton ↔ (:rankine, :brayton)
ini.tf.technology = :rebco #  ↔ (:nb3sn_iter, :rebco)
ini.tf.shape = :racetrack ↔ (:racetrack,:double_ellipse,:offset)
ini.oh.technology = :rebco # ↔ (:nb3sn_iter, :rebco)
ini.pf_active.technology = :nb3sn_iter # ↔ (:nb3sn_iter, :rebco)


### SAME AS H-MODE ###

act.ActorPedestal.model = :EPED
act.ActorFluxMatcher.rho_transport = 0.2:0.05:0.8
act.ActorTGLF.tglfnn_model ="sat0quench_em_d3d+mastu_azf+1"
act.ActorFluxMatcher.evolve_pedestal = true
ini.core_profiles.zeff = ini.core_profiles.zeff = 2.0 ↔ [1.5, 4.0]
ini.equilibrium.δ = 0.4 ↔ [0.3, 0.7]
ini.equilibrium.κ =  1.8
ini.tf.shape = :offset ↔ (:double_ellipse,:offset)
ini.requirements.lh_power_threshold_fraction = 1.

ini.core_profiles.ne_setting = :greenwald_fraction_ped
ini.core_profiles.ne_value = 0.5 ↔ [0.2, 1.0]    

if qce_regime

    act.ActorEPED.ped_factor = 0.8 # potentially scan this later
    ini.core_profiles.ne_sep_to_ped_ratio = 0.6 ↔ [0.5, 0.9]   

else
    
    act.ActorEPED.ped_factor = 1.0 # scan this 0.4, 0.6, 0.8, 1.0
    ini.core_profiles.ne_sep_to_ped_ratio = 0.25 # default

end    

# Requirements
ini.requirements.flattop_duration = 3600.
ini.requirements.log10_flattop_duration = missing # log10 versionis not needed

ini.requirements.power_electric_net = 250e6 # 250 +/- 50 MWe
ini.requirements.tritium_breeding_ratio = 1.1
ini.requirements.q95 = 3.0
ini.requirements.Psol_R = 15. 

ini.requirements.coil_j_margin = 0.1
ini.requirements.coil_stress_margin = 0.1


IMAS.update_ObjectiveFunctionsLibrary!()
IMAS.update_ConstraintFunctionsLibrary!()

OFL = deepcopy(IMAS.ObjectiveFunctionsLibrary)
CFL = deepcopy(IMAS.ConstraintFunctionsLibrary)

objective_functions = [OFL[:min_capital_cost],OFL[:max_q95]]#, OFL[:max_log10_flattop]]
transport_error_func = IMAS.ConstraintFunction(:max_transport_error, "", dd -> @ddtime(dd.transport_solver_numerics.convergence.time_step.data),<,1.1)
constraint_functions = []

if qce_regime
    constraint_functions = [
        CFL[:min_power_electric_net],
        CFL[:min_q95],
        transport_error_func,
        CFL[:min_lh_power_threshold_fraction],
        CFL[:max_Psol_R],
        CFL[:min_alpha_t],
        CFL[:max_Psol_R],
        CFL[:max_tf_coil_j],CFL[:max_oh_coil_j],
        CFL[:max_pl_stress],CFL[:max_tf_coil_stress],CFL[:max_oh_coil_stress]]
   
else
    constraint_functions = [
        CFL[:min_power_electric_net],
        CFL[:min_q95],
        transport_error_func,
        CFL[:min_lh_power_threshold_fraction],
        CFL[:max_Psol_R],
        CFL[:max_tf_coil_j],CFL[:max_oh_coil_j],
        CFL[:max_pl_stress],CFL[:max_tf_coil_stress],CFL[:max_oh_coil_stress]]
end

println("== OBJECTIVE FUNCTIONS ==")
display(objective_functions)
println()
println("== CONSTRAINT FUNCTIONS ==")
display(constraint_functions)

sty = FUSE.study_parameters(:MultiObjectiveOptimizer);

sty.server = "localhost"

sty.n_workers = N_workers = 126
sty.restart_workers_after_n_generations = 0# this is the default behavior and releases workers after running the study

sty.population_size = 500

sty.number_of_generations = 60

save_dir = "/pscratch/sd/s/slendebt/projects/ELM_free_study/results/QCE_1"

if !isdir(save_dir)
    mkdir(save_dir)
end

sty.save_folder = save_dir

println(sty.save_folder)


using Distributed
@everywhere import FUSE

study = FUSE.StudyMultiObjectiveOptimizer(sty,ini, act, constraint_functions, objective_functions); # it is possible to pass in keyword arguments to sty
