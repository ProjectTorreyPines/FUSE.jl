module FUSE

using IMAS
import IMAS: heaviside, pulse, ramp, trap, gaus, beta, sequence
import IMASutils: mirror_bound, argmin_abs, trapz
import Plots
using Plots
using HelpPlots
using Printf
using InteractiveUtils
import LinearAlgebra
using StaticArrays
import AbstractTrees: print_tree
import ProgressMeter
import Measurements: ±, Measurement

function __init__()
    # By default we disable use of threads in BLAS if using multiple Julia threads
    # BLAS threads have a positive effect on larger problem sizes and a unthreaded Julia
    # (eg. no benefit for matrices < 1000x1000 size)
    # but can have very detrimental effects when used in conjunction with Julia threads
    # https://github.com/ProjectTorreyPines/TJLF.jl/issues/8#issuecomment-1837648536
    if Threads.nthreads() > 1
        LinearAlgebra.BLAS.set_num_threads(1)
    end
end

const __FUSE__ = abspath(joinpath(@__DIR__, ".."))
const deg = pi / 180 # convert degrees to radians

#= ===== =#
#  UTILS  #
#= ===== =#
include("utils_begin.jl")

#= =================== =#
#  ABSTRACT PARAMETERS  #
#= =================== =#
include("parameters.jl")

#= ===== =#
#  CASES  #
#= ===== =#
include(joinpath("cases", "_toksys.jl"))
include(joinpath("cases", "_cases.jl"))

#= ======= =#
#  PHYSICS  #
#= ======= =#
include("physics.jl")
include("experiments.jl")

#= ====== =#
#  DDINIT  #
#= ====== =#
include(joinpath("parameters", "parameters_inits.jl"))
include(joinpath("parameters", "ini_from_ods.jl"))

include(joinpath("ddinit", "init.jl"))
include(joinpath("ddinit", "init_pulse_schedule.jl"))
include(joinpath("ddinit", "init_equilibrium.jl"))
include(joinpath("ddinit", "init_build.jl"))
include(joinpath("ddinit", "init_balance_of_plant.jl"))
include(joinpath("ddinit", "init_core_profiles.jl"))
include(joinpath("ddinit", "init_edge_profiles.jl"))
include(joinpath("ddinit", "init_hcd.jl"))
include(joinpath("ddinit", "init_core_sources.jl"))
include(joinpath("ddinit", "init_currents.jl"))
include(joinpath("ddinit", "init_pf_active.jl"))
include(joinpath("ddinit", "init_others.jl"))
include(joinpath("ddinit", "write_init_expressions.jl"))

#= ====== =#
#  ACTORS  #
#= ====== =#
# the order of include matters due to import/using statements as well as the dependency of defines structures
include("actors.jl")

include(joinpath("actors", "noop_actor.jl"))
include(joinpath("actors", "replay_actor.jl"))

include(joinpath("actors", "equilibrium", "tequila_actor.jl"))
include(joinpath("actors", "equilibrium", "fresco_actor.jl"))
include(joinpath("actors", "equilibrium", "eggo_actor.jl"))
include(joinpath("actors", "equilibrium", "chease_actor.jl"))
include(joinpath("actors", "equilibrium", "equilibrium_actor.jl"))

include(joinpath("actors", "pf", "pf_active_utils.jl"))
include(joinpath("actors", "pf", "pf_active_actor.jl"))
include(joinpath("actors", "pf", "pf_design_actor.jl"))
include(joinpath("actors", "pf", "pf_passive_actor.jl"))
include(joinpath("actors", "pf", "pf_plots.jl"))

include(joinpath("actors", "build", "oh_magnet.jl"))
include(joinpath("actors", "build", "tf_magnet.jl"))
include(joinpath("actors", "build", "stresses_actor.jl"))
include(joinpath("actors", "build", "fluxswing_actor.jl"))
include(joinpath("actors", "build", "lfs_actor.jl"))
include(joinpath("actors", "build", "hfs_actor.jl"))
include(joinpath("actors", "build", "cx_actor.jl"))

include(joinpath("actors", "nuclear", "blanket_actor.jl"))
include(joinpath("actors", "nuclear", "neutronics_actor.jl"))

include(joinpath("actors", "control", "ip_control.jl"))

include(joinpath("actors", "current", "qed_actor.jl"))
include(joinpath("actors", "current", "steadycurrent_actor.jl"))
include(joinpath("actors", "current", "current_actor.jl"))

include(joinpath("actors", "diagnostics", "fits_actor.jl"))
include(joinpath("actors", "diagnostics", "interferometer_actor.jl"))
include(joinpath("actors", "diagnostics", "magnetics_actor.jl"))

include(joinpath("actors", "hcd", "simple_common.jl"))
include(joinpath("actors", "hcd", "ec", "ec_simple_actor.jl"))
include(joinpath("actors", "hcd", "ec", "torbeam_actor.jl"))
include(joinpath("actors", "hcd", "ic_simple_actor.jl"))
include(joinpath("actors", "hcd", "lh_simple_actor.jl"))
include(joinpath("actors", "hcd", "nbi", "nb_simple_actor.jl"))
include(joinpath("actors", "hcd", "nbi", "rabbit_actor.jl"))
include(joinpath("actors", "hcd", "pl_simple_actor.jl"))
include(joinpath("actors", "hcd", "neutral_fueling_actor.jl"))
include(joinpath("actors", "hcd", "hcd_actor.jl"))

include(joinpath("actors", "pedestal", "analytic_pedestal_actor.jl"))
include(joinpath("actors", "pedestal", "EPED_actor.jl"))
include(joinpath("actors", "pedestal", "WPED_actor.jl"))
include(joinpath("actors", "pedestal", "pedestal_actor.jl"))

include(joinpath("actors", "divertors", "divertors_actor.jl"))

include(joinpath("actors", "transport", "neoclassical_actor.jl"))
include(joinpath("actors", "transport", "analytic_turbulence_actor.jl"))
include(joinpath("actors", "transport", "tglf_actor.jl"))
include(joinpath("actors", "transport", "qlgyro_actor.jl"))
include(joinpath("actors", "transport", "flux_calculator_actor.jl"))
include(joinpath("actors", "transport", "flux_matcher_actor.jl"))
include(joinpath("actors", "transport", "eped_profiles_actor.jl"))
include(joinpath("actors", "transport", "sawteeth_actor.jl"))
include(joinpath("actors", "transport", "core_transport_actor.jl"))

include(joinpath("actors", "stability", "limits_actor.jl"))
include(joinpath("actors", "stability", "limit_models.jl"))
include(joinpath("actors", "stability", "troyon_actor.jl"))
include(joinpath("actors", "stability", "vertical_actor.jl"))

include(joinpath("actors", "balance_plant", "thermal_plant_actor.jl"))
include(joinpath("actors", "balance_plant", "power_needs_actor.jl"))
include(joinpath("actors", "balance_plant", "balance_of_plant_actor.jl"))

include(joinpath("actors", "costing", "costing_utils.jl"))
include(joinpath("actors", "costing", "sheffield_costing_actor.jl"))
include(joinpath("actors", "costing", "aries_costing_actor.jl"))
include(joinpath("actors", "costing", "costing_actor.jl"))

include(joinpath("actors", "wall_loading", "particle_hf_actor.jl"))
include(joinpath("actors", "wall_loading", "corerad_hf_actor.jl"))

include(joinpath("actors", "sol", "sol_box_actor.jl"))
include(joinpath("actors", "sol", "sol_actor.jl"))

# NOTE: compound actors should be defined last
include(joinpath("actors", "compound", "stationary_plasma_actor.jl"))
include(joinpath("actors", "compound", "dynamic_plasma_actor.jl"))
include(joinpath("actors", "compound", "whole_facility_actor.jl"))

include(joinpath("parameters", "parameters_actors.jl"))

#= ============ =#
#  OPTIMIZATION  #
#= ============ =#
include("optimization.jl")

#= ======= =#
#  STUDIES  #
#= ======= =#
include(joinpath("parameters", "parameters_studies.jl"))
include("studies.jl")
include(joinpath("studies", "database_generator.jl"))
include(joinpath("studies", "multi_objective_optimization.jl"))
include(joinpath("studies", "TGLF_database.jl"))
include(joinpath("studies", "study_database.jl"))
include(joinpath("studies", "experiment_postdictive.jl"))

#= ========= =#
#  WORKFLOWS  #
#= ========= =#
include(joinpath("workflows", "optimization_workflow.jl"))
include(joinpath("workflows", "DB5_validation_workflow.jl"))

#= ======= =#
#  LOGGING  #
#= ======= =#
include("logging.jl")

#= ===== =#
#  UTILS  #
#= ===== =#
include("utils_end.jl")

#= ========== =#
#  PRECOMPILE  #
#= ========== =#
include("precompile.jl")

#= ===== =#
#  TESTS  #
#= ===== =#
include("test_cases.jl")

#= ====== =#
#= EXPORT =#
#= ====== =#
export IMAS, @ddtime, help, ±, ↔, Logging, print_tree, help_plot, help_plot!, @findall
export @checkin, @checkout
export step, pulse, ramp, trap, gaus, beta, sequence
export digest

end