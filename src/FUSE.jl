__precompile__(true)

module FUSE

using IMAS
import Plots
using Plots
using Printf
using InteractiveUtils
import LinearAlgebra
using StaticArrays

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
include(joinpath("cases", "_cases.jl"))

#= ======= =#
#  PHYSICS  #
#= ======= =#
include("physics.jl")

#= ====== =#
#  DDINIT  #
#= ====== =#
include("signal.jl")
include(joinpath("parameters", "parameters_inits.jl"))

include(joinpath("ddinit", "init.jl"))
include(joinpath("ddinit", "init_from_ods.jl"))
include(joinpath("ddinit", "init_pulse_schedule.jl"))
include(joinpath("ddinit", "init_equilibrium.jl"))
include(joinpath("ddinit", "init_build.jl"))
include(joinpath("ddinit", "init_balance_of_plant.jl"))
include(joinpath("ddinit", "init_core_profiles.jl"))
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

include(joinpath("actors", "equilibrium", "solovev_actor.jl"))
include(joinpath("actors", "equilibrium", "chease_actor.jl"))
include(joinpath("actors", "equilibrium", "tequila_actor.jl"))
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

include(joinpath("actors", "current", "qed_actor.jl"))
include(joinpath("actors", "current", "steadycurrent_actor.jl"))
include(joinpath("actors", "current", "current_actor.jl"))

include(joinpath("actors", "hcd", "simple_common.jl"))
include(joinpath("actors", "hcd", "ec_simple_actor.jl"))
include(joinpath("actors", "hcd", "ic_simple_actor.jl"))
include(joinpath("actors", "hcd", "lh_simple_actor.jl"))
include(joinpath("actors", "hcd", "nb_simple_actor.jl"))
include(joinpath("actors", "hcd", "pellet_simple_actor.jl"))
include(joinpath("actors", "hcd", "hcd_actor.jl"))

include(joinpath("actors", "pedestal", "EPED_actor.jl"))
include(joinpath("actors", "pedestal", "WPED_actor.jl"))
include(joinpath("actors", "pedestal", "pedestal_actor.jl"))

include(joinpath("actors", "divertors", "divertors_actor.jl"))

include(joinpath("actors", "transport", "neoclassical_actor.jl"))
include(joinpath("actors", "transport", "tglf_actor.jl"))
include(joinpath("actors", "transport", "qlgyro_actor.jl"))
include(joinpath("actors", "transport", "flux_calculator_actor.jl"))
include(joinpath("actors", "transport", "flux_matcher_actor.jl"))
include(joinpath("actors", "transport", "eped_profiles_actor.jl"))
include(joinpath("actors", "transport", "core_transport_actor.jl"))

include(joinpath("actors", "stability", "limits_actor.jl"))
include(joinpath("actors", "stability", "limit_models.jl"))
include(joinpath("actors", "stability", "vertical_actor.jl"))

include(joinpath("actors", "balance_plant", "thermal_plant_actor.jl"))
include(joinpath("actors", "balance_plant", "power_needs_actor.jl"))
include(joinpath("actors", "balance_plant", "balance_of_plant_actor.jl"))

include(joinpath("actors", "costing", "costing_utils.jl"))
include(joinpath("actors", "costing", "sheffield_costing_actor.jl"))
include(joinpath("actors", "costing", "aries_costing_actor.jl"))
include(joinpath("actors", "costing", "costing_actor.jl"))

include(joinpath("actors", "control", "controller_actor.jl"))

include(joinpath("actors", "wall_loading", "particle_hf_actor.jl"))
include(joinpath("actors", "wall_loading", "corerad_hf_actor.jl"))

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

#= ========= =#
#  WORKFLOWS  #
#= ========= =#
include(joinpath("workflows", "yaml_workflow.jl"))
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

#= ====== =#
#= EXPORT =#
#= ====== =#
export IMAS, @ddtime, constants, ±, ↔, Logging
export step, pulse, ramp, trap, gaus, beta

end
