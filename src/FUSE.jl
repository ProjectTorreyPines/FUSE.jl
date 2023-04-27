__precompile__(true)

module FUSE

using IMAS
import Plots
using Plots
using Printf
using InteractiveUtils
import SnoopPrecompile

#= ===== =#
#  UTILS  #
#= ===== =#
include("utils_begin.jl")

#= =================== =#
#  ABSTRACT PARAMETERS  #
#= =================== =#
include("parameters.jl")

#= ====================== =#
#  PHYSICS and TECHNOLOGY  #
#= ====================== =#
include("physics.jl")
include("technology.jl")

#= ====== =#
#  DDINIT  #
#= ====== =#
include(joinpath("ddinit", "init.jl"))
include(joinpath("ddinit", "init_equilibrium.jl"))
include(joinpath("ddinit", "init_build.jl"))
include(joinpath("ddinit", "init_core_profiles.jl"))
include(joinpath("ddinit", "init_core_sources.jl"))
include(joinpath("ddinit", "init_currents.jl"))
include(joinpath("ddinit", "init_pf_active.jl"))
include(joinpath("ddinit", "init_others.jl"))
include(joinpath("ddinit", "gasc.jl"))

#= ====== =#
#  ACTORS  #
#= ====== =#
# the order of include matters due to import/using statements as well as the dependency of defines structures
include(joinpath("actors", "abstract_actors.jl"))

include(joinpath("actors", "equilibrium", "solovev_actor.jl"))
include(joinpath("actors", "equilibrium", "chease_actor.jl"))
include(joinpath("actors", "equilibrium", "equilibrium_actor.jl"))

include(joinpath("actors", "pf", "pf_active_actor.jl"))
include(joinpath("actors", "pf", "pf_passive_actor.jl"))

include(joinpath("actors", "build", "stresses_actor.jl"))
include(joinpath("actors", "build", "fluxswing_actor.jl"))
include(joinpath("actors", "build", "lfs_actor.jl"))
include(joinpath("actors", "build", "hfs_actor.jl"))
include(joinpath("actors", "build", "cx_actor.jl"))

include(joinpath("actors", "nuclear", "blanket_actor.jl"))
include(joinpath("actors", "nuclear", "neutronics_actor.jl"))

include(joinpath("actors", "current", "qed_actor.jl"))
include(joinpath("actors", "current", "steadycurrent_actor.jl"))

include(joinpath("actors", "hcd", "simple_common.jl"))
include(joinpath("actors", "hcd", "ec_simple_actor.jl"))
include(joinpath("actors", "hcd", "ic_simple_actor.jl"))
include(joinpath("actors", "hcd", "lh_simple_actor.jl"))
include(joinpath("actors", "hcd", "nb_simple_actor.jl"))
include(joinpath("actors", "hcd", "hcd_actor.jl"))

include(joinpath("actors", "pedestal_actor.jl"))

include(joinpath("actors", "divertors_actor.jl"))

include(joinpath("actors", "transport", "tauenn_actor.jl"))
include(joinpath("actors", "transport", "neoclassical_actor.jl"))
include(joinpath("actors", "transport", "tglf_actor.jl"))
include(joinpath("actors", "transport", "flux_calculator_actor.jl"))
include(joinpath("actors", "transport", "flux_matcher_actor.jl"))
include(joinpath("actors", "transport", "core_transport_actor.jl"))

include(joinpath("actors", "stability", "limits_actor.jl"))
include(joinpath("actors", "stability", "limit_models.jl"))

include(joinpath("actors", "balance_plant", "heat_transfer_actor.jl"))
include(joinpath("actors", "balance_plant", "thermal_cycle_actor.jl"))
include(joinpath("actors", "balance_plant", "power_needs_actor.jl"))
include(joinpath("actors", "balance_plant", "balance_of_plant_actor.jl"))
include(joinpath("actors", "balance_plant", "balance_of_plant_plot.jl"))

include(joinpath("actors", "costing", "costing_actor.jl"))
include(joinpath("actors", "costing", "costing_ARIES.jl"))
include(joinpath("actors", "costing", "costing_Sheffield.jl"))
include(joinpath("actors", "costing", "costing_fuse.jl"))

# NOTE: compound actors should be defined last
include(joinpath("actors", "compound", "equilibrium_transport_actor.jl"))
include(joinpath("actors", "compound", "whole_facility_actor.jl"))

#= ========== =#
#  PARAMETERS  #
#= ========== =#
include("parameters_inits.jl")
include("parameters_actors.jl")

#= ============ =#
#  OPTIMIZATION  #
#= ============ =#
include("optimization.jl")

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

#= ====== =#
#= EXPORT =#
#= ====== =#
export IMAS, @ddtime, constants, ±, ↔, Logging

#= ======== =#
#= __INIT__ =#
#= ======== =#
function __init__()
    # Rename `julia` process as `python3` on known shared systems where FUSE is run
    if gethostname() in ["omega-a.gat.com", "omega-b.gat.com"]
        change_julia_process_name("python3")
        @warn "Running in a shared cluster environment: FUSE is obfuscating the `julia` executable name as `python3`"
    end
end

end
