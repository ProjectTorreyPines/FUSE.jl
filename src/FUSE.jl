module FUSE

__precompile__(true)

using IMAS
import Plots
using Plots
using Printf

#= ===== =#
#  UTILS  #
#= ===== =#
include("utils.jl")

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
include(joinpath("ddinit", "init_pf_active.jl"))
include(joinpath("ddinit", "init_others.jl"))
include(joinpath("ddinit", "gasc.jl"))

#= ====== =#
#  ACTORS  #
#= ====== =#
# the order of include matters due to import/using statements as well as the dependency of defines structures
include(joinpath("actors", "abstract_actors.jl"))
include(joinpath("actors", "equilibrium_actors.jl"))
include(joinpath("actors", "pf_active_actor.jl"))
include(joinpath("actors", "pf_passive_actor.jl"))
include(joinpath("actors", "stresses_actors.jl"))
include(joinpath("actors", "build_actors.jl"))
include(joinpath("actors", "blanket_actor.jl"))
include(joinpath("actors", "balance_of_plant_actor.jl"))
include(joinpath("actors", "current_actors.jl"))
include(joinpath("actors", "divertors_actor.jl"))
include(joinpath("actors", "sources_actors.jl"))
include(joinpath("actors", "tauenn_actor.jl"))
include(joinpath("actors", "costing_actor.jl"))
include(joinpath("actors", "neutronics_actor.jl"))
include(joinpath("actors", "neoclassical_actor.jl"))
include(joinpath("actors", "pedestal_actor.jl"))
include(joinpath("actors", "tglf_actor.jl"))
include(joinpath("actors", "core_transport_actor.jl"))
include(joinpath("actors", "transport_solver_actor.jl"))
include(joinpath("actors", "limits_actor.jl"))
# NOTE: compound actors should be defined last
include(joinpath("actors", "compound_actors.jl"))

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

#= ====== =#
#= EXPORT =#
#= ====== =#
export IMAS, @ddtime, constants, ±, ↔, Logging

end
