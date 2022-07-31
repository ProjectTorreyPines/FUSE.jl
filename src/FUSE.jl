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

#= ========== =#
#  PARAMETERS  #
#= ========== =#
include("parameters.jl")
include("parameters_init.jl")

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
include(joinpath("actors", "pf_active_actors.jl"))
include(joinpath("actors", "stresses_actors.jl"))
include(joinpath("actors", "build_actors.jl"))
include(joinpath("actors", "blanket_actors.jl"))
include(joinpath("actors", "balance_of_plant_actors.jl"))
include(joinpath("actors", "current_actors.jl"))
include(joinpath("actors", "divertors_actors.jl"))
include(joinpath("actors", "sources_actors.jl"))
include(joinpath("actors", "transport_actors.jl"))
include(joinpath("actors", "costing_actors.jl"))
include(joinpath("actors", "compound_actors.jl"))
include(joinpath("actors", "neutronics_actors.jl"))

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
export IMAS, @ddtime, constants, ±, ↔, global_logger

end
