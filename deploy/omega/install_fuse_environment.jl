@assert (Threads.nthreads() == 1) "Error: Installing FUSE sysimage requires running Julia with one thread"

@assert ("FUSE_ENVIRONMENT" in keys(ENV)) "Error: Must define FUSE_ENVIRONMENT environment variable"
fuse_env = ENV["FUSE_ENVIRONMENT"]
env_dir = joinpath(ENV["FUSE_HOME"], "environments", fuse_env)

import Pkg

# Setup main environment for installer
Pkg.activate()
Pkg.Registry.add(Pkg.RegistrySpec(url="https://github.com/ProjectTorreyPines/FuseRegistry.jl.git"))
Pkg.Registry.add("General")
Pkg.add("PackageCompiler")
using PackageCompiler

# Setup new environment
Pkg.activate(env_dir)
Pkg.add(["FUSE", "Plots", "IJulia", "WebIO", "Interact"])
Pkg.build("IJulia")

# Freeze Project and Manifest to read only
chmod(joinpath(env_dir, "Project.toml"),  0o444)
chmod(joinpath(env_dir, "Manifest.toml"), 0o444)

# Create precompile script
precompile_execution_file = joinpath(env_dir, "precompile_script.jl")
precompile_cmds = """
using FUSE
include(joinpath(pkgdir(FUSE), "docs", "src", "tutorial.jl"))
include(joinpath(pkgdir(FUSE), "test", "runtests.jl"))
"""
write(precompile_execution_file, precompile_cmds)

# Precompile FUSE sys image
sysimage_path = joinpath(env_dir, "sys_fuse.so")
cpu_target = ENV["JULIA_CPU_TARGET"]
create_sysimage(["FUSE"]; sysimage_path, precompile_execution_file, cpu_target)

# Create IJulia kernels (10 threads for login, 40 for worker)
import IJulia
IJulia.installkernel("Julia+FUSE (login 10-threads)",  "--sysimage=$sysimage_path"; env=Dict("JULIA_NUM_THREADS"=>"10"))
IJulia.installkernel("Julia+FUSE (worker 40-threads)", "--sysimage=$sysimage_path"; env=Dict("JULIA_NUM_THREADS"=>"40"))

# Create module file
module_file = joinpath(ENV["FUSE_HOME"], "modules", "fuse", fuse_env * ".lua")
header = """
local basedir = "/fusion/projects/codes/julia/fuse"
local fuse_env = "$fuse_env"

help([[
Module for julia with FUSE $fuse_env sysimage
Automatically created by FUSE install script:
  `julia /fusion/projects/codes/julia/fuse/install/install_fuse_environment.jl`
Maintainers: B.C. Lyons, lyonsbc@fusion.gat.com
             C.M. Clark, clarkm@fusion.gat.com
Physics Officers: O.M. Meneghini, meneghini@fusion.gat.com
                  B.C. Lyons, lyonsbc@fusion.gat.com
Known technical debt: 
The first time a custom Jupyter kernel is used, it may hang.
Restarting (sometimes twice) normally resolves the issue. 
]])

"""

base = read(joinpath(@__DIR__, "base.lua"), String)

write(module_file, header * base)
