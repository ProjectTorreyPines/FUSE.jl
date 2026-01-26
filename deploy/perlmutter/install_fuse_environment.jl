@assert (Threads.nthreads() == 1) "Error: Installing FUSE sysimage requires running Julia with one thread"

@assert ("FUSE_ENVIRONMENT" in keys(ENV)) "Error: Must define FUSE_ENVIRONMENT environment variable"
fuse_env = ENV["FUSE_ENVIRONMENT"]
env_dir = joinpath(ENV["FUSE_HOME"], "environments", fuse_env)
install_dir = ENV["FUSE_INSTALL_DIR"]
cpu_target = ENV["JULIA_CPU_TARGET"]

import Pkg

println("### Setup main environment for installer")
Pkg.activate()
Pkg.Registry.add(Pkg.RegistrySpec(url="https://github.com/ProjectTorreyPines/FuseRegistry.jl.git"))
Pkg.Registry.add("General")
Pkg.add("PackageCompiler")
Pkg.update()
using PackageCompiler

println()
println("### parse PTP packages from Makefile")
function get_packages_from_makefile()
    makefile = joinpath(@__DIR__, "..", "..", "Makefile")
    for line in eachline(makefile)
        if occursin(r"^FUSE_PACKAGES_MAKEFILE\s*:=", line)
            # Remove the variable name and the ':=' operator
            pkg_line = replace(line, r"^FUSE_PACKAGES_MAKEFILE\s*:=" => "")
            # Trim any leading/trailing whitespace
            return split(strip(pkg_line))
        end
    end
end
packages = get_packages_from_makefile()
pkgs_using = join(packages, ", ")
println("    ", packages)

println()
println("### Setup new environment")
Pkg.activate(install_dir)
Pkg.add([["FUSE", "Plots", "IJulia", "WebIO", "Interact", "EFIT", "ArgParse"]; packages])
Pkg.build("IJulia")

println()
println("### Freeze Project and Manifest to read only")
chmod(joinpath(install_dir, "Project.toml"),  0o444)
chmod(joinpath(install_dir, "Manifest.toml"), 0o444)

println()
println("### Create precompile script")
precompile_execution_file = joinpath(install_dir, "precompile_script.jl")
precompile_cmds = """
using FUSE, EFIT, $pkgs_using
include(joinpath(pkgdir(FUSE), "docs", "src", "tutorial.jl"))
include(joinpath(pkgdir(FUSE), "test", "runtests.jl"))
"""
write(precompile_execution_file, precompile_cmds)
chmod(precompile_execution_file, 0o444)

println()
println("### Precompile FUSE sys image")
sysimage_path = joinpath(install_dir, "sys_fuse.so")
create_sysimage(["FUSE"]; sysimage_path, precompile_execution_file, cpu_target)
chmod(sysimage_path, 0o555)

println()
println("### Create IJulia kernels")
module_folder = joinpath(ENV["FUSE_HOME"], "modules")
wrapper = """
#!/bin/bash
command -v module >/dev/null 2>&1 || source /usr/share/lmod/lmod/init/bash
module use $module_folder
module load fuse/$fuse_env
exec julia "\$@"
"""
println(wrapper)
wrap_file = joinpath(install_dir, "fuse_kernel.sh")
write(wrap_file, wrapper)
chmod(wrap_file, 0o555)

import IJulia
IJulia.installkernel("Julia FUSE-$fuse_env 1 Thread(s)",  "--sysimage=$sysimage_path";
                     julia = `$(joinpath(env_dir, "fuse_kernel.sh"))`,
                     env=Dict("JULIA_NUM_THREADS"=>"1"))
IJulia.installkernel("Julia FUSE-$fuse_env 8 Thread(s)", "--sysimage=$sysimage_path";
                     julia = `$(joinpath(env_dir, "fuse_kernel.sh"))`,
                     env=Dict("JULIA_NUM_THREADS"=>"8"))

println()
println("### Create fuse executable")
fuse_banner = raw"""
#!/bin/bash

# ANSI color codes
RESET="\033[0m"
BOLD="\033[1m"
BLUE="\033[34m"
RED="\033[31m"
GREEN="\033[32m"
PURPLE="\033[35m"

echo -e "  ${BOLD}${GREEN}_${RESET}  __               ${BOLD}${PURPLE}_${RESET} ${BOLD}${BLUE}_${RESET}
 ${BOLD}${GREEN}(_)${RESET}/ _|             ${BOLD}${PURPLE}(_${RESET}${BOLD}${BLUE}(_)${RESET} |  Documentation: https://fuse.help
${BOLD}${RED}(_)${RESET}| |_ _   _ ___  _${BOLD}${BLUE}(_${RESET}${BOLD}${PURPLE}(_)${RESET}  |
   |  _| | | / __|/ _ \    |  Julia REPL with a FUSE """ * "$fuse_env" * raw""" sysimage.
  ${BOLD}${RED}_${RESET}| | | |_| \__ \  __/    |  Warning: Default environment is read-only.
 ${BOLD}${RED}(_)${RESET}_|  \__,_|___/\___|${BOLD}${BLUE}_${RESET}   |           Adding packages in new environment may
${BOLD}${GREEN}(_(_)${RESET}                 ${BOLD}${BLUE}(_)${RESET}  |           cause conflicts or unexpected behavior.
"


"""
fuse_exe = "JULIA_PROJECT=$env_dir julia -i --banner=no --sysimage=$(env_dir)/sys_fuse.so \$@"
exe_file = joinpath(install_dir, "fuse")
write(exe_file, fuse_banner * fuse_exe)
chmod(exe_file, 0o555)


println()
println("### Create module file")
module_file = joinpath(install_dir, fuse_env * ".lua")
header = """
local basedir = "/global/common/software/m3739/perlmutter/fuse"
local fuse_env = "$fuse_env"

help([[
Module for julia with FUSE $fuse_env sysimage
Automatically created by FUSE install script:
  `julia <FUSE.jl git repo>/deploy/perlmutter/install_fuse_environment.jl`
Maintainers: B.C. Lyons, lyonsbc@fusion.gat.com
Physics Officers: O.M. Meneghini, meneghini@fusion.gat.com
                  B.C. Lyons, lyonsbc@fusion.gat.com
]])

"""

base = read(joinpath(@__DIR__, "base.lua"), String)

# set the cpu target to the one defined in the environment
base = replace(base, """setenv("JULIA_CPU_TARGET", "generic")""" => """setenv("JULIA_CPU_TARGET", "$cpu_target")""")

write(module_file, header * base)

