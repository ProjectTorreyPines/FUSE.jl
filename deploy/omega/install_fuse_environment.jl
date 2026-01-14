@assert ("FUSE_ENVIRONMENT" in keys(ENV)) "Error: Must define FUSE_ENVIRONMENT environment variable"
fuse_env = ENV["FUSE_ENVIRONMENT"]
env_dir = joinpath(ENV["FUSE_HOME"], "environments", fuse_env)
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
Pkg.activate(env_dir)
Pkg.add([["FUSE", "Plots", "IJulia", "WebIO", "Interact", "EFIT", "ArgParse"]; packages])
Pkg.build("IJulia")
Pkg.build("WebIO")

println()
println("### Freeze Project and Manifest to read only")
chmod(joinpath(env_dir, "Project.toml"),  0o444)
chmod(joinpath(env_dir, "Manifest.toml"), 0o444)

println()
println("### Create precompile script")
precompile_execution_file = joinpath(env_dir, "precompile_script.jl")
precompile_cmds = """
using WebIO
using FUSE, EFIT, $pkgs_using
include(joinpath(pkgdir(FUSE), "docs", "src", "tutorial.jl"))
include(joinpath(pkgdir(FUSE), "test", "runtests.jl"))
include(joinpath(pkgdir(FUSE), "deploy", "omega", "time_dependent_d3d.jl")
"""
write(precompile_execution_file, precompile_cmds)
chmod(precompile_execution_file, 0o444)

println()
println("### Precompile FUSE sys image")
sysimage_path = joinpath(env_dir, "sys_fuse.so")
create_sysimage(["FUSE", "IJulia", "WebIO", "Interact", "Plots"];
                project=env_dir,
                sysimage_path,
                precompile_execution_file,
                cpu_target)

chmod(sysimage_path, 0o555)

println()
println("### Create IJulia kernels")
import IJulia

#===
We're putting IJulia, WebIO, and Interact into the sysimage now
This causes issues with `import WebIO` not seeing jupyter
  when the sysimage is loaded
The solution is to call `WebIO.__init__()` at the beginning of the kernel
This does some fancy stuff to print warning from this to the terminal
  instead of inside the notebook where it may confuse users
===#
preload_webio_commands = """const __WEBIO_INITED__ = Ref(false)

try
    using IJulia, Logging
    IJulia.push_preexecute_hook(() -> begin
        if !__WEBIO_INITED__[]
            @info "WebIO automatically reinitialized for Julia+FUSE sysimage"
            term_logger = ConsoleLogger(IJulia.orig_stderr[], Logging.Warn)
            with_logger(term_logger) do
                # send warnings/errors to the terminal, not the notebook
                redirect_stderr(IJulia.orig_stderr[]) do
                    @eval import WebIO
                    WebIO.__init__()    # re-register provider quietly for the notebook
                end
            end
            __WEBIO_INITED__[] = true
        end
        nothing
    end)
catch e
    @warn "preload_webio failed" exception=(e, catch_backtrace())
end
"""
preload_webio_file = joinpath(env_dir, ".jupyter", "preload_webio.jl")
write(preload_webio_file, preload_webio_commands)
chmod(preload_webio_file, 0o444)
IJulia.installkernel("Julia+FUSE - single thread",
                     "--project=$env_dir",
                     "--sysimage=$sysimage_path",
                     "--load=$preload_webio_file";
                     env=Dict("JULIA_NUM_THREADS"=>"1"))
IJulia.installkernel("Julia+FUSE - 16-thread (medium)",
                     "--project=$env_dir",
                     "--sysimage=$sysimage_path",
                     "--load=$preload_webio_file";
                     env=Dict("JULIA_NUM_THREADS"=>"16"))
IJulia.installkernel("Julia+FUSE - 10-thread (long)",
                     "--project=$env_dir",
                     "--sysimage=$sysimage_path",
                     "--load=$preload_webio_file";
                     env=Dict("JULIA_NUM_THREADS"=>"10"))

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
exe_file = joinpath(env_dir, "fuse")
write(exe_file, fuse_banner * fuse_exe)
chmod(exe_file, 0o555)


println()
println("### Create module file")
module_file = joinpath(ENV["FUSE_HOME"], "modules", "fuse", fuse_env * ".lua")
header = """
local basedir = "/fusion/projects/codes/julia/fuse"
local fuse_env = "$fuse_env"

help([[
Module for julia with FUSE $fuse_env sysimage
Automatically created by FUSE install script:
  `julia <FUSE.jl git repo>/deploy/omega/install_fuse_environment.jl`
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

# set the cpu target to the one defined in the environment
base = replace(base, """setenv("JULIA_CPU_TARGET", "generic")""" => """setenv("JULIA_CPU_TARGET", "$cpu_target")""")

write(module_file, header * base)
