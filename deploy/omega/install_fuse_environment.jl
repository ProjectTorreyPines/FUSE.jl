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
Pkg.add([["FUSE", "Plots", "IJulia", "WebIO", "Interact", "EFIT", "ArgParse", "PrecompileTools"]; packages])
Pkg.build("IJulia")
Pkg.build("WebIO")

println()
println("### Disable PrecompileTools workloads")
#===
Belt and braces for the precompilation below: if a cache is invalidated between
here and `create_sysimage` the package falls back to being included from source,
and a `@compile_workload` is the most likely thing to blow up when it does. The
sysimage still gets its precompiled code from the trace further down, which
exercises the paths FUSE actually uses. Drop this if you would rather have the
workloads' coverage.

PrecompileTools must be a direct dependency for the preference to have any
effect (hence its presence in the `Pkg.add` above): `Base.get_uuid_name` maps a
UUID to a LocalPreferences.toml section by searching only the project's own
name, `deps`, `extras` and `weakdeps`, so a section naming an indirect
dependency is silently ignored.
===#
local_preferences_file = joinpath(env_dir, "LocalPreferences.toml")
write(local_preferences_file, """
[PrecompileTools]
precompile_workloads = false
""")

println()
println("### Precompile the environment the way PackageCompiler loads it")
#===
PackageCompiler runs both its tracing script and the `--output-o` process with
`--pkgimages=no` (see `PackageCompiler.get_julia_cmd`), and a cache built with
pkgimages is rejected under that flag. The `--output-o` process is additionally
forbidden from precompiling on demand, so every package that `Pkg.add` cached
but `precompile_script.jl` does not itself load gets *included from source*
there. That runs the package's load-time code in the one process where Julia
defers every module's `__init__`, so JLL libraries are never dlopened:

  - `@compile_workload` blocks that call a JLL -- IJulia starts a real Jupyter
    kernel over ZMQ, HiGHS instantiates a `HiGHS.Optimizer`:
        could not load symbol "zmq_ctx_new" / "Highs_create"
  - top-level module code that calls a JLL -- GPUCompiler builds LLVM IR at
    load time, via CUDA:
        could not load symbol "LLVMGetValueContext"

Precompiling the whole manifest with `--pkgimages=no` gives every package a
cache the sysimage build will accept, so none of that code ever runs. This must
come after LocalPreferences.toml is written, or the preference change would
invalidate what we just built.
===#
run(`$(Base.julia_cmd()) --startup-file=no --pkgimages=no --project=$env_dir -e "using Pkg; Pkg.precompile()"`)

println()
println("### Check the environment is ready for PackageCompiler")
#===
Fail here rather than minutes into `create_sysimage`, which only reports the
first package it trips over and costs a full trace run to reach.
===#
check_script, check_io = mktemp()
write(check_io, """
using Pkg
bad = String[]
for (uuid, entry) in Pkg.Types.Context().env.manifest
    pkgid = Base.PkgId(uuid, entry.name)
    Base.in_sysimage(pkgid) && continue
    Base.isprecompiled(pkgid) || push!(bad, entry.name)
end
print(join(sort(bad), " "))
""")
close(check_io)
source_loaded = readchomp(`$(Base.julia_cmd()) --startup-file=no --pkgimages=no --project=$env_dir $check_script`)
@assert isempty(source_loaded) "these packages would be included from source by create_sysimage: $source_loaded"

# The global switch is checked before the per-package one, so `workload_enabled`
# of any module is `false` exactly when the preference is being honored.
workloads_enabled = readchomp(`$(Base.julia_cmd()) --startup-file=no --pkgimages=no --project=$env_dir -e "import PrecompileTools; print(PrecompileTools.workload_enabled(Base))"`)
@assert workloads_enabled == "false" "PrecompileTools workloads still enabled ($workloads_enabled): check $local_preferences_file"

println()
println("### Freeze Project and Manifest to read only")
chmod(joinpath(env_dir, "Project.toml"),  0o444)
chmod(joinpath(env_dir, "Manifest.toml"), 0o444)
chmod(local_preferences_file, 0o444)

println()
println("### Create precompile script")
precompile_execution_file = joinpath(env_dir, "precompile_script.jl")
precompile_cmds = """
using WebIO
using FUSE, EFIT, $pkgs_using
GC.enable(false)
include(joinpath(pkgdir(FUSE), "docs", "src", "tutorial.jl"))
GC.enable(true)
include(joinpath(pkgdir(FUSE), "test", "runtests.jl"))
include(joinpath(pkgdir(FUSE), "deploy", "omega", "time_dependent_d3d.jl"))
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
mkpath(dirname(preload_webio_file))
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
Maintainers: M.G. Yoo, yoom@fusion.gat.com
             C.M. Clark, clarkm@fusion.gat.com
Physics Officers: J. McClenaghan, mcclenaghanj@fusion.gat.com
                  M.G. Yoo, yoom@fusion.gat.com
Known technical debt:
The first time a custom Jupyter kernel is used, it may hang.
Restarting (sometimes twice) normally resolves the issue.
]])

"""

base = read(joinpath(@__DIR__, "base.lua"), String)

# set the cpu target to the one defined in the environment
julia_version = "$(VERSION.major).$(VERSION.minor).$(VERSION.patch)"
base = replace(base, "JULIA_VERSION" => julia_version)
base = replace(base, """setenv("JULIA_CPU_TARGET", "generic")""" => """setenv("JULIA_CPU_TARGET", "$cpu_target")""")

write(module_file, header * base)
