@assert (Threads.nthreads() == 1) "Error: Installing FUSE sysimage requires running Julia with one thread"

# Check required environment variables
@assert ("FUSE_ENVIRONMENT" in keys(ENV)) "Error: Must define FUSE_ENVIRONMENT environment variable"
@assert ("FUSE_HOME" in keys(ENV)) "Error: Must define FUSE_HOME environment variable"

fuse_env = ENV["FUSE_ENVIRONMENT"]
base_dir = ENV["FUSE_HOME"]
env_dir = joinpath(base_dir, "environments", fuse_env)
cpu_target = get(ENV, "JULIA_CPU_TARGET", "generic")
module_dir = get(ENV, "KAIROS_MODULE_DIR", "/opt/modulefiles")

println("=== FUSE Installation Configuration ===")
println("FUSE Version: $fuse_env")
println("Base Directory: $base_dir")
println("Environment Directory: $env_dir")
println("CPU Target: $cpu_target")
println("Module Directory: $module_dir")
println("========================================")

import Pkg

println("### Setup main environment for installer")
Pkg.activate()
Pkg.Registry.add(Pkg.RegistrySpec(url="https://github.com/ProjectTorreyPines/FuseRegistry.jl.git"))
Pkg.Registry.add("General")
Pkg.add("PackageCompiler")
Pkg.update()
Pkg.instantiate()
using PackageCompiler

println("\n### Parse PTP packages from Makefile")
function get_packages_from_makefile()
    # Look for Makefile in parent directories (more robust)
    search_paths = [pwd(), joinpath(pwd(), "..", ".."), joinpath(@__DIR__, "..", "..")]
    
    for search_dir in search_paths
        makefile = joinpath(search_dir, "Makefile")
        if isfile(makefile)
            println("Found Makefile at: $makefile")
            for line in eachline(makefile)
                if occursin(r"^FUSE_PACKAGES_MAKEFILE\s*:=", line)
                    pkg_line = replace(line, r"^FUSE_PACKAGES_MAKEFILE\s*:=" => "")
                    return split(strip(pkg_line))
                end
            end
        end
    end
    
    # Error if Makefile not found
    error("Makefile not found in any of the search paths: $search_paths\n" *
          "Please ensure you are running this script from the correct directory or " *
          "that the Makefile exists in the project root.")
end
packages = get_packages_from_makefile()
pkgs_using = join(packages, ", ")
println("Packages from Makefile: ", packages)
if isempty(packages)
    error("No packages found in Makefile. Check the FUSE_PACKAGES_MAKEFILE line in Makefile.")
end

println("\n### Setup new environment")
Pkg.activate(env_dir)
Pkg.add([["FUSE", "Plots", "EFIT"]; packages])

println("\n### Freeze Project and Manifest to read only")
chmod(joinpath(env_dir, "Project.toml"),  0o444)
chmod(joinpath(env_dir, "Manifest.toml"), 0o444)

println("\n### Create precompile script")
precompile_execution_file = joinpath(env_dir, "precompile_script.jl")
precompile_cmds = """
using FUSE, EFIT, $pkgs_using
include(joinpath(pkgdir(FUSE), "docs", "src", "tutorial.jl"))
include(joinpath(pkgdir(FUSE), "test", "runtests.jl"))
"""
write(precompile_execution_file, precompile_cmds)
chmod(precompile_execution_file, 0o444)

println("\n### Precompile FUSE sys image")
sysimage_path = joinpath(env_dir, "sys_fuse.so")
create_sysimage(["FUSE"]; sysimage_path, precompile_execution_file, cpu_target)
chmod(sysimage_path, 0o555)

# println("\n### Create IJulia kernels")
# # Create kernel wrapper script for module loading
# kernel_wrapper = joinpath(env_dir, "fuse_kernel.sh")
# wrapper_content = """
# #!/bin/bash
# # FUSE kernel wrapper for Kairos
# command -v module >/dev/null 2>&1 || source /usr/share/lmod/lmod/init/bash
# module use $module_dir
# module load fuse/$fuse_env
# exec julia "\$@"
# """
# write(kernel_wrapper, wrapper_content)
# chmod(kernel_wrapper, 0o755)

# import IJulia
# IJulia.installkernel("Julia FUSE-$fuse_env 1 Thread",  "--sysimage=$sysimage_path";
#     julia=`$kernel_wrapper`, env=Dict("JULIA_NUM_THREADS"=>"1"))
# IJulia.installkernel("Julia FUSE-$fuse_env 4 Threads", "--sysimage=$sysimage_path";
#     julia=`$kernel_wrapper`, env=Dict("JULIA_NUM_THREADS"=>"4"))
# IJulia.installkernel("Julia FUSE-$fuse_env 8 Threads", "--sysimage=$sysimage_path";
#     julia=`$kernel_wrapper`, env=Dict("JULIA_NUM_THREADS"=>"8"))

println("\n### Create fuse executable")
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
   |  _| | | / __|/ _ \    |  Julia REPL with FUSE """ * "$fuse_env" * raw""" sysimage (Kairos)
  ${BOLD}${RED}_${RESET}| | | |_| \__ \  __/    |  System installation - packages read-only
 ${BOLD}${RED}(_)${RESET}_|  \__,_|___/\___|${BOLD}${BLUE}_${RESET}   |  Create user environment for custom packages
${BOLD}${GREEN}(_(_)${RESET}                 ${BOLD}${BLUE}(_)${RESET}  |"

"""
fuse_exe = "JULIA_PROJECT=$env_dir julia -i --banner=no --sysimage=$sysimage_path \$@"
exe_file = joinpath(env_dir, "fuse")
write(exe_file, fuse_banner * fuse_exe)
chmod(exe_file, 0o755)

println("\n### Create module file (optional)")
# Only create module file if module directory is writable and base.tcl template exists
base_tcl_path = joinpath(@__DIR__, "base.tcl")
module_file = joinpath(module_dir, "fuse", "$fuse_env")  # No .lua extension for TCL modules

if isfile(base_tcl_path) && isdir(dirname(module_dir))
    try
        mkpath(dirname(module_file))  # Create module directory if it doesn't exist
        
        # Read base module template and customize
        base_content = read(base_tcl_path, String)
        
        # Replace placeholders with actual values
        customized_content = replace(base_content, "__FUSE_VERSION_TAG__" => fuse_env)
        customized_content = replace(customized_content, "__FUSE_HOME_PATH__" => base_dir)
        customized_content = replace(customized_content, "generic" => cpu_target)
        
        write(module_file, customized_content)
        chmod(module_file, 0o644)
        println("Module file created: $module_file")
    catch e
        println("Warning: Could not create module file: $e")
        println("Module system will not be available, but direct execution will work")
    end
else
    println("Warning: Module file creation skipped (base.tcl not found or module directory not accessible)")
    println("Direct execution will still work: $env_dir/fuse")
end

println("\n### Installation complete!")
println("========================================")
println("FUSE $fuse_env installed successfully!")
println("")
println("To use FUSE:")
println("  module use $module_dir")
println("  module load fuse/$fuse_env")
println("  julia")
println("")
println("Or directly:")
println("  $env_dir/fuse")
println("========================================")
