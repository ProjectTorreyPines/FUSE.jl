#%Module1.0
## Base module file for FUSE on Kairos
## This file is a template that gets processed by install_fuse_environment.jl

proc ModulesHelp { } {
    puts stderr "Julia FUSE environment __FUSE_VERSION_TAG__ with precompiled systemimage"
    puts stderr "Documentation: https://fuse.help"
    puts stderr "To start: fuse"
}

module-whatis "Name: FUSE"
module-whatis "Version: __FUSE_VERSION_TAG__"  
module-whatis "Description: Julia FUSE environment with precompiled systemimage"
module-whatis "Maintainer: FUSE Development Team"

# Dependencies - load Julia if not already loaded
if { ![is-loaded julia] } {
    module load julia
}

# Set FUSE directories
set fuse_home "__FUSE_HOME_PATH__"
set fuse_env "__FUSE_VERSION_TAG__"
set env_dir "${fuse_home}/environments/${fuse_env}"
set base_depot "${env_dir}/.julia"

# Core FUSE environment variables
setenv FUSE_HOME "${fuse_home}"
setenv FUSE_ENVIRONMENT "${fuse_env}"

# Julia depot path management
# We put the user depot first so their own packages get installed there,
# then the FUSE environment's depot after so it can find packages for the
# precompiled sysimage
if { [info exists env(JULIA_DEPOT_PATH)] } {
    # Case: JULIA_DEPOT_PATH is already set
    set user_depot $env(JULIA_DEPOT_PATH)
    if { [string match "*:" $user_depot] } {
        # Case: ends with ":", append base_depot
        setenv JULIA_DEPOT_PATH "${user_depot}${base_depot}:"
    } else {
        # Case: doesn't end with ":", add separator
        setenv JULIA_DEPOT_PATH "${user_depot}:${base_depot}:"
    }
} else {
    # Case: JULIA_DEPOT_PATH is unset, use default HOME user depot
    setenv JULIA_DEPOT_PATH "$env(HOME)/.julia:${base_depot}:"
}

# The FUSE sysimage environment is the last place julia looks for packages
# when a user does `using <package>`, but this allows Julia to automatically
# find FUSE, Plots
setenv JULIA_LOAD_PATH ":${env_dir}"

# Compiler settings
setenv JULIA_CC "gcc -O3"

# CPU target (will be replaced by actual target during installation)
setenv JULIA_CPU_TARGET "generic"

# Add FUSE executable to PATH
prepend-path PATH "${env_dir}"

# Load message
if { [module-info mode load] } {
    puts stderr "FUSE ${fuse_env} environment loaded"
    puts stderr "Documentation: https://fuse.help"
    puts stderr "To start: fuse"
}
