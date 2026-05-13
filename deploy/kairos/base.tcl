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

# Dependencies - load required modules (only when loading, not unloading)
if { [module-info mode load] } {
    puts stderr "FUSE module dependency check:"

    # Check if required Julia version is already loaded
    if {[is-loaded julia/__JULIA_VERSION__]} {
        puts stderr "  julia/__JULIA_VERSION__ already loaded"
    } else {
        if {[is-loaded julia]} {
            set loaded_julia [exec sh -c "module list 2>&1 | grep 'julia/' | head -1 | awk '{print \$2}'"]
            puts stderr "  \[Switching\]: $loaded_julia --> julia/__JULIA_VERSION__"
            module unload julia
        } else {
            puts stderr "  Loading julia/__JULIA_VERSION__"
        }
        module load julia/__JULIA_VERSION__
    }

    # Check if required GCC version is already loaded
    if {[is-loaded gcc/__GCC_VERSION__]} {
        puts stderr "  gcc/__GCC_VERSION__ already loaded"
    } else {
        if {[is-loaded gcc]} {
            set loaded_gcc [exec sh -c "module list 2>&1 | grep 'gcc/' | head -1 | awk '{print \$2}'"]
            puts stderr "  \[Switching\]: $loaded_gcc --> gcc/__GCC_VERSION__"
            module unload gcc
        } else {
            puts stderr "  Loading gcc/__GCC_VERSION__"
        }
        module load gcc/__GCC_VERSION__
    }

    puts stderr "  FUSE dependencies ready!"
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
        # Case: doesn't end with ":", user explicitly excluded default depots
        puts stderr "ERROR: Cannot parse existing JULIA_DEPOT_PATH=$user_depot"
        puts stderr "It must be unset or end with ':'"
        puts stderr "A path without trailing ':' means you don't want Julia's default depots."
        puts stderr "Adding our depot would change this behavior."
        exit 1
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
