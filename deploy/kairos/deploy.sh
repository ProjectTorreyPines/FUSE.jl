#!/bin/bash

# Load required modules
module load julia

# Kairos-specific proxy configuration
export http_proxy="http://203.230.124.39:4128"
export https_proxy="http://203.230.124.39:4128"
export no_proxy="localhost,127.0.0.1,*.hpc.nfri.re.k"

# User installation directories - use environment variables with fallback defaults
basedir="${FUSE_HOME:-/scratch/yoom/fuse}"
module_dir="${FUSE_MODULE_DIR:-/scratch/yoom/fuse/modulefiles}"

# Get FUSE version - use environment variable with fallback default
fuse_env="${FUSE_ENVIRONMENT:-latest}" # Hardcoded because the api.github.com is not accessible due to network constraints on Kairos

# The following will fetch the latest release dynamically (currently not working due to network constraints).
# fuse_env=$(curl -s https://api.github.com/repos/ProjectTorreyPines/FUSE.jl/releases/latest | jq -r .name)

envdir="$basedir/environments/$fuse_env"

if [ -d "$envdir" ]; then
    echo "FUSE $fuse_env already installed in $envdir"
    echo "Exiting"
    exit 0
fi

# Environment setup
export FUSE_HOME="$basedir"
export FUSE_ENVIRONMENT="$fuse_env"
export JULIA_DEPOT_PATH="$envdir/.julia:"
export JULIA_PKG_USE_CLI_GIT=1
export JULIA_SSL_NO_VERIFY_HOSTS="*"
export JULIA_DOWNLOADS_TIMEOUT=120
# Git configuration for proxy environment
git config --global http.postBuffer 524288000
git config --global http.lowSpeedLimit 1000
git config --global http.lowSpeedTime 600
git config --global core.compression 0
git config --global http.version HTTP/1.1
git config --global http.proxy "http://203.230.124.39:4128"
git config --global https.proxy "http://203.230.124.39:4128"
export JULIA_CC="gcc -O3"
export JULIA_NUM_THREADS=1

# CPU target for Kairos (Intel Xeon Platinum 8260 - Cascade Lake)
# Use cascadelake for optimal performance, with generic fallback for compatibility
export JULIA_CPU_TARGET="generic;cascadelake,-xsaveopt,-rdrnd,clone_all"

# Export module directory for Julia script
export KAIROS_MODULE_DIR="$module_dir"

scriptdir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

julia "$scriptdir/install_fuse_environment.jl"

unset http_proxy
unset https_proxy
unset no_proxy