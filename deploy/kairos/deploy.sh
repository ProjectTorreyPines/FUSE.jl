#!/bin/bash

# Load required modules
module use /opt/nfri/glib/modulefiles
module purge
module load PrgEnv-gnu
module swap gcc/8.3.0 gcc/12.1.0
module load julia

# Check if modules are loaded and display versions
echo "Checking module availability and versions..."

if ! command -v julia &> /dev/null; then
    echo "ERROR: Julia module not loaded or not available"
    exit 1
fi

if ! command -v gcc &> /dev/null; then
    echo "ERROR: GCC module not loaded or not available"
    exit 1
fi

echo "Julia version: $(julia --version)"
echo "GCC version: $(gcc --version | head -n1)"
echo "All required modules loaded successfully"
echo ""

# Kairos-specific proxy configuration
export http_proxy="http://203.230.124.39:4128"
export https_proxy="http://203.230.124.39:4128"
export no_proxy="localhost,127.0.0.1,*.hpc.nfri.re.k"

# User installation directories - use environment variables with fallback defaults
basedir="${FUSE_HOME:-/opt/nfri/glib/fuse}"
module_dir="${FUSE_MODULE_DIR:-/opt/nfri/glib/modulefiles}"

# Get FUSE version - use environment variable with fallback to Project.toml version
scriptdir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
project_version=$(grep '^version = ' "$scriptdir/../../Project.toml" | sed 's/version = "\(.*\)"/v\1/')
fuse_env="${FUSE_ENVIRONMENT:-$project_version}" # Falls back to version from Project.toml if FUSE_ENVIRONMENT not set

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
export JULIA_DOWNLOADS_TIMEOUT=600  # 10 minutes timeout
# Git configuration for proxy environment and network instability
git config --global http.postBuffer 1048576000  # 1GB buffer
git config --global core.compression 0
git config --global http.version HTTP/1.1
git config --global http.proxy "http://203.230.124.39:4128"
git config --global https.proxy "http://203.230.124.39:4128"
git config --global pack.windowMemory 256m      # Reduce memory usage
git config --global pack.packSizeLimit 2g       # 2GB pack limit
git config --global http.maxRequestBuffer 100m  # 100MB request buffer
export JULIA_CC="gcc -O3"
export JULIA_NUM_THREADS=1

# CPU target for Kairos (Intel Xeon Platinum 8260 - Cascade Lake)
# Use cascadelake for optimal performance, with generic fallback for compatibility
export JULIA_CPU_TARGET="generic;cascadelake,-xsaveopt,-rdrnd,clone_all"

# Export module directory for Julia script
export KAIROS_MODULE_DIR="$module_dir"

# Capture loaded module versions for template substitution
echo "Capturing module versions for template..."
julia_version=$(module list julia 2>&1 | grep -oP 'julia/\K[0-9]+\.[0-9]+\.[0-9]+' | head -1)
gcc_version=$(module list gcc 2>&1 | grep -oP 'gcc/\K[0-9]+\.[0-9]+\.[0-9]+' | head -1)

echo "Julia version: $julia_version"
echo "GCC version: $gcc_version"

# Store versions for later use in template processing
export JULIA_MODULE_VERSION="$julia_version"
export GCC_MODULE_VERSION="$gcc_version"

julia "$scriptdir/install_fuse_environment.jl"

unset http_proxy
unset https_proxy
unset no_proxy