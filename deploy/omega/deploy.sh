#!/bin/bash

module load julia
module load env/gcc11.x

basedir="/fusion/projects/codes/julia/fuse"
export PATH="$basedir/miniconda3/bin:$PATH"

fuse_env=`gh release list -R ProjectTorreyPines/FUSE.jl --json name,isLatest --jq ".[] | select(.isLatest)|.name"`

envdir="$basedir/environments/$fuse_env"

if [ -d "$envdir" ]; then
    echo "FUSE $fuse_env already installed in $envdir"
    echo "Exiting"
    exit 0
fi

export FUSE_HOME="$basedir"
export FUSE_ENVIRONMENT="$fuse_env"
export JULIA_DEPOT_PATH="$envdir/.julia:"
export JULIA_PKG_USE_CLI_GIT=1
export JULIA_CC="gcc -O3"
export JULIA_NUM_THREADS=10
export JULIA_NUM_PRECOMPILE_TASKS=10

# cascadelake for Intel login; znver2 for worker; generic fallback
# remove xsaveopt since that's what julia distributions do
export JULIA_CPU_TARGET="generic;cascadelake,-xsaveopt,-rdrnd,clone_all;znver2,-xsaveopt,-rdrnd,clone_all"
export IJULIA_NODEFAULTKERNEL=1
export JUPYTER_DATA_DIR="$envdir/.jupyter"

julia $HOME/src/FUSE.jl/deploy/omega/install_fuse_environment.jl
