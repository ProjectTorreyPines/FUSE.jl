#!/bin/bash
source /usr/share/lmod/lmod/init/bash
module load julia/1.11.4

basedir="/global/common/software/m3739/perlmutter/fuse"

fuse_env=`curl -s https://api.github.com/repos/ProjectTorreyPines/FUSE.jl/releases/latest | jq -r .name`

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
export JULIA_NUM_THREADS=1

# znver3 for login and cpu & gpu worker nodes,but give generic fallback
# remove xsaveopt since that's what julia distributions do
export JULIA_CPU_TARGET="generic;znver3,-xsaveopt,-rdrnd,clone_all"
export IJULIA_NODEFAULTKERNEL=1
export JUPYTER_DATA_DIR="$envdir/.jupyter"

julia $HOME/src/FUSE.jl/deploy/perlmutter/install_fuse_environment.jl
