#!/bin/bash -v
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

installdir="$SCRATCH/fuse/environments/$fuse_env"
[ -d "$installdir" ] && rm -rf "$installdir"

export FUSE_HOME="$basedir"
export FUSE_INSTALL_DIR="$installdir"
export FUSE_ENVIRONMENT="$fuse_env"
export JULIA_DEPOT_PATH="$installdir/.julia:"
export JULIA_PKG_USE_CLI_GIT=1
export JULIA_CC="gcc -O3"
export JULIA_NUM_THREADS=1

# znver3 for login and cpu & gpu worker nodes,but give generic fallback
# remove xsaveopt since that's what julia distributions do
export JULIA_CPU_TARGET="generic;znver3,-xsaveopt,-rdrnd,clone_all"
export IJULIA_NODEFAULTKERNEL=1
export JUPYTER_DATA_DIR="$installdir/.jupyter"

scriptdir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

srun --nodes 1 --qos interactive --time 04:00:00 --constraint cpu --account m3739 julia $scriptdir/install_fuse_environment.jl

/bin/cp -r $installdir $envdir
/bin/cp $installdir/$fuse_env.lua $basedir/modules/fuse
