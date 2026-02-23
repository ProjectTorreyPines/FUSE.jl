#!/bin/bash -v
command -v module >/dev/null 2>&1 || source /usr/share/lmod/lmod/init/bash
module load julia/1.11.4

# Set up OMFIT and OMAS which is needed for the DIII-D study
export FUSE_OMFIT_HOST="localhost"
export FUSE_OMFIT_ROOT="/global/common/software/m3739/perlmutter/OMFIT-CAKE"
export FUSE_OMAS_ROOT="/global/common/software/m3739/perlmutter/FUSE_OMAS"

basedir="/global/common/software/m3739/perlmutter/fuse"

fuse_env=`curl -s https://api.github.com/repos/ProjectTorreyPines/FUSE.jl/releases/latest | jq -r .name`

envdir="$basedir/environments/$fuse_env"

logfile="$PSCRATCH/fuse_build_logs/$fuse_env.log"

# Check if log file already exists (indicates previous run)
if [[ -f "$LOG_FILE" ]]; then
    exit 1
fi

if [ -d "$envdir" ]; then
    echo "FUSE $fuse_env already installed in $envdir"
    echo "Exiting"
    exit 0
fi

send_failure_notification() {
    local exit_code=$1
    
    echo "Julia process failed with exit code: $exit_code"
    echo "Sending last 200 lines of log to $RECIPIENT"
    
    # Create email body
    {
        echo "Julia regression test/sysimage build failed with exit code: $exit_code"
        echo ""
        echo "Last 200 lines of output:"
        echo "=========================="
        tail -n 200 "$logfile"
    } | mail -s "Perlmutter FUSE build failure" "$USER"
}

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
# This can be run as a scrontab so we need to clear SLURM environment variables
unset ${!SLURM_@}

srun --nodes 1 --qos regular --time 04:00:00 --constraint cpu --account m3739 --output $logfile julia $scriptdir/install_fuse_environment.jl
JULIA_EXIT_CODE=$?

if [[ $JULIA_EXIT_CODE -ne 0 ]]; then
    send_failure_notification "$JULIA_EXIT_CODE"
    exit "$JULIA_EXIT_CODE"
else
    /bin/cp -r $installdir $envdir
    /bin/cp $installdir/$fuse_env.lua $basedir/modules/fuse
    rm $basedir/modules/fuse/default.lua
    /bin/ln -s $basedir/modules/fuse/$fuse_env.lua $basedir/modules/fuse/default.lua
    # Remove sysimage and .julia folder of laste environment
    /bin/rm -rf "$(readlink -f $basedir/environments/latest)"/.julia
    /bin/rm -f "$(readlink -f $basedir/environments/latest)"/sys_fuse.so
    # Remove old link
    /bin/rm $basedir/environments/latest
    # Recreate it
    /bin/ln -s $envdir $basedir/environments/latest
fi
