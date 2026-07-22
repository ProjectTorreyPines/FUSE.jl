#!/bin/bash
#SBATCH --job-name=fuse-container-test
#SBATCH --partition=short
#SBATCH --time=0:25:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --output=fuse_container_test_%j.out
#
# Smoke-test the FUSE SIF inside a Slurm job. Submit once per node generation
# to cover both CPU types in the omega compute fleet — all AMD: the 'amd'
# feature flag marks the newer EPYC 7513 (znver3) nodes, the unflagged ones
# are EPYC 7502 (znver2). (The sysimage's cascadelake target is kept for
# Intel login/head nodes.)
#
#   sbatch --export=ALL,SIF=/fusion/.../fuse_v2.6.0.sif -C amd deploy/omega-container/test_slurm.sbatch
#   sbatch --export=ALL,SIF=/fusion/.../fuse_v2.6.0.sif \
#          --exclude="$(sinfo -h -N -p short -o '%N %f' | awk '$2=="amd"{print $1}' | paste -sd,)" \
#          deploy/omega-container/test_slurm.sbatch
#
# Pass: prints the in-container sysimage path, `using FUSE` loads in seconds,
# Sys.CPU_NAME reports znver2 (or znver3 nodes report znver2 — the closest
# sysimage clone), ends with SLURM SMOKE OK.

set -euo pipefail

: "${SIF:?set SIF=/path/to/fuse_<version>.sif (sbatch --export=ALL,SIF=...)}"

module load singularity/3.11.3

launcher="$(dirname "$SIF")/fuse-container"
[[ -x "$launcher" ]] || launcher="$(dirname "$0")/fuse-container"

echo "node: $(hostname)"
lscpu | grep 'Model name'

# Run through fuse-container: it FUSE-mounts the SIF via squashfuse as a tracked
# child, so a Slurm time-limit / scancel (SIGKILL of the job cgroup) tears the
# mount down instead of orphaning it. The in-image `fuse` wrapper loads the
# sysimage and prepends a writable $HOME depot for the read-only rootfs.
"$launcher" "$SIF" \
    --threads="$SLURM_CPUS_PER_TASK" -e '
        println("sysimage: ", unsafe_string(Base.JLOptions().image_file))
        println("cpu: ", Sys.CPU_NAME, "  threads: ", Threads.nthreads())
        @time using FUSE
        ini, act = FUSE.case_parameters(:D3D, :L_mode)
        dd = FUSE.init(ini, act)
        println("SLURM SMOKE OK  FUSE ", pkgversion(FUSE))'
