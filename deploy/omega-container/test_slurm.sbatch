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

# squashfuse (for fast SIF FUSE-mounting) lives next to the SIF in the shared
# containers layout; fall back to extraction if absent.
sif_fuse_flag="--sif-fuse"
if [[ -x "$(dirname "$SIF")/bin/squashfuse" ]]; then
    export PATH="$(dirname "$SIF")/bin:$PATH"
elif ! command -v squashfuse >/dev/null 2>&1; then
    echo "WARNING: squashfuse not found — falling back to SIF extraction (slow)"
    sif_fuse_flag=""
fi

echo "node: $(hostname)"
lscpu | grep 'Model name'

# The in-image `fuse` wrapper runs julia with the sysimage and, on read-only
# images like a SIF, prepends a writable $HOME depot (required — package init
# crashes otherwise).
singularity exec --cleanenv $sif_fuse_flag --bind /fusion \
    "$SIF" \
    fuse --threads="$SLURM_CPUS_PER_TASK" -e '
        println("sysimage: ", unsafe_string(Base.JLOptions().image_file))
        println("cpu: ", Sys.CPU_NAME, "  threads: ", Threads.nthreads())
        @time using FUSE
        ini, act = FUSE.case_parameters(:D3D, :L_mode)
        dd = FUSE.init(ini, act)
        println("SLURM SMOKE OK  FUSE ", pkgversion(FUSE))'
