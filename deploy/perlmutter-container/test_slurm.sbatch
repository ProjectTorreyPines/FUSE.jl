#!/bin/bash
#SBATCH --job-name=fuse-container-test
#SBATCH --account=m3739
#SBATCH --constraint=cpu
#SBATCH --qos=debug
#SBATCH --time=0:25:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --output=fuse_container_test_%j.out
#
# Smoke-test the migrated FUSE podman-hpc image inside a Slurm job on a
# Perlmutter CPU node (AMD Milan, znver3 — one of the sysimage's optimized
# clone targets). Run Phase 1 of the README first (pull/build + migrate), then:
#
#   FUSE_ENVIRONMENT=v1.2.0 sbatch --export=ALL deploy/perlmutter-container/test_slurm.sbatch
#
# Optional: SQUASH_DIR=/global/cfs/cdirs/m3739/shared_images to test the shared
# image store instead of the per-user one.
#
# Pass: prints the in-container sysimage path, Sys.CPU_NAME reports znver3,
# `using FUSE` loads in seconds, the D3D L-mode flux matcher completes, ends
# with SLURM SMOKE OK.

set -euo pipefail

if [[ -n "${FUSE_ENVIRONMENT:-}" ]]; then
    version="$FUSE_ENVIRONMENT"
else
    version="$(curl -s https://api.github.com/repos/ProjectTorreyPines/FUSE.jl/releases/latest | jq -r .name)"
fi
[[ -n "$version" && "$version" != "null" ]] || { echo "ERROR: set FUSE_ENVIRONMENT" >&2; exit 1; }

podman=(podman-hpc)
[[ -n "${SQUASH_DIR:-}" ]] && podman=(podman-hpc --squash-dir "$SQUASH_DIR")

echo "node: $(hostname)"
lscpu | grep 'Model name'

# The in-image `fuse` wrapper loads the sysimage and, on the read-only squash
# rootfs, prepends a writable $HOME depot.
"${podman[@]}" run --rm -e JULIA_NUM_THREADS="$SLURM_CPUS_PER_TASK" \
    "localhost/fuse:$version" fuse -e '
        println("sysimage: ", unsafe_string(Base.JLOptions().image_file))
        println("cpu: ", Sys.CPU_NAME, "  threads: ", Threads.nthreads())
        @time using FUSE
        ini, act = FUSE.case_parameters(:D3D, :L_mode)
        dd = FUSE.init(ini, act)
        FUSE.ActorFluxMatcher(dd, act)
        println("SLURM SMOKE OK  FUSE ", pkgversion(FUSE))'
