# FUSE in a Singularity container on GA's omega cluster

This directory builds the same self-contained FUSE image as
[`../perlmutter-container`](../perlmutter-container) (FUSE.jl, all
ProjectTorreyPines packages, and a precompiled `sys_fuse.so` sysimage for
fast startup) and ships it as a **Singularity SIF** — a single read-only file
on `/fusion` that any omega node can run. The pipeline is:

1. rootless `podman build` on a worker node (storage on node-local
   `/local-scratch`, since NFS `$HOME` is slow and quota-bound for images),
2. `podman save` → `singularity build` to a SIF,
3. run everywhere with `singularity exec/run` (module `singularity/3.11.3`).

The image is built from the shared
[`Containerfile`](../perlmutter-container/Containerfile) with
`--build-arg JULIA_CPU_TARGET="generic;cascadelake,...;znver2,..."` (same
targets as the omega module deploy). The compute fleet is all AMD — EPYC 7502
(znver2, unflagged nodes) and EPYC 7513 (znver3, feature flag `amd`), both of
which run the znver2 clone — while the cascadelake target covers Intel
login/head nodes and `generic` catches anything else.

## Quick start (use the published image — no build needed)

A prebuilt image and the `squashfuse` helper live at:

```
/fusion/projects/dt/fuse_containers/            # shared with the digital_twin group
├── fuse_<version>.sif
└── bin/squashfuse
```

Run the FUSE REPL from any omega node:

```bash
module load singularity/3.11.3
export SIF_DIR=/fusion/projects/dt/fuse_containers
export PATH=$SIF_DIR/bin:$PATH                  # squashfuse, for fast SIF mounting
singularity run --cleanenv --sif-fuse --bind /fusion $SIF_DIR/fuse_<version>.sif
```

Install the Jupyter kernel for the published image:

```bash
SIF=/fusion/projects/dt/fuse_containers/fuse_<version>.sif \
  ./deploy/omega-container/install_kernel.sh
```

The rest of this README covers building a new image and the full test flow.

## Contents

| File | Purpose |
|------|---------|
| `build.sh` | podman build → SIF export (uses `../perlmutter-container/Containerfile`). |
| `install_kernel.sh` | Installs a Jupyter kernelspec wrapping `singularity exec`. |
| `kernel.json.template` | Template for the kernelspec. |
| `test_slurm.sbatch` | Smoke test as a Slurm job (submit per architecture). |
| `test_kernel_headless.py` | Headless Jupyter-kernel smoke test via `jupyter_client`. |

## Prerequisite: squashfuse (fast SIF startup)

Omega's singularity install is unprivileged (no setuid starter) and the nodes
have no `squashfuse`, so by default singularity **extracts** the SIF to a
temporary sandbox on every run — minutes of overhead for an image this size.
With `squashfuse` on `PATH` and the `--sif-fuse` flag, the SIF is FUSE-mounted
directly and startup takes seconds.

A shared build lives next to the published SIFs (`<sif-dir>/bin/squashfuse`);
all scripts in this directory pick it up automatically. To (re)build it:

```bash
git clone https://github.com/vasi/squashfuse.git && cd squashfuse
./autogen.sh && ./configure --prefix=<sif-dir> && make && make install
```

## 1. Build the image

The build runs the FUSE test suite + tutorial as the sysimage precompile
workload, so it is heavy (hours). Run it from an interactive allocation
(`short` is capped at 30 min — use `medium`), or directly on an idle worker
node. **podman storage is node-local**: run every step on the same node.

```bash
salloc -p medium -t 8:00:00 -c 16
module load singularity/3.11.3

cd <FUSE repo root>
SIF_DIR=/fusion/ga/projects/ird/ptp/$USER/fuse_containers ./deploy/omega-container/build.sh
```

`build.sh` will:

1. Resolve the image tag (`fuse:<version>`) from the latest FUSE.jl release
   (override with `FUSE_ENVIRONMENT=v2.6.0`).
2. `podman build` with omega's dual CPU target, storing images on
   `/local-scratch/$USER` (no NFS, no persistent config).
3. `podman save` → `singularity build` → `SIF_DIR/fuse_<version>.sif`
   (default `SIF_DIR` is node-local `/local-scratch/$USER`; point it at a
   `/fusion` path to make the SIF visible cluster-wide).

## 2. Run the FUSE REPL

```bash
module load singularity/3.11.3
export PATH=<sif-dir>/bin:$PATH   # squashfuse, see above

singularity run --cleanenv --sif-fuse --bind /fusion <sif-dir>/fuse_<version>.sif
```

This launches a Julia REPL with the FUSE sysimage preloaded (via the
in-image `fuse` wrapper — the container's entrypoint). Notes on flags:

- `--cleanenv` keeps host Julia variables (`JULIA_PROJECT`,
  `JULIA_DEPOT_PATH`, ... from the module-based workflow) out of the
  container. Pass runtime knobs explicitly, e.g.
  `--env JULIA_NUM_THREADS=8`.
- `--bind /fusion` exposes the shared filesystems (`$HOME`, `/tmp`, and the
  current directory are bound by default).
- The SIF is fully read-only, and Julia package init **crashes** if it cannot
  write logs/scratch to a depot. The `fuse` wrapper handles this: when the
  baked depot is not writable it prepends `$HOME/.julia_fuse_container` as a
  writable first depot. Always run through `fuse` (not bare `julia`), or set
  `--env JULIA_DEPOT_PATH=$HOME/.julia_fuse_container:/opt/fuse/.julia`
  yourself.

> Use `singularity run` only for the bare interactive REPL. For one-liners and
> scripts use `singularity exec ... <sif> fuse ...`: singularity `run`
> re-word-splits its arguments (docker-CMD emulation `eval`s them), which
> mangles quoted `-e '...'` code.

Quick smoke test:

```bash
singularity exec --cleanenv --sif-fuse --bind /fusion <sif> \
  fuse -e 'using FUSE; ini,act=FUSE.case_parameters(:D3D,:L_mode); println("FUSE OK: ", pkgversion(FUSE))'
```

Offline/self-contained test (unprivileged singularity cannot drop the network
like podman's `--network none`; `JULIA_PKG_OFFLINE` plus watching for
`Downloading` lines is the omega equivalent):

```bash
singularity exec --cleanenv --sif-fuse --bind /fusion \
  --env JULIA_PKG_OFFLINE=true <sif> \
  fuse -e 'using FUSE; ini,act=FUSE.case_parameters(:D3D,:L_mode); dd=FUSE.init(ini,act); println("OFFLINE OK")'
```

## 3. Slurm jobs

Use the same `singularity exec` line inside `sbatch`/`srun`. The provided
test job validates the image on both compute-node generations:

```bash
SIF=<sif-dir>/fuse_<version>.sif
sbatch --export=ALL,SIF=$SIF -C amd deploy/omega-container/test_slurm.sbatch
sbatch --export=ALL,SIF=$SIF \
       --exclude="$(sinfo -h -N -p short -o '%N %f' | awk '$2=="amd"{print $1}' | paste -sd,)" \
       deploy/omega-container/test_slurm.sbatch
```

Each job should print the in-container sysimage path, a `using FUSE` time of
seconds (not minutes), the node's `Sys.CPU_NAME` (`znver2` on both — the
newer EPYC 7513 nodes also select the znver2 sysimage clone), and end with
`SLURM SMOKE OK`.

## 4. Interactive use via Jupyter

```bash
module load singularity/3.11.3
SIF=<sif-dir>/fuse_<version>.sif ./deploy/omega-container/install_kernel.sh
# optionally: THREADS=8 SIF=... ./install_kernel.sh
```

This writes `$HOME/.local/share/jupyter/kernels/fuse-<version>/kernel.json`,
whose argv runs the in-container Julia via `singularity exec` (with the
absolute singularity path baked in, since Jupyter cannot `module load`).
Singularity's default `$HOME` bind makes the Jupyter connection file visible
inside the container.

Test it headless (uses `jupyter_client`, available in the omega FUSE conda
environment — see [`docs/src/install_omega.md`](../../docs/src/install_omega.md)):

```bash
python3 deploy/omega-container/test_kernel_headless.py fuse-<version>
```

Expected: `KERNEL OK <version>` and `HEADLESS KERNEL TEST PASSED`. In a real
Jupyter session (SSH-tunnel workflow from `install_omega.md`), select the
**Julia FUSE-<version> container** kernel.

> If your Jupyter writes connection files under `$XDG_RUNTIME_DIR` instead of
> `$HOME/.local/share/jupyter/runtime`, set
> `JUPYTER_RUNTIME_DIR=$HOME/.local/share/jupyter/runtime` before starting
> Jupyter so the container can see them.

## 5. Mounting data

Singularity binds `$HOME`, `/tmp`, and the current directory by default; add
`--bind /fusion` for the shared project areas. Bind anything else the same
way:

```bash
singularity exec --cleanenv --sif-fuse --bind /fusion,/local-scratch <sif> ...
```

## 6. Publishing a shared image

Images are published to the digital_twin project area
`/fusion/projects/dt/fuse_containers` (writable by the `digital_twin` group).
To publish a newly built image, copy the SIF and the `bin/squashfuse` helper
there and make them world-readable:

```bash
dest=/fusion/projects/dt/fuse_containers
mkdir -p $dest/bin
cp <sif-dir>/fuse_<version>.sif $dest/
cp <sif-dir>/bin/squashfuse     $dest/bin/
chmod a+rX -R $dest
```

Users then run with `$dest/fuse_<version>.sif` and get `squashfuse` from
`$dest/bin` automatically (all scripts here look next to the SIF).
