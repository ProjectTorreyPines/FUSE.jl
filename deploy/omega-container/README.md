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

A prebuilt image, the `squashfuse` helper, and the `fuse-container` launcher
live at:

```
/fusion/projects/dt/fuse_containers/            # shared with the digital_twin group
├── fuse_<version>.sif
├── fuse-container                              # launcher (use this to run)
└── bin/squashfuse
```

Run the FUSE REPL from any omega node:

```bash
module load singularity/3.11.3
export SIF_DIR=/fusion/projects/dt/fuse_containers
$SIF_DIR/fuse-container $SIF_DIR/fuse_<version>.sif
```

Run a script or one-liner (anything after the SIF is passed to Julia):

```bash
$SIF_DIR/fuse-container $SIF_DIR/fuse_<version>.sif \
  -e 'using FUSE; ini,act=FUSE.case_parameters(:D3D,:L_mode); println(pkgversion(FUSE))'
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
| `fuse-container` | Launcher: mounts the SIF and runs it with leak-free cleanup. |
| `install_kernel.sh` | Installs a Jupyter kernelspec that runs via `fuse-container`. |
| `kernel.json.template` | Template for the kernelspec. |
| `test_slurm.sbatch` | Smoke test as a Slurm job (submit per architecture). |
| `test_kernel_headless.py` | Headless Jupyter-kernel smoke test via `jupyter_client`. |

## How the image is run: the `fuse-container` launcher

Omega's singularity is unprivileged (no setuid starter) and the nodes have no
`squashfuse`, so by default singularity **extracts** the whole 5.8 GB SIF to a
temporary sandbox on every run — minutes of overhead. `squashfuse` avoids that
by FUSE-mounting the SIF directly (startup in seconds), but singularity's own
`--sif-fuse` **detaches** the squashfuse process to init: if the caller is then
killed abruptly (a Jupyter kernel shutdown, a Slurm time-limit or `scancel` —
anything that sends SIGKILL), the mount is orphaned and lingers on the node,
accumulating stale processes/mounts and pinning the SIF open on NFS.

`fuse-container` fixes this by mounting the SIF with `squashfuse -f`
(foreground) as a **tracked child in the same process group / cgroup**, so a
group kill reaps it and a normal or signalled exit unmounts it via a trap.
Either way nothing is left behind. Always run the image through
`fuse-container` (the kernel and the Slurm test already do); it needs only
`module load singularity/3.11.3` and finds `squashfuse` next to the SIF.

To (re)build the `squashfuse` helper for a new location:

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
<sif-dir>/fuse-container <sif-dir>/fuse_<version>.sif
```

This launches a Julia REPL with the FUSE sysimage preloaded. `fuse-container`
handles the details:

- It FUSE-mounts the SIF with tracked cleanup (see the launcher section above)
  and runs it with `--cleanenv` so host Julia variables (`JULIA_PROJECT`,
  `JULIA_DEPOT_PATH`, ... from the module workflow) don't leak in.
- `/fusion` is bound (plus `$HOME`, `/tmp`, `$PWD` by default). Bind more via
  `FUSE_BIND=/fusion,/local-scratch`.
- The SIF is fully read-only, and Julia package init **crashes** if it cannot
  write logs/scratch to a depot. The in-image `fuse` entrypoint handles this:
  it prepends `$HOME/.julia_fuse_container` as a writable first depot.

Quick smoke test (everything after the SIF is passed to Julia):

```bash
<sif-dir>/fuse-container <sif> \
  -e 'using FUSE; ini,act=FUSE.case_parameters(:D3D,:L_mode); println("FUSE OK: ", pkgversion(FUSE))'
```

Offline/self-contained test (unprivileged singularity cannot drop the network
like podman's `--network none`; `JULIA_PKG_OFFLINE` plus watching for
`Downloading` lines is the omega equivalent):

```bash
JULIA_PKG_OFFLINE=true <sif-dir>/fuse-container <sif> \
  -e 'using FUSE; ini,act=FUSE.case_parameters(:D3D,:L_mode); dd=FUSE.init(ini,act); println("OFFLINE OK")'
```

## 3. Slurm jobs

Run the image through `fuse-container` inside `sbatch`/`srun` (so a time-limit
or `scancel` tears the mount down instead of orphaning it). The provided test
job validates the image on both compute-node generations:

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
whose argv runs the in-container Julia through `fuse-container` (with the
absolute singularity bin dir baked onto PATH, since Jupyter cannot
`module load`). Because the kernel runs via the launcher, shutting the kernel
down — gracefully or by JupyterHub killing it — leaves no orphaned mount.
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

`$HOME`, `/tmp`, `$PWD`, and `/fusion` are available by default. Bind more with
`FUSE_BIND` (comma-separated), e.g. to also expose node-local scratch:

```bash
FUSE_BIND=/fusion,/local-scratch <sif-dir>/fuse-container <sif> ...
```

## 6. Publishing a shared image

Images are published to the digital_twin project area
`/fusion/projects/dt/fuse_containers` (writable by the `digital_twin` group).
To publish a newly built image, copy the SIF, the `fuse-container` launcher,
and the `bin/squashfuse` helper there and make them world-readable:

```bash
dest=/fusion/projects/dt/fuse_containers
mkdir -p $dest/bin
cp <sif-dir>/fuse_<version>.sif                     $dest/
cp deploy/omega-container/fuse-container            $dest/
cp <sif-dir>/bin/squashfuse                         $dest/bin/
chmod a+rX -R $dest
```

Users then run `$dest/fuse-container $dest/fuse_<version>.sif`; the launcher
finds `squashfuse` in `$dest/bin` automatically.
