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
[`Containerfile`](../perlmutter-container/Containerfile), whose default
`JULIA_CPU_TARGET` is a **universal** set —
`generic;cascadelake;znver2;znver3` (each with `-xsaveopt,-rdrnd,clone_all`) —
so one image is performance-optimal on omega's compute fleet (EPYC 7502 =
znver2 on the unflagged nodes, EPYC 7513 = znver3 on the `amd`-flagged ones),
on NERSC Perlmutter (znver3 Milan), on Intel cascadelake login nodes, and
falls back to `generic` on any other x86_64 machine (laptops included). The
same image is therefore published both as the omega SIF and to GHCR (see
[Publishing](#6-publishing-a-shared-image)).

## Quick start (use the published image — no build needed)

A prebuilt image, the `fuse-container` launcher, an Lmod module, and the
`squashfuse` helper live at:

```
/fusion/projects/dt/fuse_containers/            # shared with the digital_twin group
├── fuse_<version>.sif
├── fuse-container                              # launcher
├── modules/fuse-container.lua                  # Lmod module
└── bin/squashfuse
```

**Easiest — via the module** (one-time `module use`, or add it to `~/.bashrc`):

```bash
module use /fusion/projects/dt/fuse_containers/modules
module load fuse-container

fuse-container                                  # interactive FUSE REPL
fuse-container -e 'using FUSE; @show pkgversion(FUSE)'   # run a script/one-liner
```

The module puts `fuse-container` on PATH and loads singularity; the launcher
picks the newest image automatically. That's the whole thing — no SIF path, no
manual `module load singularity`.

**Without the module**, call the launcher by path (it still auto-loads
singularity and auto-selects the newest image):

```bash
/fusion/projects/dt/fuse_containers/fuse-container -e 'using FUSE; @show pkgversion(FUSE)'
```

Pin a specific version by passing the `.sif` explicitly:
`fuse-container /fusion/projects/dt/fuse_containers/fuse_v1.2.0.sif ...`.

Start a Jupyter notebook on the container (worker nodes only; installs the
kernelspec automatically, `THREADS=N` to change the kernel's thread count):

```bash
/fusion/projects/dt/fuse_containers/fuse-container notebook   # or: lab
```

Or install just the Jupyter kernel for the published image:

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
| `fuse-container.lua` | Lmod modulefile (`module load fuse-container`). |
| `install_kernel.sh` | Installs a Jupyter kernelspec that runs via `fuse-container`. |
| `kernel.json.template` | Template for the kernelspec. |
| `acceptance.sh` | Pre-publish acceptance checks against a freshly built image + SIF. |
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
   (override with `FUSE_ENVIRONMENT=v1.2.0`).
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

### Home directory over quota

The per-user depot above needs only a few MB, but if `$HOME` is at its quota
even the initial `mkdir` fails (`quota -s` to check). The launcher pre-creates
the depot and reports this clearly; without that check Julia dies at package
init with a cryptic `InitError(mkdir ... Unknown system error -122)` (EDQUOT).
Either free up home space, or redirect the writable depot to node-local
scratch (`SINGULARITYENV_` survives `--cleanenv`, plain env vars do not; the
launcher auto-binds the redirected path into the container):

```bash
mkdir -p /local-scratch/$USER/.julia_fuse_container
export SINGULARITYENV_JULIA_DEPOT_PATH="/local-scratch/$USER/.julia_fuse_container:/opt/fuse/.julia"
fuse-container
```

Note `/local-scratch` is per-node: on another Slurm node the depot starts
empty, which is fine — it only holds logs and scratch files.

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

One-step, on a worker node (refuses to run on the login node; add
`--no-browser` and an SSH tunnel for remote use):

```bash
fuse-container notebook        # or: lab; THREADS=8 fuse-container notebook
```

Or install just the kernelspec and use any Jupyter later:

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
mkdir -p $dest/bin $dest/modules
cp <sif-dir>/fuse_<version>.sif                     $dest/
cp deploy/omega-container/fuse-container            $dest/
cp deploy/omega-container/fuse-container.lua        $dest/modules/
cp <sif-dir>/bin/squashfuse                         $dest/bin/
chmod a+rX -R $dest
```

Users then `module use $dest/modules; module load fuse-container` (or call
`$dest/fuse-container` directly). The launcher finds `squashfuse` in `$dest/bin`
and selects the newest `fuse_*.sif` automatically. A site admin with write
access to the shared Lmod tree (e.g. `/fusion/usc/c8/modulefiles-git`) can drop
`fuse-container.lua` there so plain `module load fuse-container` works with no
`module use`.

### Publishing to GHCR (all sites and laptops)

`ghcr.io/projecttorreypines/fuse:<version>` and `:latest` are **manifest
lists**: one tag pointing at both per-architecture images, so users get the
right one automatically from a single command.

```
ghcr.io/projecttorreypines/fuse:<version>        manifest list
    ├── ghcr.io/projecttorreypines/fuse:<version>-amd64    built by hand on omega
    └── ghcr.io/projecttorreypines/fuse:<version>-arm64    built in CI on release
```

The two legs are built in different places, because they hit two unrelated
walls.

**amd64 is blocked by memory.** Emitting the sysimage object peaks at ~98.5 GiB
RSS, far beyond the 16 GB of standard GitHub-hosted runners (it would need the
32-core/128 GB larger-runner class), so amd64 is built by hand on omega.

**arm64 cannot have a sysimage at all** — a linker limit, not a resource one.
FUSE's sysimage embeds a ~3.2 GB blob of serialized program state, and on
aarch64 the linker must bridge it with 32-bit relative references
(`R_AARCH64_PREL32`) that max out at 2.147 GB. The build compiles for 5+ hours
successfully and then fails in the final second at the link step, identically
for every variant tried. x86_64 escapes only because the medium code model
gives it a large-data section (`.ldata`) that moves the blob out of the
critical span; aarch64 has no equivalent. The arm64 image therefore builds with
`FUSE_PRECOMPILE_WORKLOAD=none` and ships per-package pkgimages — many small
libraries instead of one giant one — so `import FUSE` takes ~30–60 s instead of
~5 s, with identical compute speed after that. Skipping the sysimage emission
is also what keeps the peak under 16 GB, so **arm64 builds automatically in CI
on every published release** ([`container`
workflow](../../.github/workflows/container.yml)).

Publishing is therefore two manual steps (arm64 takes care of itself). Both
need `gh` authenticated with `write:packages`
(`gh auth refresh --hostname github.com -s write:packages`).

**1. amd64, after `build.sh` on the same omega node** (the podman store is
node-local). Run the acceptance checks first — they are what catches a
dependency that silently vanished between releases, and they exit nonzero if
anything failed:

```bash
FUSE_ENVIRONMENT=<version> ./deploy/omega-container/acceptance.sh
FUSE_ENVIRONMENT=<version> ./deploy/omega-container/publish_ghcr.sh
```

The push only creates `:<version>-amd64`; it does not move `latest`.

**2. The manifest, once both arch tags exist and have passed their acceptance
tests.** This is the step that moves `latest` for everyone. The CI run that
built arm64 skips the manifest when amd64 is not yet published, so this is
normally run from omega right after step 1:

```bash
MANIFEST=1 FUSE_ENVIRONMENT=<version> ./deploy/omega-container/publish_ghcr.sh
docker manifest inspect ghcr.io/projecttorreypines/fuse:latest  # expect amd64 + arm64
```

Consuming the published image:

```bash
# Laptop / any machine with Docker or podman — architecture selected automatically
docker run -it ghcr.io/projecttorreypines/fuse:<version>

# Laptop JupyterLab: the image ships a Jupyter server + FUSE kernelspec, so
# one command serves JupyterLab on localhost:8888 (like `fuse-container lab`
# on omega; THREADS=N sets the kernel's Julia threads). A fixed JUPYTER_TOKEN
# gives a predictable URL — for the browser, or for VS Code's
# "Select Kernel -> Existing Jupyter Server" (bind 127.0.0.1 so the known
# token is not reachable from the LAN).
docker run -it --pull always -p 127.0.0.1:8888:8888 -e JUPYTER_TOKEN=fuse -e THREADS=8 \
    -v "$PWD":/work -w /work ghcr.io/projecttorreypines/fuse:<version> lab

# NOTE: `docker run` re-uses a locally cached tag without consulting the
# registry. If a machine pulled `latest` before the multi-arch manifest was
# published it will keep running the wrong architecture (e.g. amd64 under
# emulation on Apple Silicon, with HostCPUFeatures "CPU info is out-of-date"
# warnings) until an explicit `docker pull ghcr.io/projecttorreypines/fuse:latest`.

# NERSC Perlmutter
podman-hpc pull ghcr.io/projecttorreypines/fuse:<version>
podman-hpc migrate ghcr.io/projecttorreypines/fuse:<version>
podman-hpc run --rm -it ghcr.io/projecttorreypines/fuse:<version>

# Any cluster with singularity/apptainer (omega included)
singularity build fuse_<version>.sif docker://ghcr.io/projecttorreypines/fuse:<version>
```

The x86_64 image is built with the universal `JULIA_CPU_TARGET`
(`generic;cascadelake;znver2;znver3`), so one image is optimal on omega and
NERSC Perlmutter alike and falls back to `generic` on laptops. It carries the
precompiled sysimage; the arm64 image does not, so its first-call latency is
higher.
