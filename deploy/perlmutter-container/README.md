# FUSE in a podman-hpc container on NERSC Perlmutter

This directory builds a self-contained [podman-hpc](https://docs.nersc.gov/development/containers/podman-hpc/overview/)
image that bakes in FUSE.jl, all ProjectTorreyPines packages, and a precompiled
PackageCompiler system image (`sys_fuse.so`) for fast startup. It mirrors the
Lmod module deploy in [`../perlmutter`](../perlmutter) but ships everything in a
single OCI image instead of an environment module.

## Contents

| File | Purpose |
|------|---------|
| `Containerfile` | Image definition (`FROM docker.io/library/julia:1.11.7`), installs FUSE + builds the sysimage. |
| `install_fuse_container.jl` | Runs inside the build: adds registries, installs packages, compiles `sys_fuse.so`. |
| `build.sh` | Builds and migrates the image on Perlmutter. |
| `kernel.json.template` | Jupyter kernelspec template wrapping `podman-hpc run --jupyter`. |
| `install_kernel.sh` | Generates and installs the FUSE Jupyter kernel for the current user. |

## 1. Build the image

The build runs the FUSE test suite + tutorial as the sysimage precompile
workload, so it is heavy. Run it from an interactive CPU allocation:

```bash
salloc -N 1 -C cpu -q interactive -A m3739 -t 04:00:00

cd $PSCRATCH/.julia/dev/FUSE        # the FUSE repo root (build context)
./deploy/perlmutter-container/build.sh
```

`build.sh` will:

1. Resolve the image tag (`fuse:<version>`) from the latest FUSE.jl release
   (override with `FUSE_ENVIRONMENT=v1.1.3`).
2. `podman-hpc build -t fuse:<version> -f deploy/perlmutter-container/Containerfile .`
3. `podman-hpc migrate fuse:<version>` — converts the image to a read-only
   squashfs under `$SCRATCH` so it can be used in jobs and on any login node.
   (Set `SQUASH_DIR=/global/cfs/cdirs/m3739/shared_images` to share it across
   the project.)

The Containerfile's default `JULIA_CPU_TARGET` is a universal set
(`generic;cascadelake;znver2;znver3`) that includes Perlmutter's AMD Milan
(znver3) among its optimized targets, so the same image also runs on GA's
omega cluster and on laptops. Override with
`--build-arg JULIA_CPU_TARGET="generic;znver3,-xsaveopt,-rdrnd,clone_all"`
for a leaner Perlmutter-only build.

> Instead of building, you can pull the published universal image from the
> GitHub Container Registry (see
> [`../omega-container/README.md`](../omega-container/README.md) for how it is
> published):
>
> ```bash
> podman-hpc pull ghcr.io/projecttorreypines/fuse:<version>
> podman-hpc migrate ghcr.io/projecttorreypines/fuse:<version>
> ```

> Re-run `build.sh` (build + migrate) after any change. A migrated image is
> read-only; you must re-`migrate` to pick up a rebuild.

## Testing the shared image

The project image is published to `/global/cfs/cdirs/m3739/shared_images` and is
usable by anyone in `m3739` (no build or pull needed) by adding `--squash-dir`.
On a Perlmutter login node:

```bash
# shortcut so every command uses the shared image store
alias fuse-podman='podman-hpc --squash-dir /global/cfs/cdirs/m3739/shared_images'

# 1) confirm the image is visible
fuse-podman images | grep fuse
# expect: localhost/fuse  v1.1.3  ...  R/O

# 2) fast smoke test (loads FUSE from the baked-in sysimage, starts in seconds)
fuse-podman run --rm fuse:v1.1.3 \
  julia --sysimage=/opt/fuse/sys_fuse.so \
  -e 'using FUSE; ini,act=FUSE.case_parameters(:D3D,:L_mode); println("FUSE OK: ", pkgversion(FUSE))'

# 3) offline test (proves it is fully self-contained: no "Downloading artifact" lines)
fuse-podman run --rm --network none fuse:v1.1.3 \
  julia --sysimage=/opt/fuse/sys_fuse.so \
  -e 'using FUSE; ini,act=FUSE.case_parameters(:D3D,:L_mode); dd=FUSE.init(ini,act); println("OFFLINE OK")'

# 4) interactive REPL (optional)
fuse-podman run --rm -it fuse:v1.1.3
```

For the Jupyter equivalent, see [Interactive use via Jupyter](#3-interactive-use-via-jupyter)
(install the kernel with `SQUASH_DIR=/global/cfs/cdirs/m3739/shared_images`).

## 2. Run the FUSE REPL

```bash
podman-hpc run --rm -it fuse:<version>
```

This launches a Julia REPL with the FUSE sysimage preloaded (startup in
seconds). Quick smoke test:

```bash
podman-hpc run --rm fuse:<version> \
  julia --sysimage=/opt/fuse/sys_fuse.so \
  -e 'using FUSE; ini,act=FUSE.case_parameters(:D3D,:L_mode); println("FUSE ok")'
```

Use more threads at runtime:

```bash
podman-hpc run --rm -it -e JULIA_NUM_THREADS=8 fuse:<version>
```

## 3. Interactive use via Jupyter

```bash
FUSE_ENVIRONMENT=<version> ./deploy/perlmutter-container/install_kernel.sh
# optionally: THREADS=8 FUSE_ENVIRONMENT=<version> ./install_kernel.sh
```

This writes `$HOME/.local/share/jupyter/kernels/fuse-<version>/kernel.json`,
whose `argv` runs the in-container Julia via `podman-hpc run --rm --jupyter`.
The `--jupyter` flag bind-mounts `/tmp` and `$HOME` so the kernel can connect
and create notebooks in your home directory.

To use an image from a shared squash dir (instead of your per-user store), set
`SQUASH_DIR`; the generated kernel will pass `podman-hpc --squash-dir <dir> run`:

```bash
SQUASH_DIR=/global/cfs/cdirs/m3739/shared_images \
FUSE_ENVIRONMENT=<version> ./deploy/perlmutter-container/install_kernel.sh
```

### Test the kernel

1. Open [NERSC JupyterHub](https://jupyter.nersc.gov) and start a server on a
   Perlmutter login node.
2. Create a new notebook and select the **Julia FUSE-<version>** kernel
   (display name e.g. `Julia FUSE-v1.1.3 (1 thread(s))`).
3. Run this test cell:

   ```julia
   using FUSE
   ini, act = FUSE.case_parameters(:D3D, :L_mode)
   dd = FUSE.init(ini, act)
   @show pkgversion(FUSE)
   ```

Expected: the kernel connects within a few seconds, `using FUSE` returns almost
instantly (sysimage), the actors run, and `pkgversion(FUSE)` prints the image
version. If the kernel fails to start, check the kernelspec at
`$HOME/.local/share/jupyter/kernels/fuse-<version>/kernel.json` and confirm
`podman-hpc images` (add `--squash-dir <dir>` if shared) lists `fuse:<version>`.

## Mounting data (optional)

Containers only see `/tmp` and `$HOME` by default (plus whatever `--jupyter`
mounts). To read/write data elsewhere, bind-mount it:

```bash
podman-hpc run --rm -it \
  --volume $SCRATCH/fuse_runs:/work \
  fuse:<version>
```

For the DIII-D time-dependent study, the OMFIT/OMAS roots used by the module
deploy live on CFS and are **not** baked into the image. Mount them and set the
matching env vars if needed:

```bash
podman-hpc run --rm -it \
  --volume /global/common/software/m3739/perlmutter/OMFIT-CAKE:/omfit:ro \
  --volume /global/common/software/m3739/perlmutter/FUSE_OMAS:/omas:ro \
  -e FUSE_OMFIT_ROOT=/omfit -e FUSE_OMAS_ROOT=/omas -e FUSE_OMFIT_HOST=localhost \
  fuse:<version>
```

## Notes

- Defaults reuse the existing deploy: account `m3739` and image tag from the
  latest FUSE release. Adjust via `salloc -A ...` and `FUSE_ENVIRONMENT`.
- The image is built from the official `julia:1.11.7` base (the required Julia
  version) rather than rebuilding Julia from the NERSC module.
- See the NERSC docs: [podman-hpc overview](https://docs.nersc.gov/development/containers/podman-hpc/overview/)
  and [containers in Jupyter](https://docs.nersc.gov/services/jupyter/how-to-guides/#how-to-use-a-container-to-run-a-jupyter-kernel).
