-- Lmod modulefile for the FUSE Singularity container on omega.
--
-- Puts the `fuse-container` launcher on PATH and loads singularity, so users
-- can run the containerized FUSE with a single command:
--
--     module use /fusion/projects/dt/fuse_containers/modules
--     module load fuse-container
--     fuse-container                       # REPL (newest image)
--     fuse-container -e 'using FUSE; ...'  # run a one-liner / script
--
-- Add the `module use` line to ~/.bashrc to make it permanent. (A site admin
-- with access to the shared module tree can instead drop this file there so
-- plain `module load fuse-container` works with no `module use`.)
--
-- containers_dir must hold: the launcher `fuse-container`, the `fuse_*.sif`
-- image(s), and `bin/squashfuse`.
local containers_dir = "/fusion/projects/dt/fuse_containers"

help([[
Containerized FUSE (Singularity) on omega — self-contained image with FUSE,
all ProjectTorreyPines packages, and a precompiled sysimage. No private depot.

  fuse-container                       interactive FUSE REPL (newest image)
  fuse-container -e 'using FUSE; ...'  run a script / one-liner
  fuse-container <path.sif> ...        pin a specific image version

See deploy/omega-container/README.md in the FUSE.jl repo.
]])

whatis("Name    : fuse-container")
whatis("Version : container")

-- The launcher auto-loads singularity too, but declaring it keeps `module list`
-- honest and gives a clean error if the singularity module is unavailable.
depends_on("singularity/3.11.3")

setenv("FUSE_SIF_DIR", containers_dir)
prepend_path("PATH", containers_dir)
