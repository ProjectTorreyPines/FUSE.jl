whatis("Name    : fuse")
whatis("Version : " .. fuse_env)

depends_on("julia/1.11.3")
depends_on("env/gcc11.x")

-- FUSE environment uses its own conda install with custom jupyter kernels
conflict("python")
conflict("mamba")

local envdir = basedir .. "/environments/" .. fuse_env
local user_depot = os.getenv("JULIA_USER_DEPOT")
local base_depot = envdir .. "/.julia"

setenv("FUSE_HOME", basedir)
setenv("FUSE_ENVIRONMENT", fuse_env)

-- We put the user depot first so their own packages get installed there,
--   then the FUSE environment's depot after so it can find packages for the
--   precompiled sysimage
setenv("JULIA_DEPOT_PATH",  user_depot  .. ":" .. base_depot .. ":")

-- The FUSE sysimage enviornment is the last place julia looks for packages
--   when a user does `using <package>`, but this allows Julia to automatically
--   find FUSE, Plots, and IJulia.
setenv("JULIA_LOAD_PATH", ":" .. envdir)

setenv("JULIA_CC", "gcc -O3")

-- This lets the compiled sysimage work on login and worker nodes,
--   modeled after how the Julia binaries are built
setenv("JULIA_CPU_TARGET", "generic;sandybridge,-xsaveopt,clone_all;haswell,-rdrnd,base(1)")


prepend_path("JUPYTER_PATH", envdir .. "/.jupyter")

prepend_path("PATH", basedir .. "/miniconda3/bin")

local fuse_sysimage = envdir .. "/sys_fuse.so"
set_alias("julia", "julia --sysimage=" .. fuse_sysimage)
