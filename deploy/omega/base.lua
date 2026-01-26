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

-- Environment variables for data fetching
setenv("FUSE_OMFIT_HOST", "localhost")
setenv("FUSE_OMFIT_ROOT", "/fusion/projects/theory/fuse/d3d_data_fetching/OMFIT-source")
setenv("FUSE_OMAS_ROOT", "/fusion/projects/theory/fuse/d3d_data_fetching/omas")


-- We put the user depot first so their own packages get installed there,
--   then the FUSE environment's depot after so it can find packages for the
--   precompiled sysimage
setenv("JULIA_DEPOT_PATH",  user_depot  .. ":" .. base_depot .. ":")

-- The FUSE sysimage enviornment is the last place julia looks for packages
--   when a user does `using <package>`, but this allows Julia to automatically
--   find FUSE, Plots, and IJulia.
setenv("JULIA_LOAD_PATH", ":" .. envdir)

setenv("JULIA_CC", "gcc -O3")

-- This lets the compiled sysimage work on login and worker nodes
setenv("JULIA_CPU_TARGET", "generic")

prepend_path("JUPYTER_PATH", envdir .. "/.jupyter")

prepend_path("PATH", basedir .. "/miniconda3/bin")
prepend_path("PATH", envdir)
