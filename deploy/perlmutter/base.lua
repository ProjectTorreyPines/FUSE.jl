whatis("Name    : fuse")
whatis("Version : " .. fuse_env)

depends_on("julia/1.11.4")

local envdir = basedir .. "/environments/" .. fuse_env
local base_depot = envdir .. "/.julia"

setenv("FUSE_HOME", basedir)
setenv("FUSE_ENVIRONMENT", fuse_env)

-- We put the user depot first so their own packages get installed there,
--   then the FUSE environment's depot after so it can find packages for the
--   precompiled sysimage
setenv("JULIA_DEPOT_PATH", base_depot .. ":")

-- The FUSE sysimage enviornment is the last place julia looks for packages
--   when a user does `using <package>`, but this allows Julia to automatically
--   find FUSE, Plots, and IJulia.
setenv("JULIA_LOAD_PATH", ":" .. envdir)

setenv("JULIA_CC", "gcc -O3")

-- This lets the compiled sysimage work on login and worker nodes
setenv("JULIA_CPU_TARGET", "generic")

prepend_path("JUPYTER_PATH", envdir .. "/.jupyter")

prepend_path("PATH", envdir)