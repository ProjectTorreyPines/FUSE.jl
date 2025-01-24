whatis("Name    : fuse")
whatis("Version : " .. fuse_env)

whatis("Name    : fuse")
whatis("Version : " .. fuse_env)

purge()
depends_on("julia/11.2")
depends_on("env/gcc11.x")

local envdir = basedir .. "/environments/" .. fuse_env
local user_depot = os.getenv("JULIA_USER_DEPOT")
local base_depot = envdir .. "/.julia"

setenv("FUSE_HOME", basedir)
setenv("FUSE_ENVIRONMENT", fuse_env)

setenv("JULIA_DEPOT_PATH",  user_depot  .. ":" .. base_depot .. ":")
setenv("JULIA_LOAD_PATH", ":" .. envdir)
setenv("JULIA_CC", "gcc -O3")
setenv("JULIA_CPU_TARGET", "generic;sandybridge,-xsaveopt,clone_all;haswell,-rdrnd,base(1)")

prepend_path("JUPYTER_PATH", envdir .. "/.jupyter")

prepend_path("PATH", basedir .. "/miniconda3/bin")

local fuse_sysimage = envdir .. "/sys_fuse.so"
set_alias("julia", "julia --sysimage=" .. fuse_sysimage)
