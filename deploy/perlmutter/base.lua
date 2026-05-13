whatis("Name    : fuse")
whatis("Version : " .. fuse_env)
depends_on("julia-latest/1.0")
depends_on("julia/1.11.4")

local envdir = basedir .. "/environments/" .. fuse_env
local base_depot = envdir .. "/.julia"

setenv("FUSE_HOME", basedir)
setenv("FUSE_ENVIRONMENT", fuse_env)

-- We put the user depot first so their own packages get installed there,
--   then the FUSE environment's depot after so it can find packages for the
--   precompiled sysimage

local home       = os.getenv("HOME")                  -- user’s home directory
local user_depot = os.getenv("JULIA_DEPOT_PATH")      -- current value (may be nil)

if not user_depot then
   -- case 1:  JULIA_DEPOT_PATH is unset, so use default HOME user depot
   setenv("JULIA_DEPOT_PATH", home .. "/.julia:" .. base_depot .. ":")
elseif user_depot:sub(-1) == ":" then
   -- case 2: JULIA_DEPOT_PATH ends with “:”, so append base_depot after it
   --         Already ends in ":", so no need to append a ":" in between
   setenv("JULIA_DEPOT_PATH", user_depot .. base_depot .. ":")
else
   LmodError(
      "Cannot parse existing JULIA_DEPOT_PATH=" .. user_depot .. "\n" ..
      "It must be unset or end with ':'."
   )
end

-- Environment variables for data fetching
setenv("FUSE_OMFIT_HOST", "localhost")
setenv("FUSE_OMFIT_ROOT", "/global/common/software/m3739/perlmutter/OMFIT-CAKE")
setenv("FUSE_OMAS_ROOT", "/global/common/software/m3739/perlmutter/FUSE_OMAS")
setenv("FUSE_RESULT_ARCHIVE", pathJoin("/global/cfs/cdirs/m3739/FUSE/d3d-time-dependent", os.getenv("USER")))

-- The FUSE sysimage enviornment is the last place julia looks for packages
--   when a user does `using <package>`, but this allows Julia to automatically
--   find FUSE, Plots, and IJulia.
setenv("JULIA_LOAD_PATH", ":" .. envdir)

setenv("JULIA_CC", "gcc -O3")

-- This lets the compiled sysimage work on login and worker nodes
setenv("JULIA_CPU_TARGET", "generic")

-- symbolic link kernels
local kernsrc = envdir .. "/.jupyter/kernels"
local kerndst = home  .. "/.local/share/jupyter/kernels"
execute{
	cmd ="mkdir -p " .. kerndst .. " && ln -sf " .. kernsrc .."/* " .. kerndst,
	modeA = {"load"}
}


prepend_path("PATH", envdir)