#!/usr/bin/env julia
"""
Julia-side FUSE install steps (invoked from install_fuse_*.sh, not the REPL).

Environment:
  FUSE_INSTALL_DIR   optional fusebot install directory (NERSC: ~/.local/shared/bin)
  FUSE_SETUP_SHELL   "true" to append fusebot dir to shell rc (recommended on HPC)
"""

using Pkg

const FUSE_UUID = Base.UUID("e64856f0-3bb8-4376-b4b7-c03396503992")

function log(msg)
    println("[install_fuse] ", msg)
    flush(stdout)
end

function fuse_module()
    return Base.require(Base.PkgId(FUSE_UUID, "FUSE"))
end

function add_registries!()
    log("Adding FuseRegistry and General")
    Pkg.Registry.add(RegistrySpec(url="https://github.com/ProjectTorreyPines/FuseRegistry.jl.git"))
    Pkg.Registry.add("General")
end

function install_fuse_packages!()
    if isfile("Project.toml")
        log("Project.toml found — instantiating environment")
        Pkg.instantiate()
        log("Updating project environment")
        Pkg.update()
    end

    add_registries!()

    if haskey(Pkg.dependencies(), Base.UUID("e64856f0-3bb8-4376-b4b7-c03396503992"))
        log("FUSE already installed — resolving and updating")
    else
        log("Adding FUSE (15–30 minutes on first install)")
        Pkg.add("FUSE")
    end
    Pkg.resolve()
    log("Updating FUSE dependencies")
    Pkg.update()
    Pkg.precompile()
end

function install_revise!()
    log("Adding Revise.jl and enabling startup hook")
    Pkg.add("Revise")
    julia_dir = first(DEPOT_PATH)
    config_dir = joinpath(julia_dir, "config")
    startup = joinpath(config_dir, "startup.jl")
    mkpath(config_dir)
    existing = isfile(startup) ? read(startup, String) : ""
    if !occursin("using Revise", existing)
        open(startup, "w") do io
            println(io, "using Revise")
            print(io, existing)
        end
        log("Wrote using Revise to $startup")
    else
        log("Revise already in $startup")
    end
end

function install_fusebot!()
    FUSE = fuse_module()
    setup_shell = get(ENV, "FUSE_SETUP_SHELL", "false") == "true"
    if haskey(ENV, "FUSE_INSTALL_DIR") && !isempty(ENV["FUSE_INSTALL_DIR"])
        dir = ENV["FUSE_INSTALL_DIR"]
        log("Installing fusebot to $dir")
        Base.invokelatest(FUSE.install_fusebot, dir; setup_shell=setup_shell)
    else
        log("Installing fusebot to default location")
        Base.invokelatest(FUSE.install_fusebot; setup_shell=setup_shell)
    end
    println(Base.invokelatest(FUSE.default_fusebot_dir))
end

function smoke_test!()
    FUSE = fuse_module()
    log("Smoke test: case_parameters(:FPP) + init")
    ini, act = Base.invokelatest(FUSE.case_parameters, :FPP)
    Base.invokelatest(FUSE.init, ini, act)
    log("Smoke test passed")
end

const STEPS = Dict(
    "packages" => install_fuse_packages!,
    "revise" => install_revise!,
    "fusebot" => install_fusebot!,
    "smoke" => smoke_test!,
    "all" => function ()
        install_fuse_packages!()
        install_revise!()
        install_fusebot!()
        smoke_test!()
    end,
)

if isempty(ARGS)
    STEPS["all"]()
else
    for step in ARGS
        haskey(STEPS, step) || error("Unknown step: $step (choose: $(join(keys(STEPS), ", ")))")
        STEPS[step]()
    end
end
