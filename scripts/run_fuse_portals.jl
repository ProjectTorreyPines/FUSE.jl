#!/usr/bin/env julia
"""
Run a FUSE case whose core transport is flux-matched by PORTALS (MITIM).

FUSE initializes the case, hands the plasma state to PORTALS as an `input.gacode`
(`ActorFluxMatcherGACODE`, driven through `scripts/run_portals_gacode.py`), and loads the
flux-matched profiles + turbulent/neoclassical fluxes back into `dd`.

Usage:

    julia --project=@. scripts/run_fuse_portals.jl [flags]

Flags (all optional):

    --case=ARC              FUSE case name passed to `FUSE.case_parameters`
    --dir=<path>            run directory (default: ./portals_<case>_<timestamp>)
    --python=<path>         python with MITIM-fusion installed
    --iterations=20         maximum PORTALS (Bayesian optimization) iterations
    --rho=0.25,0.35,...     rho_tor_norm locations to flux-match
    --channels=te,ti,ne     channels to predict
    --sat=sat2              TGLF saturation rule (sat0, sat1, sat2, sat3)
    --dry-run               initialize and write `input.gacode` only, do not run PORTALS

The run directory keeps everything: the `input.gacode` FUSE handed over, PORTALS's own
`portals_run/` folder, the flux-matched `input.gacode.new`, and `dd.json` at the end.
"""

using FUSE
import IMAS
import GACODE
import Dates

#= ================= =#
#  command line flags  #
#= ================= =#

function flag(name::String, default::String)
    for arg in ARGS
        startswith(arg, "--$name=") && return split(arg, "=", limit=2)[2]
    end
    return default
end

hasflag(name::String) = "--$name" in ARGS

const case = Symbol(flag("case", "ARC"))
const python = flag("python", "python")
const portals_driver = abspath(joinpath(@__DIR__, "run_portals_gacode.py"))
const max_iterations = parse(Int, flag("iterations", "20"))
const rho_transport = parse.(Float64, split(flag("rho", "0.25,0.35,0.45,0.55,0.65,0.75,0.85"), ","))
const channels = Symbol.(split(flag("channels", "te,ti,ne"), ","))
const sat_rule = Symbol(flag("sat", "sat2"))
const dry_run = hasflag("dry-run")
const run_dir = abspath(flag("dir", "portals_$(lowercase(string(case)))_$(Dates.format(Dates.now(), "yyyymmdd_HHMMSS"))"))

#= ========= =#
#  preflight  #
#= ========= =#

# fail on a missing python/driver now rather than after the (slow) case initialization
# (`success`, not `run`: FUSE exports its own `run` for studies)
@assert Sys.which(python) !== nothing "python executable not found: $python (pass --python=<path>)"
@assert isfile(portals_driver) "PORTALS driver not found: $portals_driver"
if !dry_run
    @assert success(pipeline(`$python -c "import mitim_modules.portals"`; stdout=devnull, stderr=devnull)) "`$python` cannot import MITIM (`import mitim_modules.portals` failed)"
end

mkpath(run_dir)
@info "PORTALS run directory: $run_dir"

#= ================== =#
#  case setup and init  #
#= ================== =#

ini, act = FUSE.case_parameters(case)

# core transport: delegate the whole flux match to PORTALS
act.ActorFluxMatcher.algorithm = :external
act.ActorFluxMatcher.rho_transport = rho_transport
act.ActorFluxMatcher.max_iterations = max_iterations
act.ActorFluxMatcher.evolve_Te = :te in channels ? :flux_match : :fixed
act.ActorFluxMatcher.evolve_Ti = :ti in channels ? :flux_match : :fixed
act.ActorFluxMatcher.evolve_densities = :ne in channels ? :flux_match : :fixed
act.ActorTGLF.sat_rule = sat_rule

# how the external code is launched: `executable` is run inside `save_dir`, where it finds
# `input.gacode` and is expected to write `input.gacode.new` (+ fluxes_turb/neoc.json)
act.ActorFluxMatcherGACODE.executable = "$python $portals_driver"
act.ActorFluxMatcherGACODE.save_dir = run_dir
act.ActorFluxMatcherGACODE.clear_workdir = false  # keep the PORTALS folder for inspection
act.ActorFluxMatcherGACODE.save_fluxes = true

dd = FUSE.init(ini, act)

IMAS.imas2json(dd, joinpath(run_dir, "dd_init.json"); freeze=false, strict=false)

#= ==================== =#
#  hand off to PORTALS   #
#= ==================== =#

if dry_run
    # same file the actor would hand to PORTALS, so it can be run by hand:
    #   cd <run_dir> && python run_portals_gacode.py
    GACODE.save(GACODE.InputGACODE(dd), joinpath(run_dir, "input.gacode"))
    @info "dry run: wrote $(joinpath(run_dir, "input.gacode")) (PORTALS not launched)"
    exit(0)
end

@info "launching PORTALS: $(length(rho_transport)) radii, channels $(join(channels, ",")), $(uppercase(string(sat_rule))), up to $max_iterations iterations"
actor = FUSE.ActorFluxMatcher(dd, act)

#= ======= =#
#  results  #
#= ======= =#

IMAS.imas2json(dd, joinpath(run_dir, "dd.json"); freeze=false, strict=false)

cp1d = dd.core_profiles.profiles_1d[]
@info "flux-matched profiles: Te0 = $(round(cp1d.electrons.temperature[1] / 1e3; digits=2)) keV, " *
      "Ti0 = $(round(cp1d.t_i_average[1] / 1e3; digits=2)) keV, " *
      "ne0 = $(round(cp1d.electrons.density_thermal[1] / 1e19; digits=2)) 10¹⁹ m⁻³"
# 3-arg `show`: the 2-arg one dumps the raw ExtractFunction dict instead of the table
show(stdout, MIME"text/plain"(), IMAS.extract(dd))
println()
@info "results in $run_dir (dd.json, input.gacode, input.gacode.new, portals_run/)"
