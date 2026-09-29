#!/usr/bin/env python
"""
PORTALS flux-matcher driver for FUSE's `ActorFluxMatcherGACODE`.

Contract (matches the "executable + run dir" option of ActorFluxMatcherGACODE):
  * reads an `input.gacode` from the current working directory (or --input)
  * runs a PORTALS optimization (TGLF/NEO transport) to flux-match the core
  * writes the best (flux-matched) profiles back out as `input.gacode.new`
    (or --output) in the same directory

Typical FUSE usage sets the actor's `executable` to:
    /opt/local/bin/python /Users/mcclenaghan/.julia/dev/FUSE/scripts/run_portals_gacode.py

so the script runs with no arguments inside the actor's run directory and uses
the cwd defaults below. All knobs are overridable via CLI flags.

Environment: this must run with the MacPorts Python 3.11 interpreter that has
the editable MITIM-fusion install (`/opt/local/bin/python`). PORTALS picks up
the local TGLF/NEO + GACODE setup from MITIM's templates/config_user.json
(the `local` machine block + `neo: local`), so no extra env is needed here.
"""

import argparse
import sys
from pathlib import Path


def parse_args(argv):
    p = argparse.ArgumentParser(description="Run PORTALS from an input.gacode and write the flux-matched input.gacode")
    p.add_argument("--input", default="input.gacode", help="input.gacode to start from (default: ./input.gacode)")
    p.add_argument("--output", default="input.gacode.new", help="flux-matched input.gacode to write (default: ./input.gacode.new)")
    p.add_argument("--folder", default="portals_run", help="PORTALS working folder (default: ./portals_run)")
    p.add_argument("--predicted-rho", default="0.25,0.35,0.45,0.55,0.65,0.75,0.85", help="comma-separated rho_tor_norm locations to flux-match")
    p.add_argument("--predicted-channels", default="te,ti,ne", help="comma-separated channels (te,ti,ne,...)")
    p.add_argument("--max-iterations", type=int, default=20, help="maximum BO iterations")
    p.add_argument("--code-settings", default="SAT2", help="TGLF code settings (e.g. SAT0, SAT2, SAT3)")
    p.add_argument("--cold-start", action="store_true", help="wipe the working folder before running")
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(sys.argv[1:] if argv is None else argv)

    from mitim_tools.opt_tools import STRATEGYtools
    from mitim_modules.portals import PORTALSmain
    from mitim_modules.portals.utils import PORTALSanalysis
    from mitim_tools.gacode_tools import PROFILEStools
    from mitim_tools.misc_tools import IOtools

    input_gacode = Path(args.input).expanduser().resolve()
    output_gacode = Path(args.output).expanduser().resolve()
    folder = Path(args.folder).expanduser().resolve()

    if not input_gacode.is_file():
        raise FileNotFoundError(f"input.gacode not found: {input_gacode}")

    if args.cold_start and folder.exists():
        IOtools.shutil_rmtree(folder)
    folder.mkdir(parents=True, exist_ok=True)

    predicted_rho = [float(x) for x in args.predicted_rho.split(",") if x.strip()]
    predicted_channels = [c.strip() for c in args.predicted_channels.split(",") if c.strip()]

    # ------------------------------------------------------------------
    # Configure PORTALS
    # ------------------------------------------------------------------
    portals_fun = PORTALSmain.portals(folder)
    portals_fun.optimization_options["convergence_options"]["maximum_iterations"] = args.max_iterations
    portals_fun.portals_parameters["solution"]["predicted_rho"] = predicted_rho
    portals_fun.portals_parameters["solution"]["predicted_channels"] = predicted_channels
    portals_fun.portals_parameters["transport"]["options"]["tglf"]["run"]["code_settings"] = args.code_settings

    # ------------------------------------------------------------------
    # Load + sanitize the starting plasma state, then run
    # ------------------------------------------------------------------
    plasma_state = PROFILEStools.gacode_state(input_gacode)
    plasma_state.correct(options={
        "recalculate_ptot": True,
        "remove_fast": True,
        "quasineutrality": True,
        "enforce_same_aLn": True,
    })

    portals_fun.prep(plasma_state)

    mitim_bo = STRATEGYtools.MITIM_BO(portals_fun, cold_start=args.cold_start, askQuestions=False)
    mitim_bo.run()

    # ------------------------------------------------------------------
    # Extract the best (flux-matched) profiles and write them back out
    # ------------------------------------------------------------------
    portals_output = PORTALSanalysis.PORTALSanalyzer.from_folder(folder)
    best = portals_output.extractProfiles()  # defaults to the best evaluation (ibest)
    best.write_state(file=output_gacode)
    print(f"[run_portals_gacode] wrote flux-matched profiles to {output_gacode}")

    # ------------------------------------------------------------------
    # Also export the best-evaluation fluxes (turbulent + neoclassical, gyro-Bohm units)
    # next to the output gacode so the FUSE actor can store them in dd.core_transport.
    # ------------------------------------------------------------------
    import shutil

    ibest = portals_output.ibest
    eval_dir = folder / "Execution" / f"Evaluation.{ibest}" / "transport_simulation_folder"
    for suffix in ("turb", "neoc"):
        src = eval_dir / f"fluxes_{suffix}.json"
        if not src.exists():
            src = eval_dir / "plasma_0" / f"fluxes_{suffix}.json"  # batched single-plasma layout
        dest = output_gacode.parent / f"fluxes_{suffix}.json"
        if src.exists():
            shutil.copyfile(src, dest)
            print(f"[run_portals_gacode] copied {suffix} fluxes (eval {ibest}) to {dest}")
        else:
            print(f"[run_portals_gacode] WARNING: fluxes_{suffix}.json not found for ibest={ibest}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
