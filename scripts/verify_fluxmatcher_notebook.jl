#!/usr/bin/env julia
"""
Step 2 after install: run the first three cells of fluxmatcher.ipynb.

  Cell 0 (code):     using Revise, using Plots, using FUSE
  Cell 1 (markdown): flux-matcher introduction (nothing to execute)
  Cell 2 (code):     flux-match the DIII-D L-mode case

The wrapper scripts (verify_fluxmatcher_notebook.sh / .ps1) extract the code
from the first three cells into a temp file and point `FUSE_VERIFY_CELLS` at it,
so this runner stays faithful to whatever the notebook actually contains (the
markdown cell is skipped since there is nothing to execute).

Run standalone without `FUSE_VERIFY_CELLS` set and it falls back to cell 0
(the imports only).
"""

cells_file = get(ENV, "FUSE_VERIFY_CELLS", "")

if !isempty(cells_file) && isfile(cells_file)
    println("[verify_fluxmatcher] Running the first three cells (code cells 0 and 2; cell 1 is markdown)")
    println("[verify_fluxmatcher] NOTE: cell 2 flux-matches a DIII-D L-mode case (often ~6 minutes on one thread the first time)")
    flush(stdout)
    @time include(cells_file)
    println("[verify_fluxmatcher] First three cells passed")
else
    println("[verify_fluxmatcher] FUSE_VERIFY_CELLS not set — running cell 0 imports only")
    flush(stdout)
    @time using Revise
    @time using Plots
    @time using FUSE
    println("[verify_fluxmatcher] Cell 0 passed")
end
