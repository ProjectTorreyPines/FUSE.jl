#!/usr/bin/env julia
"""
Step 2 after install: run fluxmatcher.ipynb cell 0 (code imports).

Cell 1 is markdown and is checked by verify_fluxmatcher_notebook.sh with Python.
"""

println("[verify_fluxmatcher] Running cell 0: using Revise, Plots, FUSE")
flush(stdout)

@time using Revise
@time using Plots
@time using FUSE

println("[verify_fluxmatcher] Cell 0 passed")
