# FUSE_Design Project Memory

## User Preferences
- Never add `Co-Authored-By` lines to commit messages

## Project Structure
- Parameter definitions: `src/parameters/parameters_inits.jl` (`FUSEparameters__core_profiles` at ~line 77)
- Core profiles init: `src/ddinit/init_core_profiles.jl` (high-level wrapper + low-level implementation)
- Case files: `src/cases/*.jl` (ARC, ITER, SPARC, etc.)
- `init_core_profiles!` has two methods: one taking `(dd, ini, act)` that unpacks ini params, one taking `(cp, equil, summary; kwargs...)`