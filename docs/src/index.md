#   FUSE (**FU**sion **S**ynthesis **E**ngine)

FUSE (**FU**sion **S**ynthesis **E**ngine) is GA's proprietary framework for Fusion Power Plant (FPP) integrated design.

## FUSE objectives

* Couple physics, engineering, control, costing, and balance of plant
* Enable both stationary as well as time-dependent simulations
* Be generic and modular, supporting hierarcy of models
* Leverage parallelism and HPC systems for optimization studies
* Support sensitivity and uncertainty quantification analyses

## Basic concepts

FUSE is written completely in Julia, and is structured as follows:
1. Data is stored in the `dd` data structure, which is based on the ITER IMAS onthology
1. Physics and engineering `actors` are the fundamental building blocks of FUSE simulations
1. Actors functionality is controlled via `act` parameters
1. The data structure can be initialized starting from 0D `ini` parameters
1. FUSE comes with a series of template `use cases` for different machines (FPP, ITER, ARC, ...)
1. `workflows` perform self-contained studies/optimizations (typically running many FUSE simulations)
1. FUSE can interface with the existing GA ecosystem of modeling codes (GASC, OMFIT/OMAS) as well as IMAS

These concepts are illustrated in the diagram below:
![svg](./assets/FUSE.svg)

Here is an example, illustrating how a simple FUSE simulation can be setup and run in Julia:
```julia
# using the Revise and FUSE packages (Revise for development)
using Revise
using FUSE
global_logger(FUSE.logger);
# get `ini` and `act` for a given use-case
ini, act = FUSE.case_parameters(:FPP; version=:v1_demount, init_from=:scalars)
# initialize `dd` from 0D parameters
dd = FUSE.init(ini, act; do_plot=true)
# run an actor
FUSE.ActorEquilibriumTransport(dd, act; do_plot=true);
```