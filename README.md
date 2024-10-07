# FUSE.jl

FUSE (**FU**sion **S**ynthesis **E**ngine) is an open-source framework for the integrated design of Fusion Power Plants (FPP). Originally developed by General Atomics, FUSE is now publicly available under the [Apache 2.0 license](https://fuse.help/dev/notice.html).

## 📢 Upcoming 2024 Code Camp 📢

Join the FUSE community and help shape the future of Fusion Power Plant design:

* **Date:** Dec 9th - 13th
* **Location:** In person @ General Atomics
* **Seats:** Limited to 40 participants - **[🔥 REGISTER TODAY!](https://docs.google.com/forms/d/e/1FAIpQLSc948qTy9HmlAwj9x-kzYYTStrLvgDT4UCEkM3HwTZnovBjsw/viewform)**

## Resources

Here are some key resources for getting started with FUSE:

* 📚 **[Online documentation](https://fuse.help)**
* 📊 **[Summary presentation](https://tinyurl.com/FUSEslideDeck)**
* 📜 **[Preprint publication](https://arxiv.org/abs/2409.05894)**
* 📦 **[Package ecosystem](https://fuse.help/dev/deps.html)**
* 🆘 **[Discord community](https://discord.gg/CbjpZH9SKM)**

## Objectives

FUSE aims to achieve the following objectives:

* ⚡ Provide a highly efficient, modular framework that tightly couples models across different domains.
* 🧩 Integrate plasma physics, engineering, control, balance of plant, and costing systems.
* 🤖 Leverage machine learning to overcome the typical fidelity/speed tradeoff in simulations.
* ⏱️ Support both stationary and time-dependent simulations.
* 💻 Harness parallelism and high-performance computing (HPC) for large-scale studies.
* 🎯 Perform multi-objective constrained optimization to explore design tradeoffs.
* 🔍 Enable comprehensive sensitivity analysis and uncertainty quantification.

## Basic Concepts

FUSE is entirely written in Julia and is structured around the following core concepts:

1. **📂 Data storage**: All data is stored in the `dd` structure, which follows the ITER IMAS ontology.
2. **🧠 Actors**: The core components of FUSE simulations are physics and engineering actors.
3. **🕹️ Control**: Actor functionality is governed by `act` parameters.
4. **🚀 Initialization**: The data structure can be initialized from 0D `ini` parameters.
5. **🔧 Use cases**: FUSE includes templates for various machines (e.g., FPP, ITER, ARC).
6. **🔄 Workflows**: Self-contained studies and optimizations are conducted via workflows, typically involving multiple FUSE simulations.
7. **🌍 Interoperability**: FUSE interfaces with existing modeling tools like OMFIT/OMAS and the IMAS ecosystem.

A diagram illustrating these concepts is provided below:  
![FUSE Diagram](./docs/src/assets/FUSE.svg)

## Usage Example

Here’s a simple example of setting up and running a FUSE simulation in Julia:

```julia
using FUSE

# Obtain `ini` and `act` parameters for a specific use case
ini, act = FUSE.case_parameters(:FPP)

# Initialize the `dd` structure with 0D parameters
dd = FUSE.init(ini, act)

# Run a stationary plasma actor simulation
FUSE.ActorStationaryPlasma(dd, act)

# Get an overview of the simulation results
FUSE.digest(dd)
```

## Installation

FUSE and its related packages are available through the [FuseRegistry](https://github.com/ProjectTorreyPines/FuseRegistry.jl/). To install:

1. [Install Julia](https://github.com/JuliaLang/juliaup?tab=readme-ov-file#juliaup---julia-version-manager)

2. Add the FuseRegistry and General registries, then install FUSE:

```julia
using Pkg
Pkg.Registry.add(RegistrySpec(url="https://github.com/ProjectTorreyPines/FuseRegistry.jl.git"))
Pkg.Registry.add("General")
Pkg.add("FUSE")
```

## Citation

Please cite this work as follows:

```
@article{meneghini2024fuse,
author = {Meneghini, O. and Slendebroek, T. and Lyons, B.C. and McLaughlin, K. and McClenaghan, J. and Stagner, L. and Harvey, J. and Neiser, T.F. and Ghiozzi, A. and Dose, G. and Guterl, J. and Zalzali, A. and Cote, T. and Shi, N. and Weisberg, D. and Smith, S.P. and Grierson, B.A. and Candy, J.},
doi = {10.48550/arXiv.2409.05894},
journal = {arXiv},
title = {{FUSE (Fusion Synthesis Engine): A Next Generation Framework for Integrated Design of Fusion Pilot Plants}},
year = {2024}
}
```
