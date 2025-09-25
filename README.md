[![Ask DeepWiki](https://deepwiki.com/badge.svg)](https://deepwiki.com/ProjectTorreyPines/FUSE.jl)

# FUSE.jl

FUSE (**FU**sion **S**ynthesis **E**ngine) is an open-source framework for integrated simulations of fusion devices.
Originally developed by General Atomics, FUSE is now publicly available under the [Apache 2.0 license](https://fuse.help/dev/notice.html).

The figure below is a sneakpeak of the models implemented in FUSE:
[![FUSE capabilities](https://raw.githubusercontent.com/ProjectTorreyPines/FUSE_extra_files/refs/heads/master/carousel.jpg)](https://raw.githubusercontent.com/ProjectTorreyPines/FUSE_extra_files/refs/heads/master/carousel.jpg)

## Resources

Here are some key resources for getting started with FUSE:

* 📨 **[Sign up for more info](https://forms.gle/iQGYKWfgeNkw2fCZ7)**
* 📚 **[Online documentation](https://fuse.help)**
* 🎓 **[Intro tutorial](https://fuse.help/dev/tutorial.html)**
* 🎤 **[Recent presentation](https://github.com/ProjectTorreyPines/FUSE_extra_files/raw/master/2025_D3D/SET_mar_2025.pdf)**
* 📜 **[Preprint publication](https://arxiv.org/abs/2409.05894)**
* 🆘 **[Discord community](https://discord.gg/CbjpZH9SKM)**
* 🗒️ **[Weekly devs meeting minutes](https://github.com/ProjectTorreyPines/FUSE.jl/discussions)**

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
6. **🔄 Studies**: Studies pre-bake commonly used functionalities, typically involving multiple FUSE simulations (eg. database generation or multi objective optimizations).
7. **🌍 Interoperability**: FUSE interfaces via `dd` with existing modeling tools like OMFIT/OMAS and the IMAS ecosystem.

A diagram illustrating these concepts is provided below:
![FUSE Diagram](./docs/src/assets/FUSE.svg)

## Usage Example

Here’s a simple example of setting up and running a FUSE simulation in Julia:

```julia
using FUSE

# Obtain `ini` and `act` parameters for a specific use case
ini, act = FUSE.case_parameters(:FPP)

# Initialize the `dd` structure from 0D `ini` parameters
dd = FUSE.init(ini, act)

# Run a stationary plasma actor simulation
FUSE.ActorStationaryPlasma(dd, act)

# Get an overview of the simulation results
FUSE.digest(dd)
```

Make sure to take a look at the [introductory tutorial](https://fuse.help/dev/tutorial.html) and [examples](https://github.com/ProjectTorreyPines/FuseExamples).

## Installation

Follow these [installation instructions](https://fuse.help/dev/install.html).

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
