# FUSE.jl

https://fuse.help

## Installation

FUSE and related packages are registered at the [FuseRegistry](https://github.com/ProjectTorreyPines/FuseRegistry.jl/).

For installation:

```julia
using Pkg
Pkg.Registry.add(RegistrySpec(url="https://github.com/ProjectTorreyPines/FuseRegistry.jl.git"))
Pkg.Registry.add("General")
Pkg.add("FUSE")
```