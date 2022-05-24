# Getting started

![image](https://user-images.githubusercontent.com/1537880/167070559-aeb20212-de01-4fff-ba68-4ebe70cc2b18.png)

The FUsion Synthesis Engine (FUSE) framework is built with the following ideas in mind:
* Data is organized according to the ITER IMAS `data dictionary` onthology
* Physics and engineering `actors` are the fundamental building blocks of FUSE simulations
* `ini` parameters allow to conveniently populate the data dictionary to begin a FUSE simulation
* Actors functionality is controlled via `act` parameters

These concepts are illustrated in this simple example:
```julia
using FUSE
ini, act = FUSE.case_parameters(:FPP; version=:v1, init_from=:scalars)
dd = IMAS.dd()
FUSE.init(dd, ini, act; do_plot=true)
FUSE.ActorEquilibriumTransport(dd, act; do_plot=true)
```

In addition FUSE defines a series of **workflows** that are used to perform self-contained studies/optimizations, typically running many FUSE simulations.