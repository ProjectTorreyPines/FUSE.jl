# GA Systems Code

```@meta
CurrentModule = FUSE
```

FUSE simulations can be initialized from the results of GASC simulations, which is done in two steps:
1. The GASC simulation is first loaded as a `FUSE.GASC` data structure
2. The `FUSE.GASC` structure is converted into `ini` and `act` parameters, which can then fed to FUSE actors/workflows as usual

This is accomplished this way (see the [FPP example](https://github.com/ProjectTorreyPines/FUSE.jl/blob/master/cases/FPP.jl)):
```julia
gasc = FUSE.GASC("gasc_file.json", gasc_simulation_number)
ini, act = case_parameters(gasc)
```

!!! note
    To be read by FUSE the GASC Python pickle output files need to first converted to Json format.
    This can be done with the following piece of code in OMFIT.
    ```python
    filename = "path_to_the_gasc_output.pkl"
    casename = os.path.splitext(os.path.basename(filename))[1]
    json = OMFITjson(casename + ".json", objects_encode=False)
    json.update(OMFITpickle(filename))
    json.deploy()
    ```

Conversion of GASC ouput to FUSE `ini` and `act` parameters is done via the following functions:

```@docs
    FUSE.gasc_2_equilibrium
    FUSE.gasc_2_sources
    FUSE.gasc_2_core_profiles
    FUSE.gasc_2_build
    FUSE.gasc_2_layers
    FUSE.gasc_2_coil_technology
```
