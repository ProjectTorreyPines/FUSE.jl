# Use cases

```@meta
CurrentModule = FUSE
```

FUSE comes with a set of pre-cooked use cases.
The `case_parameters(:use_case, ...)` method returns the `ini` and `act` parameters for that specific `use_case`.
These `ini` and `act` can then be further customized before running a FUSE simulation.

To create your own case and add them to `FUSE/cases` copy one of the other cases as a template and change the ini/act parameters inside.
A handy way of generating the ini code from `your_ini` that you created in a notebook or elsewhere is to call the function `FUSE.case_parameter_creation_from_ini(your_ini)` which will return a nicely formatted code snippet.

!!! tip "Tip!"
    Click on the `Source` button of each use case to see how each is setup

```@docs
case_parameters
```