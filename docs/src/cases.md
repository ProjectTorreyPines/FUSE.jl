# Use cases

```@meta
CurrentModule = FUSE
```

FUSE comes with a set of pre-cooked use cases.
The `case_parameters(:use_case, ...)` method returns the `ini` and `act` parameters for that specific `use_case`.
These `ini` and `act` can then be further customized before running a FUSE simulation.

To create your own case and add them to `src/cases/` copy one of the other cases as a template and change the ini/act parameters inside.

!!! tip "Tip!"
    Click on the `Source` button of each use case to see how each is setup

```@docs
case_parameters
```