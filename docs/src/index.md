# Getting started

![image](https://user-images.githubusercontent.com/1537880/167070559-aeb20212-de01-4fff-ba68-4ebe70cc2b18.png)

The FUsion Synthesis Engine (FUSE) framework is built with the following ideas in mind:

* **Data** is organized according to the ITER IMAS ontology

* Physics and engineering **actors** are the fundamental building blocks of FUSE simulations

* `ini` parameters allow to conveniently populate the data dictionary to begin a FUSE simulation

* Actors functionality is controlled via `act` parameters

* Both the `ini` and `act` **parameters** can be thought of a glorified namelists, which can be:

  * Filled starting from scratch

  * Taken (and modified) from pre-defined cases (see under the `FUSE/cases/` folder)

  * Taken (and modified) from GASC outputs

* **Workflows** perform studies that start from `ini` and `act` parameters
