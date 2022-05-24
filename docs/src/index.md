# Getting started

![image](https://user-images.githubusercontent.com/1537880/167070559-aeb20212-de01-4fff-ba68-4ebe70cc2b18.png)

The FUsion Synthesis Engine (FUSE) framework is built with the following ideas in mind:

* **Data** is organized according to the ITER IMAS ontology

  * `dd = IMAS.dd()` (which stands for "data dictionary") is the root of the FUSE data structure

* Physics and engineering **actors** are the fundamental building blocks of FUSE simulations: 

  * Actors operate exclusively on the `dd` data dictionary
  
  * Actors functionality is controlled via `act` parameters

  * Actors can be combined into other actors

* To first populate the data dictionary and run any actors FUSE one can **initialize** `dd`:

  * Manually

  * By reading in an existing OMAS JSON data structure

  * Using the `FUSE.Init(dd, ini, act)` method to populates `dd` starting from 0D `ini` parameters (same spirit of OMFIT's PRO_create module) and `act` Actor parameters

* Both the `ini` and `act` **parameters** can be thought of a glorified namelists, which can be:

  * Filled starting from scratch

  * Taken (and modified) from pre-defined cases (see under the `FUSE/cases/` folder)

  * Taken (and modified) from GASC outputs

* **Workflows** perform studies that start from `ini` and `act` parameters
