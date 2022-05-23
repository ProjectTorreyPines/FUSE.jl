# Key concepts

In FUSE:

* Data is organized according to the [ITER IMAS hierarchy](https://gafusion.github.io/omas/schema.html)
  * `dd = IMAS.dd()` (which stands for "data dictionary" is the root of the data structure

* To first populate the data dictionary and run any actors FUSE one can:

  * Manually fill the `dd` data structure

  * Read in an existing OMAS JSON data structure

  * Use the `FUSE.Init(dd, ini)` method to populates `dd` starting from 0D `ini` parameters (same spirit of OMFIT's PRO_create module)

* Physics and engineering "actors" are the building blocks of simulations and operate **exclusively** on the `dd data dictionary and their functionality is controlled via `act` parameters

* Both the `ini` and `act` parameters structures can be thought of a glorified namelists, which can be:

  * Filled starting from scratch

  * Taken (and modified) from pre-defined cases (see under the `FUSE/cases/`  folder)

  * Taken (and modified) from GASC outputs

![image](https://user-images.githubusercontent.com/1537880/167070559-aeb20212-de01-4fff-ba68-4ebe70cc2b18.png)
