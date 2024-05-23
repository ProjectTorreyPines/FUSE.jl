## Getting started on the OMEGA cluster

1. Get a OMEGA account and ask to be added to the `ptp` UNIX group

1. Create a directory with your username under `/fusion/ga/projects/ird/ptp`
   ```
   mkdir /fusion/ga/projects/ird/ptp/$USER
   ```

1. Install miniconda
   ```
   cd /fusion/ga/projects/ird/ptp/$USER
   wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
   sh Miniconda3-latest-Linux-x86_64.sh
   ```
   read and accept the license, and install under `/fusion/ga/projects/ird/ptp/$USER/miniconda3`, answer questions, and restart your shell

1. install `mamba` for faster package management
   ```
   /fusion/ga/projects/ird/ptp/$USER/miniconda3/bin/conda install -c conda-forge mamba
   ```
   !!! note
       We use the full `conda` path to avoid picking up the system `conda` install. There is no system-wide `mamba` executable, so that's not necessary when running `mamba`.

1. install `jupyterlab`
   ```
   mamba install -c conda-forge jupyterlab
   ```

1. Setup your environment to run CHEASE (optional)
   ```
   export PATH=$PATH:/fusion/ga/projects/ird/ptp/chease/src-f90
   ```

1. Create a symbolic link from `/fusion/ga/projects/ird/ptp/$USER/julia/` to `~/.julia`
   ```
   mkdir -p /fusion/ga/projects/ird/ptp/$USER/julia/dev
   ln -s /fusion/ga/projects/ird/ptp/$USER/julia ~/.julia
   ```
   `~/.julia` is where the Julia will install itself by default, and this will trick it to install itself in the IR&D folder instead.

   For convenience create also a symbolic link in your `$HOME` that points to the Julia dev folder:
   ```
   ln -s /fusion/ga/projects/ird/ptp/$USER/julia/dev ~/julia_dev
   ```

1. Remove `module load defaults` from your `~/.bashrc`
   This module is used to run experimental tools like review+, efit_veiwer, etc...
   but it does not play well with the Julia executable.
   (alternatively you'll have to `module purge` or `module unload defaults`)

1. Now follow the standard Julia and FUSE installation instructions

1. Setup a multi-threaded Jupyter Julia kernel that does not take the whole login node
   ```
   cd ~/julia_dev/FUSE

   export JULIA_NUM_THREADS=10
   make IJulia

   export JULIA_NUM_THREADS=40
   make IJulia
   ```
   OMEGA login nodes are a shared resource. Each login node has 40 cores.
   This will setup a Jupyter Julia kernel with both 10 and 40 threads.
   Use 10 threads on login nodes and 40 threads on worker nodes.

## Jupyter on OMEGA cluster

1. Connect to `omega` and launch `screen`

   !!! note
       You can re-connect to an existing `screen` session with `screen -r`

1. For larger jobs especially, consider doing what follows on a login node.
   To do this (eg. a whole worker node of 40 cores for 4 days)

    srun --partition=ga-ird --nodes=1 --time=4-00:00:00 --pty bash -l

   !!! note
       Use the queue, time, cpu, and memory limits that make the most sense for your application
       see these [instructions](https://fusionga.sharepoint.com/sites/Computing/SitePages/Omega.aspx#using-slurm-to-run-interactive-tasks%E2%80%8B%E2%80%8B%E2%80%8B%E2%80%8B%E2%80%8B%E2%80%8B%E2%80%8B) for help

1. Then start the Jupyter lab server from the `screen` session (`screen` will keep `jupyter` running even when you log out)
   ```
   jupyter lab --no-browser --port 55667
   ```

   Copy the token that you see on this session it should look something like ```token=1f1e0259cbc1..................```

1. On your computer setup your `~/.ssh/config` this way (need to do this only once):
   ```
   Host cybele cybele.gat.com
   Hostname cybele.gat.com
   User meneghini
   Port 2039

   Host omegae omega.gat.com
   Hostname omega.gat.com
   User meneghini
   ProxyCommand ssh -q cybele nc %h %p

   # change XX to the worker node number you've been assigned to
   Host omegaXX omegaXX.gat.com
   Hostname omegaXX.gat.com
   User meneghini
   ProxyCommand ssh -q cybele nc %h %p
   ```

1. On your computer start a tunnel going through `cybele` to `omega`
   ```
   ssh -N -L localhost:33445:localhost:55667 omegae
   ```
   !!! note
       Keep this terminal always open. You may need to re-issue this command whenever you put your laptop to sleep.

1. On your computer open a web browser tab to `localhost:33445` to connect to the Jupyter-lab session on `omega`. Use the token when prompted.

## Using Revise on OMEGA
When working on omega it seems ones need to manually trigger revise to pick up code changes:
```
import Revise
Revise.revise()  # manual trigger
```

This is even if setting [`JULIA_REVISE_POLL=1`](https://timholy.github.io/Revise.jl/stable/config/#Polling-and-NFS-mounted-code-directories:-JULIA_REVISE_POLL)

## Using GACODE on OMEGA with Julia
Julia may be incompatible with some environments and will crash when launched.
This is the case for the GACODE environment on OMEGA.
To be able to run both GACODE and Julia on OMEGA (eg. to run NEO and TGLF) do the following:
```
module load atom
module unload gcc
module unload env
```
