# Parallel jobs

## Parallel runs on OMEGA cluster

```julia
using Distributed
using ClusterManagers

np = 500
ENV["JULIA_WORKER_TIMEOUT"] = "180"
addprocs(SlurmManager(np-nprocs()+1), partition="preemptable", ntasks_per_core=1, mem_per_cpu="4G", time="99:99:99", topology=:master_worker)
```

## Parallel Jupyter on SAGA cluster

On `cybele` via NoMachine connect to `saga`

```
jupyter lab --no-browser --port 55667
```

Copy the token that you see on this session it should look something like ```token=1f1e0259cbc1..................```

On your computer setup your `~/.ssh/config` this way (need to do this only once):
```
Host cybele cybele.gat.com
   Hostname cybele.gat.com
   User meneghini
   Port 2039

Host sagae saga.gat.com
   Hostname saga.gat.com
   User meneghini
   ProxyCommand ssh -q cybele nc %h %p
```

and start a tunnel going through `cybele` to `saga` (keep this terminal always open)
```
ssh -N -L localhost:33445:localhost:55667 sagae
```

Now, on your computer open a web browser tab to `localhost:33445` to connect to the Jupyter-lab session on `saga`.
Use the token when prompted.

And finally, in this interactive Jupyter-lab one can use `saga` nodes with:
```julia
nodes = 4
np = 30 * nodes
import Distributed
import ClusterManagers
Distributed.addprocs(ClusterManagers.SlurmManager(np), exclusive="", topology=:master_worker)
```

## Getting started on the SAGA cluster

1. Get a SAGA account and ask to have a directory created for you under `/mnt/beegfs/users`

2. Install miniconda
   ```
   cd /mnt/beegfs/users/$USER
   wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
   sh Miniconda3-latest-Linux-x86_64.sh
   ```
   read and accept the license, and install under `/mnt/beegfs/users/$USER/miniconda3`, answer questions, and restart your shell

3. install `mamba` for faster package management, and then the rest
   ```
   conda install -c conda-forge mamba
   mamba install -c conda-forge gfortran # to compile CHEASE
   mamba install -c conda-forge jupyterlab
   ```
