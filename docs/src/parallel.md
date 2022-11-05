# Parallel jobs

## Parallel runs on Omega cluster

```julia
using Distributed
using ClusterManagers

np = 500
ENV["JULIA_WORKER_TIMEOUT"] = "180"
addprocs(SlurmManager(np-nprocs()+1), partition="preemptable", ntasks_per_core=1, mem_per_cpu="4G", time="99:99:99", topology=:master_worker)
```

## Parallel Jupyter on Saga cluster

On `cybele` via NoMachine connect to `saga`

```
jupyter lab --no-browser --port 55667
```

Copy the token that you see on this session it should look something like ```token=1f1e0259cbc1..................```


On your computer setup your `~/.ssh/config` this way:
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

and start a tunnel going through `cybele` to `saga` (keep this open)
```
ssh -N -L localhost:33445:localhost:55667 sagae
```

Use the token that you coppied from your nomachine session and use it in the http://localhost:33445 login page

Now, on your computer open a web browser tab to `localhost:33445` to connect to the Jupyter-lab session on `saga`.

And finally, in this interactive Jupyter-lab one can use `saga` nodes with:
```julia
nodes = 4
np = 30 * nodes
import Distributed
import ClusterManagers
Distributed.addprocs(ClusterManagers.SlurmManager(np), exclusive="", topology=:master_worker)
```