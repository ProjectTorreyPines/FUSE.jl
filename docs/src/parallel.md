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

On your computer setup your `~/.ssh.config` this way:
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

Then on your computer open a web browser tab to `localhost:33445`
<img width="212" alt="image" src="https://user-images.githubusercontent.com/1537880/166413141-a3db3746-d603-4895-96c1-4e33d4491881.png">

Now in an interactive Jupyter-lab one can use `saga` nodes with:
```julia
nodes = 4
np = 30 * nodes
import Distributed
import ClusterManagers
Distributed.addprocs(ClusterManagers.SlurmManager(np), exclusive="")
```