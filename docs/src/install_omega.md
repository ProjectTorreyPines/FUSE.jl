# On GA's OMEGA cluster

!!! note
    FUSE should only be run on `somega` or on a worker node, never on the `omega` login nodes.

## FUSE module

If you only intend to use FUSE and don't plan to develop the source code, release versions of FUSE
have been installed on the Omega cluster. All available versions can be found with `module avail fuse`.
To load the latest version, do:
```
module load fuse
```

The FUSE module does several things:

1. It exposes the `fuse` executable, which runs Julia in a read-only environment with the FUSE
   codebase precompiled into a
   [sysimage](https://julialang.github.io/PackageCompiler.jl/dev/sysimages.html).
   This greatly reduces the time-to-first-execution (TTFX) for many functions in the FUSE code suite,
   at the expense of "locking" those packages and functions to the versions with which they
   were compiled.

1. The Julia module is also loaded, so development work can be done by running the `julia` command.
   This gives you access to a public Julia installation with your own private "depot" in which each
   user can add or develop their own packages. The location of this private depot is given by the
   environment variable `JULIA_USER_DEPOT`.

1. FUSE is already available when you launch `julia` or `fuse`, so there's no need to do
   `Pkg.add("FUSE")`. You can simply do `using FUSE` and begin working.

1. A custom conda installation is made available to you that has Jupyterlab with
   precompiled Julia kernels that include the FUSE sysimage. You can just do `jupyter lab` to
   start a Jupyter session and select the desired kernels. There is a kernel with 1 thread meant
   for the login nodes and 10 & 16 threads for the worker nodes.
   !!! warning
       **Potential Problem**: There's a bug that may occur when a new user first launches one of these
       Julia + FUSE Jupyter kernels.
       In your terminal, you may see output about precompiling IJulia, which is expected.
       Once the precompilation is done, it will report `Starting kernel event loops` but then the
       kernel may hang and your notebook may not work. It is unclear why this happens, but it is
       only the first time for each user.

       **Solution**: The latest installations on Omega seemed to have fixed this issue, but if it happens,
       restart the kernel. Occasionally this needs to be done twice, perhaps if you
       restart too quickly and the precompilation was not finished. In any case, if the problem
       does not resolve after restarting the kernel twice, reach out to the FUSE developers.

## Connecting to a Jupyter-lab server running on OMEGA from your laptop

1. Connect to `omega` and launch `screen`

   !!! note
       You can re-connect to an existing `screen` session with `screen -r`

1. **If (and only if) you want to run jupyter-lab on a worker node** do as follows:

    `srun --partition=ga-ird --nodes=1 --time=4-00:00:00 --pty bash -l`

   !!! note
       Use the queue, time, CPU, and memory limits that make the most sense for your application
       see these [instructions](https://fusionga.sharepoint.com/sites/Computing/SitePages/Omega.aspx#using-slurm-to-run-interactive-tasks%E2%80%8B%E2%80%8B%E2%80%8B%E2%80%8B%E2%80%8B%E2%80%8B%E2%80%8B) for help

1. Then start the Jupyter lab server from the `screen` session (`screen` will keep `jupyter`
   running even when you log out)
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
       Keep this terminal always open. You may need to re-issue this command whenever you put your
       laptop to sleep.

1. On your computer open a web browser tab to `localhost:33445` to connect to the Jupyter-lab
   session on `omega`. Use the token when prompted.

## Three ways to run parallel jobs

Keep in mind that each worker node on Omega has 128 CPUs

1. Screen + Jupyter on the login node, workers on the worker nodes

   OK when the master process will not be doing a lot of work, and we need multiple nodes

   Here we will use the `FUSE.parallel_environment("omega", ...)` call.

1. Screen on the login node, Jupyter and workers on one worker node

   OK when the master process will be doing a lot of work, and we don't need more than one node

   Here we will use the `FUSE.parallel_environment("localhost", ...)` call.

1. Screen on the login node, Jupyter on a worker node, workers on different worker nodes

   OK when the master process will be doing a lot of work, and we need multiple nodes

   This is more complex, and finicky. Avoid if possible.

   Here we will use the `FUSE.parallel_environment("omega", ...)` call.

## Running Julia/FUSE via VScode

You can use [VScode and remote SHH connection](https://code.visualstudio.com/docs/remote/ssh-tutorial) to run Julia and FUSE on Omega.

To do this:

* Open a remote connection to Omega

* Install the `Julia` and the `Jupyter` VScode extensions in on Omega

* Open the `Code > Settings... > Settings` menu

  * Select the `Remote [SSH: omega]` tab

  * Search for `julia executable` in the search bar

  * Edit the `julia: Executable Path` to `/fusion/projects/codes/julia/fuse/julia_with_fuse`

Now Julia scripts and notebooks can be run directly from this remote VScode session.

## Using Revise on Omega
When working on omega it seems ones need to manually trigger revise to pick up code changes:
```
import Revise
Revise.revise()  # manual trigger
```

This is even if setting [`JULIA_REVISE_POLL=1`](https://timholy.github.io/Revise.jl/stable/config/#Polling-and-NFS-mounted-code-directories:-JULIA_REVISE_POLL)

## Using GACODE on Omega with Julia
Julia may be incompatible with some environments and will crash when launched.
This is the case for the GACODE environment on Omega.
To be able to run both GACODE and Julia on Omega (eg. to run NEO and TGLF) do the following:
```
module load atom
module unload gcc
module unload env
```
