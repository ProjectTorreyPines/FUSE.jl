## Getting started on the OMEGA cluster

1. Get a OMEGA account and ask to be added to the `ptp` UNIX group

2. Create a directory with your username under `/fusion/ga/projects/ird/ptp`
   ```
   mkdir /fusion/ga/projects/ird/ptp/$USER
   ```

3. Install miniconda
   ```
   cd /fusion/ga/projects/ird/ptp/$USER
   wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
   sh Miniconda3-latest-Linux-x86_64.sh
   ```
   read and accept the license, and install under `/fusion/ga/projects/ird/ptp/$USER/miniconda3`, answer questions, and restart your shell

4. install `mamba` for faster package management
   ```
   /fusion/ga/projects/ird/ptp/$USER/miniconda3/bin/conda install -c conda-forge mamba
   ```
   !!! note
       We use the full `conda` path to avoid picking up the system `conda` install. There is no system-wide `mamba` executable, so that's not necessary when running `mamba`.

5. install `jupyterlab`
   ```
   mamba install -c conda-forge jupyterlab
   ```

6. Setup your environment to run CHEASE (optional)
   ```
   export PATH=$PATH:/fusion/ga/projects/ird/ptp/chease/src-f90
   ```

7. Create a symbolic link from `/fusion/ga/projects/ird/ptp/$USER/julia/` to `~/.julia`
   ```
   mkdir -p /fusion/ga/projects/ird/ptp/$USER/julia/dev
   ln -s /fusion/ga/projects/ird/ptp/$USER/julia ~/.julia
   ```
   `~/.julia` is where the Julia will install itself by default, and this will trick it to install itself in the IR&D folder instead.

   For convenience create also a symbolic link in your `$HOME` that points to the Julia dev folder:
   ```
   ln -s /fusion/ga/projects/ird/ptp/$USER/julia/dev ~/julia_dev
   ```

8. Now follow the standard Julia and FUSE installation instructions

9. Hide your Julia processes on OMEGA
   ```
   cd ~/julia_dev/FUSE
   make hide
   ```
   OMEGA is a shared resource used by non-GA employees. This will rename all your `julia` session into a much more discrete `python3`.

10. Setup a multi-threaded Jupyter Julia kernel that does not take the whole login node
   ```
   export JULIA_NUM_THREADS=10
   cd ~/julia_dev/FUSE
   make IJulia
   ```
   OMEGA is a shared resource. Each login node has 40 cores. This will setup a Jupyter Julia kernel with 10 threads.

## Jupyter on OMEGA cluster

1. Connect to `omega` and launch `screen`

   !!! note
       You can re-connect to an existing `screen` session with `screen -r`

2. Then start the Jupyter lab server from the `screen` session (`screen` will keep `jupyter` running even when you log out)
   ```
   jupyter lab --no-browser --port 55667
   ```

   Copy the token that you see on this session it should look something like ```token=1f1e0259cbc1..................```

3. On your computer setup your `~/.ssh/config` this way (need to do this only once):
   ```
   Host cybele cybele.gat.com
   Hostname cybele.gat.com
   User meneghini
   Port 2039

   Host omegae omega.gat.com
   Hostname omega.gat.com
   User meneghini
   ProxyCommand ssh -q cybele nc %h %p
   ```

4. On your computer start a tunnel going through `cybele` to `omega`
   ```
   ssh -N -L localhost:33445:localhost:55667 omegae
   ```
   !!! note
       Keep this terminal always open. You may need to re-issue this command whenever you put your laptop to sleep.

5. On your computer open a web browser tab to `localhost:33445` to connect to the Jupyter-lab session on `omega`. Use the token when prompted.
