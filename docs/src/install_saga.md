# On SAGA cluster

## Getting started on the SAGA cluster

1. Get a SAGA account and ask to have a directory created for you under `/mnt/beegfs/users/$USER`

2. Install miniconda
   ```
   cd /mnt/beegfs/users/$USER
   wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
   sh Miniconda3-latest-Linux-x86_64.sh
   ```
   read and accept the license, and install under `/mnt/beegfs/users/$USER/miniconda3`, answer questions, and restart your shell

3. install `mamba` for faster package management
   ```
   conda install -c conda-forge mamba
   ```

4. install `jupyterlab`
   ```
   mamba install -c conda-forge jupyterlab
   ```
5. Setup your environment to run CHEASE (optional)
   ```
   export PATH=$PATH:/mnt/beegfs/users/meneghini/chease/src-f90
   ```

6. Create a symbolic link from `/mnt/beegfs/users/$USER/julia/` to `~/.julia`
   ```
   mkdir -p /mnt/beegfs/users/$USER/julia/dev
   ln -s /mnt/beegfs/users/$USER/julia ~/.julia
   ```
   `~/.julia` is where the Julia will install itself by default, and this will trick it to install itself in the IR&D folder instead.

   For convenience create also a symbolic link in your `$HOME` that points to the Julia `dev` folder:
   ```
   ln -s /mnt/beegfs/users/$USER/julia/dev ~/julia_dev
   ```

8. Now follow the standard Julia and FUSE installation instructions

## Jupyter on SAGA cluster

1. Connect to `saga` and launch `screen`

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

   Host sagae saga.gat.com
   Hostname saga.gat.com
   User meneghini
   ProxyCommand ssh -q cybele nc %h %p
   ```

4. On your computer start a tunnel going through `cybele` to `saga`
   ```
   ssh -N -L localhost:33445:localhost:55667 sagae
   ```
   !!! note
       Keep this terminal always open. You may need to re-issue this command whenever you put your laptop to sleep.

5. On your computer open a web browser tab to `localhost:33445` to connect to the Jupyter-lab session on `saga`. Use the token when prompted.
