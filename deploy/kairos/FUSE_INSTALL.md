# FUSE Installation Guide for Kairos HPC System

## Overview
This guide provides step-by-step instructions for installing the FUSE package and compiling it into a system image on the Kairos HPC server.

## Requirements
- Julia module available
- GCC module available
- Network access through proxy (203.230.124.39:4128)
- Sufficient disk space (~10 GB per installation)

## Environment Variables

Configure the following environment variables before installation:

| Variable | Description | Default Value |
|----------|-------------|---------------|
| `FUSE_HOME` | Base directory for FUSE installation | `/opt/nfri/glib/fuse` |
| `FUSE_MODULE_DIR` | Directory for module files | `/opt/nfri/glib/modulefiles` |
| `FUSE_ENVIRONMENT` | FUSE version/environment name (optional) | Auto-detected from Project.toml (latest release, currently v0.8.10) |

## Installation Steps

### 1. Clone or Copy FUSE Repository
**Note:** Currently, you need to checkout the `deploy/kairos` branch. Once this branch is merged into the main repository, you can clone without specifying the branch.

```bash
git clone https://github.com/ProjectTorreyPines/FUSE.jl.git

cd FUSE.jl/deploy/kairos
```

### 2. Set Environment Variables (Optional)
Only set if you want to override the defaults:
```bash
export FUSE_HOME="/path/to/your/fuse/installation"
export FUSE_MODULE_DIR="/path/to/your/modulefiles"
# export FUSE_ENVIRONMENT="v0.8.10"  # Only set to install a specific version instead of latest
```

### 3. Run Installation Script
**Note:** The installation process takes approximately 2 hours to complete.

```bash
./deploy.sh
```

The script will:
- Load required modules (Julia, GCC)
- Configure proxy settings for Kairos network
- Install FUSE and dependencies
- Create precompiled system image
- Generate module file for easy loading
- Create `fuse` executable

### 4. Verify Installation
Once complete, you'll see:
```
FUSE [version] installed successfully!
```

## Usage

```bash
# 'module use' is only needed when using non-default modulefiles
module use /path/to/your/modulefiles
module load fuse
fuse  # Julia REPL will start in a few seconds (~10s)
```


### Example Session
The FUSE-optimized Julia REPL runs most commands with minimal compilation time thanks to the precompiled system image:

```julia
yoom@kairos1:~/FUSE> fuse
  _  __               _ _
 (_)/ _|             (_(_) |  Documentation: https://fuse.help
(_)| |_ _   _ ___  _(_(_)  |
   |  _| | | / __|/ _ \    |  Julia REPL with FUSE v0.8.10 sysimage (Kairos)
  _| | | |_| \__ \  __/    |  System installation - packages read-only
 (_)_|  \__,_|___/\___|_   |  Create user environment for custom packages
(_(_)                 (_)  |

julia> @time using FUSE
  0.052444 seconds (173 allocations: 12.648 KiB)

julia> @time ini, act = FUSE.case_parameters(:KDEMO);
  0.805015 seconds (18.23 k allocations: 814.219 KiB, 0.16% compilation time)

julia> @time dd = FUSE.init(ini, act);
actors: Equilibrium
actors:  TEQUILA
actors: CXbuild
actors: HCD
actors:  SimpleEC
actors:  SimpleIC
actors:  NeutralFueling
actors: Current
actors:  QED
actors: PassiveStructures
 23.941435 seconds (134.85 M allocations: 5.063 GiB, 9.07% gc time, 0.08% compilation time)
```


## Notes
- Each version is installed in a separate environment
- Package files are read-only after installation to ensure consistency
- The system image is optimized for Intel Cascade Lake CPUs (Xeon Platinum 8260)