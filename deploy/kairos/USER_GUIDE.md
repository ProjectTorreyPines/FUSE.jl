# FUSE User Guide for Kairos HPC

This guide explains how to use the FUSE environment on Kairos HPC system.

## Module Usage

FUSE is now officially installed system-wide on Kairos HPC system.

### Standard Usage

To use FUSE:

```bash
# Add the module path
module use /opt/nfri/glib/modulefiles

# Load FUSE module
module load fuse

# Start FUSE-optimized Julia REPL
fuse
```

The `fuse` command will start a Julia REPL with FUSE preloaded. Initial startup takes ~10 seconds.

## Working with FUSE

### Quick Start Example

The FUSE environment comes with a precompiled system image, providing fast execution with minimal compilation time:

```julia
yoom@kairos1:~> module use /opt/nfri/glib/modulefiles
yoom@kairos1:~> module load fuse
FUSE v0.8.10 environment loaded
Documentation: https://fuse.help
To start: fuse

yoom@kairos1:~> fuse
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

Note the minimal compilation time (0.08-0.16%) - most functions are already compiled in the system image.


### GUI Usage

For graphical output on Kairos, use TigerVNC for remote desktop access. 
Note: Jupyter Notebook setting for FUSE is currently unavailable due to HTTPS proxy restrictions on Kairos.

```julia
julia> using Plots

# Plots are not displayed by default (headless mode)
julia> p1 = plot(dd.equilibrium)
Plot{Plots.GRBackend() n=48}
Captured extra kwargs:
  Series{47}:
    coordinate: psi_norm

# Save plot to file
julia> savefig(p1, "test.png")
"/home/yoom/test.png"

# Display plot in GUI window (requires VNC session)
julia> plot(dd.equilibrium, show=true);

# Set default to always show plots
julia> default(show=true)

julia> plot(dd.equilibrium);  # Now displays automatically

```

### Available Resources

- **Documentation**: https://fuse.help
- **Precompiled Packages**: FUSE, EFIT, Plots, and all dependencies
- **Performance**: Optimized for Intel Cascade Lake CPUs (Xeon Platinum 8260)

### Tips for Efficient Usage

1. **Fast Startup**: The `fuse` command provides a Julia REPL with precompiled FUSE, eliminating first-run compilation delays

2. **Thread Configuration**: Default configuration uses optimal thread settings for Kairos nodes

3. **Custom Packages**: While the system FUSE installation is read-only, you can still add personal packages to your user Julia depot

## Module Information

To check available FUSE versions:
```bash
module avail fuse
```

To see module details:
```bash
module show fuse
```

## Support

For issues or questions about FUSE usage on Kairos, contact the FUSE development team or refer to the documentation at https://fuse.help