# How to Use FUSE - Knowledge Base

## Overview

FUSE is a comprehensive fusion energy modeling framework that follows the ITER IMAS ontology. This document captures key concepts, workflows, and best practices learned from studying FUSE tutorials and examples.

## Core Concepts

### 1. Data Storage (`dd`) - IMAS Data Dictionary
- All data is stored in the `dd` structure following ITER IMAS ontology
- Create empty data dictionary: `dd = IMAS.dd()`
- The `dd` contains various Interface Data Structures (IDSs): equilibrium, core_profiles, build, etc.
- View available IDSs: `keys(dd)`
- Get help for specific IDS: `help(dd.equilibrium; maxdepth=2)`

#### IMAS Data Fundamentals
- **IDSs (Interface Data Structures)**: Standardized data structures for tokamak data
- **Data Dictionary (`dd`)**: Root structure containing all IDSs
- **Coordinates**: Must be assigned before data values: `IMAS.coordinates(dd.equilibrium.vacuum_toroidal_field, :b0)`  
- **Units**: Check units with: `IMAS.units(dd.equilibrium.vacuum_toroidal_field, :b0)`
- **Time coordinate**: Set time before setting time-dependent arrays: `dd.equilibrium.time = [1.0]`

### 2. Actors (🧠)
- Core components that perform physics and engineering calculations
- Execute actors by passing `dd` and `act`: `FUSE.ActorName(dd, act)`
- Actors can call other actors to create workflows
- Some actors have custom plotting methods: `plot(actor)`

### 3. Control Parameters (`act`)
- Govern how actors behave and what models they use
- Can be modified to change actor behavior: `act.ActorName.parameter = value`
- Parameters can also be passed directly to actors: `FUSE.Actor(dd, act; parameter=value)`

### 4. Initialization Parameters (`ini`)
- 0D parameters used to bootstrap the `dd` with plausible initial data
- Modify before initialization: `ini.section.parameter = value`
- Initialize data structure: `FUSE.init(dd, ini, act)`
- **Important**: `init()` provides a starting point, NOT a self-consistent solution

### 5. Use Cases
- FUSE includes predefined templates for various machines
- Available cases: `methods(FUSE.case_parameters)`
- Get parameters: `ini, act = FUSE.case_parameters(:CASE_NAME)`
- Examples: `:KDEMO`, `:ITER`, `:FPP`, `:ARC`

## Basic Workflow

### 1. Setup and Initialization
```julia
using Plots
using FUSE  # Also imports IMAS

# Get predefined case parameters
ini, act = FUSE.case_parameters(:KDEMO)

# Modify parameters as needed
ini.equilibrium.B0 = 7.8
ini.equilibrium.R0 = 6.5
act.ActorCoreTransport.model = :FluxMatcher

# Initialize data structure
dd = IMAS.dd()
FUSE.init(dd, ini, act)
```

### 2. Checkpoint System
- Save progress: `@checkin :tag_name dd ini act`
- Restore state: `@checkout :tag_name dd ini act`
- Useful for iterative development and avoiding restarts

### 3. Individual Actor Execution
```julia
# PF coil design
FUSE.ActorPFdesign(dd, act; do_plot=true)

# Self-consistent plasma solution
FUSE.ActorStationaryPlasma(dd, act)

# Build and sizing
FUSE.ActorHFSsizing(dd, act)  # High field side sizing
FUSE.ActorLFSsizing(dd, act)  # Low field side sizing
FUSE.ActorCXbuild(dd, act)    # 2D cross-section build

# Other physics actors
FUSE.ActorPassiveStructures(dd, act)
FUSE.ActorVerticalStability(dd, act)
FUSE.ActorNeutronics(dd, act)
FUSE.ActorBlanket(dd, act)
FUSE.ActorDivertors(dd, act)
FUSE.ActorBalanceOfPlant(dd, act)
FUSE.ActorCosting(dd, act)
```

### 4. Whole Facility Design
```julia
# Complete integrated design workflow
FUSE.ActorWholeFacility(dd, act)
```

## Key Actors and Their Functions

### Core Physics
- **ActorStationaryPlasma**: Iterates between transport, pedestal, equilibrium, and sources for self-consistent plasma solution
- **ActorPFdesign**: Positions PF coils for desired plasma shape
- **ActorVerticalStability**: Evaluates plasma vertical stability
- **ActorFluxMatcher**: Transport model for flux matching (configurable algorithm: `:anderson`)

### Build and Engineering  
- **ActorHFSsizing**: Adjusts OH and TF layer thickness on high field side for current/stress constraints
- **ActorLFSsizing**: Positions outer TF leg to meet ripple requirements  
- **ActorCXbuild**: Builds 2D cross-sections after sizing changes
- **ActorPassiveStructures**: Generates passive structures (vacuum vessel)

### Neutronics and Power
- **ActorNeutronics**: Calculates heat flux on first wall
- **ActorBlanket**: Optimizes first wall, breeder, shield thickness and Li6 enrichment for target TBR
- **ActorDivertors**: Calculates divertor heat flux
- **ActorBalanceOfPlant**: Optimizes cooling flow rates and electricity conversion efficiency

### Economics
- **ActorCosting**: Breaks down capital and operational costs

## Transport and Flux Matching

### Flux Matcher Overview
FUSE includes a sophisticated flux-matcher developed at General Atomics with algorithms based on TGYRO flux-matcher and advanced optimization. The flux matcher is used to find self-consistent plasma profiles by matching experimental heat and particle fluxes.

### Basic Flux Matching Setup
```julia
# Set up flux matching transport model
act.ActorCoreTransport.model = :FluxMatcher

# Configure transport model (TGLF Neural Network)
act.ActorTGLF.model = :TGLFNN                          # Options: :TGLFNN, :TJLF, :TGLF, :QLGYRO
act.ActorTGLF.tglfnn_model = "sat3_em_d3d_azf-1"       # Specific NN model
act.ActorTGLF.warn_nn_train_bounds = false             # Suppress warnings about NN bounds
act.ActorTGLF.electromagnetic = true                   # EM effects (not relevant for NNs)
act.ActorTGLF.onnx_model = false                       # Use custom PyTorch/ONNX models

# Flux matcher configuration
act.ActorFluxMatcher.rho_transport = 0.1:0.05:0.85     # Radial points for flux matching
act.ActorFluxMatcher.max_iterations = 300              # Maximum iterations
act.ActorFluxMatcher.algorithm = :custom               # Algorithm choice
act.ActorFluxMatcher.custom_algorithm = FUSE.NonlinearSolve.SimpleDFSane()
act.ActorFluxMatcher.step_size = 1.0                   # Step size for iterations
act.ActorFluxMatcher.verbose = true                    # Show progress
act.ActorFluxMatcher.evolve_rotation = :flux_match     # Options: :flux_match, :fixed
act.ActorFluxMatcher.evolve_pedestal = false           # Use experimental pedestal
act.ActorFluxMatcher.evolve_densities = :flux_match    # Options: :flux_match, :fixed
```

### Transport Models Available

#### TGLF Neural Network (TGLFNN)
- **Model**: `act.ActorTGLF.model = :TGLFNN`
- **Available models**: Various NN models trained on TGLF data
- **Example model**: `"sat3_em_d3d_azf-1"` (SAT3 electromagnetic D3D model)
- **Performance**: Fastest option, good for rapid iteration
- **Custom models**: Place ONNX models in `~/.julia/dev/TGLFNN/models/`

#### TJLF (Julia TGLF)
- **Model**: `act.ActorTGLF.model = :TJLF`
- **Description**: Julia implementation of TGLF transport model
- **Configuration**: 
  ```julia
  act.ActorTGLF.sat_rule = :sat3                       # SAT rule selection
  act.ActorTGLF.lump_ions = true                      # Ion lumping
  act.ActorTGLF.electromagnetic = true                # Include EM effects
  ```
- **Performance**: More accurate than NN, slower than NN

#### Other Transport Models
- **TGLF**: Original Fortran TGLF (requires GACODE installation)
- **QLGYRO**: Quasilinear GYRO (requires GACODE installation)
- **GKNN**: Gyrokinetic Neural Network (experimental)

### Flux Matching Algorithms
- **:custom**: Use custom NonlinearSolve algorithm
- **:simple_dfsane**: Simple DFSane algorithm
- **:basic_polyalg**: Alternative algorithm (better convergence for difficult cases)

### Flux Matching Workflow
```julia
# Initialize case
ini, act = FUSE.case_parameters(:D3D, :L_mode)  # or :H_mode
dd = IMAS.dd()
FUSE.init(dd, ini, act)

# Configure flux matcher (see configuration above)
act.ActorCoreTransport.model = :FluxMatcher
# ... configure parameters ...

# Run flux matcher
dd_result = deepcopy(dd)
actor_transport = FUSE.ActorFluxMatcher(dd_result, act)

# Compare results
plot(dd.core_profiles, label=" Experiment")
plot!(dd_result.core_profiles, label=" Flux Matched")
plot(dd_result.core_transport)  # Show transport coefficients
```

### Using Custom ODS Data
```julia
# Load experimental data from ODS/IMAS JSON
ods = IMAS.json2imas("/path/to/ods_shot_time.json"; show_warnings=false)
ods.core_profiles.profiles_1d[1].time = 4.6
ods.core_profiles.time = [4.6]

# Initialize with custom data
ini, act = FUSE.case_parameters(:D3D, ods)

# Disable heating models if not needed
act.ActorHCD.ec_model = :none      # Electron cyclotron
act.ActorHCD.ic_model = :none      # Ion cyclotron  
act.ActorHCD.lh_model = :none      # Lower hybrid
act.ActorHCD.nb_model = :none      # Neutral beam
act.ActorHCD.neutral_model = :none # Neutral particles

# Continue with flux matching workflow...
```

### Advanced Flux Matching Techniques

#### Using TGLF-NN as Preconditioner
```julia
# First pass with fast NN model
act.ActorTGLF.model = :TGLFNN
act.ActorTGLF.tglfnn_model = "sat3_em_d3d_azf-1"
dd_nn = deepcopy(dd)
FUSE.ActorFluxMatcher(dd_nn, act)

# Second pass with more accurate TJLF
act.ActorTGLF.model = :TJLF
act.ActorFluxMatcher.max_iterations = 50  # Fewer iterations needed
dd_tjlf = deepcopy(dd_nn)  # Start from NN solution
FUSE.ActorFluxMatcher(dd_tjlf, act)
```

#### Algorithm Selection Tips
- **For L-mode and negative triangularity**: Can extend `rho_transport` closer to edge
- **For difficult convergence**: Try `:basic_polyalg` instead of `:simple_dfsane`
- **For H-mode**: Standard `rho_transport = 0.1:0.05:0.85` works well

### Model Discovery
```julia
# To see available NN models, trigger error message:
act.ActorTGLF.tglfnn_model = ""  # Will show all available models in error
```

### Flux Matching Best Practices
1. **Start with standard cases** (`:D3D :L_mode`, `:D3D :H_mode`) before custom data
2. **Use NN models for rapid iteration** and parameter studies
3. **Use TJLF for final results** when accuracy is critical
4. **Monitor convergence** with `verbose = true`
5. **Adjust radial range** based on plasma regime (L-mode vs H-mode)
6. **Use preconditioning** (NN → TJLF) for difficult cases
7. **Check NN training bounds** by setting `warn_nn_train_bounds = true`

## IMAS Physics Functions

IMAS.jl provides [extensive physics functions](https://projecttorreypines.github.io/IMAS.jl/dev/api/) that operate directly on the `dd` structure. These cover the full breadth of tokamak plasma modeling:

### Core Physics Functions  
- **Field line tracing**: `IMAS.trace_field_line(eqt, r, z; max_turns=100, step_size=0.1)`
- **Flux surface calculations**: Automatic flux surface tracing and analysis
- **Transport analysis**: Advanced transport coefficient calculations
- **Heating and current drive**: NBI, ECRH, ICRH, and other auxiliary heating models
- **Nuclear fusion reactions**: Reaction rate calculations and neutron production
- **Radiation physics**: Bremsstrahlung, synchrotron, and line radiation
- **Plasma-material interactions**: Wall erosion, particle deposition, heat flux

### Benefits
- **Validated implementations**: High-performance, tested physics algorithms
- **Common foundation**: Avoid reinventing fundamental calculations
- **Consistency**: Ensures consistent results across different simulation workflows
- **Community development**: Accelerates fusion research through shared tools

## Plotting and Visualization

### Built-in Plotting with Plots.jl Recipes
FUSE uses Plots.jl recipes for visualizing IMAS data structures. Different plots are shown automatically based on the data type.

```julia
using Plots

# Plot various data structures  
plot(dd.build)           # Tokamak cross-section
plot(dd.equilibrium)     # Plasma equilibrium
plot(dd.core_profiles)   # Temperature, density profiles  
plot(dd.core_sources)    # Heating/current drive sources
plot(dd.core_transport)  # Transport coefficients
plot(dd.costing)         # Cost breakdown

# Some actors have custom plots
actor = FUSE.ActorPFdesign(dd, act)
plot(actor)
```

### Advanced Plotting Features
```julia
# Get plot help and options
help_plot(dd.equilibrium; core_profiles_overlay=true, levels_in=21, levels_out=5, 
          show_secondary_separatrix=true, coordinate=:rho_tor_norm)

# Compose plots using plot!()
plot(dd.equilibrium; color=:gray, cx=true)
plot!(dd.build.layer)
plot!(dd.pf_active)
plot!(dd.pf_passive)

# Plot individual fields with coordinates and units
plot(dd.core_profiles.profiles_1d[], :pressure_thermal)

# Customize plot attributes
plot(dd.core_profiles.profiles_1d[], :pressure_thermal; 
     label="", linewidth=2, color=:red, labelfontsize=25)

# Plot multiple fields using search
plot(findall(dd.equilibrium.time_slice[], r"\.psi"))
```

### Time-Dependent Plotting
```julia
# Interactive time-dependent plots with Interact.jl
using Interact
@manipulate for time0 in dd.equilibrium.time
    plot(dd.equilibrium; cx=true, time0)
    plot!(dd.wall)
    plot!(dd.pf_active; time0)
    plot!(dd.pf_passive; time0)
    plot!(size=(600,600))
end
```

### Field Line Tracing Example
```julia
# Physics function example: field line tracing  
eqt = dd.equilibrium.time_slice[]
r, z = 2.3457, 0.1
t1 = IMAS.trace_field_line(eqt, r, z; max_turns=100, step_size=0.1)
t2 = IMAS.trace_field_line(eqt, r, z; max_turns=100, step_size=-0.1)

plot(eqt; cx=true, primary=false)
plot!(t1.r, t1.z; label="Upwind")
plot!(t2.r, t2.z; label="Downwind")
scatter!([r], [z]; label="Starting point", color=:black)
```

### Comparison Plots
```julia
# Compare before/after
peq = plot(dd.equilibrium; label="before")
# ... run actor ...
plot!(peq, dd.equilibrium; label="after")
```

## IMAS Data Manipulation

### Working with Arrays of Structures
- **Resize before use**: `resize!(dd.equilibrium.time_slice, 1)`
- **Resize with conditions**: `source1 = resize!(dd.core_sources.source, "identifier.index" => 1)`
- **Preserve content**: Use `wipe=false` to avoid clearing data: `resize!(dd.core_sources.source, "identifier.index" => 1; wipe=false)`
- **Array management**: Use Julia array functions: `empty!()`, `pop!()`, `popat!()`, `deleteat!()`, `push!()`, `pushfirst!()`

### Time-Dependent Data Management  
#### Global Time
- Set current time: `dd.global_time = 1.0`
- Used throughout FUSE/IMAS for time-dependent operations

#### Time-Dependent Arrays of Structures
- **Add time slice at specific time**: `resize!(dd.equilibrium.time_slice, 1.0)`
- **Add at global time**: `resize!(dd.equilibrium.time_slice)` (recommended)
- **Access by index**: `eqt = dd.equilibrium.time_slice[2]`
- **Access by time**: `eqt = dd.equilibrium.time_slice[1.0]`  
- **Access at global time**: `eqt = dd.equilibrium.time_slice[]` (recommended)
- **Note**: Cannot add timeslices "in the past", uses causal nearest neighbor

#### Time-Dependent Arrays of Data
- **@ddtime macro**: Manipulate time-dependent arrays at `dd.global_time`
- **Get value**: `@ddtime(dd.equilibrium.vacuum_toroidal_field.b0)`
- **Set value**: `@ddtime(dd.equilibrium.vacuum_toroidal_field.b0=4.0)`
- **Causal interpolation**: Uses nearest neighbor, cannot interpolate back in time
- **Forward extrapolation**: Uses last available value for future times

### Dynamic Expressions
- Some fields are expressions (functions) that calculate values dynamically
- **View structure**: `print_tree(dd.core_profiles.profiles_1d[]; maxdepth=1)` 
- **Auto-evaluation**: Accessing a dynamic expression automatically evaluates it
- **Freeze expressions**: `IMAS.freeze(ids)` - evaluates all expressions in tree
- **Freeze single field**: `IMAS.freeze(ids, field::Symbol)`
- **Force re-evaluation**: `IMAS.refreeze!(ids, field)`  
- **Revert to expression**: `IMAS.empty!(ids, field::Symbol)`

### Loading and Saving Data
- **Load from JSON**: `dd = IMAS.json2imas("filename.json")`
- **Other formats available**: See [IMAS.jl IO documentation](https://projecttorreypines.github.io/IMASdd.jl/stable/api/#IO)

### Data Exploration and Search
- **Print tree structure**: `print_tree(dd.equilibrium.time_slice[1].boundary; maxdepth=1)`
- **Search fields**: `findall(dd.equilibrium.time_slice[], r"\.psi")` (uses regex)
- **Compare IDSs**: `IMAS.diff(dd1.equilibrium, dd2.equilibrium; verbose=true)`
- **Summary**: `IMAS.extract(dd)` - snapshot of 0D quantities at `dd.global_time`

## Data Structure Access

### Key IMAS Structures
- `dd.build`: Tokamak geometry and build
- `dd.build.layer`: Layer summary (has custom `show()` method)
- `dd.equilibrium`: Plasma equilibrium
- `dd.core_profiles`: Plasma profiles
- `dd.core_sources`: Heating and current drive
- `dd.core_transport`: Transport coefficients
- `dd.solid_mechanics.center_stack.stress`: Center stack stresses
- `dd.neutronics`: Neutron transport and wall loading
- `dd.blanket`: Blanket design and tritium breeding
- `dd.divertors`: Divertor design and heat removal
- `dd.balance_of_plant`: Power conversion systems
- `dd.costing`: Economic analysis

### Useful Functions
- `FUSE.digest(dd)`: Summary of key 0D quantities
- `IMAS.freeze()`: Freeze data structure state  
- `print_tree()`: Print hierarchical structure
- `IMAS.coordinates(ids, field)`: Get coordinate information
- `IMAS.units(ids, field)`: Get units for a field
- `findall(ids, regex)`: Search fields using regular expressions
- `IMAS.diff(ids1, ids2)`: Compare two IDSs and show differences
- `IMAS.extract(dd)`: Get snapshot of 0D quantities at global time

## Best Practices

### Development Workflow
1. Start with a predefined use case
2. Modify `ini` and `act` parameters as needed
3. Initialize with `FUSE.init()`
4. Use checkpoints (`@checkin`/`@checkout`) for iterative development
5. Run individual actors to understand their effects
6. Use `ActorWholeFacility` for complete integrated design

### Performance Notes
- Julia uses JIT compilation - first execution is slower, subsequent calls are fast
- Use checkpoints to avoid recomputation
- Consider using RemoteREPL.jl for packages with long startup times

### Parameter Configuration
- Always check available parameters in documentation
- Some actors support runtime parameter overrides: `Actor(dd, act; param=value)`
- Transport models are configurable: `act.ActorCoreTransport.model = :FluxMatcher`
- Algorithm selection: `act.ActorFluxMatcher.algorithm = :anderson`

### Integration Notes
- Some actors depend on others (e.g., sizing actors require subsequent `ActorCXbuild`)
- Passive structures must exist before running vertical stability
- PF design should be run after build changes
- FUSE interfaces with OMFIT/OMAS and the broader IMAS ecosystem

### IMAS Data Best Practices
- **Time coordinate precedence**: Always set time coordinates before data values
- **Array resizing**: Use `resize!()` before working with arrays of structures
- **Global time usage**: Use `dd.global_time` and `[]` indexing for current time operations
- **Expression handling**: Be aware of dynamic expressions vs. frozen data
- **Search and exploration**: Use `findall()` with regex for finding fields efficiently
- **Data comparison**: Use `IMAS.diff()` to identify changes between datasets

## Time-Dependent Simulations

### Overview
FUSE supports time-dependent plasma simulations, allowing users to model plasma evolution through different phases like rampup, flattop, and rampdown. Time-dependent parameters are defined as functions of time in the `ini` structure.

### Experimental Data-Based Time-Dependent Simulations

FUSE can initialize time-dependent simulations directly from experimental data, allowing for post-diction studies that replay experimental conditions.

#### D3D Shot Analysis Example
```julia
using Plots
using FUSE
using Interact
FUSE.ProgressMeter.ijulia_behavior(:clear)

# Load D3D shot data with automatic profile fitting
use_local_cache = false
shot = 168830
@time ini, act = FUSE.case_parameters(:D3D, shot; fit_profiles=true, use_local_cache)

# Initialize from experimental data - this populates pulse_schedule for post-diction
ini.time.simulation_start = ini.general.dd.equilibrium.time_slice[2].time
dd = IMAS.dd()
@time FUSE.init!(dd, ini, act)
```

The key insight is that `FUSE.init!()` with experimental data automatically populates `dd.pulse_schedule` to allow post-diction (reproducing experimental conditions). For "what if" scenarios, users would modify the pulse schedule after initialization.

### Time-Dependent Parameter Configuration

#### Pulse Schedule Time Basis
```julia
# Define time grid for simulation
ini.time.pulse_shedule_time_basis = range(0, 300; step=1.0)  # 0-300s with 1s steps
ini.time.simulation_start = 50.0  # Start simulation at t=50s
```

#### Time-Dependent Functions
Parameters can be defined as functions of time `t` instead of scalar values:

```julia
# Define time-dependent plasma current
rampup_ends = 10.0
ini.equilibrium.ip = t -> ramp(t / rampup_ends) * 14E6 + ramp((t - 100) / 100) * 1E6

# Time-dependent core pressure
ini.equilibrium.pressure_core = t -> ramp(t / rampup_ends) .^ 2 * 0.643e6

# Time-dependent auxiliary heating power
ini.nb_unit[1].power_launched = t -> (1 .+ ramp((t - 100) / 100)) * 16.7e6
ini.ec_launcher[1].power_launched = t -> (1 .+ ramp((t - 100) / 100)) * 10E6
ini.ic_antenna[1].power_launched = t -> (1 .+ ramp((t - 100) / 100)) * 12E6
ini.lh_antenna[1].power_launched = t -> (1 .+ ramp((t - 100) / 100)) * 5E6

# Time-dependent pellet injection frequency
ini.pellet_launcher[1].frequency = t -> (1 .+ ramp((t - 100) / 100)) * 0.01  # Hz
```

#### Pulse Shaping Functions
FUSE provides convenient pulse shaping functions for common time dependencies:

- **`ramp(t)`**: Linear ramp from 0 to 1
- **`step(t)`**: Step function
- **`pulse(t)`**: Pulse shape
- **`trap(t, ramp_fraction)`**: Trapezoidal pulse
- **`gaus(t, order)`**: Gaussian pulse
- **`beta(t, mode)`**: Beta function pulse

#### Sequence Function
For complex time profiles, use `sequence(t, t_y_tuple_sequence)`:
```julia
# Define plasma current with specific time points and values
ini.equilibrium.ip = t -> sequence(t, [(0.0, 0.0), (10.0, 13.0E6), (100.0, 13.0E6), (200.0, 15.0E6)])
```

### Rampup Configuration
```julia
ini.rampup.side = :lfs                      # Low field side rampup
ini.rampup.ends_at = 10.0                   # Rampup ends at t=10s
ini.rampup.diverted_at = rampup_ends * 0.8  # Divertion at 80% of rampup time
```

### Time-Dependent Workflow

#### 1. Initialize with Time-Dependent Parameters
```julia
# Get time-dependent case parameters
ini, _ = FUSE.case_parameters(:ITER; init_from=:scalars, time_dependent=true)

# Configure time grid and start time
ini.time.pulse_shedule_time_basis = range(0, 300; step=1.0)
ini.time.simulation_start = 50.0

# Define time-dependent parameters (see examples above)
# ...

# Initialize hardware from ODS if needed
ini_hw, act = FUSE.case_parameters(:ITER; init_from=:ods)
dd = IMAS.dd()
FUSE.init(dd, ini_hw, act)  # Initialize hardware

# Initialize pulse schedule
FUSE.init(dd, ini, act; initialize_hardware=false)
```

#### 2. Achieve Self-Consistent Initial State
```julia
# Configure stationary plasma solver
act.ActorStationaryPlasma.convergence_error = 2E-2
act.ActorStationaryPlasma.max_iterations = 5
act.ActorSteadyStateCurrent.current_relaxation_radius = 0.2
act.ActorFluxMatcher.verbose = true
act.ActorFluxMatcher.relax = 0.5

# Find self-consistent initial state
FUSE.ActorStationaryPlasma(dd, act; verbose=true)
```

#### 3. L-H Transition Analysis and Configuration

For experimental post-diction, FUSE can analyze L-H transitions from experimental data:

```julia
# Analyze L-H transitions from experimental data
experiment_LH = FUSE.LH_analysis(dd; do_plot=true)

# Configure dynamic pedestal model based on experimental analysis
act.ActorPedestal.model = :dynamic
act.ActorPedestal.tau_n = experiment_LH.tau_n              # Density decay time
act.ActorPedestal.tau_t = experiment_LH.tau_t              # Temperature decay time
act.ActorWPED.ped_to_core_fraction = experiment_LH.W_ped_to_core_fraction
act.ActorEPED.ped_factor = 1.0
act.ActorPedestal.T_ratio_pedestal = 1.0                   # Ti/Te in pedestal

# Option 1: Use experimental density and Zeff profiles as-is
act.ActorPedestal.density_ratio_L_over_H = 1.0
act.ActorPedestal.zeff_ratio_L_over_H = 1.0

# Option 2: Allow density to transition at different times
# act.ActorPedestal.density_ratio_L_over_H = experiment_LH.ne_L_over_H
# act.ActorPedestal.zeff_ratio_L_over_H = experiment_LH.zeff_L_over_H
# dd.pulse_schedule.density_control.n_e_line.reference = experiment_LH.ne_H
# dd.pulse_schedule.density_control.zeff_pedestal.reference = experiment_LH.zeff_H

# L-H transition timing options
# Option 1: From L-H scaling law
# act.ActorPedestal.mode_transitions = missing

# Option 2: User-defined transition times (typical for post-diction)
act.ActorPedestal.mode_transitions = experiment_LH.mode_transitions
act.ActorPedestal.mode_transitions[5.2] = :L_mode          # Force L-mode at t=5.2s
```

#### 4. Advanced Actor Configuration for Time-Dependent Simulations

```julia
# Equilibrium solver configuration
act.ActorEquilibrium.model = :FRESCO    # Options: :EGGO, :FRESCO
act.ActorFRESCO.nR = 65                 # Radial grid points
act.ActorFRESCO.nZ = 65                 # Vertical grid points

# Neutral fueling configuration
act.ActorNeutralFueling.τp_over_τe = 0.25  # Particle confinement time ratio

# Flux matcher for time-dependent transport
act.ActorFluxMatcher.evolve_plasma_sources = false
act.ActorFluxMatcher.algorithm = :simple
act.ActorFluxMatcher.max_iterations = -10          # Negative to suppress warnings
act.ActorFluxMatcher.evolve_pedestal = false
act.ActorFluxMatcher.relax = 0.5

# TGLF configuration for time-dependent transport
act.ActorTGLF.tglfnn_model = "sat1_em_d3d"

# Time step configuration
δt = 0.05                              # Time step size (s)
dd.global_time = ini.general.dd.equilibrium.time_slice[2].time
final_time = ini.general.dd.equilibrium.time[end]
act.ActorDynamicPlasma.Nt = Int(ceil((final_time - dd.global_time) / δt))
act.ActorDynamicPlasma.Δt = final_time - dd.global_time
```

#### 5. Run Time-Dependent Simulation
```julia
# Configure dynamic plasma evolution
act.ActorDynamicPlasma.evolve_current = true          # Evolve current profile
act.ActorDynamicPlasma.evolve_equilibrium = true      # Evolve equilibrium
act.ActorDynamicPlasma.evolve_transport = true        # Evolve transport
act.ActorDynamicPlasma.evolve_hcd = true              # Evolve heating/current drive
act.ActorDynamicPlasma.evolve_pf_active = false       # Keep PF coils fixed
act.ActorDynamicPlasma.evolve_pedestal = true         # Evolve pedestal

# Optional replay modes (for debugging or comparison)
# act.ActorCurrent.model = :replay
# act.ActorEquilibrium.model = :replay  
# act.ActorCoreTransport.model = :replay
# act.ActorPedestal.model = :replay
# act.ActorHCD.ec_model = :replay
# act.ActorHCD.ic_model = :replay
# act.ActorHCD.lh_model = :replay
# act.ActorHCD.nb_model = :replay
# act.ActorHCD.pellet_model = :replay
# act.ActorHCD.neutral_model = :none

# Execute time-dependent simulation
@time actor = FUSE.ActorDynamicPlasma(dd, act; verbose=true)
```

#### 4. Continue Simulation
```julia
# Simulation can be continued for additional steps
# NOTE: actor has its own internal copy of dd
dd = actor.dd  # Get updated data from actor
FUSE.finalize(FUSE.step(actor; n_steps=10))  # Continue for 10 more steps
```

### Visualization of Time-Dependent Results

#### Interactive Time Slider with Profile Visualization
```julia
using Interact

# Simple profile time slider with experimental data overlay
@manipulate for time0 in dd.equilibrium.time
    plot(dd.core_profiles.profiles_1d[time0]; thomson_scattering=true, charge_exchange=true)
end
```

#### Plasma Overview Visualization
```julia
# Interactive plasma overview with comparison to initial state
@manipulate for time0 in slider(dd.equilibrium.time, value=dd.global_time/2.0, label="time")
    try
        # Plot comprehensive plasma overview with reference to initial state
        FUSE.plot_plasma_overview(dd, Float64(time0); dd1, aggregate_hcd=true)
    catch e
        plot()  # Fallback if plotting fails
    end
end
```

#### Animated GIF Generation
```julia
using Printf

# Create animation with detailed plasma overview plots
a = @animate for (k, time0) in enumerate(dd.equilibrium.time[2:2:end])  # Skip every other time point
    try
        FUSE.plot_plasma_overview(dd, Float64(time0); dd1, aggregate_hcd=true)
        
        # Set consistent y-axis limits for better animation
        IMAS.ylim(Dict{Int,Float64}(
            -2 => 0.0,                          # Lower bound for temperature
            -3 => 0.0, 3 => 2.0,               # Pressure limits
            -4 => 0.0, 4 => 4.0,               # Current density limits  
            -5 => 0.0, 5 => 4.0,               # Safety factor limits
            6 => 1.1E20,                        # Density upper bound
            
            -8 => -0.2, 8 => 2.0,              # Transport coefficients
            -9 => -0.25, 9 => .5,              # Heat flux
            -10 => -0.25, 10 => .5,            # Particle flux
            -11 => -1E20, 11 => 1.0E20,        # Sources/sinks
            
            -14 => 0.0, 14 => 0.101,           # Radial bounds
            -15 => 0.0, 15 => 0.101,           # Normalized radial bounds
            -16 => -2.0E19, 16 => 2.0E19))     # Power density bounds
        
        # Save individual frames for GIF creation
        savefig("D3D_$(shot)/D3D_$(shot)___$(@sprintf("%04d", k)).png")
        
    catch e
        plot()  # Fallback plot if errors occur
    end
end

# Create GIF from animation
gif(a, "D3D_$(shot)/D3D_$(shot).gif", fps=12)
```

#### Terminal Command for GIF Creation
The notebook also shows how to create optimized GIFs using ImageMagick:
```bash
# Use ImageMagick to create optimized GIF from PNG frames
magick -delay 2 -loop 0 D3D_168830___*.png -layers Optimize D3D_168830.gif
```

#### Plot Pulse Schedule
```julia
# Plot the configured pulse schedule
plot(ini)

# Interactive pulse schedule visualization
using Interact
@manipulate for time0 in ini.time.pulse_shedule_time_basis
    plot(dd.pulse_schedule; time0, ini.time.simulation_start)
end
```

### Saving and Extracting Results

#### Save Time-Dependent Simulation
```julia
IMAS.imas2json(dd, "ITER_time_dep.json"; freeze=true, strict=true)
```

#### Extract Specific Time Slice
```julia
# Extract data at current global time
dd0 = IMAS.get_timeslice(dd, dd.global_time)
FUSE.plot_plasma_overview(dd0)

# Extract summary of 0D quantities
IMAS.extract(dd0)
```

### Key Actors for Time-Dependent Simulations

- **ActorDynamicPlasma**: Main actor for time-dependent plasma evolution
- **ActorStationaryPlasma**: Establishes self-consistent initial conditions
- **ActorSteadyStateCurrent**: Handles current profile relaxation during initialization

### Time-Dependent Simulation Best Practices

1. **Start with hardware initialization** from ODS data when available
2. **Define realistic time dependencies** using pulse shaping functions
3. **Establish self-consistent initial state** with ActorStationaryPlasma before evolution
4. **Configure evolution flags** based on physics of interest
5. **Use appropriate time step sizes** for numerical stability
6. **Monitor simulation progress** with verbose output
7. **Save checkpoints** at key phases for restart capability
8. **Extract specific time slices** for detailed analysis

### External Control Integration

FUSE supports co-simulation with external controllers through the FuseExchangeProtocol:
```julia
# Connect to external control system (commented in example)
# IMAS.fxp_connect(dd)
```

## Studies and Distributed Computing

### FUSE Studies Framework

FUSE includes a studies framework for managing complex distributed computing workflows. Studies are similar to actors but designed for parameter sweeps, database generation, and large-scale computations.

#### Study Types
- **StudyTGLFdb**: Generate TGLF model databases with parameter sweeps
- **StudyDatabaseGenerator**: Generate databases with parameter distributions for uncertainty quantification and sensitivity analysis

#### Study Parameters Structure (`sty`)

Studies use a `sty` (study parameters) structure similar to `act` (actor parameters) but with study-specific configuration:

```julia
# Get study parameters - similar to case_parameters but for studies
sty, act = FUSE.study_parameters(:TGLFdb)

# Common study parameters
sty.server = "localhost"           # Distributed computing server  
sty.n_workers = 1                  # Number of parallel workers
sty.database_folder = "/path/to/db" # Where to store database
sty.save_folder = "/path/to/save"   # Results save location
sty.file_save_mode = :safe_write    # File writing mode
sty.release_workers_after_run = true # Clean up workers after completion

# Study-specific parameters for TGLFdb
sty.sat_rules = missing                        # Use custom models, or [:sat1, :sat2, :sat3]
sty.custom_tglf_models = ["sat3_em_d3d_azf-1_withnegD"]  # Custom TGLF models to use
```

#### Study Workflow

```julia
using Distributed

# 1. Setup study parameters
sty, act = FUSE.study_parameters(:TGLFdb)
sty.server = "localhost"
sty.n_workers = 4
sty.database_folder = "/path/to/database"
sty.save_folder = "/path/to/results"
sty.custom_tglf_models = ["sat3_em_d3d_azf-1_withnegD"]

# 2. Create and setup study
study = FUSE.StudyTGLFdb(sty, act)  # Automatically calls FUSE.setup()
# OR manually setup: study = FUSE.setup(study)

# 3. Import FUSE on all workers (required for distributed computing)
@everywhere import FUSE

# 4. Configure actor parameters for the study
study.act.ActorFluxMatcher.evolve_rotation = :flux_match
study.act.ActorFluxMatcher.rho_transport = 0.1:0.05:0.75
study.act.ActorFluxMatcher.max_iterations = 300
study.act.ActorFluxMatcher.algorithm = :simple_dfsane
study.act.ActorFluxMatcher.step_size = 1.0

# 5. Run the study
FUSE.run(study)  # Executes distributed computation and saves dataframes

# 6. Access results
study.dataframes_dict["sat3_em_d3d_azf-1_withnegD"]  # DataFrame with results

# 7. Analyze study results
FUSE.analyze(study)
```

#### Study Data Management

- **Results storage**: Study results are stored in `study.dataframes_dict`
- **Key-value access**: Results indexed by model name (e.g., custom TGLF model names)
- **DataFrames**: Results stored as DataFrames for easy analysis and plotting
- **Automatic saving**: Study results are automatically saved to specified folders

#### Distributed Computing Notes

1. **Worker management**: Studies automatically manage distributed workers
2. **FUSE import**: Must import FUSE on all workers with `@everywhere import FUSE`
3. **Worker cleanup**: Workers are released after study completion (configurable)
4. **Resource scaling**: Configure `n_workers` based on available computational resources

### Database Generator Study Framework

The Database Generator study provides a systematic approach to generate comprehensive databases of FUSE simulations with parameter variations. This is particularly useful for uncertainty quantification, sensitivity analysis, and machine learning training data.

#### Basic Database Generator Setup

```julia
# Set up case parameters with distributions
ini, act = FUSE.case_parameters(:ITER; init_from=:scalars)

# Configure actor parameters
act.ActorPedestal.density_match = :ne_line

# Set up core profiles with base values and distributions
ini.core_profiles.ne_setting = :greenwald_fraction
ini.core_profiles.ne_value = 0.2 ↔ [0.2, 1.0]  # Uniform distribution

# Parameter variations using ↔ operator
ini.core_profiles.impurity = :Kr ↔ (:Kr, :Ne, :Xe)      # Categorical options
ini.core_profiles.plasma_mode = :H_mode ↔ (:H_mode, :L_mode)  # Categorical options  
ini.core_profiles.zeff = 1.1 ↔ [1.1, 10.]               # Uniform range
```

#### Advanced Parameter Distributions

FUSE supports sophisticated probability distributions via Distributions.jl:

```julia
using FUSE.SimulationParameters.Distributions

# 1. Uniform distribution (default with ↔ operator)
uniform_zeff = 1.1 ↔ [1.1, 10.]

# 2. Truncated Normal distribution
truncated_normal = 5.0 ↔ truncated(Normal(5.0, 1.5), lower=2.0, upper=8.0)

# 3. Mixed distribution (bimodal example)
low_zeff = truncated(Normal(2, 1.0), lower=1.0, upper=Inf)
high_zeff = truncated(Normal(7.0, 1.0), lower=0.5, upper=Inf)
mixed_zeff = 2.0 ↔ MixtureModel([low_zeff, high_zeff], [0.3, 0.7])

# Apply distribution to parameter
ini.core_profiles.zeff = mixed_zeff
```

#### Distribution Sampling and Testing

```julia
# Test parameter sampling before running study
Nsample = 1e4

# Generate samples from distribution
samples = [rand(getfield(ini.core_profiles, :zeff)) for _ in 1:Nsample]

# Visualize distribution
using Plots
histogram(samples; normalize=:pdf, nbins=50, alpha=0.4)
xlabel!("Zeff")
ylabel!("Probability")
```

#### Database Generator Study Configuration

```julia
# Get study parameters for database generation
sty, _ = FUSE.study_parameters(:DatabaseGenerator)

# Configure distributed computing
sty.server = "localhost"
sty.n_workers = 4
sty.file_save_mode = :safe_write
sty.save_folder = mktempdir()  # Or specify permanent directory
sty.n_simulations = 10         # Number of random samples to generate
```

#### Database Generator Workflow

```julia
# 1. Create study with configured parameters
study = FUSE.StudyDatabaseGenerator(sty, ini, act)

# 2. Set up distributed computing
using Distributed
@everywhere import FUSE
@everywhere ProgressMeter = FUSE.ProgressMeter

# 3. Define custom workflow function
@everywhere function workflow_DatabaseGenerator(dd::FUSE.IMAS.dd, ini::FUSE.ParametersAllInits, act::FUSE.ParametersAllActors)
    FUSE.init(dd, ini, act)
    # Add additional physics calculations if needed:
    # FUSE.ActorFluxMatcher(dd, act)
    # FUSE.ActorStationaryPlasma(dd, act)
    return nothing
end

# 4. Assign workflow to study
study.workflow = workflow_DatabaseGenerator

# 5. Execute study
FUSE.run(study)

# 6. Access results
study.dataframe  # Contains all simulation results and parameters
```

#### Analyzing Database Results

```julia
# Plot correlations between parameters and outputs
using Plots
scatter(study.dataframe.zeff_ped, abs.(study.dataframe.Prad_tot), 
        xlabel="Zeff", ylabel="Total radiation [MW]", 
        yscale=:log10, label=nothing)

# Load specific cases for detailed analysis
dd, ini, act = FUSE.load(joinpath(splitpath(study.dataframe.dir[end])[1:end-1]))
FUSE.digest(dd)

# Run additional analysis on individual cases
FUSE.ActorFluxMatcher(dd, act; do_plot=true)
```

#### Parametric Scanning Mode

Instead of random sampling, you can specify exact parameter combinations:

```julia
# Create specific parameter combinations
inis = [deepcopy(ini), deepcopy(ini)]
inis[1].core_profiles.ne_value = 0.8  # First case
inis[2].core_profiles.ne_value = 0.5  # Second case

# Use parametric list instead of random sampling
sty.save_folder = mktempdir()
study = FUSE.StudyDatabaseGenerator(sty, inis, act)  # Pass array of inis
study.workflow = workflow_DatabaseGenerator

FUSE.run(study)
```

#### Study vs Actor Differences

| Aspect | Actors | Studies |
|--------|---------|---------|
| Purpose | Single physics calculation | Parameter sweeps, databases, optimization |
| Parameters | `act` structure | `sty` + `act` structures |
| Execution | Single-threaded | Distributed computing |
| Results | Modifies `dd` in-place | Returns DataFrames |
| Setup | Manual actor calls | Automatic workflow management |
| Parameter Variation | Manual modification | Built-in distribution sampling |

### Multi-Objective Optimization Studies

FUSE provides a sophisticated multi-objective optimization framework through the `StudyMultiObjectiveOptimizer` study type. This enables systematic optimization of fusion reactor designs subject to engineering and physics constraints.

#### Multi-Objective Optimization Setup

```julia
# Setup IMAS objective and constraint function libraries
IMAS.update_ObjectiveFunctionsLibrary!()
IMAS.update_ConstraintFunctionsLibrary!()

# Create local copies for modification
OFL = deepcopy(IMAS.ObjectiveFunctionsLibrary)
CFL = deepcopy(IMAS.ConstraintFunctionsLibrary)

# Define objective functions (what to minimize/maximize)
objective_functions = [OFL[:min_capital_cost]]

# Define constraint functions (engineering/physics limits)
constraint_functions = [
    CFL[:power_electric_net],           # Net electrical power constraint
    CFL[:min_lh_power_threshold_fraction], # L-H threshold constraint
    CFL[:max_tf_coil_j],                # TF coil current density limit
    CFL[:max_oh_coil_j],                # OH coil current density limit  
    CFL[:max_pl_stress],                # Plasma-facing component stress
    CFL[:max_tf_coil_stress],           # TF coil stress limit
    CFL[:max_oh_coil_stress]            # OH coil stress limit
]
```

#### Optimization-Specific Actor Configuration

```julia
# Get base case parameters
ini, act = FUSE.case_parameters(:KDEMO)

# Configure actors for optimization
act.ActorEquilibrium.model = :TEQUILA
act.ActorCoreTransport.model = :FluxMatcher
act.ActorPFdesign.model = :optimal              # Enable optimal PF coil design
act.ActorFluxSwing.operate_oh_at_j_crit = true  # Maximize flattop duration
act.ActorWholeFacility.update_plasma = true     # Enable plasma updates during optimization

# Disable error throwing for optimization (handled by constraints)
act.ActorHFSsizing.error_on_performance = false
act.ActorHFSsizing.error_on_technology = false

# Configure transport and pedestal models
act.ActorPedestal.model = :EPED
act.ActorFluxMatcher.evolve_densities = :flux_match
act.ActorFluxMatcher.max_iterations = 500
act.ActorFluxMatcher.rho_transport = 0.2:0.05:0.8
act.ActorTGLF.tglfnn_model = "sat0quench_em_d3d+mastu_azf+1"
act.ActorFluxMatcher.algorithm = :simple
act.ActorFluxMatcher.evolve_pedestal = true
```

#### Parameter Bounds Definition

Define optimization variables and their bounds using the `↔` operator:

```julia
# Define optimization variables with bounds [min, max]
ini.equilibrium.B0 = ini.equilibrium.B0 ↔ [3.0, 15.0]    # Toroidal field
ini.equilibrium.ip = ini.equilibrium.ip ↔ [4.0e6, 22e6]   # Plasma current  
ini.equilibrium.R0 = ini.equilibrium.R0 ↔ [4.0, 8.0]      # Major radius

# Remove parameters from optimization by setting to missing
ini.equilibrium.pressure_core = missing
ini.requirements.flattop_duration = missing
ini.requirements.log10_flattop_duration = missing

# Set fixed requirement targets
ini.requirements.power_electric_net = 50e6  # Target 50 MWe
ini.requirements.tritium_breeding_ratio = 1.1
ini.requirements.coil_j_margin = 0.1        # 10% current margin
ini.requirements.coil_stress_margin = 0.1   # 10% stress margin
```

#### Study Parameters Configuration

```julia
# Get study parameters for multi-objective optimization
sty, _ = FUSE.study_parameters(:MultiObjectiveOptimizer)

# Configure distributed computing
sty.server = "localhost"  # Can be "saga", "omega", or cluster name
sty.n_workers = 2

# Set up save directory
dirr = joinpath(pwd(), "test_moop")
if !isdir(dirr)
    sty.save_folder = mkdir(dirr)
else
    sty.save_folder = dirr
end

# Configure genetic algorithm parameters
sty.restart_workers_after_n_generations = 5  # Worker cleanup frequency
sty.population_size = 20      # For testing (use 300 for realistic studies)
sty.number_of_generations = 5 # For testing (use 100 for realistic studies)
```

#### Multi-Objective Optimization Workflow

```julia
# 1. Create optimization study
study = FUSE.StudyMultiObjectiveOptimizer(sty, ini, act, constraint_functions, objective_functions)

# 2. Set up distributed computing
using Distributed
@everywhere import FUSE
@everywhere import IJulia  # If running in Jupyter

# 3. Verify workers are available
Distributed.workers()

# 4. Execute optimization study
FUSE.run(study)  # Runs genetic algorithm optimization

# Note: Workers are automatically released after completion
```

#### Available Objective Functions

Common objective functions in the ObjectiveFunctionsLibrary:
- **`:min_capital_cost`**: Minimize total capital cost
- **`:max_power_electric_net`**: Maximize net electrical power output
- **`:min_tritium_breeding_ratio`**: Minimize tritium breeding ratio (or maximize when negative)
- **`:min_major_radius`**: Minimize reactor size
- **`:max_fusion_power`**: Maximize fusion power

#### Available Constraint Functions

Common constraint functions in the ConstraintFunctionsLibrary:
- **Power constraints**: `:power_electric_net`, `:min_lh_power_threshold_fraction`
- **Magnetic coil limits**: `:max_tf_coil_j`, `:max_oh_coil_j`, `:max_cs_coil_j`
- **Stress constraints**: `:max_pl_stress`, `:max_tf_coil_stress`, `:max_oh_coil_stress`
- **Physics limits**: `:max_beta_n`, `:min_q95`, `:max_neutron_wall_loading`
- **Engineering bounds**: `:max_heat_flux_divertor`, `:min_tritium_breeding_ratio`

#### Multi-Objective Optimization Best Practices

1. **Start with small populations** (20-50) for testing, scale to 300+ for production
2. **Use realistic generation counts** (100+ for convergence)
3. **Define appropriate constraints** to ensure feasible designs
4. **Monitor worker performance** and adjust `restart_workers_after_n_generations`
5. **Set up adequate computing resources** - optimization is computationally intensive
6. **Validate constraint functions** match your design requirements
7. **Use bounds wisely** - too narrow limits solution space, too wide causes convergence issues

#### Optimization Results Analysis

Results from multi-objective optimization studies are stored in the study object and can be analyzed for Pareto-optimal solutions, parameter sensitivity, and design trade-offs.

## Tutorial Files Studied

1. **tutorial.ipynb** - Basic FUSE concepts, actor execution, and workflows
   - Core concepts and data structures
   - Individual actor execution examples  
   - Whole facility design workflow
   - Checkpoint system usage
   - Plotting and visualization examples

2. **tutorial_imas.ipynb** - IMAS data structures and manipulation
   - IMAS/IMASdd.jl fundamentals and basic usage
   - Working with arrays of structures and time-dependent data
   - Time series manipulation with `@ddtime` macro and global time
   - Dynamic expressions and data freezing
   - Advanced plotting with Plots.jl recipes
   - Physics functions and field line tracing examples
   - Data exploration, search, and comparison tools

3. **fluxmatcher.ipynb** - Transport modeling and flux matching
   - FUSE flux matcher configuration and usage
   - Transport model options (TGLFNN, TJLF, TGLF, QLGYRO, GKNN)
   - Flux matching algorithms and convergence strategies
   - Working with D3D L-mode and H-mode cases
   - Using custom ODS/IMAS experimental data
   - Neural network model selection and configuration
   - Advanced techniques: preconditioning and algorithm selection
   - Heating and current drive model configuration

4. **study_TGLFdb.ipynb** - FUSE studies framework and distributed computing
   - Study parameters structure (`sty`) vs actor parameters (`act`)
   - StudyTGLFdb for TGLF model database generation
   - Distributed computing setup with `@everywhere import FUSE`
   - Study workflow: setup, configuration, execution, and analysis
   - Results management with DataFrames and automatic saving
   - Custom TGLF models and saturation rules configuration
   - Worker management and resource scaling

5. **time_dependent_iter.ipynb** - Time-dependent plasma simulations
   - Time-dependent parameter configuration using functions of time
   - Pulse shaping functions: ramp(), step(), pulse(), trap(), gaus(), beta()
   - Pulse schedule time basis and simulation start time configuration
   - Rampup configuration and divertion timing
   - ActorDynamicPlasma for time evolution with configurable physics flags
   - ActorStationaryPlasma for self-consistent initial conditions
   - Visualization of time-dependent results with interactive sliders and animations
   - Simulation continuation and step management
   - External control integration with FuseExchangeProtocol
   - Time slice extraction and result saving strategies

6. **study_database_generator.ipynb** - Database generation with parameter distributions
   - StudyDatabaseGenerator for systematic parameter variation studies
   - Parameter distribution specification using ↔ operator with ranges and categorical options
   - Advanced probability distributions via Distributions.jl (uniform, truncated normal, mixed)
   - Distribution sampling and testing before full study execution
   - Custom workflow functions for database generation studies  
   - Distributed computing setup with @everywhere import FUSE
   - Two modes: random sampling from distributions vs parametric scanning with specific parameter lists
   - Database result analysis with correlation plots and individual case loading
   - Integration with FUSE physics actors (ActorFluxMatcher, etc.) in workflow functions

7. **study_multi_objective_optimizer.ipynb** - Multi-objective optimization framework
   - StudyMultiObjectiveOptimizer for systematic reactor design optimization
   - IMAS ObjectiveFunctionsLibrary and ConstraintFunctionsLibrary integration
   - Parameter bounds specification using ↔ operator for optimization variables
   - Optimization-specific actor configuration (ActorPFdesign.model=:optimal, etc.)
   - Genetic algorithm parameters: population size, generations, worker management
   - Objective functions: capital cost minimization, power maximization, etc.
   - Constraint functions: coil current/stress limits, power requirements, physics bounds
   - Distributed computing setup for computationally intensive optimization
   - Integration with whole facility design workflow (ActorWholeFacility.update_plasma=true)

8. **time_dependent_d3d.ipynb** - Experimental post-diction time-dependent simulations
   - Loading D3D shot data with automatic profile fitting (fit_profiles=true)
   - Post-diction workflow: experimental data automatically populates pulse_schedule
   - L-H transition analysis using FUSE.LH_analysis() for experimental parameter extraction
   - Dynamic pedestal modeling with ActorPedestal.model=:dynamic configuration
   - Advanced time-dependent actor configuration (FRESCO equilibrium, neutral fueling)
   - Replay mode options for debugging and component-by-component analysis
   - Comprehensive visualization with FUSE.plot_plasma_overview() and experimental overlay
   - Animation creation with consistent axis limits for smooth time evolution
   - Profile visualization with experimental data overlay (thomson_scattering, charge_exchange)
   - ImageMagick integration for optimized GIF creation from PNG frames