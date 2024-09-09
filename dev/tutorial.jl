# # FUSE Introductory Tutorial

# Import the necessary packages
using Plots
using FUSE

# ## Starting from a use-case
# Here, we start by picking from one of the existing test cases.
# Note: Some test cases are for non-nuclear experiments.
# In those cases, certain Actors like Blankets or BalanceOfPlant will not perform any actions.
FUSE.test_cases

# Get initial parameters (`ini`) and actions (`act`) for a given use-case
ini, act = FUSE.case_parameters(:KDEMO);

# Modifying `ini` parameters
ini.equilibrium.B0 = 7.8
ini.equilibrium.R0 = 6.5
ini.equilibrium

# Modifying `act` parameters
act.ActorCoreTransport.model = :FluxMatcher

# Initialize the data dictionary (`dd`) using the 0D parameters
dd = FUSE.init(ini, act);

## Using checkpoints to save and restore states (we'll use this later)
chk = FUSE.Checkpoint()
chk[:init] = dd, ini, act;

# ## Exploring the data dictionary
# * FUSE stores data following the IMAS data schema.
# * The root of the data structure is `dd`, which stands for "Data Dictionary".
# * More details are available in the [documentation](https://fuse.help/stable/dd.html).

# Display part of the equilibrium data in `dd`
dd.equilibrium.time_slice[2].boundary

# Understanding the IMAS structure
# `dd` and its fields, such as `dd.equilibrium`, are instances of Julia structs defined in the `IMAS.jl` package.
IMAS.equilibrium__time_slice___boundary

# ## Plotting data from `dd`
# FUSE provides Plots.jl recipes for visualizing data from `dd`, this means different plots are shown by calling the same `plot()` function on different items in the data structure.
# Learn more about Plots.jl [here](https://docs.juliaplots.org)

# For example plotting the equilibrium...
plot(dd.equilibrium)

# ...or the core profiles
plot(dd.core_profiles)

# These plots can be composed by calling `plot!()` instead of `plot()`
plot(dd.equilibrium; color=:gray, cx=true)
plot!(dd.build; equilibrium=false, pf_active=false)
plot!(dd.pf_active)

# Plotting an array...
plot(dd.core_profiles.profiles_1d[1].pressure_thermal)

# ...is different from plotting a field from the IDS
plot(dd.core_profiles.profiles_1d[1], :pressure_thermal)

# Customizing plot attributes:
plot(dd.core_profiles.profiles_1d[1], :pressure_thermal; label="", linewidth=2, color=:red, labelfontsize=25)

# ## Working with time series

# The IMAS data structure supports time-dependent data, and IMAS.jl provides ways to handle time data efficiently.

# Each `dd` has a global time defined. Most actors operate at the time defined by the `dd.global_time`
dd.global_time

# Here we see that equilibrium has mulitiple time_slices
dd.equilibrium.time

# Accessing time-dependent arrays of structures, via integer index
eqt = dd.equilibrium.time_slice[2]
eqt.time

# At a given time, by passing the time as a floating point number (in seconds)
eqt = dd.equilibrium.time_slice[0.0]
eqt.time

# At the global time, leaving the square brackets empty
eqt = dd.equilibrium.time_slice[]
eqt.time

# Using the `@ddtime` macro to access and modify time-dependent arrays at `dd.global_time`:
dd.equilibrium.vacuum_toroidal_field.b0

# Accessing data at `dd.global_time`
my_b0 = @ddtime(dd.equilibrium.vacuum_toroidal_field.b0)

# Writin data at `dd.global_time`
@ddtime(dd.equilibrium.vacuum_toroidal_field.b0 = my_b0 + 1)

#
dd.equilibrium.vacuum_toroidal_field.b0

# ## Expressions and consistency in IMAS data

# Some fields in the data dictionary are expressions (ie. Functions), that are dynamically evaluated upon access
print_tree(dd.core_profiles.profiles_1d[1]; maxdepth=1)

# for example, pressure is dynamically calculated as the product of thermal densities and temperature with addition of fast ions contributions
dd.core_profiles.profiles_1d[1].electrons.pressure

# Expressions in the tree can be evaluated using `IMAS.freeze()`
print_tree(IMAS.freeze(dd.core_profiles.profiles_1d[1]); maxdepth=1)

# ## Whole facility design

# Restore init checkpoint
dd, ini, act = chk[:init];

# Run ActorWholeFacility to get a self-consistent stationary whole facility design
FUSE.ActorWholeFacility(dd, act);

# Checkpoint results
chk[:awf] = dd, ini, act;

# ## Running a custom workflow

# Let's runs a series of actors similar to what ActorWholeFacility does

dd, ini, act = chk[:init];
plot(dd.equilibrium; label="before")
FUSE.ActorEquilibrium(dd, act; ip_from=:equilibrium);
plot!(dd.equilibrium; label="after")

# The stationary plasma actor iterates between plasma transport, pedestal, equilibrium and sources to return a self-consistent plasma solution
plot(dd.core_profiles; color=:gray)
FUSE.ActorStationaryPlasma(dd, act);
plot!(dd.core_profiles)

#
plot(dd.core_sources)

#
plot(dd.core_transport)

# HFS sizing actor changes the thickness of the OH and TF layers on the high field side to satisfy current and stresses constraints
plot(dd.build)
FUSE.ActorHFSsizing(dd, act);
plot!(dd.build; cx=false)

#
plot(dd.solid_mechanics.center_stack.stress)

# LFS sizing actors change location of the outer TF leg to meet ripple requirements
plot(dd.build)
FUSE.ActorLFSsizing(dd, act);
plot!(dd.build; cx=false)

# ActorHFSsizing and ActorLFSsizing change the layer's thicknesses
dd.build.layer

# We then need to trigger a build of the 2D cross-sections
FUSE.ActorCXbuild(dd, act);
plot(dd.build)

# We can then position the PF coils to best match the desired plasma shape
actor = FUSE.ActorPFdesign(dd, act; update_equilibrium=true);
plot(actor)

# We need to update the wall shape, since the equilibrium may have changed based on the new PF coils locations
FUSE.ActorCXbuild(dd, act; rebuild_wall=true);
plot!(dd.build)

# The Neutronics actor calculates the heat flux on the first wall
FUSE.ActorNeutronics(dd, act);
p = plot(; layout=2, size=(900, 350))
plot!(p, dd.neutronics.time_slice[].wall_loading, subplot=1)
plot!(p, FUSE.define_neutrons(dd, 100000)[1], dd.equilibrium.time_slice[]; subplot=1, colorbar_entry=false)
plot!(p, dd.neutronics.time_slice[].wall_loading; cx=false, subplot=2, ylabel="")

# The BlanketActor will change the thickess of the first wall, breeder, shield, and Li6 enrichment to achieve target TBR
FUSE.ActorBlanket(dd, act);
print_tree(IMAS.freeze(dd.blanket); maxdepth=5)

# The Divertors actor calculates the divertors heat flux
FUSE.ActorDivertors(dd, act);
print_tree(IMAS.freeze(dd.divertors); maxdepth=4)

# Here we calculate the optimal cooling flowrates for the specific sources of heat and electricity conversion cycle
actor = FUSE.ActorBalanceOfPlant(dd, act);
plot(actor)

# Capital and operational cost breakdown
FUSE.ActorCosting(dd, act)
plot(dd.costing)

# Checkpoint results
chk[:manual] = dd, ini, act;

# ## Saving and loading data
tutorial_temp_dir = tempdir()
filename = joinpath(tutorial_temp_dir, "$(ini.general.casename).json")

# When saving data to be shared outside of FUSE, one can set `freeze=true` so that all expressions in the dd are evaluated and saved to file.
IMAS.imas2json(dd, filename; freeze=false, strict=false);

# Load from JSON
dd1 = IMAS.json2imas(filename);

# ## Comparing two IDSs
# We can introduce a change in the `dd1` and spot it with the `diff` function
dd1.equilibrium.time_slice[1].time = -100.0
IMAS.diff(dd.equilibrium, dd1.equilibrium)

# ## Summary
# Snapshot of `dd` in 0D quantities (evaluated at `dd.global_time`)
FUSE.extract(dd)

# Extract + plots saved to PDF (or printed to screen it `filename` is omitted)
filename = joinpath(tutorial_temp_dir, "$(ini.general.casename).pdf")
FUSE.digest(dd, filename)
