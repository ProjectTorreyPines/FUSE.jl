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
# FUSE provides Plots.jl recipes for visualizing data from `dd`. These can be customized as needed.

# 
plot(dd.equilibrium)

# 
plot(dd.core_profiles)

#
plot(dd.core_sources)

# Composing multiple plots:
plot(dd.equilibrium; color=:gray, cx=true)
plot!(dd.build; equilibrium=false, pf_active=false)
plot!(dd.pf_active)

# Plotting an array...
plot(dd.core_profiles.profiles_1d[1].pressure_thermal)

# ...is different from plotting a field from the IDS. See how
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

# 
my_b0 = @ddtime(dd.equilibrium.vacuum_toroidal_field.b0)

#
@ddtime(dd.equilibrium.vacuum_toroidal_field.b0 = my_b0 + 1)

#
@ddtime(dd.equilibrium.vacuum_toroidal_field.b0)

# ## Expressions and consistency in IMAS data

# Some fields in the IMAS structure are calculated dynamically to ensure consistency.
dd.core_profiles.profiles_1d[1].electrons.pressure

# Evaluate all expressions using `IMAS.freeze()`
IMAS.freeze(dd.core_profiles.profiles_1d[1])

# ## Running ActorWholeFacility to get a self-consistent stationary whole facility design

# Restore init checkpoint
dd, ini, act = chk[:init];

# Run ActorWholeFacility
FUSE.ActorWholeFacility(dd, act);

# Checkpoint results
chk[:awf] = dd, ini, act;

# ## Running a custom workflow

# Let's runs a series of actors similar to what ActorWholeFacility does

dd, ini, act = chk[:init];
FUSE.ActorEquilibrium(dd, act; ip_from=:equilibrium, do_plot=true);

#
FUSE.ActorHFSsizing(dd, act; do_plot=true);

#
FUSE.ActorLFSsizing(dd, act; do_plot=true);

#
FUSE.ActorCXbuild(dd, act; do_plot=true);

#
FUSE.ActorPFactive(dd, act; do_plot=true, update_equilibrium=true);

# Update wall shape after equilibrium update
FUSE.ActorCXbuild(dd, act; rebuild_wall=true, do_plot=true);

#
FUSE.ActorNeutronics(dd, act; do_plot=true);

#
FUSE.ActorBlanket(dd, act);

#
FUSE.ActorDivertors(dd, act);

#
FUSE.ActorBalanceOfPlant(dd, act);

#
FUSE.ActorCosting(dd, act)
chk[:manual] = dd, ini, act;

# ## Saving and loading data to file
tutorial_temp_dir = tempdir()
filename = joinpath(tutorial_temp_dir, "$(ini.general.casename).json")

# When saving data to be shared outside of FUSE, one can set `freeze=true` so that all expressions in the dd are evaluated and saved to file.
IMAS.imas2json(dd, filename; freeze=false, strict=false)

# Load from JSON
dd1 = IMAS.json2imas(filename);

# ## Comparing two IDSs
dd1.equilibrium.time_slice[1].time = -100.0 # Introducing a change
IMAS.diff(dd.equilibrium, dd1.equilibrium)

# ## Summary
# Snapshot of `dd` in 0D quantities (evaluated at `dd.global_time`)
FUSE.extract(dd)

# Extract + plots
FUSE.digest(dd)

# Extract + plots saved to PDF
filename = joinpath(tutorial_temp_dir, "$(ini.general.casename).pdf")
FUSE.digest(dd, filename)
