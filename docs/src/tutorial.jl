# # FUSE Introductory Tutorial

# Download this tutorial from the [FuseExamples repository](https://github.com/ProjectTorreyPines/FuseExamples/blob/master/tutorial.ipynb)

# Import the necessary packages

using Plots # for plotting
using FUSE # this will also import IMAS in the current namespace

# ## Starting from a use-case
# FUSE comes with some predefined [use-cases](https://fuse.help/stable/cases.html), some of which are used for regression testing.
# Note that some use cases are for non-nuclear experiments and certain Actors like Blankets or BalanceOfPlant will not perform any actions.

FUSE.test_cases

# Get initial parameters (`ini`) and actions (`act`) for a given use-case

ini, act = FUSE.case_parameters(:KDEMO);

# Modifying [`ini` parameters](https://fuse.help/stable/ini.html).

ini.equilibrium.B0 = 7.8
ini.equilibrium.R0 = 6.5;

# Modifying [`act` parameters](https://fuse.help/stable/act.html).

act.ActorCoreTransport.model = :FluxMatcher;

# Initialize the data dictionary (`dd`) using the 0D parameters

dd = FUSE.init(ini, act);

# Using checkpoints to save and restore states (we'll use this later)

chk = FUSE.Checkpoint()
@checkin chk :init dd ini act

# ## Exploring the data dictionary
# * FUSE stores data following the IMAS data schema.
# * The root of the data structure is `dd`, which stands for "Data Dictionary".
# * More details are available in the [documentation](https://fuse.help/stable/dd.html).

# Display part of the equilibrium data in `dd`

dd.equilibrium.time_slice[2].boundary

# this can be done up to a certain depth with `print_tree`

print_tree(dd.equilibrium.time_slice[2].boundary; maxdepth=1)

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

# Each `dd` has a `global_time` attribute, and actors operate at such time

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

dd.equilibrium.vacuum_toroidal_field.b0

# ## Expressions in `dd`

# Some fields in the data dictionary are expressions (ie. Functions).
# For example `dd.core_profiles.profiles_1d[].pressure` is dynamically calculated as the product of thermal densities and temperature with addition of fast ions contributions

print_tree(dd.core_profiles.profiles_1d[1]; maxdepth=1)

# accessing a dynamic expression, automatically evaluates it (in the `pressure` example, we get an array with data)

dd.core_profiles.profiles_1d[1].electrons.pressure

# In addition to evaluating expressions by accessing them, expressions in the tree can be evaluated using `IMAS.freeze()`

print_tree(IMAS.freeze(dd.core_profiles.profiles_1d[1]); maxdepth=1)

# ## Whole facility design

# Here we restore the `:init` checkpoint that we had previously stored. Resetting any changes to `dd`, `ini`, and `act` that we did in the meantime.

@checkout chk :init dd ini act

# Actors in FUSE can be executed by passing two arguments to them: `dd` and `act`.
# Internally, actors can call other actors, creating workflows.
# For example, the `ActorWholeFacility` can be used to to get a self-consistent stationary whole facility design.
# The `actors:` print statements with their nested output tell us what actors are calling other actors.

FUSE.ActorWholeFacility(dd, act);

# Like before we can checkpoint results for later use

@checkin chk :awf dd ini act

# ## Running a custom workflow

# Let's now run a series of actors similar to what `ActorWholeFacility` does
# and play around with plotting to get a sense of what each individual actor does.

# Let's start again from after the initialization stage

@checkout chk :init dd ini act

# Let's start by positioning the PF coils, so that we stand a chance to reproduce the desired plasma shape.
# This will be important to ensure the stability of the `ActorStationaryPlasma` that we are going to run next.

actor = FUSE.ActorPFdesign(dd, act);

# The `ActorStationaryPlasma` iterates between plasma transport, pedestal, equilibrium and sources to return a self-consistent plasma solution

peq = plot(dd.equilibrium; label="before")
pcp = plot(dd.core_profiles; color=:gray, label="before")
FUSE.ActorStationaryPlasma(dd, act);

# we can compare equilibrium before and after the self-consistency loop

plot!(peq, dd.equilibrium; label="after")

# we can compare core_profiles before and after the self-consistency loop

plot!(pcp, dd.core_profiles; label="after")

# here are the sources

plot(dd.core_sources)

# and the flux-matched transport

plot(dd.core_transport)

# HFS sizing actor changes the thickness of the OH and TF layers on the high field side to satisfy current and stresses constraints

plot(dd.build)
FUSE.ActorHFSsizing(dd, act);
plot!(dd.build; cx=false)

# The stresses on the center stack are stored in the `solid_mechanics` IDS

plot(dd.solid_mechanics.center_stack.stress)

# LFS sizing actors change location of the outer TF leg to meet ripple requirements

plot(dd.build)
FUSE.ActorLFSsizing(dd, act);
plot!(dd.build; cx=false)

# A custom `show()` method is defined to print the summary of `dd.build.layer`

dd.build.layer

# ActorHFSsizing and ActorLFSsizing only change the layer's thicknesses, so we then need to trigger a build of the 2D cross-sections after them:

FUSE.ActorCXbuild(dd, act);
plot(dd.build)

# Generate passive structures information (for now the vacuum vessel)

FUSE.ActorPassiveStructures(dd, act)
plot(dd.pf_passive)

# We can now give the PF coils their final position given the new build

actor = FUSE.ActorPFdesign(dd, act);
plot(actor)

# With information about both pf_active and pf_passive we can now evaluate vertical stability

ActorVerticalStability(dd, act)
IMAS.freeze(dd.mhd_linear)

# The `ActorNeutronics` calculates the heat flux on the first wall

FUSE.ActorNeutronics(dd, act);
p = plot(; layout=2, size=(900, 350))
plot!(p, dd.neutronics.time_slice[].wall_loading, subplot=1)
plot!(p, FUSE.define_neutrons(dd, 100000)[1], dd.equilibrium.time_slice[]; subplot=1, colorbar_entry=false)
plot!(p, dd.neutronics.time_slice[].wall_loading; cx=false, subplot=2, ylabel="")

# The `ActorBlanket` will change the thickess of the first wall, breeder, shield, and Li6 enrichment to achieve target TBR

FUSE.ActorBlanket(dd, act);
print_tree(IMAS.freeze(dd.blanket); maxdepth=5)

# The `ActorDivertors` actor calculates the divertors heat flux

FUSE.ActorDivertors(dd, act);
print_tree(IMAS.freeze(dd.divertors); maxdepth=4)

# The `ActorBalanceOfPlant` calculates the optimal cooling flow rates for the heat sources (breeder, divertor, and wall) and get an efficiency for the electricity conversion cycle

actor = FUSE.ActorBalanceOfPlant(dd, act);
plot(actor)

# `ActorCosting` will break down the capital and operational costs

FUSE.ActorCosting(dd, act)
plot(dd.costing)

# Let's checkpoint our results

@checkin chk :manual dd ini act

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
FUSE.digest(dd)#, filename)


