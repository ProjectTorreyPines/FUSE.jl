# # FUSE Introductory Tutorial

# Download this tutorial from the [FuseExamples repository](https://github.com/ProjectTorreyPines/FuseExamples/blob/master/tutorial.ipynb)

# ## Basic concepts
# 
# To make sense of this tutorial, you'll need to know the following organization concepts of FUSE:
# 
# 1. **📂 Data storage**: All data is stored in the `dd` structure, which follows the ITER IMAS ontology.
# 2. **🧠 Actors**: The core components of FUSE simulations are physics and engineering actors.
# 3. **🕹️ Control**: Actor functionality is governed by `act` parameters.
# 4. **🚀 Initialization**: The data structure can be initialized from 0D `ini` parameters.
# 5. **🔧 Use cases**: FUSE includes templates for various machines (e.g., FPP, ITER, ARC).
# 6. **🔄 Workflows**: Self-contained studies and optimizations are conducted via workflows, typically involving multiple FUSE simulations.
# 7. **🌍 Interoperability**: FUSE interfaces with existing modeling tools like OMFIT/OMAS and the IMAS ecosystem.
# 
# A diagram illustrating these concepts is provided below:
# ![image.png](attachment:5ff4862c-b856-485b-9649-02fec5bdd0d7.png)

# # Let's get started!
# -----------
# 
# **NOTE**: Julia is a Just In Time (JIT) programming language. The first time something is executed it will take longer because of the compilation process. Subsequent calls the the same code will be blazingly fast.
# 
# -----------
# 
# Import the necessary packages

using Plots # for plotting
using FUSE # this will also import IMAS in the current namespace

# ## Starting from a use-case
# 
# FUSE comes with some predefined [use-cases](https://fuse.help/stable/cases.html), some of which are used for regression testing.
# 
# Note that some use cases are for non-nuclear experiments and certain actors like Blankets or BalanceOfPlant will not perform any actions.

# Get initial parameters (`ini`) and actions (`act`) for a given use-case, let's use KDEMO for example

ini, act = FUSE.case_parameters(:KDEMO);

# The `ini` data structure contains 0D parameters that will be used to bootstrap the `dd` with plausible data.
# 
# The [`ini` parameters](https://fuse.help/stable/ini.html) can be modified.

ini.equilibrium.B0 = 7.8
ini.equilibrium.R0 = 6.5;

# The `act` data structure contains parameters that define how the actors (ie the models) will behave.
# 
# The [`act` parameters](https://fuse.help/stable/act.html) can also be modified.

act.ActorCoreTransport.model = :FluxMatcher;

# `ini` and `act` can now be used to initialize the data dictionary (`dd`) using the 0D parameters.
# 
# **NOTE:** `init()` does not return a self-consistent solution, just a plausible starting point to initialize our simulations!

dd = IMAS.dd() # an empty dd
FUSE.init(dd, ini, act);

# Let's see what we got

plot(dd.build)

plot(dd.equilibrium)

plot(dd.core_profiles)

plot(dd.core_sources)

# We can `@checkin` and `@checkout` variables with an associated tag.
# 
# This is handy to save and restore (checkpoint) our progress without having to always start from scratch (we'll use this later).

@checkin :init dd ini act

# ## Running Actors

# Let's now run a series of actors and play around with plotting to get a sense of what each individual actor does.

# Here's how we can restore things back to after the initialization stage (in case we did anything else in between)

@checkout :init dd ini act

# Actors in FUSE can be executed by passing two arguments to them: `dd` and `act`.
# 
# Let's start by positioning the PF coils, so that we stand a chance to reproduce the desired plasma shape.
# This will be important to ensure the stability of the `ActorStationaryPlasma` that we are going to run next.

FUSE.ActorPFdesign(dd, act; do_plot=true); # instead of setting `act.ActorPFdesign.do_plot=true` we can just pass `do_plot=true` as argument without chaning `act`

# The `ActorStationaryPlasma` iterates between plasma transport, pedestal, equilibrium and sources to return a self-consistent plasma solution

peq = plot(dd.equilibrium; label="before")
pcp = plot(dd.core_profiles; color=:gray, label="before")
#act.ActorFluxMatcher.verbose = true
act.ActorFluxMatcher.algorithm = :anderson
#act.ActorFluxMatcher.step_size = 0.1
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
plot(actor) # some actors define their own plot

# With information about both pf_active and pf_passive we can now evaluate vertical stability

FUSE.ActorVerticalStability(dd, act)
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

FUSE.ActorBalanceOfPlant(dd, act);
IMAS.freeze(dd.balance_of_plant)

# `ActorCosting` will break down the capital and operational costs

FUSE.ActorCosting(dd, act)
plot(dd.costing)

# Let's checkpoint our results

@checkin :manual dd ini act

# ## Whole facility design

# Here we restore the `:init` checkpoint that we had previously stored. Resetting any changes to `dd`, `ini`, and `act` that we did in the meantime.

@checkout :init dd ini act

# Actors can call other actors, creating workflows.
# For example, the `ActorWholeFacility` can be used to to get a self-consistent stationary whole facility design.

FUSE.ActorWholeFacility(dd, act);

# Let's check what we got at a glance with the `FUSE.digest(dd)` function:

FUSE.digest(dd)

# Like before we can checkpoint results for later use

@checkin :awf dd ini act

# # Working with the IMAS data structure

# See the IMAS tutorial https://github.com/ProjectTorreyPines/FuseExamples/blob/master/tutorial_imas.ipynb
# 
# Understanding how to work with the IMAS data structure is a must for working within the FUSE ecosystem!

# ## Summary
# Snapshot of `dd` in 0D quantities (evaluated at `dd.global_time`).
# 
# Extract + plots saved to PDF (printed to screen if `filename` is omitted). NOTE: For PDF creation to work, one may need to install of `DejaVu Sans Mono` font.

tutorial_temp_dir = tempdir()
filename = joinpath(tutorial_temp_dir, "$(ini.general.casename).pdf")
display(filename)
FUSE.digest(dd)#, filename)

