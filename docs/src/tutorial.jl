# # FUSE introductory tutorial

# Import Plots and FUSE packages in current namespace
using Plots
using FUSE

# ## Select a use-case

# To start we can pick from one of the exiting test cases
# NOTE: Some of the test cases are for non-nuclear experiments
#       For these some of the Actors will not do anything (eg. Blankets, BalanceOfPlant)
FUSE.test_cases

# Get `ini` and `act` for a given use-case
ini, act = FUSE.case_parameters(:FPP);

# ## Initialize `dd` from 0D parameters

# 
dd = FUSE.init(ini, act);

# ## Run our first Actor
FUSE.ActorEquilibrium(dd, act; ip_from=:equilibrium);

# ## Playing around with `dd`
# * FUSE stores data according to the IMAS data schema
# * The root of the data structure where FUSE stores data is generally referred to as `dd` (which stands for `Data Dictionary`)
# * See the [`IMAS.dd()` documentation](https://fuse.help/stable/dd.html)
# 
# ### Show `dd` content

# We take a look at only one portion of dd
dd.equilibrium.time_slice[2].boundary

# ### Show IMAS stucture
# 
# * `dd`, `dd.equilibrium`, ... are instances of Julia `strut`s that are defined in the `IMASDD.jl` package
#   * `dd = IMAS.dd()`
#   * `dd.equilibrium` is of type `IMAS.equilibrium`
#   * `dd.equilibrium.time_slice` is of type `IMAS.equilibrium__time_slice`
#   * `dd.equilibrium.time_slice[1].boundary` is of type `IMAS.equilibrium__time_slice___boundary`
#   
# * The curious reader can take a look at the [IMASdd/src/dd.jl](https://github.com/ProjectTorreyPines/IMASdd.jl/blob/master/src/dd.jl) file to see those definitions

# Whenever things start with `IMAS.`, then this is a Julia type. This will print the IMAS documentation for that type.
# the convention is that `.` get replaced by a double dash `__` and `[]` get replaced by a triple dash `___`
IMAS.equilibrium__time_slice___boundary

# ### Plot data in `dd`
# 
# * There are [Plots.jl recipies](https://docs.juliaplots.org/latest/recipes/) defined for different IMASdd.jl types of structs
# * These recipies are defined in [IMAS/src/plot.jl](https://github.com/ProjectTorreyPines/IMASdd.jl/blob/master/src/plot.jl)
# * Plots can be [customized](https://docs.juliaplots.org/latest/generated/attributes_series)
# * NOTE: use `display()` is used to force the plot to show when `plot` is not called at the end of the cell

plot(dd.equilibrium)

#

plot(dd.core_profiles)

#

plot(dd.core_sources)

# plots can be composed
plot(dd.equilibrium; color=:gray, cx=true)
plot!(dd.build; equilibrium=false, pf_active=false)
plot!(dd.pf_active)

# plot of an array ...
plot(dd.core_profiles.profiles_1d[1].pressure_thermal)

# ... is different than plotting of a field in an IDS
plot(dd.core_profiles.profiles_1d[1], :pressure_thermal)

# plots can be [customized](https://docs.juliaplots.org/latest/generated/attributes_series):
plot(dd.core_profiles.profiles_1d[1], :pressure_thermal; label="", linewidth=2, color=:red, labelfontsize=25)

# ### Working with time series
# 
# * The IMAS data structure can accomodate time-dependent data
#   
# * Manually handling time in IMAS is tedius and error-prone
#   
# * IMAS.jl provides convenient ways to handle time
#   * `dd.global_time`® sets the "current working time" throughout all of the `dd`
#   * different IDSs can have time arrays of different lengths
#   * data returned based on nearest neighbour

@show dd.global_time

#

@show dd.equilibrium.time;

# To access a time-dependent array of structures at time slice use **Integer** index
eqt = dd.equilibrium.time_slice[2]
@show eqt.time;

# To access a time-dependent array of structures at `dd.globaltime` use []
eqt = dd.equilibrium.time_slice[]
@show eqt.time;

# To access a time-dependent array of structures at a given time use **Float** time
eqt = dd.equilibrium.time_slice[0.0]
@show eqt.time;

# To access a time-dependent array we use @ddtime macro
@show dd.equilibrium.vacuum_toroidal_field.b0; # this is a time dependent data array

# GET data of time-dependent array at `dd.globaltime` (use `IMAS.get_time_array()` to access at other times)
my_b0 = @ddtime(dd.equilibrium.vacuum_toroidal_field.b0)

# SET data of time-dependent array at `dd.globaltime` (use `IMAS.set_time_array()` to access at other times)
my_b0 = round(my_b0; digits=3)
@ddtime(dd.equilibrium.vacuum_toroidal_field.b0 = my_b0)

# ### Expressions
# 
# * Many fields in the IMAS data structure are related to one another. Given some fields others can be calculated. This leads to an issue of consistency, where if some field is updated, all the fields that relies on that data should be updated too.
# 
# * IMAS.jl solves this problem by assigning expressions to certain fields. The value of these expressions are then evaluated on the fly, when the field is requested.
# 
# * Expressions are defined in [IMAS/src/expressions/dynamic.jl](../../../IMAS/src/expressions/dynamic.jl) and [IMAS/src/expressions/onetime.jl](../../../IMAS/src/expressions/onetime.jl)

dd.core_profiles.profiles_1d[1]

# accessing the field evaluates the expression
dd.core_profiles.profiles_1d[1].electrons.pressure

# `IMAS.freeze()` evaluates all expressions
IMAS.freeze(dd.core_profiles.profiles_1d[1])

# ### Saving/loading data
# 
# IMAS.jl can dump data in JSON format as well as HDF5 hierarchical (`imas2hdf` and `hdf2imas`) or tensorized (`imas2h5i` and `h5i2imas`).

# Dump data to JSON
# * `freeze` says whether we want to evaluate all expressions when saving data to file
# * `strict` says whether we want the data to adhere strictly to ITER IMAS dd, or FUSE extensions should also be saved
tutorial_temp_dir = tempdir()
save_name = "FPP_v1_starting_point.json"
@show filename = joinpath(tutorial_temp_dir, save_name)
IMAS.imas2json(dd, filename; freeze=false, strict=false);

# Check for differences between IDSs
dd1 = IMAS.json2imas(filename);
dd1.equilibrium.time_slice[1].time = -100.0 # let's change something
IMAS.diff(dd.equilibrium, dd1.equilibrium)

# ## Playing around with `ini`
# * See the [`ParametersInit()` documentation](https://fuse.help/stable/ini.html)
# * Organized by topical areas
# * Can have multiple layers of nesting
# * Most field defaulted to `missing`

new_ini = FUSE.ParametersInits()

# access a sub-tree
new_ini.equilibrium

# access a leaf
new_ini.equilibrium.Z0

# ### Access detailed descriptions [online](https://fuse.help/stable/ini_details.html#ini.equilibrium.Z0)
# * Notice **units**, **descriptons**
# * Some fields only allow a limited set of **options**
# * Each field stores three things:
#     * **default**: default value when `ini = ParametersInit()` is first called
#     * **base**: value when `FUSE.set_new_base!(ini)` is called
#     * **value**: current value

# Error if trying to access something that is not initialized
try
    new_ini.equilibrium.ip
catch e
    Base.showerror(stderr, e)
end

# Some ini fields only allow a limited set of options
try
    new_ini.tf.technology.material = "something"
catch e
    Base.showerror(stderr, e)
end

# Making changes highlight things in *red* 
new_ini.equilibrium.B0 = 6.12
new_ini.equilibrium.R0 = 5.1
new_ini.equilibrium

# setting a new `base`  (this is done at the end of use-cases definitions), highlights entries in *blue*
FUSE.set_new_base!(new_ini)
new_ini.equilibrium

# Making changes highlight things in *red* 
new_ini.equilibrium.B0 = 6.123
new_ini.equilibrium

# ### Usecases return pre-filled `ini` and `act` parameters
# 
# See for example:
# * [FUSE/src/cases/FPP.jl](../../../FUSE/src/cases/FPP.jl)
# * [FUSE/src/cases/ITER.jl](../../../FUSE/src/cases/ITER.jl)
# * [FUSE/src/cases/ARC.jl](../../../FUSE/src/cases/ARC.jl)

ini, act = FUSE.case_parameters(:FPP)
ini.equilibrium

# ## Playing around with `act`
# * [`FUSE.ParametersActor()` documentation](https://fuse.help/stable/act.html)

# ### Show overall organization
# * Organized by Actor
# * Works the same was as `ini` Parameters

act.ActorEquilibrium

# ## Running actors 
# * [Actors documentation](https://fuse.help/actors.html)
# * We'll manually step through what **Actor** do:
#   * **ActorEquilibrium**: equilibrium (Solovev, TEQUILA, CHEASE) 
#   * **ActorHFSsizing**: Sizes the High Field Side of radial build (plug, OH, TF) superconductors and stresses
#   * **ActorLFSsizing**: Sizes the Low Field Side of radial build
#   * **ActorCXbuild**: Generate the 2D cross section of the build
#   * **ActorPFactive**: Find optimal PF coil locations and currents to match equilibrium boundary shape and a field-null region
#   * **ActorNeutronics**: Calculate neutron loading on the wall
#   * **ActorBlanket**: Blanket tritium breeding ration and heating
#   * **ActorDivertors**: Divertor heat flux
#   * **ActorBalanceOfPlant**: Calculate the net electrical power output
#   * **ActorCosting**: Calculate the cost of the fusion power plant

# Start from scratch
ini, act = FUSE.case_parameters(:FPP);
ini.pf_active.n_coils_outside = 6;
ini.core_profiles.zeff = 2.0;
dd = FUSE.init(ini, act; do_plot=false);

# Checkpoint system allows us to save the state of `dd`, `ini` and `act`
chk = FUSE.Checkpoint()
chk[:my_checkpoint_name] = dd, ini, act;

FUSE.ActorEquilibrium(dd, act; ip_from=:equilibrium, do_plot=true);

# We can use checkpoint system to restore previous conditions to go back and re-run the actor.
# This is useful when working in a Jupyter notebook. One can store a checkpoint at the end of a cell, and restore it at the beginning of another allowing to experiment with parameters and conditions in the second cell without having to re-run the first cell.
dd, ini, act = chk[:my_checkpoint_name]
FUSE.ActorEquilibrium(dd, act; ip_from=:equilibrium);

#

FUSE.ActorHFSsizing(dd, act; do_plot=true);
display(plot(dd.solid_mechanics.center_stack.stress));
IMAS.freeze(dd.build.oh)

#

FUSE.ActorLFSsizing(dd, act; do_plot=true);
IMAS.freeze(dd.build.tf)

#

FUSE.ActorCXbuild(dd, act; do_plot=true);

# weight_currents: make smaller to force currents to fit in the limits
# update_equilibrium: overwrite the target equilibrium with what the coils can actually generate
FUSE.ActorPFactive(dd, act; do_plot=true, update_equilibrium=true);

# Let's run the ActorCXbuild to update the first wall shape based on the new equilibrium
FUSE.ActorCXbuild(dd, act, rebuild_wall=true; do_plot=false);
plot(dd.build; legend=false)

#

FUSE.ActorNeutronics(dd, act; do_plot=true);

#

FUSE.ActorBlanket(dd, act)#; do_plot=true);
IMAS.freeze(dd.blanket)

#

FUSE.ActorDivertors(dd, act)#; do_plot=true);
IMAS.freeze(dd.divertors)

#

FUSE.ActorBalanceOfPlant(dd, act)#; do_plot=true);
IMAS.freeze(dd.balance_of_plant)

#

FUSE.ActorCosting(dd, act)#; do_plot=true);
IMAS.freeze(dd.costing)

# Save `dd` to JSON, so we can share our results with others
@show filename = joinpath(tutorial_temp_dir, "$(ini.general.casename).json")
IMAS.imas2json(dd, filename; freeze=true, strict=false);

# NOTE: can be loaded in OMFIT this way
# OMFIT["FUSE"] = ODS().load(".../fpp_v1_FUSE.json", consistency_check=False)
