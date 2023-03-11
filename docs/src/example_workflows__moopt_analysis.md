# Multi objective optimization post-processing analysis


```@julia
import Revise
import Interact
using FUSE
using Plots
```


```@julia
extract = Dict(
    :Paux => (dd, ini, act) -> ini.ec_launchers.power_launched / 1E6,
    :zeff => (dd, ini, act) -> ini.core_profiles.zeff,
    :κ => (dd, ini, act) -> ini.equilibrium.κ,
    :δ => (dd, ini, act) -> ini.equilibrium.δ,
    :ζ => (dd, ini, act) -> ini.equilibrium.ζ,
    :B0 => (dd, ini, act) -> ini.equilibrium.B0,
    :ip => (dd, ini, act) -> ini.equilibrium.ip / 1E6,
    :R0 => (dd, ini, act) -> ini.equilibrium.R0,
    :beta_n => (dd, ini, act) -> dd.equilibrium.time_slice[].global_quantities.beta_normal,
    :Te0 => (dd, ini, act) -> dd.core_profiles.profiles_1d[].electrons.temperature[1],
    :levelized_CoE => (dd, ini, act) -> dd.costing.levelized_CoE,
    :Pfusion => (dd, ini, act) -> IMAS.fusion_power(dd.core_profiles.profiles_1d[]) / 1E6,
    :Pelectric => (dd, ini, act) -> @ddtime(dd.balance_of_plant.power_electric_net) / 1E6,
    :log10_flattop => (dd, ini, act) -> log10(dd.build.oh.flattop_duration / 3600.0)
)

path= "/mnt/beegfs/users/meneghini/optimization_runs_supersimple_sms"
path= "/mnt/beegfs/users/meneghini/optimization_runs_supersimple_NSGA2"
path= "/mnt/beegfs/users/meneghini/optimization_runs_simple_NSGA2"
dirs = filter(isdir,sort(readdir(path; join=true)))
println(length(dirs))
dirs = filter(x->!isfile(joinpath(x,"error.txt")) && isfile(joinpath(x,"dd.json")),dirs)
println(length(dirs))
outputs=FUSE.load(dirs, [extract], filter_invalid=true)[1];
```


```@julia
scatter(outputs[:,"Pelectric"])
display(hline!([200]))

scatter(outputs[:,"log10_flattop"])
display(hline!([1]))
```


```@julia
plot([histogram(outputs[:,name],title=name,label="") for name in names(outputs)]...)
```


```@julia
# x axis
xname="levelized_CoE"
x=outputs[:,xname]

# y axis
yname="beta_n"
y=outputs[:,yname]

# color
cname="Pelectric"
c=outputs[:,cname]
clim=(-Inf,Inf)

# subselection criterion
s=outputs[:,"Pelectric"]
n=1000.0
index=findall(cc->cc>200-n && cc<200+n, s)#[200:end]
#index=(1:length(x))[100:end]

# plotting
x,y,c,s=x[index], y[index], c[index], s[index]
scatter(x, y, marker_z=c ,xlabel=xname, ylabel=yname, colorbar_title=cname, marker=:circle,
    markersize=5, markerstrokewidth=0, label="", clim=clim,alpha=0.5,xlim=(0,1))
```


```@julia
n=sortperm(c)[1]
println(dirs[n])
dd,ini,act=FUSE.load(dirs[n])
plot(dd.equilibrium,cx=true)
plot!(dd.build)
#plot(dd.solid_mechanics.center_stack.stress)
#dd.blanket
#ini
#plot(dd.core_profiles)
#plot(dd.equilibrium)
#(dd.costing)
#dd.costing
```
