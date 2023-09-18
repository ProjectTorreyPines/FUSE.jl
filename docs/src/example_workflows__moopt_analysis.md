# Multi objective optimization post-processing analysis


```@julia
#import WarmupFUSE
using FUSE
using Plots;
gr();
```


```@julia
# load results from cache file alone
cache_path = "/Users/meneghini/Downloads/";
path = nothing;
outputs = FUSE.extract(path; filter_invalid=:cols, cache=joinpath(cache_path, "extract.csv"), read_cache=true, write_cache=true);
```


```@julia
# on SAGA
path = "/mnt/beegfs/users/meneghini/optimization_run_STEP_nodelta"
path = "/mnt/beegfs/users/meneghini/optimization_run_STEP_nodelta_zerohm"
path = "/mnt/beegfs/users/meneghini/optimization_run_STEP_nodelta_zerohm_expressions"
path = "/mnt/beegfs/users/meneghini/optimization_run_STEP_nodelta_zerohm_expressions_jfix"
path = "/mnt/beegfs/users/meneghini/optimization_run_STEP_nodelta_zerohm_expressions_jfix_noonetime"
path = "/mnt/beegfs/users/meneghini/optimization_run_STEP_nodelta_zerohm_expressions_jfix_slvv"
path = "/mnt/beegfs/users/meneghini/optimization_run_STEP_nodelta_zerohm_expressions_jfix_neo"
path = "/mnt/beegfs/users/meneghini/optimization_run_STEP_nodelta_zerohm_expressions_jfix_neo_SMSEMOA"
path = "/mnt/beegfs/users/meneghini/optimization_run_STEP_nodelta_zerohm_expressions_jfix_neo_SPEA2"
path = "/mnt/beegfs/users/meneghini/optimization_run_STEP_nodelta_zerohm_expressions_jfix_neo_SPEA2_HCD"

# on OMEGA
path = "optimization_run_LTSorHTS_ohtf"
path = "optimization_run_LTSorHTS_ohtf_CCMO"
path = "optimization_run_LTSorHTS_ohtf_fixHFS2"
path = "optimization_run_LTSorHTS_ohtf_fixHFS_req"
path = "optimization_run_LTSorHTS_ohtf_fixHFS_req_dens"
path = "optimization_run_LTSorHTS_ohtf_fixHFS_req_dens_flattop"
path = "optimization_run_LTSorHTS_ohtf_fixHFS_req_dens_flattop_fixHFSagain"
path = "optimization_run_LTSorHTS_ohtf_fixHFS_req_dens_flattop_fixHFSagain_no0ohm"
path = "optimization_run_LTSorHTS_ohtf_fixHFS_req_dens_flattop_fixHFSagain_no0ohm_explore"
path = "optimization_run_LTSorHTS_ohtf_fixHFS_req_dens_flattop_fixHFSagain_no0ohm_explore2"
path = "optimization_run_LTSorHTS_ohtf_fixHFS_req_dens_flattop_fixHFSagain_no0ohm_explore2_minopt"
path = "optimization_run_LTSorHTS_ohtf_fixHFS_req_dens_flattop_fixHFSagain_no0ohm_explore2_minopt_0ohm"
path = "optimization_run_LTSorHTS_ohtf_fixHFS_req_dens_flattop_fixHFSagain_no0ohm_explore2_minopt_maxflat"
path = "optimization_run_LTSorHTS_ohtf_fixHFS_req_dens_flattop_fixHFSagain_no0ohm_explore2_minopt_maxflat_aratio"
cache_path = path
```


```@julia
# load results from directories of optimization run and generate/update cache file

all_dirs = filter(isdir, sort(readdir(path; join=true)))
println(length(all_dirs))

dirs = sort(filter(x -> !isfile(joinpath(x, "error.txt")) && isfile(joinpath(x, "dd.json")), all_dirs))
println(length(dirs))

IMAS.update_ExtractFunctionsLibrary!() # to pick up any ongoing development to extract function library
outputs = FUSE.extract(dirs; filter_invalid=:cols, cache=joinpath(cache_path, "extract.csv"), read_cache=true, write_cache=true);
```


```@julia
# error analysis
errors = FUSE.categorize_errors(all_dirs; do_plot=true, show_first_line=true)
#println(read("optimization_run_LTSorHTS_ohtf/2023-04-28T11:34:32.300__3823851/error.txt",String))
```


```@julia
# this is just to list all the fields that can be queried
IMAS.ExtractFunctionsLibrary
```


```@julia
# plot evolution of key scalar quantities
N = length(outputs[:, "Pelectric_net"])
X = LinRange(1, N, 50)

gr()

import Statistics
function y_auto_range(y; σ=5, N=50)
    y_nonan = y[@. !isnan.(y)]
    m = Statistics.median(y_nonan)
    s = Statistics.median(Statistics.median(abs.(y_nonan .- m))) * σ
    Y = LinRange(max(m - s, minimum(y_nonan)), min(m + s, maximum(y_nonan)), N)
end

y = outputs[:, "Pelectric_net"]
yname = "Pelectric net [MW]"
Y = y_auto_range(y)
p = histogram2d(y, bins=(X, Y), ylabel=yname, xlabel="Indiviuals")
#scatter!(y, ylim=(minimum(Y),maximum(Y)), label="")
hline!([200.0], ls=:dash, label="")
display(p)

y = outputs[:, "βn"]
yname = "βn net"
Y = y_auto_range(y)
p = histogram2d(y, bins=(X, Y), ylabel=yname, xlabel="Indiviuals")
#scatter!(y, ylim=(minimum(Y),maximum(Y)), label="")
display(p)

y = outputs[:, "capital_cost"]
yname = "capital_cost [\$B]"
Y = y_auto_range(y)
p = histogram2d(y, bins=(X, Y), ylabel=yname, xlabel="Indiviuals")
#scatter!(y, ylim=(minimum(Y),maximum(Y)), label="")
display(p)

y = log10.(outputs[:, "flattop"])
yname = "Log10(flattop)"
Y = y_auto_range(y)
p = histogram2d(y, bins=(X, Y), ylabel=yname, xlabel="Indiviuals")
#scatter!(y, ylim=(minimum(Y),maximum(Y)), label="")
display(p)

y = outputs[:, "ip_ohm"]
yname = "Ip ohmic [MA]"
Y = y_auto_range(y)
p = histogram2d(y, bins=(X, Y), ylabel=yname, xlabel="Indiviuals")
hline!([0.0], ls=:dash, label="")
#scatter!(y, ylim=(minimum(Y),maximum(Y)), label="")
display(p)

y = (outputs[:, "TF_material"] .== "ReBCO") .+ (outputs[:, "OH_material"] .== "ReBCO") .* 2.0;
yname = "HTS";
Y = LinRange(0, 3.00001, 50)
p = histogram2d(y, bins=(X, Y), ylabel=yname, xlabel="Indiviuals")
#scatter!(y, ylim=(minimum(Y),maximum(Y)), label="")
display(p)
```


```@julia
# plot evolution of all scalar quantities extracted from the dataset

gr()

for name in names(outputs)
    if !(typeof(outputs[1, name]) <: Number) || all(isnan.(outputs[:, name]))
        continue
    end
    y = outputs[:, name]
    y_nonan = y[@. !isnan.(y)]
    m = Statistics.median(y_nonan)
    σ = Statistics.median(Statistics.median(abs.(y_nonan .- m))) * 10
    Y = LinRange(max(m - σ, minimum(y_nonan)), min(m + σ, maximum(y_nonan)), 100)
    if σ != 0.0
        display(histogram2d(y, bins=(X, Y), ylabel=name))
    end
end
```


```@julia
# plot in optimization space

# x axis
xname = "capital_cost";
x = outputs[:, xname];
#xname="βpol"; x=outputs[:,xname]
#xname="Pec [MW]"; x=outputs[:,"Pec"]

# y axis
yname = "βn";
y = outputs[:, yname];
#yname="ip_bs/ip_bs_aux_ohm"; y=outputs[:,"ip_bs"]./outputs[:,"ip_bs_aux_ohm"];
#yname="Pec=R0*ne*ip_aux(5+zeff)/(Te*0.09)"; y= @. outputs[:,"<ne>"]/1E20*outputs[:,"R0"]*outputs[:,"ip_aux"]*(5.0+outputs[:,"<zeff>"])/(0.09*outputs[:,"<Te>"])

# z axis
zname = "log10(flattop)";
z = log10.(outputs[:, "flattop"]);

# color
#cname="individual"; c=1:length(outputs[:,xname]);clim=(1,Inf);
#cname="log10(flattop)"; c=log10.(outputs[:,"flattop"]);clim=(-Inf,Inf);
#cname="Pelectric_net"; c=outputs[:,cname];clim=(-Inf,Inf);
#cname="Ip aux [MA]"; c=outputs[:,"ip_aux"];clim=(-Inf,Inf);
#cname="fGW"; c=outputs[:,"fGW"];clim=(-Inf,Inf);
#cname="<zeff>"; c=outputs[:,"<zeff>"];clim=(-Inf,Inf);
#cname="R0 [m]"; c=outputs[:,"R0"];clim=(-Inf,Inf);
#cname="B0 [T]"; c=outputs[:,"B0"];clim=(-Inf,Inf);
#cname="B0*R0 [T*m]"; c=outputs[:,"B0"].*outputs[:,"R0"];clim=(-Inf,Inf);
#cname="ip_ohm"; c=outputs[:,cname];clim=(-Inf,Inf);
#cname="ip"; c=outputs[:,cname];clim=(-Inf,Inf);
#cname = "q95"; c = outputs[:, cname]; clim = (-Inf, Inf);
cname = "1/ϵ"; c = outputs[:, "R0"]./outputs[:, "a"]; clim = (-Inf, Inf);
#cname="Pec [MW]"; c=outputs[:,"Pec"];clim=(-Inf,Inf);
#cname="ip_bs_aux_ohm/ip"; c=outputs[:,"ip_bs_aux_ohm"]./outputs[:,"ip"];clim=(-Inf,Inf);
#cname="ip_ohm"; c=abs.(outputs[:,"ip_ohm"]);clim=(-Inf,Inf);
#cname="ip_bs_aux_ohm/ip"; c=outputs[:,"ip_bs_aux_ohm"]./outputs[:,"ip"];clim=(-Inf,Inf);
#cname="Ip aux [MA]"; c=outputs[:,"ip_aux"];clim=(-Inf,Inf);

#cname="HTS"; c=(outputs[:,"TF_material"].=="ReBCO").+(outputs[:,"OH_material"].=="ReBCO").*2.0; clim=(0,4)

# selection
sname0 = "Pelectric_net";
s0 = outputs[:, sname0];
sname1 = "ip_ohm";
s1 = outputs[:, sname1];
sname2 = "flattop";
s2 = outputs[:, sname2];

# Selection criterion
min_Pelectric = -10000
min_flattop = 00

n = length(x)
index = 1:n
index = findall((s0 .> min_Pelectric) .&& (s2 .> min_flattop) .&& (index .> 0000))
println("$(length(index)) points")

gr()

# plotting
annot = map(x -> (x, :center, 3, "courier"), index)
P = scatter(x[index], y[index], marker_z=c[index], xlabel=xname, ylabel=yname, colorbar_title=cname, marker=:circle,
    markersize=5, markerstrokewidth=0, label="", clim=clim, alpha=0.5,
    #    series_annotations=annot,
    xlim=(0, 20), ylim=(0, 3.5))

# Pareto front
pindex = index[FUSE.pareto_front([[x[index[k]], y[index[k]]] for k in 1:length(index)])]
sort!(pindex, by=i -> y[i])
sort!(pindex, by=i -> x[i])
println(length(pindex))
println(pindex)
pannot = map(x -> ("\n$x", :right, 3, "courier", :red), pindex)
plot!(P, x[pindex], y[pindex], series_annotations=pannot, color=:blue, label="pareto front", lw=2)
```


```@julia
# Interactive 2D plotting

using Interact
@manipulate for min_Pelectric in LinRange(-1000, 1000, 101), min_flattop in LinRange(0, 1000, 101)

    n = length(x)
    index = 1:n
    index = findall((s0 .> min_Pelectric) .&& (s2 .> min_flattop) .&& (index .> 0000))
    #println("$(length(index)) points")

    # plotting
    annot = map(x -> (x, :center, 3, "courier"), index)
    P = scatter(x[index], y[index], marker_z=c[index], xlabel=xname, ylabel=yname, colorbar_title=cname, marker=:circle,
        markersize=5, markerstrokewidth=0, label="", clim=clim, alpha=0.5,
        xlim=(0, 20), ylim=(0, 3.5))

    # Pareto front
    pindex = index[FUSE.pareto_front([[x[index[k]], y[index[k]]] for k in 1:length(index)])]
    sort!(pindex, by=i -> y[i])
    sort!(pindex, by=i -> x[i])
    pannot = map(x -> ("\n$x", :right, 3, "courier", :red), pindex)
    plot!(P, x[pindex], y[pindex], series_annotations=pannot, color=:blue, label="pareto front", lw=2)
end
```


```@julia
# Interactive 3D plotting

plotly()
min_Pelectric = 200
index = findall((s0 .> min_Pelectric) .&& (x .< 15))
pindex = index[FUSE.pareto_front([[x[index[k]], y[index[k]], z[index[k]]] for k in 1:length(index)])]
println("$(length(pindex)) points")
pannot = map(x -> ("\n$x", :right, 3, "courier", :red), pindex)
plot(x[pindex], y[pindex], z[pindex], marker_z=c[pindex],
    xlabel=xname, ylabel=yname, zlabel=zname, colorbar_title=cname, label=cname,
    marker=:dot, xlim=(0, Inf), ylim=(0, Inf), zlim=(0, Inf))
```


```@julia
n = 17134
dir = outputs[n, "dir"]
#n=length(dirs)
println(dir)
dd, ini, act = FUSE.load(dir)
FUSE.digest(dd)
```
