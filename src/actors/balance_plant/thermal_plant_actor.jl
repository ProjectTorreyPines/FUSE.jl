#= ================= =#
#  ActorThermalPlant
#= ================= =#
# ACTOR FOR THE INTERMEDIATE HEAT TRANSFER SYSTEM
import ModelingToolkit as MTK, DifferentialEquations
import ThermalSystem_Models
TSMD    = ThermalSystem_Models.Dynamics
MTK.@variables t

Base.@kwdef mutable struct FUSEparameters__ActorThermalPlant{T} <: ParametersActor where {T<:Real}
    _parent::WeakRef = WeakRef(Nothing)
    _name::Symbol   = :not_set
    model::Switch{Symbol}      = Switch{Symbol}([:brayton,:rankine], "-", "user can only select one of these"; default = :rankine) # for future additions
    heat_load_from::Switch{Symbol} = Switch{Symbol}([:dd, :actor],"-",""; default = :dd)
    do_plot::Entry{Bool} = act_common_parameters(do_plot=false)
    verbose::Entry{Bool} = act_common_parameters(verbose=false)
end

mutable struct ActorThermalPlant{D,P} <: FacilityAbstractActor{D,P}
    dd::IMAS.dd{D}
    par::FUSEparameters__ActorThermalPlant{P}       # Actors must carry with them the parameters they are run with
    components::Vector{ModelingToolkit.ODESystem}   # Vector of type ODESystem
    connections::Vector{ModelingToolkit.Equation}   # Connection equations
    odeparams::Vector{ModelingToolkit.Num}          # Circuit Parameters 
    odedict::Dict{MTK.Symbol,MTK.ODESystem}         # Dictionary where symbol name => symbol
    buildstatus::Bool
    fullbuild                                       #::MTK.ODESystem - high level ODESystem
    plant       #::MTK.ODESystem                        # simplified ODESystem
    prob        #::MTK.ODEProblem
    G
    gplot
    sym2var::Dict{}
    var2val::Dict{}
    optpar::Vector{}
    x
    u
    # Inner constructor for the actor starting from `dd` and `par` (we generally refer to `par` as `act.ThermalPlant`)
    # NOTE: Computation should not happen here since in workflows it is normal to instantiate
    #       an actor once `ThermalPlant(dd, act.ThermalPlant)` and then call `finalize(step(actor))` several times as needed.
    function ActorThermalPlant(dd::IMAS.dd{D}, par::FUSEparameters__ActorThermalPlant{P}; kw...) where {D<:Real,P<:Real}
        logging_actor_init(ActorThermalPlant)
        par = par(kw...)
        dd.balance_of_plant.power_plant.power_cycle_type = lowercase(string(par.model))
        return new{D,P}(dd, 
                        par, 
                        ModelingToolkit.ODESystem[],
                        ModelingToolkit.Equation[],
                        ModelingToolkit.Num[],
                        Dict{ModelingToolkit.Symbol,ModelingToolkit.ODESystem}(),
                        false,
                        nothing,
                        nothing,
                        nothing,
                        nothing,
                        nothing,
                        Dict{}(),
                        Dict{}(),
                        Symbol[],
                        nothing,
                        nothing);
    end
end

"""
    ActorThermalPlant(dd::IMAS.dd, act::ParametersAllActors; kw...)

!!! note 
    Stores data in `dd.balance_of_plant`
"""
function ActorThermalPlant(dd::IMAS.dd, act::ParametersAllActors; stepkw, kw...)
    par = act.ActorThermalPlant(kw...)
    actor = ActorThermalPlant(dd, par)
    dd.balance_of_plant.power_plant.power_cycle_type = lowercase(string(par.model))
    actor = step(actor; stepkw...)
    finalize(actor)
    return actor
end

function _step(actor::ActorThermalPlant)

    
    # if use_actor_u is true then the actor will use the loading values in Actor.u instead of from dd
    # println("Actor run - use u $(use_actor_u)")
    dd  = actor.dd
    par = actor.par

    use_actor_u = par.heat_load_from == :dd ? false : true


    breeder_heat_load   = (use_actor_u == false) ? (isempty(dd.blanket.module) ? 0.0 : sum(bmod.time_slice[].power_thermal_extracted for bmod in dd.blanket.module)) : actor.u[1];
    divertor_heat_load  = (use_actor_u == false) ? (isempty(dd.divertors.divertor) ? 0.0 : sum((@ddtime(div.power_incident.data)) for div in dd.divertors.divertor)) : actor.u[2];
    wall_heat_load      = (use_actor_u == false) ?  abs.(IMAS.radiation_losses(dd.core_sources)) : actor.u[3];
    
    actor.u             =  [breeder_heat_load, divertor_heat_load, wall_heat_load]
    # Buidling the TSM System
    if actor.buildstatus == false
        @info "Rebuilding ActorThermalPlant"

        # This is for graph construction which currently relies on the heat flow trajectory through the full BOP 
        # In one of the steps to to plot the graph, (TSMD.create_plot_graph()) there is a routine which 
        # has to take a directed cyclic graph and convert it into a directed acyclic graph, to do this correctly TSM
        # relies on the heat flow (solution) trajectory through the entire plant in order to find the correct edge direction ( within the MetaGraph object), 
        # since balance equations are directionless without context, 
        # This only matters for the first step since that is when TSMD builds the graph object, afterwards 0 values are not an issue
        if 0 ∈ actor.u 
            @warn "Invalid initial thermal loading is [Qbreeder, Qdivertor, Qwall] = $(actor.u), \n Setting load to default for construction step, resetting to [500e6,100e6,100e6]), \n Rerun FUSE.step(act.ActorThermalPlant) to update loading to the correct value"
            breeder_heat_load = 500e6;
            divertor_heat_load = 100e6;
            wall_heat_load = 100e6;
        end

        # default parameters
        tspan = (0.0,10)
        Tmax_wall    = 950;     # Maximum cooling temperature for first wall
        Tmin_wall    = 350;     # Minimum cooling temperature for first wall
        Tmax_div     = 1000;    # Maximum cooling temperature for Divertor
        Tmin_div     = 350;     # Minimum cooling temperature for Divertor
        Tmax_breeder = par.model == :rankine ? 1136 : 1300; # Maximum cooling temperature for Breeder blanket
        Tmin_breeder = par.model == :rankine ? 674.7 : 900; # Minimum cooling temperature for Breeder blanket
        Nhx         = 4;                  # Nhx = the number of heat exchangers to add to the loop, 4: 3 to connect to cooling loops, 1 to connect to primary power cycle
        flowrate    = 300;                # mass flow rate of the intermediate loop (kg/s)
        Tmin_interloop =  350; # minimum temperature for inter loop (Kelvin)

        
        energy_sys, sts, edict = TSMD.default_energy_sys();
        η_cycle, η_bop = sts;   # add them to current namespace 
        
        # Wall circuit, Helium
        wall_sys,     wall_connections,     wparams, wdict  = TSMD.wall_circuit(; load = wall_heat_load, Tmin = Tmin_wall, Tmax = Tmax_wall);

        # Divertor circuit, Helium
        divertor_sys, divertor_connections, dparams, ddict  = TSMD.divertor_circuit(; load = divertor_heat_load, Tmin = Tmin_div, Tmax = Tmax_div);

        # Breeder Circuit (PbLi
        breeder_sys,  breeder_connections,  bparams, bdict  = TSMD.breeder_circuit(; load = breeder_heat_load, Tmin = Tmin_breeder, Tmax = Tmax_breeder);
        
        # intermediate loop
        inter_loop_sys, inter_loop_connections, iparams, idict = TSMD.intermediate_loop(; Nhx = Nhx, flowrate = flowrate, Tmin = Tmin_interloop);
        
        if par.model == :rankine
            cycle_flowrate = 250;      # kg/s
            ηpump          = 0.7;      # isentropic effeciency of the pump
            ηturbine       = 0.95;    # Isentropic effeciency of the turbine

            # Plant
            steam_systems, steam_connections, sparams, sdict = TSMD.feedwater_rankine(; flowrate = cycle_flowrate, ηpump = ηpump, ηturbine = ηturbine);
            
            # Create heat exchangers which will couple the indepentent loops
            MTK.@named hx1 = TSMD.Gen_HeatExchanger(
                B = idict[:inter_loop_hx1],
                A = wdict[:wall_hx],
                returnmode = :eq,
            );

            MTK.@named hx2 = TSMD.Gen_HeatExchanger(
                B = idict[:inter_loop_hx2],
                A = ddict[:divertor_hx],
                returnmode = :eq,
            );

            MTK.@named hx3 = TSMD.Gen_HeatExchanger(
                B = idict[:inter_loop_hx3],
                A = bdict[:breeder_hx],
                returnmode = :eq,
            );

            MTK.@named boilhx = TSMD.S2G_HeatExchanger(
                A = sdict[:steam_boiler],
                B = idict[:inter_loop_hx4],
                returnmode = :eq,
            );

            # Connect all energy reservoirs to external interfacing components
            energy_connections = vcat(
                TSMD.work_connect(
                    edict[:Electric],
                    wdict[:wall_circulator].w,
                    ddict[:divertor_circulator].w,
                    bdict[:breeder_circulator].w,
                    idict[:inter_loop_circulator].w,
                    sdict[:steam_hp_pump].w,
                    sdict[:steam_lp_pump].w,
                    sdict[:steam_turbine].hp.w,
                    sdict[:steam_turbine].lp.w,
                ),
                TSMD.heat_connect(
                    edict[:HotUtility],
                    wdict[:wall_heat].q,
                    ddict[:divertor_heat].q,
                    bdict[:breeder_heat].q,
                ),
                TSMD.heat_connect(
                    edict[:ColdUtility],
                    wdict[:wall_relief].q,
                    ddict[:divertor_relief].q,
                    bdict[:breeder_relief].q,
                    idict[:inter_loop_relief].q,
                    sdict[:steam_condensor].q,
                ),
                η_cycle ~ 1 - abs(sdict[:steam_condensor].q.Q̇ / sdict[:steam_boiler].q.Q̇),
                η_bop ~ 1 - abs(edict[:ColdUtility].Q̇ / edict[:HotUtility].Q̇),
            );

            
            # Create vector of all parameters
            plant_params = vcat(wparams, dparams, bparams, iparams, sparams);

            # Create total vector for all connecting equations for flow connected components and the energy_connections
            plant_connections = vcat(
                steam_connections,
                inter_loop_connections,
                wall_connections,
                divertor_connections,
                breeder_connections,
                energy_connections,
            );

            # Create total vector for all components
            plant_systems = vcat(steam_systems, inter_loop_sys, wall_sys, divertor_sys, breeder_sys, energy_sys);

            # add heat exchanger equations to plant_connections
            push!(plant_connections, hx1...);
            push!(plant_connections, hx2...);
            push!(plant_connections, hx3...);
            push!(plant_connections, boilhx...);

            # Create total ODESystem for the plant
            MTK.@named sys = MTK.ODESystem(plant_connections, t, sts, plant_params; systems = plant_systems);

            # Check DOF and problem size
            # TSMD.system_details(sys);

            # Simplify using ModelingToolkit's model reduction methods
            simple_sys = MTK.structural_simplify(sys);


            actor.components  = plant_systems
            actor.connections = plant_connections
            actor.odeparams   = plant_params
            actor.odedict     = Dict(edict..., wdict..., ddict..., bdict..., idict..., sdict...);
            actor.buildstatus = true
            actor.fullbuild   = sys;
            actor.plant       = simple_sys;

            actor.optpar = [:steam_ṁ,
                        :inter_loop_ṁ,
                        :inter_loop_supply₊T,
                        :wall_supply₊T,
                        :wall_heat₊Tout,
                        :divertor_supply₊T,
                        :divertor_heat₊Tout,
                        :breeder_supply₊T,
                        :breeder_heat₊Tout];


            actor.prob = MTK.ODEProblem(simple_sys, [], tspan);

            ode_sol  = DifferentialEquations.solve(actor.prob, DifferentialEquations.ImplicitEuler());
            soln(v) = ode_sol[v][end]

            utility_vector = [:HotUtility, :ColdUtility, :Electric];
            actor.G = TSMD.system2metagraph(sys, utility_vector; soln = soln, verbose = false);
            
            para_vars = MTK.parameters(simple_sys);
            para_syms = TSMD.variable2symbol(MTK.parameters(simple_sys))
            para_vals = actor.prob.p;

            actor.sym2var = Dict(para_syms[i] => para_vars[i] for i =1:length(para_vars));
            actor.var2val = Dict(para_vars[i] => para_vals[i] for i =1:length(para_vars));

            
            gcopy = TSMD.create_plot_graph(actor.G; toignore = [:steam_condensor], verbose = false);
            xLayReqs, vSortReqs, xs, ys, paths, lay2node = TSMD.layers_to_force!(gcopy);
            TSMD.initialize_plot_props!(gcopy, lay2node,xs,ys,paths);
            TSMD.add_plot_elments!(gcopy);
            TSMD.set_default_node_prop!(gcopy, :height, 1.0);
            xLayReqs, vSortReqs, xs, ys, paths, lay2node = TSMD.layers_to_force!(gcopy);
            x, y = TSMD.setVerticalSpacing!(gcopy; vspan = 40.0);
            TSMD.setLayerWidth!(gcopy; pad = 2.5, verbose = false);
            xLayReqs, vSortReqs, xs, ys, paths, lay2node = TSMD.layers_to_force!(gcopy);
            TSMD.edgeroute_nodes(gcopy; voff = 0.1);
            TSMD.set_plot_props!(gcopy);
            
            actor.gplot = gcopy;

        elseif par.model == :brayton
            cyclesys, cconnections, cparams, cdict = TSMD.brayton_regenerator(;flowrate = 300)
            energy_con = vcat(
            TSMD.work_connect(
                edict[:Electric],
                wdict[:wall_circulator].w,
                ddict[:divertor_circulator].w,
                bdict[:breeder_circulator].w,
                idict[:inter_loop_circulator].w,
                cdict[:cycle_compressor_hp].w,
                cdict[:cycle_compressor_mp].w,
                cdict[:cycle_compressor_lp].w,
                cdict[:cycle_turbine].w,
            ),
            TSMD.heat_connect(
                edict[:HotUtility],
                wdict[:wall_heat].q,
                ddict[:divertor_heat].q,
                bdict[:breeder_heat].q,
            ),
            TSMD.heat_connect(
                edict[:ColdUtility],
                wdict[:wall_relief].q,
                ddict[:divertor_relief].q,
                bdict[:breeder_relief].q,
                idict[:inter_loop_relief].q,
                cdict[:cycle_cooler].q,
                cdict[:cycle_intercooler_1],
                cdict[:cycle_intercooler_2],
            ),
            η_cycle ~ 1 - abs((cdict[:cycle_cooler].q.Q̇+cdict[:cycle_intercooler_1].Q̇ + cdict[:cycle_intercooler_2].Q̇) / cdict[:cycle_heat].q.Q̇),
            η_bop ~ 1 - abs(edict[:ColdUtility].Q̇ / edict[:HotUtility].Q̇))
            
            plant_params = vcat(wparams, dparams, bparams, iparams, cparams)
            plant_connections = vcat(
                cconnections,
                inter_loop_connections,
                wall_connections,
                divertor_connections,
                breeder_connections,
                energy_con,
            );

            plant_systems = vcat(cyclesys, inter_loop_sys, wall_sys, divertor_sys, breeder_sys, energy_sys)

            MTK.@named hx1 = TSMD.Gen_HeatExchanger(
                B = idict[:inter_loop_hx1],
                A = wdict[:wall_hx],
                returnmode = :eq,
            )
            MTK.@named hx2 = TSMD.Gen_HeatExchanger(
                B = idict[:inter_loop_hx2],
                A = ddict[:divertor_hx],
                returnmode = :eq,
            )
            MTK.@named hx3 = TSMD.Gen_HeatExchanger(
                B = idict[:inter_loop_hx3],
                A = bdict[:breeder_hx],
                returnmode = :eq,
            )
            MTK.@named hx4 = TSMD.Gen_HeatExchanger(
                A = cdict[:cycle_heat],
                B = idict[:inter_loop_hx4],
                returnmode = :eq,
            )

            push!(plant_connections, hx1...)
            push!(plant_connections, hx2...)
            push!(plant_connections, hx3...)
            push!(plant_connections, hx4...)

            # Create total ODESystem for the plant
            MTK.@named sys = MTK.ODESystem(plant_connections, t, sts, plant_params; systems = plant_systems);

            # Check DOF and problem size
            # TSMD.system_details(sys);

            # Simplify using ModelingToolkit's model reduction methods
            simple_sys = MTK.structural_simplify(sys);


            actor.components  = plant_systems;
            actor.connections = plant_connections;
            actor.odeparams   = plant_params
            actor.odedict     = Dict(edict..., wdict..., ddict..., bdict..., idict..., cdict...);
            actor.buildstatus = true
            actor.fullbuild = sys;
            actor.plant = simple_sys;

            actor.optpar = [:cycle_ṁ,
                            :inter_loop_ṁ,
                            :inter_loop_supply₊T,
                            :wall_supply₊T,
                            :wall_heat₊Tout,
                            :divertor_supply₊T,
                            :divertor_heat₊Tout,
                            :breeder_supply₊T,
                            :breeder_heat₊Tout];

            tspan    = (0.0, 10)
            actor.prob = MTK.ODEProblem(simple_sys, [], tspan);

            ode_sol  = DifferentialEquations.solve(actor.prob, DifferentialEquations.ImplicitEuler());
            sol(v) = ode_sol[v][end]

            utility_vector = [:HotUtility, :ColdUtility, :Electric];
            actor.G = TSMD.system2metagraph(sys, utility_vector; soln = sol, verbose = false);
            
            para_vars = MTK.parameters(simple_sys);
            para_syms = TSMD.variable2symbol(MTK.parameters(simple_sys))
            para_vals = actor.prob.p;

            actor.sym2var = Dict(para_syms[i] => para_vars[i] for i =1:length(para_vars));
            actor.var2val = Dict(para_vars[i] => para_vals[i] for i =1:length(para_vars));

            gcopy = TSMD.create_plot_graph(actor.G; toignore = [:cycle_cooler], verbose = false);
            xLayReqs, vSortReqs, xs, ys, paths, lay2node = TSMD.layers_to_force!(gcopy);
            TSMD.initialize_plot_props!(gcopy, lay2node,xs,ys,paths);
            TSMD.add_plot_elments!(gcopy);
            TSMD.set_default_node_prop!(gcopy, :height, 1.0);
            xLayReqs, vSortReqs, xs, ys, paths, lay2node = TSMD.layers_to_force!(gcopy);
            x, y = TSMD.setVerticalSpacing!(gcopy; vspan = 40.0);
            TSMD.setLayerWidth!(gcopy; pad = 2.5, verbose = false);
            xLayReqs, vSortReqs, xs, ys, paths, lay2node = TSMD.layers_to_force!(gcopy);
            TSMD.edgeroute_nodes(gcopy; voff = 0.1);
            TSMD.set_plot_props!(gcopy);
            actor.gplot = gcopy;
        end
        actor.x = [getval(a,actor) for a in actor.optpar]
    end
    
    soln = plant_wrapper(actor)
    TSMD.updateGraphSoln(actor.G,soln);
    TSMD.updateGraphSoln(actor.gplot,soln);

    # write to dd if ddwrite = true
    initddbop(actor; soln = soln)

    if par.do_plot == true
        sysnamedict = Dict([
                "cycle_" => "Brayton Helium",
                "steam_" => "Feedwater Rankine",
                "divertor_" => "Divertor Helium",
                "breeder_" => "Breeder PbLi",
                "inter_loop_" => "Inter. loop Helium",
                "wall_" => "first wall helium",
                "ColdUtility" => "Sink",
                "HotUtility" => "Fusion Core",
            ])
        p = TSMD.plotplant(
            actor.gplot;
            numbering = false,
            mode = :path,
            nsize = 2.0,
            compnamesubs = (
                "enfw" => "en\nfw",
                "_" => "\n",
                "circulator" => "pump",
                "hotutility" => "Fusion\nCore",
                "coldutility" => "Sink",
                "relief" => "Trim\ncooler",
                "condensor" => "Cool\nHX",
            ),
            compnameattr = (:right,  4),
            compnamerot = 90,
            sysnamedict = sysnamedict,
            legpad = 0.5,
            legwid =13,
            legheight = 1.5,
            legoffset = 2.0,
            pathattr = (
                linewidth = 1,
                marker = false,
                markersize = 0.0,
                markercolor = :red,
                alpha = 0.7,
                legend = false,
            ),
            figattr = (
                grid = false,
                aspect_ratio = :equal,
                showaxis = false,
                xlim = [-15, 75],
                ylim = [-21, 21],
                xticks = [0, 1, 2, 3, 4, 5, 6, 7],
                plot_title = "Balance Of Plant Circuit: $(titlecase(string(par.model)))",
                plot_titlefonthalign = :hcenter,
                plot_titlefontvalign = :bottom,
                dpi = 200,
                plot_titlevspan = .0001,
            ), #aspect_ratio = :equal
        )
        # display(p)
    end
    return actor
end

"""
    getval(v,act::ActorThermalPlant)

Returns the value of a system parameter identified by symbol v.
Inputs:
    v::Symbol corresponding to model parameters
    act::ActorThermalPlant
Outputs:
    The actors default value for that parameter
"""
function getval(v,act::ActorThermalPlant)
    return act.var2val[act.sym2var[v]]
end

"""
    getvar(v,act::ActorThermalPlant)

Returns the variable object of a system parameter identified by symbol v.
Inputs:
    v::Symbol corresponding to model parameters
    act::ActorThermalPlant
Outputs:
    The variable object associated with symbol v 
Example:
    getvar(:η_bop, actor) === actor.plant.η_bop
    getvar(:Electric₊Ẇ) === actor.plant.Electric.Ẇ ===  actor.odedict[:Electric].Ẇ
"""  
function getvar(v,act::ActorThermalPlant)
    return act.sym2var[v]
end

"""
    getsol(act::ActorThermalPlant)

Returns the current solution object of a the ThermalPlantActor. The returned function will output the solution point for any VARIABLE object within the plant.
"""
function getsol(act::ActorThermalPlant)
    return act.gplot.gprops[:soln]
end


"""
    plant_wrapper(x, u, simple_sys, keypara, var2val, sym2var; tspan = (0,10))

Evaluates the system described by simple_sys::ODESystem.
Inputs:
    u = heat loading [Q_breeder, Q_divertor, Q_wall] (Watts)
    x = parameters values to use during evaluation, these are identified in keypara, 
    keypara = variabels associated with the data in indices of x
    var2val and sym2var are dicts for getting the actual variable objects, they are the same as in the actor structure
Ouput: 
    soln = solution object
"""
function plant_wrapper(x, u, simple_sys, keypara, var2val, sym2var; tspan = (0,10))
    # x are parameters 
    # u are heat loads 
    #   u[1] = Qbreeder
    #   u[2] = Qdivertor
    #   u[3] = Qwall
    # y are output variables
    # 

    # new parameters dict
    pwrapped = var2val;

    pwrapped[sym2var[:Qbreeder]]      = u[1];
    pwrapped[sym2var[:Qdivertor]]     = u[2];
    pwrapped[sym2var[:Qwall]]         = u[3];

    # x are parameters
    for (i,xi) ∈ enumerate(keypara)
        pwrapped[sym2var[xi]]   =   x[i]
    end

    node_prob = MTK.ODEProblem(simple_sys, [], tspan, pwrapped);
    node_sol  = DifferentialEquations.solve(node_prob, DifferentialEquations.ImplicitEuler());
    soln(v) = node_sol[v][end]
    return soln
end

"""
    plant_wrapper(x, u, ``yvars``, simple_sys, keypara, var2val, sym2var; tspan = (0,10))

Evaluates the system described by simple_sys::ODESystem.
Inputs:
    yvars = vector of output variable objects
Ouput: 
    soln.(yvars)
"""
function plant_wrapper(x, u, yvars, simple_sys, keypara, var2val, sym2var; tspan = (0,10))
    # x are parameters 
    # u are heat loads 
    #   u[1] = Qbreeder
    #   u[2] = Qdivertor
    #   u[3] = Qwall
    # y are output variables
    # 

    # new parameters dict
    pwrapped = var2val;

    pwrapped[sym2var[:Qbreeder]]      = u[1];
    pwrapped[sym2var[:Qdivertor]]     = u[2];
    pwrapped[sym2var[:Qwall]]         = u[3];

    # x are parameters
    for (i,xi) ∈ enumerate(keypara)
        pwrapped[sym2var[xi]]   =   x[i]
    end

    node_prob = MTK.ODEProblem(simple_sys, [], tspan, pwrapped);
    node_sol  = DifferentialEquations.solve(node_prob,DifferentialEquations.ImplicitEuler());
    soln(v) = node_sol[v][end]
    return soln.(yvars)
end

"""
    plant_wrapper(act::ActorThermalPlant)

Evaluates the system described by act.plant
Ouput: 
    soln, the updated solution object
"""
function plant_wrapper(act::ActorThermalPlant)
    return plant_wrapper(act.x, act.u, act.plant, act.optpar, act.var2val, act.sym2var)
end

"""
    plant_wrapper(act::ActorThermalPlant, yvars)

Evaluates the system described by act.plant
Ouput: 
    soln.(yvars) at the updated solution for the vars stored in yvars
"""
function plant_wrapper(act::ActorThermalPlant,yvars)
    return plant_wrapper(act.x, act.u, yvars, act.plant, act.optpar, act.var2val, act.sym2var)
end

"""
    plant_wrapper(act::ActorThermalPlant,yvars,yfunc)

Evaluates the system described by act.plant and returns objective function value of yfunc. This function can be used during optimization trials.
Ouput: 
    yfunc(soln.(yvars))
"""
function plant_wrapper(act::ActorThermalPlant,yvars,yfunc)
    y = plant_wrapper(act,yvars);
    return yfunc(y)
end

"""
    initddbop(act::ActorThermalPlant; soln = nothing)

Maps data stored in the TSM objects and metagraph to dd. By default the function will use 
the internal solution value in the actor, which is updated during every step call and plant_wrapper call.
If you want to write to dd based off a different solution object, it can be passed in the kwargs 
"""
function initddbop(act::ActorThermalPlant; soln = nothing)
    gcopy = act.gplot;
    dd = act.dd;
    # syslabs = TSMD.get_prop(gcopy, :system_labels)
    # empty!(dd.balance_of_plant);
    # begin
    bop       = dd.balance_of_plant;
    bop_plant   = bop.power_plant;
    bopsys      = bop_plant.system; 

    compnamesubs = Dict(
        "enfw" => "en fw",
        "_" => " ",
        "circulator" => "pump",
        "hotutility" => "Fusion_Core",
        "coldutility" => "Sink",
        "relief" => "Trim cooler",
        "condensor" => "Cool HX",
    )

    gp = getproperty(gcopy,:gprops);
    np = getproperty(gcopy,:vprops);
    nv_g = maximum(collect(keys(np)))
    soln = (isnothing(soln) ? gp[:soln] : soln)

    # names of the internal subgraph objects
    syslabs = [titlecase(replace(lowercase(string(sl)),compnamesubs...)) for sl in TSMD.get_prop(gcopy, :system_labels)]
    format_name(x) = titlecase(replace(lowercase(string(x)),compnamesubs...))

    # initializing the 1st level of dd
    if length(bop_plant.system) != length(syslabs)
        empty!(dd.balance_of_plant.power_plant.system)
        resize!(bop_plant.system,length(syslabs))
        for i =1:length(syslabs)
            bop_plant.system[i].name  = syslabs[i] 
        end
    end
    bopsys = bop_plant.system; 
    bops_dict = Dict(bopsys[i].name => i for i =1:length(syslabs)); # dict where name => index
    valid_s   = collect(keys(bops_dict)); 


    nparent_dict = TSMD.node_propdict(gcopy,:parent);
    nname_dict   = TSMD.node_propdict(gcopy,:name);
    isreal_dict  = TSMD.node_propdict(gcopy,:nodeType);
    sysdict      = TSMD.node_propdict(gcopy,:sys);

    for i = 1:nv_g
    
        if isreal_dict[i] == :fake
            continue
        end
        parent_ = format_name(string(nparent_dict[i]))

        if parent_ ∈ valid_s
            # parent index in dd
            dd_sys_idx = bops_dict[parent_]
            
            component_names = [c.name for c in bopsys[dd_sys_idx].component[:]]
            comp        = sysdict[i];
            compname = format_name(nname_dict[i]) != parent_ ?  format_name(nname_dict[i]) : titlecase(string(nname_dict[i]))
            
            # if it is not already within the system
            if !(compname ∈ component_names)
                resize!(bopsys[dd_sys_idx].component,length(bopsys[dd_sys_idx].component)+1);
                bopsys[dd_sys_idx].component[end].name = compname
            end

            component_names = [c.name for c in bopsys[dd_sys_idx].component[:]]
    
            # index of component in parent system
            idx     = findfirst(x -> x == compname, component_names)
            bopcomp =  bopsys[dd_sys_idx].component[idx];
        
            # comp        = sysdict[i];
            # compname = string(nname_dict[i]);
            # bopsys[dd_sys_idx].component[end].name = format_name(nname_dict[i]) != parent_ ?  format_name(nname_dict[i]) : titlecase(string(nname_dict[i]))

            pps         = propertynames(comp);
            hasnext     = [hasproperty(getproperty(comp,p),:ṁ) for p in pps];    # has a fluid port
            toadd       = pps[findall(x -> x == true, hasnext)];                 # fluid port names
            # bopcomp     = bopsys[dd_sys_idx].component[end];
            flow2add    = length(toadd);

            # if length(bopcomp.port) < flow2add
            #     empty!(bopcomp.port)
            #     resize!(bopcomp.port,flow2add)
            # end

            for j =1:flow2add
                if length(bopcomp.port) < j
                    resize!(bopcomp.port,flow2add)
                    bopcomp.port[j].name = string(toadd[j])
                end

                sysp = getproperty(comp, toadd[j])
                sysT = soln(getproperty(sysp, :T))
                sysP = soln(getproperty(sysp, :P))
                sysm = soln(getproperty(sysp, :ṁ))
                @ddtime(bopcomp.port[j].temperature = sysT-273.15);
                @ddtime(bopcomp.port[j].pressure = sysP);
                @ddtime(bopcomp.port[j].massflow = sysm);
                # @ddtime(bopcomp.port[j].pressure = soln(getproperty(getproperty(comp,toadd[j]),P)))
            end

            if hasproperty(comp,:w) 
                if length(bopcomp.port) < flow2add+1
                    resize!(bopcomp.port,flow2add+1) 
                    bopcomp.port[end].name = "Pdv_Conserving"
                end
                @ddtime(bopcomp.port[end].mechanicalPower = soln(getproperty(getproperty(comp,:w),:Ẇ)));
                # println("PORT $(compname) . w")
            elseif hasproperty(comp,:Ẇ)
                if length(bopcomp.port) < flow2add+1
                    resize!(bopcomp.port,flow2add+1) 
                    bopcomp.port[end].name = "Pdv_Conserving"
                end
                @ddtime(bopcomp.port[end].mechanicalPower = soln(getproperty(comp,:Ẇ)));
            end

            
            if hasproperty(comp,:q) 
                if length(bopcomp.port) < flow2add+1
                    resize!(bopcomp.port,length(bopcomp.port)+1) 
                    bopcomp.port[end].name = "Tds_Conserving"
                end
                # println("-- $(compname) --")
                @ddtime(bopcomp.port[end].thermalPower = soln(getproperty(getproperty(comp,:q),:Q̇)));
                # println("--> @ddtime set to: $(soln(getproperty(getproperty(comp,:q),:Q̇))))")
                # println("--> @ddtime post check: $(@ddtime(bopcomp.port[end].thermalPower))")
                # println("--> @ddtime pre-check $(@ddtime(bopcomp.port[end].thermalPower))")
            elseif hasproperty(comp,:Q̇)
                if length(bopcomp.port) < flow2add+1
                    resize!(bopcomp.port,length(bopcomp.port)+1) 
                    bopcomp.port[end].name = "Tds_Conserving"
                end
                @ddtime(bopcomp.port[end].thermalPower = soln(getproperty(comp,:Q̇)))
                # println("--> $(compname) Q = $(soln(getproperty(comp,:Q̇)))")
                # @ddtime(bopcomp.port[end].mechanicalPower = soln(getproperty(comp,:Ẇ)/10^6));
            end
        #     # hasproperty(comp,:q) ? resize!(bopcomp.port,flowadd+1)
        else
            display("Failed to find parent system for $(nname_dict[i]), invalid parent $(parent_)")
        end
    end

    @ddtime(bop_plant.power_electric_generated = soln(act.odedict[:Electric].Ẇ))
    @ddtime(bop_plant.total_heat_rejected =  soln(act.odedict[:ColdUtility].Q̇))
    @ddtime(bop_plant.total_heat_supplied = -soln(act.odedict[:HotUtility].Q̇))
    @ddtime(bop.thermal_efficiency_plant = soln(:η_bop))
    @ddtime(bop.thermal_efficiency_cycle = soln(:η_cycle))
end

"""
    setxATP!(x,actorATP::ActorThermalPlant)

Setter for actor.x
"""
function setxATP!(x,actorATP::ActorThermalPlant)

    # @assert length(x)==length(actorATP.x) "x incorrect length for actor.x"

    for i =1:length(x)
        actorATP.x[i] = x[i]
        actorATP.var2val[actorATP.sym2var[actorATP.optpar[i]]] = x[i]
    end
    return actorATP
end