import JSON

#= ================== =#
#  StudyPredictiveRT  #
#= ================== =#

"""
    study_parameters(::Val{:PredictiveRT})

Runs predictive simulations using FINN transport model for specified device and shots
"""
function study_parameters(::Val{:PredictiveRT})
    return FUSEparameters__ParametersStudyPredictiveRT{Real}()
end

Base.@kwdef mutable struct FUSEparameters__ParametersStudyPredictiveRT{T<:Real} <: ParametersStudy{T}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :StudyPredictiveRT
    server::Switch{String} = study_common_parameters(; server="localhost")
    n_workers::Entry{Int} = study_common_parameters(; n_workers=missing)
    release_workers_after_run::Entry{Bool} = study_common_parameters(; release_workers_after_run=true)
    save_folder::Entry{String} = Entry{String}("-", "Folder to save the predictive runs into")
    kw_case_parameters::Entry{Dict{Symbol,Any}} = Entry{Dict{Symbol,Any}}("-", "Keyword arguments passed to case_parameters"; default=Dict{Symbol,Any}())
    redirect_output::Entry{Bool} = Entry{Bool}("-", "Redirect stdout and stderr to log.txt file"; default=true)
    verbose::Entry{Bool} = Entry{Bool}("-", "Turn on verbose progress output"; default=true)

    # Predictive-specific parameters
    device::Entry{Symbol} = Entry{Symbol}("-", "Device to run predictive simulations for")
    shots::Entry{Vector{Int}} = Entry{Vector{Int}}("-", "List of shot numbers")
    reconstruction::Entry{Bool} = Entry{Bool}("-", "Run prediction in reconstruction mode")
    start_time::Entry{Float64} = Entry{Float64}("-", "Override simulation start time [s]"; default=missing)
    end_time::Entry{Float64} = Entry{Float64}("-", "Override simulation end time [s]"; default=missing)

    # Pedestal predictor (run the study once per model to compare them, see `compare_pedestal_nn_runs`)
    nn_model::Switch{Symbol} = Switch{Symbol}([:fuse29, :fuse_efit], "-", "Pedestal predictor driving the edge: `:fuse29` (coil currents) or `:fuse_efit` (Ip, boundary shape and loop voltage)"; default=:fuse29)
    nn_state_source::Switch{Symbol} = Switch{Symbol}([:simulation, :replay], "-", "fuse_efit only: take Ip, boundary shape and loop voltage from FUSE's own equilibrium (`:simulation`) or from the experimental one (`:replay`)"; default=:simulation)
    nn_efit_form::Switch{Symbol} = Switch{Symbol}([:single, :averaged], "-", "fuse_efit only: one network (`:single`) or the average of eight (`:averaged`)"; default=:single)
end

mutable struct StudyPredictiveRT{T<:Real} <: AbstractStudy
    sty::OverrideParameters{T,FUSEparameters__ParametersStudyPredictiveRT{T}}
    act::ParametersActors
end

function StudyPredictiveRT(sty::ParametersStudy; kw...)
    sty = OverrideParameters(sty; kw...)
    study = StudyPredictiveRT(sty, ParametersActors())
    parallel_environment(sty.server, sty.n_workers)
    return study
end

"""
    _run(study::StudyPredictiveRT)

Runs the PredictiveRT study with sty settings in parallel on designated cluster
"""
function _run(study::StudyPredictiveRT)
    sty = study.sty

    @assert (sty.n_workers == 0 || sty.n_workers == length(Distributed.workers())) "The number of workers = $(length(Distributed.workers())) isn't the number of workers you requested = $(sty.n_workers)"

    # parallel run
    println("running $(length(sty.shots)) predictive RT simulations with $(sty.n_workers) workers on $(sty.server)")

    ProgressMeter.@showprogress map(shot -> run_predictive_rt_case(study, shot; sty.kw_case_parameters), sty.shots)

    # Release workers after run
    if sty.release_workers_after_run
        Distributed.rmprocs(Distributed.workers())
        @info "released workers"
    end

    return study
end

"""
    run_predictive_rt_case(study::StudyPredictiveRT, shot::Int; kw_case_parameters::Dict{Symbol,Any})

Run a single predictive case using FINN transport for a given device and shot
"""
function run_predictive_rt_case(study::StudyPredictiveRT, shot::Int; kw_case_parameters::Dict{Symbol,Any})
    sty = study.sty
    device = sty.device

    savedir = abspath(joinpath(sty.save_folder, "$(device)_$(shot)__$(Dates.now())__$(getpid())__$(sty.nn_model)"))
    @info savedir
    if !isdir(savedir)
        mkdir(savedir)
    end
    SimulationParameters.par2json(sty, joinpath(savedir, "sty.json"))

    io_out = sty.redirect_output ? joinpath(savedir, "log.txt") : nothing
    io_err = sty.redirect_output ? joinpath(savedir, "log.txt") : nothing
    redirect_stdio(stdout=io_out, stderr=io_err) do
        cd(savedir) do
            run_predictive_rt_case(device, shot; savedir, sty.reconstruction, sty.verbose, sty.start_time, sty.end_time,
                sty.nn_model, sty.nn_state_source, sty.nn_efit_form, kw_case_parameters)
        end
    end
end

function run_predictive_rt_case(device::Symbol, shot::Int; kw_case_parameters::Dict{Symbol,Any}, kw...)
    dd = IMAS.dd()
    dd_exp = IMAS.dd()
    run_predictive_rt_case!(dd, dd_exp, device, shot; kw_case_parameters, kw...)
    return (dd=dd, dd_exp=dd_exp)
end

function run_predictive_rt_case!(
    dd::IMAS.DD,
    dd_exp::IMAS.DD,
    device::Symbol,
    shot::Int;
    savedir::AbstractString=abspath("."),
    save_gif::Bool=false,
    reconstruction::Bool,
    verbose::Bool,
    start_time::Union{Float64,Missing}=missing,
    end_time::Union{Float64,Missing}=missing,
    nn_model::Symbol=:fuse29,
    nn_state_source::Symbol=:simulation,
    nn_efit_form::Symbol=:single,
    kw_case_parameters::Dict{Symbol,Any}
)

    isterminal = isa(stdout, Base.TTY) && (get(ENV, "CI", nothing) != "true")
    if !isterminal && verbose
        @info "StudyPredictiveRT: verbose output requested but stdout is not a terminal; setting verbose=false"
        verbose = false
    end

    # Get case parameters
    @info "StudyPredictiveRT: case_parameters($(repr(device)), $shot; $(repr(kw_case_parameters))...)"
    ini, act = FUSE.case_parameters(device, shot; kw_case_parameters...)

    # Override start/end time if specified
    if !ismissing(start_time)
        ini.time.simulation_start = start_time
    end

    ini.core_profiles.ne_setting    = :ne_ped
    act.ActorPedestal.density_match = :ne_ped
    act.ActorPedestal.nn_warmup_time = 0.0 # replay actuators from shot start so the recurrent state is warm at t_start

    # init
    @info "ini.time.simulation_start = $(ini.time.simulation_start)"
    FUSE.init!(dd, ini, act)

    # keep aside the dd with experimental data
    IMAS.fill!(dd_exp, dd)
    # identify LH transitions
    @info "LH_analysis"
    #experiment_LH = FUSE.LH_analysis(dd; do_plot=false)

    # LH-transition at user-defined times
    #act.ActorPedestal.mode_transitions = experiment_LH.mode_transitions

    act.ActorPedestal.model = :dynamic
    act.ActorPedestal.tau_n = 0.3 #experiment_LH.tau_n # 0.3
    act.ActorPedestal.tau_t = 0.15 #experiment_LH.tau_t # 0.15
    act.ActorEPED.ped_factor = 1.
    act.ActorPedestal.T_ratio_pedestal = 1.0 # Ti/Te in the pedestal 
    act.ActorWPED.ped_to_core_fraction = missing

    # ZMQ coupling to GSLite/GSEvolve
    act.ActorZMQ.enabled = false # true
    act.ActorPedestal.fpe_source = :dd #:zmq # dd
    act.ActorPedestal.ne_from = :nn_predictor  #:pulse_schedule 
    act.ActorPedestal.nn_model = nn_model # :fuse29 or :fuse_efit
    act.ActorPedestal.nn_state_source = nn_state_source # fuse_efit: Ip, shape and loop voltage from :simulation or :replay
    act.ActorPedestal.nn_efit_form = nn_efit_form # fuse_efit: :single or :averaged
    act.ActorPedestal.nn_ped_quantities = :all # :ne_lh (ne_ped + L/H only, Te/Ti from EPED/WPED) or :all (also apply NN te_ped, ti_ped, rotation)

    # density and Zeff from experiment
    act.ActorPedestal.density_ratio_L_over_H = 1.0
    act.ActorPedestal.zeff_ratio_L_over_H = 1.0

    act.ActorEquilibrium.model = :FRESCO
    act.ActorFRESCO.nR = 65
    act.ActorFRESCO.nZ = 65

    act.ActorNeutralFueling.τp_over_τe = 0.25

    act.ActorFluxMatcher.evolve_plasma_sources = true
    act.ActorFluxMatcher.algorithm = :default
    act.ActorFluxMatcher.max_iterations = -100 # negative to avoid print of warnings
    act.ActorFluxMatcher.evolve_pedestal = false
    act.ActorFluxMatcher.evolve_Te = :flux_match
    act.ActorFluxMatcher.evolve_Ti = :flux_match
    act.ActorFluxMatcher.evolve_densities = :flux_match
    act.ActorFluxMatcher.evolve_rotation = :replay
    act.ActorFluxMatcher.relax = 1.0
    act.ActorFINN.evolve_rotation = :replay


    # FINN transport — replaces FluxMatcher + TGLF
    #act.ActorCoreTransport.model = :FINN
    #act.ActorFINN.finn_model = "finn_sat3_d3d_withnegD"
    #act.ActorFINN.rho_transport = 0.1:0.025:0.85
    #act.ActorFINN.warn_nn_train_bounds = false

    act.ActorTGLF.model = :TGLFNN
    act.ActorTGLF.tglfnn_model = "sat3_em_d3d_azf-1_withnegD"

    act.ActorPedestal.rotation_model = :replay

    act.ActorSawteethSource.flat_factor = 1.0
    act.ActorSawteethSource.period = 0.25 # turn off flattening after 0.25s of no sawteeth events

    # time
    δt = 0.05
    dd.global_time = ini.time.simulation_start
    final_time = ismissing(end_time) ? ini.general.dd.equilibrium.time[end] : end_time
    act.ActorDynamicPlasma.Nt = Int(ceil((final_time - dd.global_time) / δt))
    act.ActorDynamicPlasma.Δt = final_time - dd.global_time
    # If GSLite does not provide start/end time, uncomment below and let GSLite
    # control the simulation duration via done signal or timeout.
    # Δt = 20.0
    # act.ActorDynamicPlasma.Nt = Int(ceil(Δt / δt))
    # act.ActorDynamicPlasma.Δt = Δt

    # choose what to evolve
    act.ActorDynamicPlasma.evolve_current = true
    act.ActorDynamicPlasma.evolve_equilibrium = true
    act.ActorDynamicPlasma.evolve_transport = true
    act.ActorDynamicPlasma.evolve_sources = true
    act.ActorDynamicPlasma.evolve_pf_active = false
    act.ActorDynamicPlasma.evolve_pedestal = true
    act.ActorDynamicPlasma.evolve_sawteeth = true

    if reconstruction
        act.ActorCoreTransport.model = :replay
        act.ActorPedestal.model = :replay

        act.ActorPFactive.boundary_weight = 1.0
        act.ActorPFactive.magnetic_probe_weight = 1.0
        act.ActorPFactive.flux_loop_weight = 1.0
        act.ActorPFactive.strike_points_weight = 0.0
        act.ActorPFactive.x_points_weight = 1.0
    else
        act.ActorPFactive.boundary_weight = 1.0
        act.ActorPFactive.magnetic_probe_weight = 0.0
        act.ActorPFactive.flux_loop_weight = 0.0
        act.ActorPFactive.strike_points_weight = 0.0
        act.ActorPFactive.x_points_weight = 1.0
    end

    # Run the simulation
    try
        @info "StudyPredictiveRT: ActorDynamicPlasma(dd, act)"
        FUSE.ActorDynamicPlasma(dd, act; verbose)
    catch e
        if isa(e, InterruptException)
            rethrow(e)
        else
            @error repr(e)
        end
    end

    Nt_OK = max(2, round(Int, (dd.global_time - ini.time.simulation_start) / δt))
    times = range(ini.time.simulation_start, dd.global_time, Nt_OK)

    # Validation residuals
    @info "StudyPredictiveRT: computing validation residuals"
    validation = compute_validation_residuals(dd, dd_exp)

    @info "IMAS.benchmark(dd, dd_exp, dd.core_profiles.time);"
    bnch = IMAS.benchmark(dd, dd_exp, dd.core_profiles.time);

    # save simulation data to directory
    if !isempty(savedir)
        @info "saving simulation results to: $(savedir)"
        @info "save act.json"
        tmp = act.ActorReplay.replay_dd
        act.ActorReplay.replay_dd = IMAS.dd()
        SimulationParameters.par2json(act, joinpath(savedir, "act.json"))
        act.ActorReplay.replay_dd = tmp

        #@info "StudyPostdictive: save dd_sim.json"
        IMAS.imas2json(dd, joinpath(savedir, "dd_sim.json"))

        #@info "StudyPostdictive: save dd_exp.json"
        IMAS.imas2json(dd_exp, joinpath(savedir, "dd_exp.json"))

        @info "save dd_benchmark.json"
        IMAS.imas2json(bnch.dd, joinpath(savedir, "dd_benchmark.json"))

        @info "save traces.json"
        traces = extract_traces(dd, dd_exp)
        open(joinpath(savedir, "traces.json"), "w") do f
            JSON.print(f, traces, 2)
        end

        @info "save validation_residuals.json"
        open(joinpath(savedir, "validation_residuals.json"), "w") do f
            JSON.print(f, validation, 2)
        end

        pedestal_nn = pedestal_nn_validation(dd, dd_exp)
        if pedestal_nn !== nothing
            @info "save pedestal_nn.json"
            open(joinpath(savedir, "pedestal_nn.json"), "w") do f
                JSON.print(f, pedestal_nn, 2)
            end
        end

        if save_gif
            @info "save animated gif"
            mkdir(joinpath(savedir, "gif"))
            animated_plasma_overview(dd, joinpath(savedir, "gif"), dd_exp)
        end
    end

    return nothing
end

"""
    extract_traces(dd::IMAS.DD, dd_exp::IMAS.DD)

Extract time traces of li and beta_pol from FUSE simulation and EFIT,
both on FUSE's time grid (EFIT is interpolated to match).
"""
function extract_traces(dd::IMAS.DD, dd_exp::IMAS.DD)
    times_sim = dd.equilibrium.time
    times_exp = dd_exp.equilibrium.time

    li_sim = [dd.equilibrium.time_slice[t].global_quantities.li_1 for t in times_sim]
    bp_sim = [dd.equilibrium.time_slice[t].global_quantities.beta_pol for t in times_sim]

    li_exp = IMAS.interp1d(times_exp, [dd_exp.equilibrium.time_slice[t].global_quantities.li_1 for t in times_exp]).(times_sim)
    bp_exp = IMAS.interp1d(times_exp, [dd_exp.equilibrium.time_slice[t].global_quantities.beta_pol for t in times_exp]).(times_sim)

    return Dict(
        "time"   => times_sim,
        "li_sim" => li_sim,
        "bp_sim" => bp_sim,
        "li_exp" => li_exp,
        "bp_exp" => bp_exp,
    )
end

function _median(x::AbstractVector)
    s = sort(x)
    n = length(s)
    return isodd(n) ? s[(n+1)÷2] : (s[n÷2] + s[n÷2+1]) / 2
end

"""
    compute_validation_residuals(dd::IMAS.DD, dd_exp::IMAS.DD)

Compute validation residuals comparing time derivatives of FUSE vs EFIT.
Residual = dq_FUSE/dt - dq_EFIT/dt at each consecutive FUSE time step.
EFIT quantities are interpolated onto FUSE's time grid before differencing.
Returns max_abs and median for: dli_dt, dbetap_dt. Rp is TODO.
"""
function compute_validation_residuals(dd::IMAS.DD, dd_exp::IMAS.DD)
    times_sim = dd.equilibrium.time
    times_exp = dd_exp.equilibrium.time

    if length(times_sim) < 2
        @warn "compute_validation_residuals: not enough time slices"
        return Dict{String,Any}()
    end

    # FUSE traces on FUSE time grid
    li_sim = [dd.equilibrium.time_slice[t].global_quantities.li_1 for t in times_sim]
    bp_sim = [dd.equilibrium.time_slice[t].global_quantities.beta_pol for t in times_sim]

    # EFIT traces interpolated onto FUSE time grid
    li_exp = IMAS.interp1d(times_exp, [dd_exp.equilibrium.time_slice[t].global_quantities.li_1 for t in times_exp]).(times_sim)
    bp_exp = IMAS.interp1d(times_exp, [dd_exp.equilibrium.time_slice[t].global_quantities.beta_pol for t in times_exp]).(times_sim)

    # dq/dt residuals using consecutive time steps
    res_li = Float64[]
    res_bp = Float64[]
    for k in 2:length(times_sim)
        dt = times_sim[k] - times_sim[k-1]
        push!(res_li, (li_sim[k] - li_sim[k-1]) / dt - (li_exp[k] - li_exp[k-1]) / dt)
        push!(res_bp, (bp_sim[k] - bp_sim[k-1]) / dt - (bp_exp[k] - bp_exp[k-1]) / dt)
    end
    # TODO: Rp = d(ψ_plasma)/dt / Ip

    summary = Dict{String,Any}()
    for (name, res) in [("dli_dt", res_li), ("dbetap_dt", res_bp)]
        if !isempty(res)
            summary[name] = Dict("max_abs" => maximum(abs.(res)), "median" => _median(res))
        end
    end

    return summary
end


#= ============================================== =#
#  Pedestal predictor validation (fuse29/fuse_efit) #
#= ============================================== =#

# value of a core_profiles quantity at rho_tor_norm = rho; `nothing` when the slice does not carry it
function _cp1d_at_rho(cp1d, rho::Float64, getter::Function)
    try
        y = getter(cp1d)
        v = IMAS.interp1d(cp1d.grid.rho_tor_norm, y).(rho)
        return isfinite(v) ? Float64(v) : nothing
    catch e
        e isa InterruptException && rethrow()
        return nothing
    end
end

"""
    pedestal_nn_validation(dd::IMAS.DD, dd_exp::IMAS.DD; rho=0.85, max_dt=0.05)

Everything needed to judge a pedestal predictor on one run and to compare two
runs of the same shot made with different models. `nothing` when no predictor
ran. Returns a JSON-ready dictionary:

- the per-tick model record of [`pedestal_nn_history`](@ref) (`model`, `time`,
  `inputs`, `indicators`, `predictions`, `is_h_mode`, `n_ticks`);
- `experiment`: `ne` [1e19 m⁻³], `te_ped`, `ti_ped` [keV] and `t_rot_ped`
  [krad/s] of the experimental profiles at `rho`, on the model's time base
  (`nothing` where no experimental profile lies within `max_dt` seconds);
- `rmse` and `n`: root-mean-square error of each of those four heads against
  the experiment, and the number of ticks it is computed on.

The model heads are the predictions the model made in the run; with
`nn_state_source = :simulation` they depend on FUSE's own equilibrium, which is
the point of comparing against `:replay` and against fuse29.
"""
function pedestal_nn_validation(dd::IMAS.DD, dd_exp::IMAS.DD; rho::Float64=0.85, max_dt::Float64=0.05)
    out = pedestal_nn_history(dd)
    out === nothing && return nothing
    times = out["time"]

    getters = Dict{String,Function}(
        "ne" => cp1d -> cp1d.electrons.density_thermal ./ 1e19,
        "te_ped" => cp1d -> cp1d.electrons.temperature ./ 1e3,
        "ti_ped" => cp1d -> cp1d.t_i_average ./ 1e3,
        "t_rot_ped" => cp1d -> cp1d.ion[1].rotation_frequency_tor ./ 1e3)
    exp_times = [cp1d.time for cp1d in dd_exp.core_profiles.profiles_1d]
    experiment = Dict{String,Any}(name => Any[] for name in keys(getters))
    for t in times
        k = isempty(exp_times) ? nothing : argmin(abs.(exp_times .- t))
        for (name, getter) in getters
            v = (k === nothing || abs(exp_times[k] - t) > max_dt) ? nothing : _cp1d_at_rho(dd_exp.core_profiles.profiles_1d[k], rho, getter)
            push!(experiment[name], v)
        end
    end
    out["rho"] = rho
    out["experiment"] = experiment

    rmse = Dict{String,Any}()
    count = Dict{String,Any}()
    for name in keys(getters)
        pred = out["predictions"][name]
        err2 = Float64[(p - e)^2 for (p, e) in zip(pred, experiment[name]) if p !== nothing && e !== nothing]
        count[name] = length(err2)
        rmse[name] = isempty(err2) ? nothing : sqrt(sum(err2) / length(err2))
    end
    out["rmse"] = rmse
    out["n"] = count
    return out
end

"""
    compare_pedestal_nn_runs(dir_a::AbstractString, dir_b::AbstractString; verbose=true)

Side-by-side comparison of two `StudyPredictiveRT` runs of the same shot made
with different pedestal predictors (typically `nn_model = :fuse29` and
`nn_model = :fuse_efit`), from the `pedestal_nn.json` each run saved.

Returns a dictionary with, per run, the model name and the RMSE of `ne`,
`te_ped`, `ti_ped`, `t_rot_ped` at rho = 0.85 against the experiment — on each
run's own ticks (`rmse`) and on the ticks the two runs share (`rmse_common`) —
plus `hmode_agreement`, the fraction of shared ticks on which the two models'
H-mode gates agree, and `n_common`.
"""
function compare_pedestal_nn_runs(dir_a::AbstractString, dir_b::AbstractString; verbose::Bool=true)
    function load(dir)
        path = isdir(dir) ? joinpath(dir, "pedestal_nn.json") : dir
        isfile(path) || error("compare_pedestal_nn_runs: $path not found (was the run made with ne_from=:nn_predictor?)")
        return JSON.parsefile(path)
    end
    a, b = load(dir_a), load(dir_b)
    heads = ("ne", "te_ped", "ti_ped", "t_rot_ped")

    # ticks present in both runs (same model time within 1 ms)
    pairs = Tuple{Int,Int}[]
    for (i, t) in enumerate(a["time"])
        j = findfirst(τ -> abs(τ - t) < 1e-3, b["time"])
        j === nothing || push!(pairs, (i, j))
    end

    function rmse_on(run, idx, head)
        err2 = Float64[]
        for i in idx
            p, e = run["predictions"][head][i], run["experiment"][head][i]
            (p === nothing || e === nothing) || push!(err2, (p - e)^2)
        end
        return isempty(err2) ? nothing : sqrt(sum(err2) / length(err2))
    end

    result = Dict{String,Any}("n_common" => length(pairs))
    for (key, run, idx) in (("a", a, first.(pairs)), ("b", b, last.(pairs)))
        result[key] = Dict{String,Any}(
            "model" => run["model"],
            "rmse" => Dict{String,Any}(h => run["rmse"][h] for h in heads),
            "rmse_common" => Dict{String,Any}(h => rmse_on(run, idx, h) for h in heads))
    end
    result["hmode_agreement"] = isempty(pairs) ? nothing :
                                count(a["is_h_mode"][i] == b["is_h_mode"][j] for (i, j) in pairs) / length(pairs)

    if verbose
        fmt(v) = v === nothing ? "      n/a" : @sprintf("%9.4f", v)
        units = Dict("ne" => "1e19 m^-3", "te_ped" => "keV", "ti_ped" => "keV", "t_rot_ped" => "krad/s")
        println("RMSE against the experiment at rho = $(a["rho"]) on the $(length(pairs)) ticks both runs share")
        @printf("%-10s %-10s %12s %12s\n", "head", "units", result["a"]["model"], result["b"]["model"])
        for h in heads
            println(@sprintf("%-10s %-10s", h, units[h]), "    ", fmt(result["a"]["rmse_common"][h]), "    ", fmt(result["b"]["rmse_common"][h]))
        end
        result["hmode_agreement"] === nothing || @printf("H-mode gates agree on %.1f %% of shared ticks\n", 100 * result["hmode_agreement"])
    end
    return result
end

"""
    plot_pedestal_nn_comparison(dir_a::AbstractString, dir_b::AbstractString; savefig=nothing, title="")

Figure for two `StudyPredictiveRT` runs of the same shot made with different
pedestal predictors, from the `pedestal_nn.json` each run saved: one panel per
head that has an experimental counterpart (`ne`, `te_ped`, `ti_ped`,
`t_rot_ped` at rho = 0.85; the experiment as points, each model as a line, its
RMSE in the panel title), the H-mode probability of each model, and the state of
each model's H-mode gate.

`savefig` is a file name (`.png`, `.pdf`, ...) to write the figure to. Returns
the plot.
"""
function plot_pedestal_nn_comparison(dir_a::AbstractString, dir_b::AbstractString;
                                     savefig::Union{Nothing,AbstractString}=nothing, title::AbstractString="")
    function load(dir)
        path = isdir(dir) ? joinpath(dir, "pedestal_nn.json") : dir
        isfile(path) || error("plot_pedestal_nn_comparison: $path not found (was the run made with ne_from=:nn_predictor?)")
        return JSON.parsefile(path)
    end
    runs = (load(dir_a), load(dir_b))
    # a model keeps its colour whichever argument it is passed as
    model_colors = Dict("fuse29" => "#2a78d6", "fuse_efit" => "#eb6834")
    colors = [get(model_colors, run["model"], "") for run in runs]
    (any(isempty, colors) || colors[1] == colors[2]) && (colors = ["#2a78d6", "#eb6834"])
    styles = (:solid, :dash)   # the second run dashed, so that coinciding traces both show
    ink = "#3d3d3a"
    num(v) = v === nothing ? NaN : Float64(v)
    reference = runs[argmax([length(run["time"]) for run in runs])]
    rho = reference["rho"]
    names = [run["model"] for run in runs]
    names[1] == names[2] && (names = ["$(names[1]) (first run)", "$(names[2]) (second run)"])

    common = (; xlabel="time [s]", gridalpha=0.15, legend=false, titlefontsize=9, guidefontsize=9, tickfontsize=8, titlelocation=:left)

    # one legend for the whole figure, in a strip of its own above the panels
    strip = Plots.plot(; framestyle=:none, legend=:top, legend_columns=-1, foreground_color_legend=nothing, background_color_legend=nothing,
                       legendfontsize=9, xlims=(0, 1), ylims=(0, 1))
    Plots.scatter!(strip, [NaN], [NaN]; color=ink, markersize=4, markerstrokewidth=0, label="experiment")
    for (name, color, linestyle) in zip(names, colors, styles)
        Plots.plot!(strip, [NaN], [NaN]; color, linestyle, linewidth=2, label=name)
    end

    panels = (("ne", "ne [1e19 m^-3]"), ("te_ped", "Te [keV]"), ("ti_ped", "Ti [keV]"), ("t_rot_ped", "rotation [krad/s]"))
    plots = []
    for (head, ylabel) in panels
        rmse = join([@sprintf("%s %.3g", name, run["rmse"][head]) for (name, run) in zip(names, runs) if run["rmse"][head] !== nothing], ", ")
        p = Plots.plot(; ylabel, title="$head at rho = $rho" * (isempty(rmse) ? "" : "    RMSE: $rmse"), common...)
        Plots.scatter!(p, reference["time"], num.(reference["experiment"][head]); color=ink, markersize=3, markerstrokewidth=0)
        for (run, color, linestyle) in zip(runs, colors, styles)
            Plots.plot!(p, run["time"], num.(run["predictions"][head]); color, linestyle, linewidth=2)
        end
        push!(plots, p)
    end

    p = Plots.plot(; ylabel="probability", title="H-mode probability", ylims=(-0.03, 1.03), common...)
    for (run, color, linestyle) in zip(runs, colors, styles)
        Plots.plot!(p, run["time"], num.(run["predictions"]["hmode"]); color, linestyle, linewidth=2)
    end
    push!(plots, p)

    p = Plots.plot(; ylabel="gate", title="H-mode gate (0 = L-mode, 1 = H-mode)", ylims=(-0.1, 1.1), yticks=[0, 1], common...)
    for (run, color, linestyle) in zip(runs, colors, styles)
        Plots.plot!(p, run["time"], Float64.(run["is_h_mode"]); seriestype=:steppost, color, linestyle, linewidth=2)
    end
    push!(plots, p)

    fig = Plots.plot(strip, plots...; layout=Plots.@layout([a{0.04h}; Plots.grid(3, 2)]), size=(1100, 980), plot_title=title, plot_titlefontsize=11,
                     left_margin=6 * Plots.mm, bottom_margin=3 * Plots.mm, background_color="#fcfcfb")
    savefig === nothing || Plots.savefig(fig, savefig)
    return fig
end
