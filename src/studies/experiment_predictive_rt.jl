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

    savedir = abspath(joinpath(sty.save_folder, "$(device)_$(shot)__$(Dates.now())__$(getpid())"))
    @info savedir
    if !isdir(savedir)
        mkdir(savedir)
    end
    SimulationParameters.par2json(sty, joinpath(savedir, "sty.json"))

    io_out = sty.redirect_output ? joinpath(savedir, "log.txt") : nothing
    io_err = sty.redirect_output ? joinpath(savedir, "log.txt") : nothing
    redirect_stdio(stdout=io_out, stderr=io_err) do
        cd(savedir) do
            run_predictive_rt_case(device, shot; savedir, sty.reconstruction, sty.verbose, kw_case_parameters)
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
    dd::IMAS.dd,
    dd_exp::IMAS.dd,
    device::Symbol,
    shot::Int;
    savedir::AbstractString=abspath("."),
    save_gif::Bool=false,
    reconstruction::Bool,
    verbose::Bool,
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

    # init
    @info "ini.time.simulation_start = $(ini.time.simulation_start)"
    FUSE.init!(dd, ini, act)

    # keep aside the dd with experimental data
    IMAS.fill!(dd_exp, dd)

    # identify LH transitions
    @info "LH_analysis"
    experiment_LH = FUSE.LH_analysis(dd; do_plot=false)

    act.ActorPedestal.model = :dynamic
    act.ActorPedestal.tau_n = experiment_LH.tau_n
    act.ActorPedestal.tau_t = experiment_LH.tau_t
    act.ActorEPED.ped_factor = 0.8
    act.ActorPedestal.T_ratio_pedestal = 1.0 # Ti/Te in the pedestal
    act.ActorWPED.ped_to_core_fraction = missing

    # density and Zeff from experiment
    act.ActorPedestal.density_ratio_L_over_H = 1.0
    act.ActorPedestal.zeff_ratio_L_over_H = 1.0

    # LH-transition at user-defined times
    act.ActorPedestal.mode_transitions = experiment_LH.mode_transitions

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
    act.ActorFluxMatcher.z_max = (core=20.0, edge=100.0, rho_transition=0.80) 


    # FINN transport — replaces FluxMatcher + TGLF
    act.ActorCoreTransport.model = :FINN
    act.ActorFINN.finn_model = "finn_sat3_d3d_withnegD"
    act.ActorFINN.rho_transport = 0.1:0.025:0.85
    act.ActorFINN.warn_nn_train_bounds = false

    act.ActorPedestal.rotation_model = :replay

    act.ActorSawteethSource.flat_factor = 1.0
    act.ActorSawteethSource.period = 0.25 # turn off flattening after 0.25s of no sawteeth events

    # time
    δt = 0.05
    dd.global_time = ini.time.simulation_start
    final_time = ini.general.dd.equilibrium.time[end]
    act.ActorDynamicPlasma.Nt = Int(ceil((final_time - dd.global_time) / δt))
    act.ActorDynamicPlasma.Δt = final_time - dd.global_time

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

    Nt_OK = round(Int, (dd.global_time - ini.time.simulation_start) / δt)
    times = range(ini.time.simulation_start, dd.global_time, Nt_OK)

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

        @info "save dd_sim.json"
        IMAS.imas2json(dd, joinpath(savedir, "dd_sim.json"))

        @info "save dd_exp.json"
        IMAS.imas2json(dd_exp, joinpath(savedir, "dd_exp.json"))

        @info "save dd_benchmark.json"
        IMAS.imas2json(bnch.dd, joinpath(savedir, "dd_benchmark.json"))

        if save_gif
            @info "save animated gif"
            mkdir(joinpath(savedir, "gif"))
            animated_plasma_overview(dd, joinpath(savedir, "gif"), dd_exp)
        end
    end

    return nothing
end
