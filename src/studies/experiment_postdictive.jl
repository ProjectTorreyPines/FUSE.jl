#= ================ =#
#  StudyPostdictive  #
#= ================ =#

"""
    study_parameters(::Val{:Postdictive})

Runs postdictive simulations for specified device and shots with configurable parameters
"""
function study_parameters(::Val{:Postdictive})
    return FUSEparameters__ParametersStudyPostdictive{Real}()
end

Base.@kwdef mutable struct FUSEparameters__ParametersStudyPostdictive{T<:Real} <: ParametersStudy{T}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :StudyPostdictive
    server::Switch{String} = study_common_parameters(; server="localhost")
    n_workers::Entry{Int} = study_common_parameters(; n_workers=missing)
    release_workers_after_run::Entry{Bool} = study_common_parameters(; release_workers_after_run=true)
    save_folder::Entry{String} = Entry{String}("-", "Folder to save the postdictive runs into")

    # Postdictive-specific parameters
    device::Entry{Symbol} = Entry{Symbol}("-", "Device to run postdictive simulations for")
    shots::Entry{Vector{Int}} = Entry{Vector{Int}}("-", "List of shot numbers")
    fit_profiles::Entry{Bool} = Entry{Bool}("-", "Whether to fit profiles in case_parameters"; default=true)
    reconstruction::Entry{Bool} = Entry{Bool}("-", "Run postdiction in reconstruction mode")
    use_local_cache::Entry{Bool} = Entry{Bool}("-", "Whether to use local cache in case_parameters"; default=false)
end

mutable struct StudyPostdictive{T<:Real} <: AbstractStudy
    sty::OverrideParameters{T,FUSEparameters__ParametersStudyPostdictive{T}}
    act::ParametersActors
end

function StudyPostdictive(sty::ParametersStudy; kw...)
    sty = OverrideParameters(sty; kw...)
    study = StudyPostdictive(sty, ParametersActors())
    parallel_environment(sty.server, sty.n_workers)
    return study
end

"""
    _run(study::StudyPostdictive)

Runs the Postdictive study with sty settings in parallel on designated cluster
"""
function _run(study::StudyPostdictive)
    sty = study.sty

    @assert sty.n_workers == length(Distributed.workers()) "The number of workers = $(length(Distributed.workers())) isn't the number of workers you requested = $(sty.n_workers)"

    # parallel run
    println("running $(length(sty.shots)) postdictive simulations with $(sty.n_workers) workers on $(sty.server)")

    ProgressMeter.@showprogress map(shot -> run_postdictive_case(study, shot), sty.shots)

    # Release workers after run
    if sty.release_workers_after_run
        Distributed.rmprocs(Distributed.workers())
        @info "released workers"
    end

    return study
end

"""
    run_postdictive_case(study::StudyPostdictive, shot::Int)

Run a single postdictive case for a given device and shot
"""
function run_postdictive_case(study::StudyPostdictive, shot::Int)
    sty = study.sty
    device = sty.device

    original_dir = pwd()

    # Redirect stdout and stderr to the file
    original_stdout = stdout
    original_stderr = stderr

    savedir = abspath(joinpath(sty.save_folder, "$(device)_$(shot)__$(Dates.now())__$(getpid())"))
    @info savedir
    if !isdir(savedir)
        mkdir(savedir)
    end
    SimulationParameters.par2json(sty, joinpath(savedir, "sty.json"))
    file_log = open(joinpath(savedir, "log.txt"), "w")

    try
        redirect_stdout(file_log)
        redirect_stderr(file_log)
        cd(savedir)

        run_postdictive_case(device, shot; user_act=study.act, sty.fit_profiles, sty.use_local_cache, savedir, sty.reconstruction)

        # catch e
        #     if isa(e, InterruptException)
        #         rethrow(e)
        #     end
        #     @error "Error in postdictive case for $(device) shot $(shot): $e"
    finally
        redirect_stdout(original_stdout)
        redirect_stderr(original_stderr)
        cd(original_dir)
        close(file_log)
    end
end

function run_postdictive_case(device::Symbol, shot::Int; kw...)
    dd = IMAS.dd()
    dd_exp = IMAS.dd()
    run_postdictive_case!(dd, dd_exp, device, shot; kw...)
    return (dd=dd, dd_exp=dd_exp)
end

function run_postdictive_case!(
    dd::IMAS.dd,
    dd_exp::IMAS.dd,
    device::Symbol,
    shot::Int;
    user_act::ParametersActors,
    fit_profiles::Bool,
    use_local_cache::Bool,
    savedir::AbstractString=abspath("."),
    reconstruction::Bool
)

    # Get case parameters
    @info "case_parameters($(repr(device)), $shot; fit_profiles=$fit_profiles, use_local_cache=$use_local_cache)"
    ini, act = FUSE.case_parameters(device, shot; fit_profiles, use_local_cache)

    # Override act with user-specific actor parameters
    #merge!(act, user_act)

    # init
    ini.time.simulation_start = ini.general.dd.equilibrium.time_slice[2].time
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
    act.ActorWPED.ped_to_core_fraction = experiment_LH.W_ped_to_core_fraction
    act.ActorEPED.ped_factor = 1.0
    act.ActorPedestal.T_ratio_pedestal = 1.0 # Ti/Te in the pedestal

    if true
        # density and Zeff from experiment
        act.ActorPedestal.density_ratio_L_over_H = 1.0
        act.ActorPedestal.zeff_ratio_L_over_H = 1.0
    else
        # density can go from L to H mode at a different time
        act.ActorPedestal.density_ratio_L_over_H = experiment_LH.ne_L_over_H
        act.ActorPedestal.zeff_ratio_L_over_H = experiment_LH.zeff_L_over_H
        dd.pulse_schedule.density_control.n_e_line.reference = experiment_LH.ne_H
        dd.pulse_schedule.density_control.zeff_pedestal.reference = experiment_LH.zeff_H
    end

    if true
        # LH-transition at user-defined times
        act.ActorPedestal.mode_transitions = experiment_LH.mode_transitions
    else
        # LH-transition from LH scaling law
        act.ActorPedestal.mode_transitions = missing
    end

    act.ActorEquilibrium.model = :FRESCO
    act.ActorFRESCO.nR = 65
    act.ActorFRESCO.nZ = 65

    act.ActorNeutralFueling.τp_over_τe = 0.25

    act.ActorFluxMatcher.evolve_plasma_sources = false
    act.ActorFluxMatcher.algorithm = :simple
    act.ActorFluxMatcher.max_iterations = -10 # negative to avoid print of warnings
    act.ActorFluxMatcher.evolve_pedestal = false
    act.ActorFluxMatcher.evolve_Te = :flux_match
    act.ActorFluxMatcher.evolve_Ti = :flux_match
    act.ActorFluxMatcher.evolve_densities = :flux_match
    act.ActorFluxMatcher.evolve_rotation = :flux_match
    act.ActorPedestal.rotation_model = :replay

    act.ActorFluxMatcher.relax = 0.5
    act.ActorEGGO.timeslice_average = 4

    #act.ActorTGLF.tglfnn_model = "sat1_em_d3d"

    # time
    δt = 0.025
    dd.global_time = ini.time.simulation_start # start_time should be early in the shot, when otherwise ohmic current will be wrong
    final_time = ini.general.dd.equilibrium.time[end]
    act.ActorDynamicPlasma.Nt = Int(ceil((final_time - dd.global_time) / δt))
    act.ActorDynamicPlasma.Δt = final_time - dd.global_time

    act.ActorDynamicPlasma.evolve_current = true
    act.ActorDynamicPlasma.evolve_equilibrium = true
    act.ActorDynamicPlasma.evolve_transport = true
    act.ActorDynamicPlasma.evolve_hcd = true
    act.ActorDynamicPlasma.evolve_pf_active = false
    act.ActorDynamicPlasma.evolve_pedestal = true

    # act.ActorCurrent.model = :replay
    # act.ActorEquilibrium.model = :replay
    if reconstruction
        act.ActorCoreTransport.model = :replay
        act.ActorPedestal.model = :replay
    end
    # act.ActorHCD.ec_model = :replay
    # act.ActorHCD.ic_model = :replay
    # act.ActorHCD.lh_model = :replay
    # act.ActorHCD.nb_model = :replay
    # act.ActorHCD.pellet_model = :replay
    # act.ActorHCD.neutral_model = :none

    # Run the simulation
    try
        @info "ActorDynamicPlasma(dd, act)"
        FUSE.ActorDynamicPlasma(dd, act; verbose=true)
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
        SimulationParameters.par2json(act,joinpath(savedir, "act.json"))
        act.ActorReplay.replay_dd = tmp

        @info "save dd_sim.json"
        IMAS.imas2json(dd, joinpath(savedir, "dd_sim.json"))

        @info "save dd_exp.json"
        IMAS.imas2json(dd_exp, joinpath(savedir, "dd_exp.json"))

        @info "save dd_benchmark.json"
        IMAS.imas2json(bnch.dd, joinpath(savedir, "dd_benchmark.json"))

        @info "save animated gif"
        mkdir(joinpath(savedir, "gif"))
        animated_plasma_overview(dd, joinpath(savedir, "gif"), dd_exp)
    end

    return nothing
end

function animated_plasma_overview(dd::IMAS.dd, dir::AbstractString, dd1::Union{IMAS.dd,Nothing}=nothing; aggregate_hcd::Bool=true, fps::Int=12)
    fulldir = abspath(dir)
    @assert isdir(fulldir) "$fulldir directory does not exist"
    #a = Interact.@animate
    for (k, time0) in enumerate(dd.equilibrium.time)
        try
            FUSE.plot_plasma_overview(dd, Float64(time0); dd1, aggregate_hcd)
            savefig(abspath(joinpath(fulldir, "$(@sprintf("%04d", k)).png")))
        catch e
            plot()
        end
    end
    #Interact.gif(a, abspath(joinpath(fulldir, "dd.gif")); fps)
end