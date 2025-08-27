#= ================ =#
#  StudyPostdictive  #
#= ================ =#

"""
    study_parameters(::Val{:Postdictive})::Tuple{FUSEparameters__ParametersStudyPostdictive,ParametersAllActors}

Runs postdictive simulations for specified device and shots with configurable parameters
"""
function study_parameters(::Val{:Postdictive})::Tuple{FUSEparameters__ParametersStudyPostdictive,ParametersAllActors}
    return FUSEparameters__ParametersStudyPostdictive{Real}()
end

Base.@kwdef mutable struct FUSEparameters__ParametersStudyPostdictive{T<:Real} <: ParametersStudy{T}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :StudyPostdictive
    server::Switch{String} = study_common_parameters(; server="localhost")
    n_workers::Entry{Int} = study_common_parameters(; n_workers=1)
    file_save_mode::Switch{Symbol} = study_common_parameters(; file_save_mode=:safe_write)
    release_workers_after_run::Entry{Bool} = study_common_parameters(; release_workers_after_run=true)
    save_folder::Entry{String} = Entry{String}("-", "Folder to save the postdictive runs into")

    # Postdictive-specific parameters
    device::Entry{Symbol} = Entry{Symbol}("-", "Device to run postdictive simulations for")
    shots::Entry{Vector{Int}} = Entry{Vector{Int}}("-", "List of shot numbers")
    fit_profiles::Entry{Bool} = Entry{Bool}("-", "Whether to fit profiles in case_parameters"; default=true)
    use_local_cache::Entry{Bool} = Entry{Bool}("-", "Whether to use local cache in case_parameters"; default=false)
end

mutable struct StudyPostdictive{T<:Real} <: AbstractStudy
    sty::OverrideParameters{T,FUSEparameters__ParametersStudyPostdictive{T}}
    ini::Union{ParametersAllInits,Missing}
    act::Union{ParametersAllActors,Missing}
    dataframe::Union{DataFrame,Missing}
    workflow::Union{Function,Missing}
end

function StudyPostdictive(sty::ParametersStudy; kw...)
    sty = OverrideParameters(sty; kw...)
    study = StudyPostdictive(sty, missing, missing, missing, missing)

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
    n_simulations = length(sty.shots)

    # parallel run
    println("running $(n_simulations) postdictive simulations with $(sty.n_workers) workers on $(sty.server)")

    if !isdir(sty.save_folder)
        mkdir(sty.save_folder)
    end

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
    file_log = open("log.txt", "w")

    try
        redirect_stdout(file_log)
        redirect_stderr(file_log)

        # Get default parameters for the study
        _, act = study_parameters(:Postdictive)

        run_postdictive_case(device, shot; sty.fit_profiles, sty.use_local_cache, sty.save_folder)

        # catch error
        #     if isa(error, InterruptException)
        #         rethrow(error)
        #     end
        #     @error "Error in postdictive case for $(device) shot $(shot): $error"
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

function run_postdictive_case!(dd::IMAS.dd, dd_exp::IMAS.dd, device::Symbol, shot::Int; fit_profiles::Bool, use_local_cache::Bool, save_folder::AbstractString=abspath("."))
    savedir = abspath(joinpath(save_folder, "$(device)_$(shot)__$(Dates.now())__$(getpid())"))
    if !isdir(savedir)
        mkdir(savedir)
    end

    # Get case parameters
    ini, act = FUSE.case_parameters(device, shot; fit_profiles, use_local_cache)

    # Override act with study-specific actor parameters
    #merge!(act, act)

    # init
    ini.time.simulation_start = ini.general.dd.equilibrium.time_slice[2].time
    FUSE.init!(dd, ini, act)

    # keep aside the dd with experimental data
    IMAS.fill!(dd_exp, dd)

    # identify LH transitions
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

    if false
        # LH-transition from LH scaling law
        act.ActorPedestal.mode_transitions = missing
    else
        # LH-transition at user-defined times
        act.ActorPedestal.mode_transitions = experiment_LH.mode_transitions
        act.ActorPedestal.mode_transitions[5.2] = :L_mode
    end

    act.ActorEquilibrium.model = :FRESCO #:EGGO or FRESCO
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

    #act.ActorTGLF.tglfnn_model = "sat1_em_d3d"

    # time
    δt = 0.025
    dd.global_time = ini.general.dd.equilibrium.time_slice[2].time # start_time should be early in the shot, when otherwise ohmic current will be wrong
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
    # act.ActorCoreTransport.model = :replay
    # act.ActorPedestal.model = :replay
    # act.ActorHCD.ec_model = :replay
    # act.ActorHCD.ic_model = :replay
    # act.ActorHCD.lh_model = :replay
    # act.ActorHCD.nb_model = :replay
    # act.ActorHCD.pellet_model = :replay
    # act.ActorHCD.neutral_model = :none

    # Run the simulation
    FUSE.ActorDynamicPlasma(dd, act; verbose=true)

    # save simulation data to directory
    if !ismissing(save_folder)
        IMAS.imas2json(dd, joinpath(savedir, "dd_sim.json"))
        IMAS.imas2json(dd_exp, joinpath(savedir, "dd_exp.json"))
        mkdir(joinpath(savedir, "gif"))
        animated_plasma_overview(dd, joinpath(savedir, "gif"), dd_exp)
    end

    return nothing
end

function animated_plasma_overview(dd::IMAS.dd, dir::AbstractString, dd1::Union{IMAS.dd,Nothing}=nothing; aggregate_hcd::Bool=true, fps::Int=12)
    #a = Interact.@animate 
    fulldir = abspath(dir)
    @assert isdir(fulldir) "$fulldir directory does not exist"
    for (k, time0) in enumerate(dd.equilibrium.time)
        try
            FUSE.plot_plasma_overview(dd, Float64(time0); dd1, aggregate_hcd)
            savefig(abspath(joinpath(fulldir, "$(@sprintf("%04d", k)).png")))
            # magick -delay 2 -loop 0 D3D_168830___*.png -layers Optimize D3D_168830.gif
        catch e
            rethrow(e)
            plot()
        end
    end
    #Interact.gif(a, "dir/dd.gif"; fps)
end
