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
    kw_case_parameters::Entry{Dict{Symbol,Any}} = Entry{Dict{Symbol,Any}}("-", "Keyword arguments passed to case_parameters"; default=Dict{Symbol,Any}())
    redirect_output::Entry{Bool} = Entry{Bool}("-", "Redirect stdout and stderr to log.txt file"; default=true)
    verbose::Entry{Bool} = Entry{Bool}("-", "Turn on verbose progress output"; default=true)

    # Postdictive-specific parameters
    device::Entry{Symbol} = Entry{Symbol}("-", "Device to run postdictive simulations for")
    shots::Entry{Vector{Int}} = Entry{Vector{Int}}("-", "List of shot numbers")
    reconstruction::Entry{Bool} = Entry{Bool}("-", "Run postdiction in reconstruction mode")
    start_time::Entry{Float64} = Entry{Float64}("s", "Override simulation start time"; default=missing)
    end_time::Entry{Float64} = Entry{Float64}("s", "Simulation end time (missing: end of the equilibrium time basis; later times are clipped to it)"; default=missing)
    δt::Entry{Float64} = Entry{Float64}("s", "Macro time step of ActorDynamicPlasma (with δt_exact = false, the largest step ≤ δt that ends at end_time)"; default=0.05, check=x -> @assert x > 0 "δt must be > 0")
    δt_exact::Entry{Bool} = Entry{Bool}("-", "true: steps of exactly δt (Nt = floor(Δt / δt)), ending up to one step before end_time; false: Nt = ceil(Δt / δt) steps that end at end_time"; default=false)
    qed_δt::Entry{Float64} = Entry{Float64}("s", "ActorQED internal time step (sets ActorQED.Nt per macro step; missing: ActorQED.Nt as is)"; default=missing, check=x -> @assert x > 0 "qed_δt must be > 0")
    qmin_desired::Entry{Float64} = Entry{Float64}("-", "ActorQED and ActorSawteethSource qmin_desired"; default=1.0, check=x -> @assert x >= 0 "qmin_desired must be >= 0")
    sawteeth_period::Entry{Float64} = Entry{Float64}("s", "ActorSawteethSource.period: keep flattening this long after the last sawtooth"; default=0.25)
end

mutable struct StudyPostdictive{T<:Real} <: AbstractStudy
    sty::OverrideParameters{T,FUSEparameters__ParametersStudyPostdictive{T}}
    act::ParametersActors # actor parameters forced on every case; unset (missing) ones keep the case and study settings
end

function StudyPostdictive(sty::ParametersStudy; kw...)
    sty = OverrideParameters(sty; kw...)
    study = StudyPostdictive(sty, unset_parameters!(ParametersActors()))
    parallel_environment(sty.server, sty.n_workers)
    return study
end

"""
    _run(study::StudyPostdictive)

Runs the Postdictive study with sty settings in parallel on designated cluster
"""
function _run(study::StudyPostdictive)
    sty = study.sty

    @assert (sty.n_workers == 0 || sty.n_workers == length(Distributed.workers())) "The number of workers = $(length(Distributed.workers())) isn't the number of workers you requested = $(sty.n_workers)"

    # parallel run
    @info "StudyPostdictive: running $(length(sty.shots)) postdictive simulations with $(sty.n_workers) workers on $(sty.server)"

    ProgressMeter.@showprogress map(shot -> run_postdictive_case(study, shot; sty.kw_case_parameters), sty.shots)

    # Release workers after run
    if sty.release_workers_after_run
        Distributed.rmprocs(Distributed.workers())
        @info "StudyPostdictive: released workers"
    end

    return study
end

"""
    run_postdictive_case(study::StudyPostdictive, shot::Int; kw_case_parameters::Dict{Symbol,Any})

Run a single postdictive case for a given device and shot
"""
function run_postdictive_case(study::StudyPostdictive, shot::Int; kw_case_parameters::Dict{Symbol,Any})
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
            # these may be unset (missing), which reading them as sty.<name> refuses
            start_time, end_time, qed_δt = (getproperty(sty, name, missing) for name in (:start_time, :end_time, :qed_δt))
            run_postdictive_case(device, shot; savedir, sty.reconstruction, sty.verbose, start_time, end_time, sty.δt, sty.δt_exact, qed_δt,
                sty.qmin_desired, sty.sawteeth_period, act_overrides=study.act, kw_case_parameters)
        end
    end
end

function run_postdictive_case(device::Symbol, shot::Int; kw_case_parameters::Dict{Symbol,Any}, kw...)
    dd = IMAS.dd()
    dd_exp = IMAS.dd()
    run_postdictive_case!(dd, dd_exp, device, shot; kw_case_parameters, kw...)
    return (dd=dd, dd_exp=dd_exp)
end

function run_postdictive_case!(
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
    δt::Float64=0.05,
    δt_exact::Bool=false,
    qed_δt::Union{Float64,Missing}=missing,
    qmin_desired::Float64=1.0,
    sawteeth_period::Float64=0.25,
    act_overrides::Union{Nothing,ParametersAllActors}=nothing,
    kw_case_parameters::Dict{Symbol,Any}
)

    isterminal = isa(stdout, Base.TTY) && (get(ENV, "CI", nothing) != "true")
    if !isterminal && verbose
        @info "StudyPostdictive: verbose output requested but stdout is not a terminal; setting verbose=false"
        verbose = false
    end

    # Get case parameters
    @info "StudyPostdictive: case_parameters($(repr(device)), $shot; $(repr(kw_case_parameters))...)"
    ini, act = FUSE.case_parameters(device, shot; kw_case_parameters...)
    if !ismissing(start_time)
        ini.time.simulation_start = start_time
    end

    # time window
    t0 = ini.time.simulation_start # start_time should be early in the shot, when otherwise ohmic current will be wrong
    final_time = ini.general.dd.equilibrium.time[end]
    if !ismissing(end_time)
        final_time = min(end_time, final_time)
    end
    Nt = δt_exact ? floor(Int, (final_time - t0) / δt) : ceil(Int, (final_time - t0) / δt)
    @assert Nt >= 1 "empty time window: simulation start $t0 s, end $final_time s, δt $δt s"
    Δt = δt_exact ? Nt * δt : final_time - t0

    # init
    @info "StudyPostdictive: ini.time.simulation_start = $(ini.time.simulation_start)"
    FUSE.init!(dd, ini, act)

    # keep aside the dd with experimental data
    IMAS.fill!(dd_exp, dd)

    # identify LH transitions
    @info "StudyPostdictive: LH_analysis"
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
    #act.ActorEGGO.timeslice_average = 4

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
    act.ActorFluxMatcher.rho_transport = 0.15:0.1:0.85

    act.ActorPedestal.rotation_model = :replay

    act.ActorSawteethSource.flat_factor = 1.0
    act.ActorSawteethSource.period = sawteeth_period # turn off flattening after this long with no sawteeth events
    act.ActorQED.qmin_desired = qmin_desired
    act.ActorSawteethSource.qmin_desired = qmin_desired
    act.ActorTGLF.model = :GKNN
    act.ActorTGLF.tglfnn_model = "sat3_em_d3d_azf-1_withnegD"

    # time
    dd.global_time = t0
    act.ActorDynamicPlasma.Nt = Nt
    act.ActorDynamicPlasma.Δt = Δt
    if !ismissing(qed_δt)
        act.ActorQED.Nt = round(Int, Δt / Nt / qed_δt)
    end

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

    if act_overrides !== nothing
        apply_act_overrides!(act, act_overrides)
    end

    if !reconstruction
        # it's best to start from a transport simulation that is flux-matched
        ped_mod = act.ActorPedestal.model
        act.ActorPedestal.model = :WPED
        FUSE.ActorFluxMatcher(dd, act; verbose, evolve_plasma_sources=true, evolve_pedestal=true, max_iterations=1000, do_plot=false);
        act.ActorPedestal.model = ped_mod
    end

    # Run the simulation
    try
        @info "StudyPostdictive: ActorDynamicPlasma(dd, act)"
        FUSE.ActorDynamicPlasma(dd, act; verbose)
    catch e
        if isa(e, InterruptException)
            rethrow(e)
        else
            @error repr(e)
        end
    end

    # save simulation data to directory
    if !isempty(savedir)
        @info "StudyPostdictive: IMAS.benchmark(dd, dd_exp, dd.core_profiles.time);"
        bnch = IMAS.benchmark(dd, dd_exp, dd.core_profiles.time);

        @info "StudyPostdictive: saving simulation results to: $(savedir)"
        @info "StudyPostdictive: save act.json"
        tmp = act.ActorReplay.replay_dd
        act.ActorReplay.replay_dd = IMAS.dd()
        SimulationParameters.par2json(act,joinpath(savedir, "act.json"))
        act.ActorReplay.replay_dd = tmp

        @info "StudyPostdictive: save dd_sim.json"
        IMAS.imas2json(dd, joinpath(savedir, "dd_sim.json"))

        @info "StudyPostdictive: save dd_exp.json"
        IMAS.imas2json(dd_exp, joinpath(savedir, "dd_exp.json"))

        @info "StudyPostdictive: save dd_benchmark.json"
        IMAS.imas2json(bnch.dd, joinpath(savedir, "dd_benchmark.json"))

        open(joinpath(savedir, "step_provenance.json"), "w") do io
            JSON.print(io, Dict("start" => t0, "δt" => Δt / Nt, "Nt" => Nt, "end" => t0 + Δt, "reached" => dd.global_time), 1)
        end

        if save_gif
            @info "StudyPostdictive: save animated gif"
            mkdir(joinpath(savedir, "gif"))
            animated_plasma_overview(dd, joinpath(savedir, "gif"), dd_exp)
        end
    end

    return (; act, t0, Δt, Nt)
end

function animated_plasma_overview(dd::IMAS.DD, dir::AbstractString, dd1::Union{IMAS.dd,Nothing}=nothing; aggregate_hcd::Bool=true, fps::Int=12)
    fulldir = abspath(dir)
    @assert isdir(fulldir) "$fulldir directory does not exist"
    a = Plots.@animate for (k, time0) in enumerate(dd.equilibrium.time)
        try
            FUSE.plot_plasma_overview(dd, Float64(time0); dd1, aggregate_hcd)
            savefig(abspath(joinpath(fulldir, "$(@sprintf("%04d", k)).png")))
        catch e
            plot()
        end
    end
    Plots.gif(a, abspath(joinpath(fulldir, "dd.gif")); fps)
end