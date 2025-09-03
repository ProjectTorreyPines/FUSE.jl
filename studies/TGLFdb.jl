using LaTeXStrings
import JSON
using TGLFNN

#= ================= =#
#  StudyTGLFdb  #
#= ================= =#

"""
study_parameters(::Type{Val{:TGLFdb}})::Tuple{FUSEparameters__ParametersStudyTGLFdb,ParametersAllActors}
"""
function study_parameters(::Type{Val{:TGLFdb}})::Tuple{FUSEparameters__ParametersStudyTGLFdb,ParametersAllActors}

    sty = FUSEparameters__ParametersStudyTGLFdb{Real}()
    act = ParametersActors()

    # Change act for the default TGLFdb run
    act.ActorCoreTransport.model = :FluxMatcher
    act.ActorFluxMatcher.evolve_pedestal = false
    act.ActorTGLF.warn_nn_train_bounds = false

    # finalize 
    set_new_base!(sty)
    set_new_base!(act)

    return sty, act
end

Base.@kwdef mutable struct FUSEparameters__ParametersStudyTGLFdb{T<:Real} <: ParametersStudy{T}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :StudyTGLFdb
    server::Switch{String} = study_common_parameters(; server="localhost")
    n_workers::Entry{Int} = study_common_parameters(; n_workers=missing)
    file_save_mode::Switch{Symbol} = study_common_parameters(; file_save_mode=:safe_write)
    release_workers_after_run::Entry{Bool} = study_common_parameters(; release_workers_after_run=true)
    save_dd::Entry{Bool} = study_common_parameters(; save_dd=true)
    sat_rules::Entry{Vector{Symbol}} = Entry{Vector{Symbol}}("-", "TGLF saturation rules to run")
    lump_ions::Entry{Bool} = Entry{Bool}("-", "Lumps the fuel species (D,T) as well as the impurities together"; default=true)
    custom_tglf_models::Entry{Vector{String}} = Entry{Vector{String}}("-", "This will run custom TGLFNN models stored in TGLFNN/models")
    save_folder::Entry{String} = Entry{String}("-", "Folder to save the database runs into")
    database_folder::Entry{String} = Entry{String}("-", "Folder with input database")
end

mutable struct StudyTGLFdb <: AbstractStudy
    sty::FUSEparameters__ParametersStudyTGLFdb
    act::ParametersAllActors
    dataframes_dict::Union{Dict{String,DataFrame},Missing}
    iterator::Union{Vector{Union{String,Symbol}},Missing}
end

function TGLF_dataframe()
    return DataFrame(;
        shot=Int[], time=Int[], ne0=Float64[],
        Te0=Float64[], Ti0=Float64[], ne0_exp=Float64[],
        Te0_exp=Float64[], Ti0_exp=Float64[], WTH_exp=Float64[],
        rot0_exp=Float64[], WTH=Float64[], #, t98=Float64[], ptot=Float64[], tau=Float64[],
        rot0=Float64[],timef=Int[], rho=Vector{Float64}[],
        Qe_target=Vector{Float64}[], Qe_TGLF=Vector{Float64}[], Qe_neoc=Vector{Float64}[],
        Qi_target=Vector{Float64}[], Qi_TGLF=Vector{Float64}[], Qi_neoc=Vector{Float64}[],
        particle_target=Vector{Float64}[], particle_TGLF=Vector{Float64}[], particle_neoc=Vector{Float64}[],
        momentum_target=Vector{Float64}[], momentum_TGLF=Vector{Float64}[],
        Q_GB=Vector{Float64}[], particle_GB=Vector{Float64}[], momentum_GB=Vector{Float64}[])
end

function StudyTGLFdb(sty, act; kw...)
    sty = sty(kw...)
    study = StudyTGLFdb(sty, act, missing, missing)
    return setup(study)
end

function _setup(study::StudyTGLFdb)
    sty = study.sty

    @assert !ismissing(getproperty(sty, :database_folder, missing)) "Specify the database_folder in sty"
    @assert !ismissing(readdir(sty.database_folder)) "There are no input cases in $(sty.database_folder)"

    check_and_create_file_save_mode(sty)

    preparse_input(sty.database_folder)

    parallel_environment(sty.server, sty.n_workers)

    return study
end

function _run(study::StudyTGLFdb)
    sty = study.sty
    act = study.act

    @assert sty.n_workers == length(Distributed.workers()) "The number of workers =  $(length(Distributed.workers())) isn't the number of workers you requested = $(sty.n_workers)"
    @assert ismissing(getproperty(sty, :sat_rules, missing)) ⊻ ismissing(getproperty(sty, :custom_tglf_models, missing)) "Specify either sat_rules or custom_tglf_models"

    cases_files = [
        joinpath(sty.database_folder, "fuse_prepared_inputs", item) for
        item in readdir(joinpath(sty.database_folder, "fuse_prepared_inputs")) if endswith(item, ".json")
    ]

    if !ismissing(getproperty(sty, :sat_rules, missing))
        iterator = sty.sat_rules
    elseif !ismissing(sty.custom_tglf_models)
        iterator = sty.custom_tglf_models
    end
    study.iterator = iterator
    study.dataframes_dict = Dict(string(name) => TGLF_dataframe() for name in study.iterator)

    # loop serially through saturation rules
    println("Running StudyTGLFdb on $(length(cases_files)) cases on $(iterator) with $(sty.n_workers) workers on $(sty.server)")
    for item in iterator
        if !ismissing(getproperty(sty, :sat_rules, missing))
            act.ActorTGLF.sat_rule = item
        else
            act.ActorTGLF.user_specified_model = item
        end
        act.ActorTGLF.lump_ions = sty.lump_ions
        if item == "wrapped_model.onnx"
            act.ActorTGLF.onnx_model=true
            act.ActorTGLF.user_specified_model = item #"/mnt/beegfs/users/neisert/.julia/dev/TGLFNN/models/wrapped_model.onnx"
        end

        # paraller run
        results = pmap(filename -> run_case(filename, study, item), cases_files)

        # populate DataFrame
        for row in results
            if row isa NamedTuple || row isa AbstractArray || row isa DataFrameRow || row isa AbstractDict
                push!(study.dataframes_dict[string(item)], row)
            else
                @warn "Invalid row type encountered: $row"
            end
        end

        # Save JSON to a file
        json_data = JSON.json(study.dataframes_dict[string(item)])
        open("$(sty.save_folder)/data_frame_$(item).json", "w") do f
            return write(f, json_data)
        end

    end

    # Release workers after run
    if sty.release_workers_after_run
        Distributed.rmprocs(Distributed.workers())
        @info "released workers"
    end
    return study
end

function _analyze(study::StudyTGLFdb)
    plot_xy_wth_hist2d(study; quantity=:WTH, save_fig=false, save_path="")
    return study
end

function preprocess_dd(filename)
    dd = IMAS.json2imas(filename; verbose=false)

    dd.summary.local.pedestal.n_e.value = [IMAS.pedestal_finder(dd.core_profiles.profiles_1d[].electrons.density_thermal, dd.core_profiles.profiles_1d[].grid.psi_norm)[1]]
    dd.summary.local.pedestal.zeff.value = [2.2] # [IMAS.pedestal_finder(dd.core_profiles.profiles_1d[].zeff, dd.core_profiles.profiles_1d[].grid.psi_norm)[1]]
    dd.pulse_schedule.tf.time = dd.summary.time
    dd.pulse_schedule.tf.b_field_tor_vacuum_r.reference = dd.equilibrium.vacuum_toroidal_field.b0

    return dd
end

function run_case(filename, study, item)
    act = study.act
    sty = study.sty

    dd = preprocess_dd(filename)

    #act.ActorFluxMatcher.evolve_densities = FUSE.setup_density_evolution_fixed(dd)
    act.ActorFluxMatcher.evolve_densities = FUSE.setup_density_evolution_electron_flux_match_impurities_fixed(dd)

    # find time from filename
    timefn = match(r"(\d+)\.\w+$", filename)
    if timefn !== nothing
        timef = parse(Int,timefn.captures[1])
    else
        timef = 0
    end

    cp1d = dd.core_profiles.profiles_1d[]
    exp_values = [
        cp1d.electrons.density_thermal[1], cp1d.electrons.temperature[1],
        cp1d.ion[1].temperature[1], @ddtime(dd.summary.global_quantities.energy_thermal.value),
        cp1d.rotation_frequency_tor_sonic[1], timef]

    name = split(splitpath(filename)[end], ".")[1]
    output_case = joinpath(sty.save_folder, name)

    try
        if !isdir(output_case)
            mkdir(output_case)
        end

        actor_transport = workflow_actor(dd, act)
        empty!(dd.summary.global_quantities)

        if sty.save_dd
            IMAS.imas2json(dd, joinpath(output_case, "result_dd_$(item).json"))
            save_inputtglfs(actor_transport, output_case, name, item)
            if parse(Bool, get(ENV, "FUSE_MEMTRACE", "false"))
                save(FUSE.memtrace, joinpath(output_case, "memtrace.txt"))
            end
        end

        return create_data_frame_row(dd, exp_values)
    catch e
        open("$output_case/error.txt", "w") do file
            return showerror(file, e, catch_backtrace())
        end
    end
end

function workflow_actor(dd, act)
    # Actors to run on the input dd

    actor_transport = ActorCoreTransport(dd, act)

    # actor_equilibrium = ActorEquilibrium(dd,act)

    return actor_transport
end

function save_inputtglfs(actor_transport, output_dir, name, item)
    rho_transport = actor_transport.tr_actor.actor_ct.par.rho_transport
    for k in 1:length(rho_transport)
        TGLFNN.save(actor_transport.tr_actor.actor_ct.actor_turb.input_tglfs[k], joinpath(output_dir, "input.tglf_$(name)_$(k)_$(item)"))
    end
end

function create_data_frame_row(dd::IMAS.dd, exp_values::AbstractArray)
    cp1d = dd.core_profiles.profiles_1d[]
    eqt  = dd.equilibrium.time_slice[]

    rho_transport = dd.core_transport.model[1].profiles_1d[].grid_flux.rho_tor_norm

    # baseline TGLF and IMAS total‐flux objects (we still use these for the “TGLF” & “neoc” columns)
    ct1d_tglf   = dd.core_transport.model[1].profiles_1d[]

    ini, act = FUSE.case_parameters(:ITER; init_from = :scalars)
    act.ActorFluxMatcher.rho_transport   = rho_transport
    act.ActorFluxMatcher.evolve_rotation = :flux_match

    nr    = length(rho_transport)
    q_mat = reshape(
      FUSE.flux_match_targets(dd, act.ActorFluxMatcher, nothing),
      nr, 4
    )
    qi, qe, qp, qg = eachcol(q_mat)
    qybro_bohms = [
      IMAS.gyrobohm_energy_flux(cp1d, eqt),
      IMAS.gyrobohm_particle_flux(cp1d, eqt),
      IMAS.gyrobohm_momentum_flux(cp1d, eqt),
    ]
    rho_cp = cp1d.grid.rho_tor_norm

    IMAS.interp1d(rho_cp, qybro_bohms[1]).(rho_transport)

    return (
      shot            = dd.dataset_description.data_entry.pulse,
      time            = Int(dd.summary.time[1] * 1000),
      ne0             = cp1d.electrons.density_thermal[1],
      Te0             = cp1d.electrons.temperature[1],
      Ti0             = cp1d.ion[1].temperature[1],
      WTH             = IMAS.@ddtime(dd.summary.global_quantities.energy_thermal.value),
      rot0            = cp1d.rotation_frequency_tor_sonic[1],

      ne0_exp         = exp_values[1],
      Te0_exp         = exp_values[2],
      Ti0_exp         = exp_values[3],
      WTH_exp         = exp_values[4],
      rot0_exp        = exp_values[5],
      timef           = exp_values[6],

      rho             = rho_transport,
      Qe_target       = qe,
      Qe_TGLF         = ct1d_tglf.electrons.energy.flux,
      Qe_neoc         = dd.core_transport.model[2].profiles_1d[].electrons.energy.flux,

      Qi_target       = qi,
      Qi_TGLF         = ct1d_tglf.total_ion_energy.flux,
      Qi_neoc         = dd.core_transport.model[2].profiles_1d[].total_ion_energy.flux,

      particle_target = qg,
      particle_TGLF   = ct1d_tglf.electrons.particles.flux,
      particle_neoc   = dd.core_transport.model[2].profiles_1d[].electrons.particles.flux,

      momentum_target = qp,
      momentum_TGLF   = ct1d_tglf.momentum_tor.flux,

      Q_GB            = IMAS.interp1d(rho_cp, qybro_bohms[1]).(rho_transport),
      particle_GB     = IMAS.interp1d(rho_cp, qybro_bohms[2]).(rho_transport),
      momentum_GB     = IMAS.interp1d(rho_cp, qybro_bohms[3]).(rho_transport),
    )
end

function preparse_input(database_folder)
    if isdir(joinpath(database_folder, "fuse_prepared_inputs"))
        return
    end

    fuse_prepared_inputs = joinpath(database_folder, "fuse_prepared_inputs")
    mkdir(fuse_prepared_inputs)

    machine = "DIII-D"

    for item in readdir(database_folder)
        if occursin("_mod", item) || occursin(".jl", item) || isdir(joinpath(database_folder, item))
            continue
        end
        filename = split("$(database_folder)/$item", ".")[1]

        pulse = parse(Int, split(splitpath(filename)[end], "_")[2])

        json_data = JSON.parsefile("$filename.json")

        delete!(json_data["equilibrium"], "code")

        time = json_data["equilibrium"]["time"][1]
        timea = [time]

        json_data["core_profiles"]["profiles_1d"][1]["time"] = time
        json_data["core_profiles"]["time"] = timea
        json_data["dataset_description"] = Dict("data_entry" => Dict("pulse" => pulse, "machine" => machine))

        for source in keys(json_data["core_sources"]["source"])
            json_data["core_sources"]["source"][source]["profiles_1d"][1]["time"] = time
        end


        json_data["summary"]["time"] = [time]
        json_string = JSON.json(json_data)

        open("$(joinpath(fuse_prepared_inputs,item))", "w") do file
            return write(file, json_string)
        end
    end
end


function plot_xy_wth_hist2d(study; quantity=:WTH, save_fig=false, save_path="")

    if study.act.ActorTGLF.electromagnetic
        EM_contribution = :EM
    else
        EM_contribution = :ES
    end

    for item in study.iterator
        plot_xy_wth_hist2d(study.dataframes_dict[string(item)], string(item), EM_contribution, quantity, save_fig, save_path)
    end
end

function plot_xy_wth_hist2d(df::DataFrame, name::String, EM_contribution::Symbol, quantity::Symbol, save_fig::Bool, save_path::String)

    x = df[!, "$(quantity)_exp"]
    y = df[!, "$(quantity)"]

    MRE = round(100 * mean_relative_error(x, y); digits=2)

    bins = 10 .^ (4:0.05:7)

    ticks = ([10^4, 10^5, 10^6, 10^7, 10^8], [L"10^4", L"10^5", L"10^6", L"10^7", L"10^8"])
    xy_lim = [bins[1], bins[end]]

    p = histogram2d(x, y; xscale=:log10, yscale=:log10, bins=(bins, bins), xticks=ticks, yticks=ticks,
        color=cgrad(:magma; rev=true), colorbar=true, show_empty_bins=true,
        ylabel="Thermal stored energy predicted [J]", xlabel="Thermal stored energy experiment [J]",
        ylim=xy_lim, xlim=xy_lim, tickfont=font(12, "Computer Modern"), fontfamily="Computer Modern",
        xguidefontsize=15, yguidefontsize=14)

    plot!(xy_lim, xy_lim; linestyle=:dash, color=:black, label=nothing, title="$(replace(uppercase(name), "_" => " "))")
    println("MRE $(replace(uppercase(name), "_" => " ")) W_thermal = $MRE % with N = $(length(x))")

    if save_fig
        savefig(p, save_path)
    else
        display(p)
    end
end
