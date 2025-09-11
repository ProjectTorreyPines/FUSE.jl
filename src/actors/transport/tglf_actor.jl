import TurbulentTransport
import TurbulentTransport: InputTGLF, InputQLGYRO, InputTJLF, InputCGYRO
import TJLF
import GACODE

#= ========= =#
#  ActorTGLF  #
#= ========= =#
@actor_parameters_struct ActorTGLF{T} begin
    model::Switch{Symbol} = Switch{Symbol}([:TGLF, :TGLFNN, :GKNN, :TJLF], "-", "Implementation of TGLF"; default=:TGLFNN)
    onnx_model::Entry{Bool} = Entry{Bool}("-", "use onnx model"; default=false)
    sat_rule::Switch{Symbol} = Switch{Symbol}([:sat0, :sat0quench, :sat1, :sat1geo, :sat2, :sat3], "-", "Saturation rule"; default=:sat1)
    electromagnetic::Entry{Bool} = Entry{Bool}("-", "Electromagnetic or electrostatic"; default=true)
    tglfnn_model::Entry{String} = Entry{String}(
        "-",
        "Use a user specified TGLF-NN model stored in TGLFNN/models";
        check=x -> @assert x in TurbulentTransport.available_models() "ActorTGLF.tglfnn_model must be one of:\n  \"$(join(TurbulentTransport.available_models(),"\"\n  \""))\""
    )
    rho_transport::Entry{AbstractVector{T}} = Entry{AbstractVector{T}}("-", "rho_tor_norm values to compute tglf fluxes on"; default=0.25:0.1:0.85)
    warn_nn_train_bounds::Entry{Bool} = Entry{Bool}("-", "Raise warnings if querying cases that are certainly outside of the training range"; default=false)
    custom_input_files::Entry{Union{Vector{<:InputTGLF},Vector{<:InputTJLF}}} =
        Entry{Union{Vector{<:InputTGLF},Vector{<:InputTJLF}}}("-", "Sets up the input file that will be run with the custom input file as a mask")
    lump_ions::Entry{Bool} = Entry{Bool}("-", "Lumps the fuel species (D,T) as well as the impurities together"; default=true)
    save_input_tglfs_to_folder::Entry{String} = Entry{String}("-", "Save the intput.tglf files in designated folder"; default="")
end

mutable struct ActorTGLF{D,P} <: SingleAbstractActor{D,P}
    dd::IMAS.dd{D}
    par::OverrideParameters{P,FUSEparameters__ActorTGLF{P}}
    input_tglfs::Union{Vector{InputTGLF{D}},Vector{InputTJLF{D}}}
    flux_solutions::Vector{GACODE.FluxSolution{D}}
end

"""
    ActorTGLF(dd::IMAS.dd, act::ParametersAllActors; kw...)

Evaluates turbulent transport fluxes using TGLF-based models including TGLF, TGLFNN, GKNN, and TJLF.

Supported models:
- `:TGLF`: Full TGLF gyrokinetic transport model
- `:TGLFNN`: Neural network surrogate of TGLF for fast evaluation
- `:GKNN`: Neural network surrogate with gyrokinetic fidelity
- `:TJLF`: Weiland model for ITG/TEM turbulence with adaptive spectral widths

The actor creates input files for the selected model based on local plasma parameters,
runs the transport calculation, and stores results as FluxSolution objects containing
electron/ion energy fluxes, particle fluxes, and momentum flux. Results are normalized
in gyro-Bohm units and written to `dd.core_transport` as anomalous transport.
"""
function ActorTGLF(dd::IMAS.dd, act::ParametersAllActors; kw...)
    actor = ActorTGLF(dd, act.ActorTGLF; kw...)
    step(actor)
    finalize(actor)
    return actor
end

function ActorTGLF(dd::IMAS.dd{D}, par::FUSEparameters__ActorTGLF; kw...) where {D<:Real}
    logging_actor_init(ActorTGLF)
    par = OverrideParameters(par; kw...)
    if par.model ∈ [:TGLF, :TGLFNN, :GKNN]
        input_tglfs = Vector{InputTGLF{D}}(undef, length(par.rho_transport))
    elseif par.model == :TJLF
        input_tglfs = Vector{InputTJLF{D}}(undef, length(par.rho_transport))
    end
    return ActorTGLF(dd, par, input_tglfs, GACODE.FluxSolution{D}[])
end

"""
    _step(actor::ActorTGLF)

Runs the selected TGLF-based transport model to evaluate turbulent fluxes on radial grid points.

Creates input files for each transport grid point containing local plasma parameters,
gradients, and geometry. For TJLF, maintains memory of spectral widths between iterations.
Supports custom input file parameters and can save input files to disk for debugging.
Runs the transport model (TGLF/TGLFNN/GKNN via TurbulentTransport.jl or TJLF) and
stores the resulting fluxes as FluxSolution objects.
"""
function _step(actor::ActorTGLF{D,P}) where {D<:Real, P<:Real}
    par = actor.par
    dd = actor.dd

    input_tglfs = InputTGLF(dd, par.rho_transport, par.sat_rule, par.electromagnetic, par.lump_ions)
    for k in eachindex(par.rho_transport)
        input_tglf = input_tglfs[k]
        if par.model ∈ [:TGLF, :TGLFNN, :GKNN]
            actor.input_tglfs[k] = input_tglf
        elseif par.model == :TJLF
            if !isassigned(actor.input_tglfs, k) # this is done to keep memory of the widths
                actor.input_tglfs[k] = InputTJLF{D}(input_tglf)
            else
                TJLF.update_input_tjlf!(actor.input_tglfs[k], input_tglf)
            end
        end

        # Overwrite TGLF / TJLF parameters with the custom parameters mask
        if !ismissing(par, :custom_input_files)
            for field_name in fieldnames(typeof(actor.input_tglfs[k]))
                if !ismissing(getproperty(par.custom_input_files[k], field_name))
                    setproperty!(actor.input_tglfs[k], field_name, getproperty(par.custom_input_files[k], field_name))
                end
            end
        end
        if isdir(par.save_input_tglfs_to_folder)
            name = lowercase(string(par.model))
            save(actor.input_tglfs[k] ,joinpath(par.save_input_tglfs_to_folder, "input.$(name)_$(Dates.format(Dates.now(), "yyyymmddHHMMSS"))_$(par.rho_transport[k])"))
        end

    end

    if par.model ∈ [:TGLFNN, :GKNN]
        if !par.onnx_model
            uncertain = eltype(dd) <: IMAS.Measurements.Measurement
            actor.flux_solutions = TurbulentTransport.run_tglfnn(actor.input_tglfs; uncertain, par.warn_nn_train_bounds, model_filename=model_filename(par), fidelity=par.model)

        elseif par.onnx_model
            sols32 = TurbulentTransport.run_tglfnn_onnx(actor.input_tglfs, par.tglfnn_model,
                [
                    "RLTS_3",
                    "KAPPA_LOC",
                    "ZETA_LOC",
                    "TAUS_3",
                    "VPAR_1",
                    "Q_LOC",
                    "RLNS_1",
                    "TAUS_2",
                    "Q_PRIME_LOC",
                    "P_PRIME_LOC",
                    "ZMAJ_LOC",
                    "VPAR_SHEAR_1",
                    "RLTS_2",
                    "S_DELTA_LOC",
                    "RLTS_1",
                    "RMIN_LOC",
                    "DRMAJDX_LOC",
                    "AS_3",
                    "RLNS_3",
                    "DZMAJDX_LOC",
                    "DELTA_LOC",
                    "S_KAPPA_LOC",
                    "ZEFF",
                    "VEXB_SHEAR",
                    "RMAJ_LOC",
                    "AS_2",
                    "RLNS_2",
                    "S_ZETA_LOC",
                    "BETAE_log10",
                    "XNUE_log10",
                    "DEBYE_log10"
                ], [
                    "OUT_G_elec",
                    "OUT_Q_elec",
                    "OUT_Q_ions",
                    "OUT_P_ions"
                ]; intra_threads=1, inter_threads=1)
                n = length(sols32)
                fs64 = Vector{GACODE.FluxSolution{Float64}}(undef, n)
                @inbounds for i in 1:n
                    s = sols32[i]
                    f1 = Float64(getfield(s, 1))
                    f2 = Float64(getfield(s, 2))
                    f3 = Float64(getfield(s, 3))
                    f4raw = getfield(s, 4)
                    f4 = isa(f4raw, AbstractArray) ? Float64.(f4raw) : Float64(f4raw)
                    f5 = Float64(getfield(s, 5))
                    fs64[i] = GACODE.FluxSolution{Float64}(f1, f2, f3, f4, f5)
                end
                actor.flux_solutions = fs64
        end

    elseif par.model == :TGLF
        actor.flux_solutions = TurbulentTransport.run_tglf(actor.input_tglfs)

    elseif par.model == :TJLF
        QL_fluxes_out = TJLF.run_tjlf(actor.input_tglfs)
        actor.flux_solutions =
            [GACODE.FluxSolution{D}(TJLF.Qe(QL_flux_out), TJLF.Qi(QL_flux_out), TJLF.Γe(QL_flux_out), TJLF.Γi(QL_flux_out), TJLF.Πi(QL_flux_out)) for QL_flux_out in QL_fluxes_out]
    end

    return actor
end

"""
    _finalize(actor::ActorTGLF)

Writes calculated transport fluxes to `dd.core_transport.model[:anomalous]`.

Sets the model identifier to include the model name, configuration details (saturation rule,
electromagnetic vs electrostatic for TGLF), and neural network model filename for TGLFNN/GKNN.
Converts flux results from GACODE format to IMAS format for electron/ion energy, particle,
and momentum fluxes using `GACODE.flux_gacode_to_imas`.
"""
function _finalize(actor::ActorTGLF)
    dd = actor.dd
    par = actor.par
    cp1d = dd.core_profiles.profiles_1d[]
    eqt = dd.equilibrium.time_slice[]

    model = resize!(dd.core_transport.model, :anomalous; wipe=false)
    model.identifier.name = string(par.model) * " " * model_filename(par)
    if par.onnx_model
        model.identifier.name *= " onnx"
    end
    m1d = resize!(model.profiles_1d)
    m1d.grid_flux.rho_tor_norm = par.rho_transport

    GACODE.flux_gacode_to_imas((:electron_energy_flux, :ion_energy_flux, :electron_particle_flux, :ion_particle_flux, :momentum_flux), actor.flux_solutions, m1d, eqt, cp1d)

    return actor
end

function model_filename(par::OverrideParameters{P,FUSEparameters__ActorTGLF{P}}) where {P<:Real}
    if par.model ∈ [:TGLFNN, :GKNN]
        filename = par.tglfnn_model
    else
        filename = string(par.sat_rule) * "_" * (par.electromagnetic ? "em" : "es")
    end
    return filename
end

function Base.show(io::IO, ::MIME"text/plain", input::Union{InputTGLF,InputTJLF})
    for field_name in fieldnames(typeof(input))
        println(io, " $field_name = $(getfield(input,field_name))")
    end
end

"""
    save(input::Union{InputCGYRO, InputQLGYRO,  InputTGLF}, filename::String)

Common save method for all the various input types
"""
function save(input::Union{InputCGYRO, InputQLGYRO,  InputTGLF}, filename::String)
    return TurbulentTransport.save(input, filename)
end

"""
    save(input::InputTJLF, filename::String)

Common save method for all the various input types
"""
function save(input::InputTJLF, filename::String)
    return TJLF.save(input, filename)
end