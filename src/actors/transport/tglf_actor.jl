import TGLFNN
import TGLFNN: InputTGLF
import TJLF
import TJLF: InputTJLF

#= ========= =#
#  ActorTGLF  #
#= ========= =#
Base.@kwdef mutable struct FUSEparameters__ActorTGLF{T<:Real} <: ParametersActor{T}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    _time::Float64 = NaN
    model::Switch{Symbol} = Switch{Symbol}([:TGLF, :TGLFNN, :TJLF], "-", "Implementation of TGLF"; default=:TGLFNN)
    sat_rule::Switch{Symbol} = Switch{Symbol}([:sat0, :sat0quench, :sat1, :sat1geo, :sat2, :sat3], "-", "Saturation rule"; default=:sat1)
    electromagnetic::Entry{Bool} = Entry{Bool}("-", "Electromagnetic or electrostatic"; default=true)
    tglfnn_model::Entry{String} = Entry{String}(
        "-",
        "Use a user specified TGLF-NN model stored in TGLFNN/models";
        check=x -> @assert x in TGLFNN.available_models() "ActorTGLF.tglfnn_model must be one of:\n  \"$(join(TGLFNN.available_models(),"\"\n  \""))\""
    )
    rho_transport::Entry{AbstractVector{T}} = Entry{AbstractVector{T}}("-", "rho_tor_norm values to compute tglf fluxes on"; default=0.25:0.1:0.85)
    warn_nn_train_bounds::Entry{Bool} = Entry{Bool}("-", "Raise warnings if querying cases that are certainly outside of the training range"; default=false)
    custom_input_files::Entry{Union{Vector{<:InputTGLF},Vector{<:InputTJLF}}} =
        Entry{Union{Vector{<:InputTGLF},Vector{<:InputTJLF}}}("-", "Sets up the input file that will be run with the custom input file as a mask")
    lump_ions::Entry{Bool} = Entry{Bool}("-", "Lumps the fuel species (D,T) as well as the impurities together"; default=true)
end

mutable struct ActorTGLF{D,P} <: SingleAbstractActor{D,P}
    dd::IMAS.dd{D}
    par::FUSEparameters__ActorTGLF{P}
    input_tglfs::Union{Vector{<:InputTGLF},Vector{<:InputTJLF}}
    flux_solutions::Union{Vector{<:IMAS.FluxSolution},Any}
end

"""
    ActorTGLF(dd::IMAS.dd, act::ParametersAllActors; kw...)

Evaluates the TGLF predicted turbulence
"""
function ActorTGLF(dd::IMAS.dd, act::ParametersAllActors; kw...)
    actor = ActorTGLF(dd, act.ActorTGLF; kw...)
    step(actor)
    finalize(actor)
    return actor
end

function ActorTGLF(dd::IMAS.dd, par::FUSEparameters__ActorTGLF; kw...)
    logging_actor_init(ActorTGLF)
    par = par(kw...)
    if par.model ∈ [:TGLF, :TGLFNN]
        input_tglfs = Vector{InputTGLF}(undef, length(par.rho_transport))
    elseif par.model == :TJLF
        input_tglfs = Vector{InputTJLF}(undef, length(par.rho_transport))
    end
    return ActorTGLF(dd, par, input_tglfs, IMAS.FluxSolution[])
end

"""
    _step(actor::ActorTGLF)

Runs TGLF actor to evaluate the turbulence flux on a vector of gridpoints
"""
function _step(actor::ActorTGLF)
    par = actor.par
    dd = actor.dd

    input_tglfs = InputTGLF(dd, par.rho_transport, par.sat_rule, par.electromagnetic, par.lump_ions)
    for k in eachindex(par.rho_transport)
        input_tglf = input_tglfs[k]
        if par.model ∈ [:TGLF, :TGLFNN]
            actor.input_tglfs[k] = input_tglf
        elseif par.model == :TJLF
            if !isassigned(actor.input_tglfs, k) # this is done to keep memory of the widths
                nky = TJLF.get_ky_spectrum_size(input_tglf.NKY, input_tglf.KYGRID_MODEL)
                actor.input_tglfs[k] = InputTJLF{Float64}(input_tglf.NS, nky)
                actor.input_tglfs[k].WIDTH_SPECTRUM .= 1.65
                actor.input_tglfs[k].FIND_WIDTH = true # first case should find the widths
            end
            update_input_tjlf!(actor.input_tglfs[k], input_tglf)
        end

        # Overwrite TGLF / TJLF parameters with the custom parameters mask
        if !ismissing(par, :custom_input_files)
            for field_name in fieldnames(typeof(actor.input_tglfs[k]))
                if !ismissing(getproperty(par.custom_input_files[k], field_name))
                    setproperty!(actor.input_tglfs[k], field_name, getproperty(par.custom_input_files[k], field_name))
                end
            end
        end
    end

    if par.model == :TGLFNN
        actor.flux_solutions = TGLFNN.run_tglfnn(actor.input_tglfs; par.warn_nn_train_bounds, model_filename=model_filename(par))

    elseif par.model == :TGLF
        actor.flux_solutions = TGLFNN.run_tglf(actor.input_tglfs)

    elseif par.model == :TJLF
        QL_fluxes_out = TJLF.run_tjlf(actor.input_tglfs)
        actor.flux_solutions =
            [IMAS.FluxSolution(TJLF.Qe(QL_flux_out), TJLF.Qi(QL_flux_out), TJLF.Γe(QL_flux_out), TJLF.Γi(QL_flux_out), TJLF.Πi(QL_flux_out)) for QL_flux_out in QL_fluxes_out]
    end

    return actor
end

"""
    _finalize(actor::ActorTGLF)

Writes results to dd.core_transport
"""
function _finalize(actor::ActorTGLF)
    dd = actor.dd
    par = actor.par
    cp1d = dd.core_profiles.profiles_1d[]
    eqt = dd.equilibrium.time_slice[]

    model = resize!(dd.core_transport.model, :anomalous; wipe=false)
    model.identifier.name = string(par.model) * " " * model_filename(par)
    m1d = resize!(model.profiles_1d)
    m1d.grid_flux.rho_tor_norm = par.rho_transport

    IMAS.flux_gacode_to_fuse((:electron_energy_flux, :ion_energy_flux, :electron_particle_flux, :ion_particle_flux, :momentum_flux), actor.flux_solutions, m1d, eqt, cp1d)

    return actor
end

function model_filename(par::FUSEparameters__ActorTGLF)
    if par.model == :TGLFNN
        filename = par.tglfnn_model
    else
        filename = string(par.sat_rule) * "_" * (par.electromagnetic ? "em" : "es")
    end
    return filename
end

"""
    update_input_tjlf!(input_tglf::InputTGLF)

Modifies an InputTJLF from a InputTGLF
"""
function update_input_tjlf!(input_tjlf::InputTJLF, input_tglf::InputTGLF)
    input_tjlf.NWIDTH = 21

    for fieldname in fieldnames(typeof(input_tglf))
        if occursin(r"\d", String(fieldname)) || fieldname == :_Qgb # species parameter
            continue
        end
        setfield!(input_tjlf, fieldname, getfield(input_tglf, fieldname))
    end

    for i in 1:input_tglf.NS
        input_tjlf.ZS[i] = getfield(input_tglf, Symbol("ZS_", i))
        input_tjlf.AS[i] = getfield(input_tglf, Symbol("AS_", i))
        input_tjlf.MASS[i] = getfield(input_tglf, Symbol("MASS_", i))
        input_tjlf.RLNS[i] = getfield(input_tglf, Symbol("RLNS_", i))
        input_tjlf.RLTS[i] = getfield(input_tglf, Symbol("RLTS_", i))
        input_tjlf.TAUS[i] = getfield(input_tglf, Symbol("TAUS_", i))
        input_tjlf.VPAR[i] = getfield(input_tglf, Symbol("VPAR_", i))
        input_tjlf.VPAR_SHEAR[i] = getfield(input_tglf, Symbol("VPAR_SHEAR_", i))
    end

    # Defaults
    input_tjlf.KY = 0.3
    input_tjlf.ALPHA_E = 1.0
    input_tjlf.ALPHA_P = 1.0
    input_tjlf.XNU_FACTOR = 1.0
    input_tjlf.DEBYE_FACTOR = 1.0
    input_tjlf.RLNP_CUTOFF = 18.0
    input_tjlf.WIDTH = 1.65
    input_tjlf.WIDTH_MIN = 0.3
    input_tjlf.BETA_LOC = 0.0
    input_tjlf.KX0_LOC = 1.0
    input_tjlf.PARK = 1.0
    input_tjlf.GHAT = 1.0
    input_tjlf.GCHAT = 1.0
    input_tjlf.WD_ZERO = 0.1
    input_tjlf.LINSKER_FACTOR = 0.0
    input_tjlf.GRADB_FACTOR = 0.0
    input_tjlf.FILTER = 2.0
    input_tjlf.THETA_TRAPPED = 0.7
    input_tjlf.ETG_FACTOR = 1.25
    input_tjlf.DAMP_PSI = 0.0
    input_tjlf.DAMP_SIG = 0.0

    input_tjlf.FIND_EIGEN = true
    input_tjlf.NXGRID = 16

    input_tjlf.ADIABATIC_ELEC = false
    input_tjlf.VPAR_MODEL = 0
    input_tjlf.NEW_EIKONAL = true
    input_tjlf.USE_BISECTION = true
    input_tjlf.USE_INBOARD_DETRAPPED = false
    input_tjlf.IFLUX = true
    input_tjlf.IBRANCH = -1
    input_tjlf.KX0_LOC = 0.0

    # for now settings
    input_tjlf.ALPHA_ZF = -1  # smooth   

    # check converison
    TJLF.checkInput(input_tjlf)

    return input_tjlf
end

function Base.show(io::IO, ::MIME"text/plain", input::Union{InputTGLF,InputTJLF})
    for field_name in fieldnames(typeof(input))
        println(io, " $field_name = $(getfield(input,field_name))")
    end
end
