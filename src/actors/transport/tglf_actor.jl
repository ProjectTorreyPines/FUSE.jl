import TGLFNN
import TGLFNN: InputTGLF
import TJLF
import TJLF: InputTJLF

#= ========= =#
#  ActorTGLF  #
#= ========= =#
Base.@kwdef mutable struct FUSEparameters__ActorTGLF{T<:Real} <: ParametersActorPlasma{T}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    _time::Float64 = NaN
    model::Switch{Symbol} = Switch{Symbol}([:TGLF, :TGLFNN, :TJLF], "-", "Implementation of TGLF"; default=:TGLFNN)
    sat_rule::Switch{Symbol} = Switch{Symbol}([:sat0, :sat0quench, :sat1, :sat1geo, :sat2, :sat3], "-", "Saturation rule"; default=:sat1)
    electromagnetic::Entry{Bool} = Entry{Bool}("-", "Electromagnetic or electrostatic"; default=true)
    user_specified_model::Entry{String} = Entry{String}("-", "Use a user specified TGLF-NN model stored in TGLFNN/models"; default="")
    rho_transport::Entry{AbstractVector{T}} = Entry{AbstractVector{T}}("-", "rho_tor_norm values to compute tglf fluxes on"; default=0.25:0.1:0.85)
    warn_nn_train_bounds::Entry{Bool} = Entry{Bool}("-", "Raise warnings if querying cases that are certainly outside of the training range"; default=false)
    custom_input_files::Entry{Union{Vector{<:InputTGLF},Vector{<:InputTJLF}}} =
        Entry{Union{Vector{<:InputTGLF},Vector{<:InputTJLF}}}("-", "Sets up the input file that will be run with the custom input file as a mask")
    lump_ions::Entry{Bool} = Entry{Bool}("-", "Lumps the fuel species (D,T) as well as the impurities together"; default=true)
    onnx_model::Entry{Bool} = Entry{Bool}("-", "use onnx model"; default=false)
end

mutable struct ActorTGLF{D,P} <: SingleAbstractActor{D,P}
    dd::IMAS.dd{D}
    par::FUSEparameters__ActorTGLF{P}
    input_tglfs::Union{Vector{<:InputTGLF},Vector{<:InputTJLF}}
    flux_solutions::Union{Vector{<:IMAS.flux_solution},Any}
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
    par = par(kw...)
    if par.model ∈ [:TGLF, :TGLFNN]
        input_tglfs = Vector{InputTGLF}(undef, length(par.rho_transport))
    elseif par.model == :TJLF
        input_tglfs = Vector{InputTJLF}(undef, length(par.rho_transport))
    end
    return ActorTGLF(dd, par, input_tglfs, IMAS.flux_solution[])
end

"""
    _step(actor::ActorTGLF)

Runs TGLF actor to evaluate the turbulence flux on a vector of gridpoints
"""
function _step(actor::ActorTGLF)
    par = actor.par
    dd = actor.dd

    eq1d = dd.equilibrium.time_slice[].profiles_1d
    cp1d = dd.core_profiles.profiles_1d[]
    ix_eq = [argmin(abs.(eq1d.rho_tor_norm .- rho)) for rho in par.rho_transport]
    ix_cp = [argmin(abs.(cp1d.grid.rho_tor_norm .- rho)) for rho in par.rho_transport]


    ϵ_st40 = 1 / 1.9
    ϵ_D3D = 0.67 / 1.67

    ϵ = dd.equilibrium.time_slice[].boundary.minor_radius / dd.equilibrium.time_slice[].boundary.geometric_axis.r
    theta_0 = (0.7 - 0.2) / (ϵ_D3D - ϵ_st40) * (ϵ - ϵ_st40) + 0.2
    theta_1 = (0.7 - 0.8) / (ϵ_D3D - ϵ_st40) * (ϵ - ϵ_st40) + 0.8
    theta_trapped = range(theta_0, theta_1, length(cp1d.grid.rho_tor_norm))

    for (k, gridpoint_cp) in enumerate(ix_cp)
        input_tglf = InputTGLF(dd, gridpoint_cp, par.sat_rule, par.electromagnetic, par.lump_ions)
        if par.model ∈ [:TGLF, :TGLFNN]
            actor.input_tglfs[k] = input_tglf
            #if par.model == :TGLFNN
            #    # TGLF-NN has some difficulty with the sign of rotation / shear
            #    actor.input_tglfs[k].VPAR_SHEAR_1 = abs(actor.input_tglfs[k].VPAR_SHEAR_1)
            #    actor.input_tglfs[k].VPAR_1 = abs(actor.input_tglfs[k].VPAR_1)
            #end
        elseif par.model == :TJLF
            if !isassigned(actor.input_tglfs, k)
                nky = TJLF.get_ky_spectrum_size(input_tglf.NKY, input_tglf.KYGRID_MODEL)
                actor.input_tglfs[k] = InputTJLF{Float64}(input_tglf.NS, nky)
                actor.input_tglfs[k].WIDTH_SPECTRUM .= 1.65
                actor.input_tglfs[k].FIND_WIDTH = true # first case should find the widths
            end
            update_input_tjlf!(actor.input_tglfs[k], input_tglf)
        end

        #if ϵ > ϵ_D3D
        #    actor.input_tglfs[k].THETA_TRAPPED = theta_trapped[gridpoint_cp]
        #end

        # Setting up the TJLF / TGLF run with the custom parameter mask (this overwrites all the above)
        if !ismissing(par, :custom_input_files)
            for field_name in fieldnames(typeof(actor.input_tglfs[k]))
                if !ismissing(getproperty(par.custom_input_files[k], field_name))
                    setproperty!(actor.input_tglfs[k], field_name, getproperty(par.custom_input_files[k], field_name))
                end
            end
        end
    end

    if par.model == :TGLFNN
        #actor.flux_solutions = TGLFNN.run_tglfnn(actor.input_tglfs; par.warn_nn_train_bounds, model_filename=model_filename(par))
        if par.onnx_model == false
            actor.flux_solutions = TGLFNN.run_tglfnn(actor.input_tglfs; par.warn_nn_train_bounds, model_filename=model_filename(par))
        elseif par.onnx_model == true
            actor.flux_solutions = TGLFNN.run_tglfnn_onnx(actor.input_tglfs, par.user_specified_model, [
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
            ];)
        end
    elseif par.model == :TGLF
        actor.flux_solutions = TGLFNN.run_tglf(actor.input_tglfs)
    elseif par.model == :TJLF
        QL_fluxes_out = TJLF.run_tjlf(actor.input_tglfs)
        actor.flux_solutions = [IMAS.flux_solution(TJLF.Qe(QL_flux_out), TJLF.Qi(QL_flux_out), TJLF.Γe(QL_flux_out), TJLF.Γi(QL_flux_out), TJLF.Πi(QL_flux_out)) for QL_flux_out in QL_fluxes_out]
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
    #model.identifier.name = string(par.model) * " " * model_filename(par)
    if par.onnx_model == false
        model.identifier.name = string(par.model) * " " * model_filename(par)
    elseif par.onnx_model == true
        model.identifier.name = string(par.model) * " " * "onnx model"
    end
    m1d = resize!(model.profiles_1d)
    m1d.grid_flux.rho_tor_norm = par.rho_transport

    IMAS.flux_gacode_to_fuse([:ion_energy_flux, :electron_energy_flux, :electron_particle_flux, :momentum_flux], actor.flux_solutions, m1d, eqt, cp1d)

    return actor
end

function model_filename(par::FUSEparameters__ActorTGLF)
    if !isempty(par.user_specified_model)
        filename = par.user_specified_model
    else
        filename = string(par.sat_rule) * "_" * (par.electromagnetic ? "em" : "es")
        filename *= "_d3d" # will be changed to FPP soon
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

    TJLF.checkInput(input_tjlf)
    # check converison
    field_names = fieldnames(InputTJLF)
    for field_name in field_names
        field_value = getfield(input_tjlf, field_name)

        if typeof(field_value) <: Missing || typeof(field_value) <: Real
            @assert !ismissing(field_value) || !isnan(field_value) "Did not properly populate input_tjlf for $field_name"
        end

        if typeof(field_value) <: Vector && field_name != :KY_SPECTRUM && field_name != :EIGEN_SPECTRUM && field_name != :EIGEN_SPECTRUM2
            for val in field_value
                @assert !isnan(val) "Did not properly populate input_tjlf for array $field_name"
            end
        end
    end

    return input_tjlf
end


function Base.show(input::Union{InputTGLF,InputTJLF})
    for field_name in fieldnames(typeof(input))
        println(" $field_name = $(getfield(input,field_name))")
    end
end