import TGLFNN: InputTGLF, run_tglf, run_tglfnn
import TJLF: InputTJLF, run_tjlf, get_ky_spectrum_size, checkInput
#= ========= =#
#  ActorTGLF  #
#= ========= =#
Base.@kwdef mutable struct FUSEparameters__ActorTGLF{T} <: ParametersActor where {T<:Real}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    model::Switch{Symbol} = Switch{Symbol}([:TGLF, :TGLFNN, :TJLF], "-", "Implementation of TGLF"; default=:TGLFNN)
    sat_rule::Switch{Symbol} = Switch{Symbol}([:sat0, :sat0quench, :sat1, :sat1geo, :sat2], "-", "Saturation rule"; default=:sat1)
    electromagnetic::Entry{Bool} = Entry{Bool}("-", "Electromagnetic or electrostatic"; default=true)
    user_specified_model::Entry{String} = Entry{String}("-", "Use a user specified TGLF-NN model stored in TGLFNN/models"; default="")
    rho_transport::Entry{AbstractVector{T}} = Entry{AbstractVector{T}}("-", "rho_tor_norm values to compute tglf fluxes on"; default=0.25:0.1:0.85)
    warn_nn_train_bounds::Entry{Bool} = Entry{Bool}("-", "Raise warnings if querying cases that are certainly outside of the training range"; default=false)
    custom_input_files::Entry{Union{Vector{InputTGLF}, Vector{InputTJLF{Float64}}, Bool}}  = Entry{ Union{Vector{InputTGLF}, Vector{InputTJLF{Float64}}, Bool}}("-", "Sets up the input file that will be run with the custom input file as a mask";  default=false)
end

mutable struct ActorTGLF{D,P} <: PlasmaAbstractActor{D,P}
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
    theta_0 = (0.7-0.2)/(ϵ_D3D- ϵ_st40) * (ϵ - ϵ_st40) + 0.2
    theta_1 = (0.7-0.8)/(ϵ_D3D- ϵ_st40) * (ϵ - ϵ_st40) + 0.8
    theta_trapped = range(theta_0,theta_1,length(cp1d.grid.rho_tor_norm))

    for (k, (gridpoint_eq, gridpoint_cp)) in enumerate(zip(ix_eq, ix_cp))
        if par.model ∈ [:TGLF, :TGLFNN]
            actor.input_tglfs[k] = InputTGLF(dd, gridpoint_eq, gridpoint_cp, par.sat_rule, par.electromagnetic)
            if par.model == :TGLFNN
                # TGLF-NN has some difficulty with the sign of rotation / shear
                actor.input_tglfs[k].VPAR_SHEAR_1 = abs(actor.input_tglfs[k].VPAR_SHEAR_1)
                actor.input_tglfs[k].VPAR_1 = abs(actor.input_tglfs[k].VPAR_1)
            end
        elseif par.model == :TJLF
            actor.input_tglfs[k] = create_input_tjlf(InputTGLF(dd, gridpoint_eq, gridpoint_cp, par.sat_rule, par.electromagnetic))
        end

        if ϵ > ϵ_D3D
            actor.input_tglfs[k].THETA_TRAPPED = theta_trapped[gridpoint_cp]
        end

        # Setting up the TJLF / TGLF run with the custom parameter mask (this overwrites all the above)
        if par.custom_input_files != false
            for k in 1:length(actor.input_tglfs)
                for field_name in fieldnames(typeof(actor.input_tglfs[k]))
                    if !ismissing(getproperty(par.custom_input_files[k], field_name))
                        setproperty!(actor.input_tglfs[k], field_name, getproperty(par.custom_input_files[k], field_name))
                    end
                end
            end
        end
    end

    if par.model == :TGLFNN
        actor.flux_solutions = run_tglfnn(actor.input_tglfs; par.warn_nn_train_bounds, model_filename=model_filename(par))
    elseif par.model == :TGLF
        actor.flux_solutions = run_tglf(actor.input_tglfs)
    elseif par.model == :TJLF
        QL_fluxes_out = run_tjlf(actor.input_tglfs) 
        actor.flux_solutions = [IMAS.flux_solution(single[1,1,1],single[1,2,3],single[1,1,2],single[1,2,2]) for single in QL_fluxes_out]
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
    create_input_tjlf(inputTGLF::InputTGLF)

Creates an InputTJLF from a InputTGLF
"""
function create_input_tjlf(inputTGLF::InputTGLF)
    nky = get_ky_spectrum_size(inputTGLF.NKY, inputTGLF.KYGRID_MODEL)
    inputTJLF = InputTJLF{Float64}(inputTGLF.NS, nky)
    inputTJLF.NWIDTH = 21
    
    for fieldname in fieldnames(typeof(inputTGLF))
        if occursin(r"\d",String(fieldname)) || fieldname==:_Qgb # species parameter
            continue
        end
        setfield!(inputTJLF,fieldname,getfield(inputTGLF,fieldname))
    end
    for i in 1:inputTGLF.NS
        inputTJLF.ZS[i] = getfield(inputTGLF,Symbol("ZS_",i))
        inputTJLF.AS[i] = getfield(inputTGLF,Symbol("AS_",i))
        inputTJLF.MASS[i] = getfield(inputTGLF,Symbol("MASS_",i))
        inputTJLF.RLNS[i] = getfield(inputTGLF,Symbol("RLNS_",i))
        inputTJLF.RLTS[i] = getfield(inputTGLF,Symbol("RLTS_",i))
        inputTJLF.TAUS[i] = getfield(inputTGLF,Symbol("TAUS_",i))
        inputTJLF.VPAR[i] = getfield(inputTGLF,Symbol("VPAR_",i))
        inputTJLF.VPAR_SHEAR[i] = getfield(inputTGLF,Symbol("VPAR_SHEAR_",i))
    end
    inputTJLF.WIDTH_SPECTRUM .= 0.0

    # Defaults
    inputTJLF.KY= 0.3
    inputTJLF.ALPHA_E = 1.0
    inputTJLF.ALPHA_P= 1.0
    inputTJLF.XNU_FACTOR = 1.0
    inputTJLF.DEBYE_FACTOR = 1.0
    inputTJLF.RLNP_CUTOFF = 18.0
    inputTJLF.WIDTH = 1.65
    inputTJLF.WIDTH_MIN = 0.3
    inputTJLF.BETA_LOC = 0.0
    inputTJLF.KX0_LOC = 1.0
    inputTJLF.PARK = 1.0
    inputTJLF.GHAT = 1.0
    inputTJLF.GCHAT = 1.0
    inputTJLF.WD_ZERO = 0.1
    inputTJLF.LINSKER_FACTOR = 0.0
    inputTJLF.GRADB_FACTOR = 0.0
    inputTJLF.FILTER = 2.0
    inputTJLF.THETA_TRAPPED = 0.7
    inputTJLF.ETG_FACTOR = 1.25
    inputTJLF.DAMP_PSI = 0.0
    inputTJLF.DAMP_SIG = 0.0
   
    
    inputTJLF.FIND_WIDTH = true # first case should find the widths
    inputTJLF.FIND_EIGEN = true
    inputTJLF.NXGRID = 16

    inputTJLF.ADIABATIC_ELEC = false
    inputTJLF.VPAR_MODEL = 0
    inputTJLF.NEW_EIKONAL = true
    inputTJLF.USE_BISECTION= true
    inputTJLF.USE_INBOARD_DETRAPPED = false
    inputTJLF.IFLUX = true
    inputTJLF.IBRANCH = -1 
    inputTJLF.WIDTH_SPECTRUM .= 0.0
    inputTJLF.KX0_LOC = 0.0
    
    # for now settings
    inputTJLF.ALPHA_ZF = -1  # smooth   

    checkInput(inputTJLF)
    # check converison
    field_names = fieldnames(InputTJLF)
    for field_name in field_names
        field_value = getfield(inputTJLF, field_name)
        
         if typeof(field_value)<:Missing || typeof(field_value)<:Real
             @assert !ismissing(field_value) || !isnan(field_value) "Did not properly populate inputTJLF for $field_name"
         end

         if typeof(field_value) <: Vector && field_name != :KY_SPECTRUM && field_name != :EIGEN_SPECTRUM && field_name != :EIGEN_SPECTRUM2
            for val in field_value
                @assert !isnan(val) "Did not properly populate inputTJLF for array $field_name"
            end
        end
    end

    return inputTJLF
end


function Base.show(input::Union{InputTGLF,InputTJLF})
    for field_name in fieldnames(typeof(input))
        println(" $field_name = $(getfield(input,field_name))")
    end
end