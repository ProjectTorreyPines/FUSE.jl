import TGLFNN
import TGLFNN: InputQLGYRO, InputCGYRO

#= ========= =#
#  ActorQLGYRO  #
#= ========= =#
Base.@kwdef mutable struct FUSEparameters__ActorQLGYRO{T<:Real} <: ParametersActorPlasma{T}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    _time::Float64 = NaN
    sat_rule::Switch{Symbol} = Switch{Symbol}([:sat1, :sat1geo, :sat2, :sat3], "-", "Saturation rule"; default=:sat1)
    nfields::Entry{Int} = Entry{Int}("-", "Electromagnetic or electrostatic"; default=1)
    dt:::Entry{Int} = Entry{Float64}("-", "step size "; default=0.005)
    rho_transport::Entry{AbstractVector{T}} = Entry{AbstractVector{T}}("-", "rho_tor_norm values to compute QLGYRO fluxes on"; default=0.25:0.1:0.85)
    lump_ions::Entry{Bool} = Entry{Bool}("-", "Lumps the fuel species (D,T) as well as the impurities together"; default=true)
end

mutable struct ActorQLGYRO{D,P} <: SingleAbstractActor{D,P}
    dd::IMAS.dd{D}
    par::FUSEparameters__ActorQLGYRO{P}
    input_qlgyros::Vector{<:InputQLGYRO}
    input_cgyros::Vector{<:InputCGYRO}
    flux_solutions::Union{Vector{<:IMAS.flux_solution},Any}
end

"""
    ActorQLGYRO(dd::IMAS.dd, act::ParametersAllActors; kw...)

Evaluates the QLGYRO predicted turbulence
"""
function ActorQLGYRO(dd::IMAS.dd, act::ParametersAllActors; kw...)
    actor = ActorQLGYRO(dd, act.ActorQLGYRO; kw...)
    step(actor)
    finalize(actor)
    return actor
end

function ActorQLGYRO(dd::IMAS.dd, par::FUSEparameters__ActorQLGYRO; kw...)
    par = par(kw...)
    if par.model ∈ [:QLGYRO]
        input_QLGYROs = Vector{InputQLGYRO}(undef, length(par.rho_transport))

    return ActorQLGYRO(dd, par, input_QLGYROs, IMAS.flux_solution[])
end

"""
    _step(actor::ActorQLGYRO)

Runs QLGYRO actor to evaluate the turbulence flux on a vector of gridpoints
"""
function _step(actor::ActorQLGYRO)
    par = actor.par
    dd = actor.dd

    eq1d = dd.equilibrium.time_slice[].profiles_1d
    cp1d = dd.core_profiles.profiles_1d[]
    ix_eq = [argmin(abs.(eq1d.rho_tor_norm .- rho)) for rho in par.rho_transport]
    ix_cp = [argmin(abs.(cp1d.grid.rho_tor_norm .- rho)) for rho in par.rho_transport]


    for (k, (gridpoint_eq, gridpoint_cp)) in enumerate(zip(ix_eq, ix_cp))
        input_qlgyro = TGLFNN.InputQLGYRO()

        input_qlgyro.N_PARALLEL = par.nky
        input_qlgyro.N_RUNS = 1
        input_qlgyro.GAMMA_E = 0.0
        input_qlgyro.CODE= -1
        input_qlgyro.NKY = par.nky * par.cpu_per_ky
        input_qlgyro.KYGRID_MODEL= par.kygrid_model
        input_qlgyro.KY= par.ky
        input_qlgyro.SAT_RULE = par.sat_rule
        input_QLGYRO = InputQLGYRO(dd, gridpoint_cp, par.electromagnetic, par.lump_ions)
        input_CGYRO = InputCGYRO(dd)

        # Setting up the QLGYRO run with the custom parameter mask (this overwrites all the above)
        if !ismissing(par, :custom_input_files)
            for field_name in fieldnames(typeof(actor.input_QLGYROs[k]))
                if !ismissing(getproperty(par.custom_input_files[k], field_name))
                    setproperty!(actor.input_QLGYROs[k], field_name, getproperty(par.custom_input_files[k], field_name))
                end
            end
        end
    end

    actor.flux_solutions = TGLFNN.run_QLGYRO(actor.input_QLGYROs,actor.input_CGYROs)

    return actor
end

"""
    _finalize(actor::ActorQLGYRO)

Writes results to dd.core_transport
"""
function _finalize(actor::ActorQLGYRO)
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

function model_filename(par::FUSEparameters__ActorQLGYRO)
    if !isempty(par.user_specified_model)
        filename = par.user_specified_model
    else
        filename = string(par.sat_rule) * "_" * (par.electromagnetic ? "em" : "es")
        filename *= "_d3d" # will be changed to FPP soon
    end
    return filename
end

function Base.show(input::InputQLGYRO)
    for field_name in fieldnames(typeof(input))
        println(" $field_name = $(getfield(input,field_name))")
    end
end