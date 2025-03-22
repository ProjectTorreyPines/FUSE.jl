import TGLFNN
import TGLFNN: InputQLGYRO, InputCGYRO
import GACODE

#= =========== =#
#  ActorQLGYRO  #
#= =========== =#
Base.@kwdef mutable struct FUSEparameters__ActorQLGYRO{T<:Real} <: ParametersActor{T}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    _time::Float64 = NaN
    model::Switch{Symbol} = Switch{Symbol}([:QLGYRO], "-", "Implementation of QLGYRO"; default=:QLGYRO)
    ky::Entry{Float64} = Entry{Float64}("-", "Max ky"; default=1.6)
    nky::Entry{Int} = Entry{Int}("-", "Number of ky modes"; default=16)
    cpu_per_ky::Entry{Int} = Entry{Int}("-", "Number of cpus per ky"; default=1)
    kygrid_model::Entry{Int} = Entry{Int}("-", "TGLF ky grid model"; default=0)
    sat_rule::Switch{Symbol} = Switch{Symbol}([:sat1, :sat2, :sat3], "-", "Saturation rule"; default=:sat1)
    n_field::Entry{Int} = Entry{Int}( "-", "1:phi, 2:phi+apar, 3:phi+apar+bpar"; default=1, check=x -> @assert 1 <= x <= 3 "n_fields must be either 1,2 or 3")
    delta_t::Entry{Float64} = Entry{Float64}("-", "CGYRO step size "; default=0.005)
    max_time::Entry{Float64} = Entry{Float64}("-", "Max simulation time (a/cs)"; default=100.0)
    rho_transport::Entry{AbstractVector{T}} = Entry{AbstractVector{T}}("-", "rho_tor_norm values to compute QLGYRO fluxes on"; default=0.25:0.1:0.85)
    lump_ions::Entry{Bool} = Entry{Bool}("-", "Lumps the fuel species (D,T) as well as the impurities together"; default=true)
end

mutable struct ActorQLGYRO{D,P} <: SingleAbstractActor{D,P}
    dd::IMAS.dd{D}
    par::FUSEparameters__ActorQLGYRO{P}
    input_qlgyros::Vector{<:InputQLGYRO}
    input_cgyros::Vector{<:InputCGYRO}
    flux_solutions::Union{Vector{<:GACODE.FluxSolution},Any}
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
    logging_actor_init(ActorQLGYRO)
    par = par(kw...)
    input_QLGYROs = Vector{InputQLGYRO}(undef, length(par.rho_transport))
    input_CGYROs = Vector{InputCGYRO}(undef, length(par.rho_transport))
    return ActorQLGYRO(dd, par, input_QLGYROs, input_CGYROs, GACODE.FluxSolution[])
end

"""
    _step(actor::ActorQLGYRO)

Runs QLGYRO actor to evaluate the turbulence flux on a vector of gridpoints
"""
function _step(actor::ActorQLGYRO)
    par = actor.par
    dd = actor.dd

    cp1d = dd.core_profiles.profiles_1d[]
    ix_cp = [argmin(abs.(cp1d.grid.rho_tor_norm .- rho)) for rho in par.rho_transport]

    for (k, gridpoint_cp) in enumerate(ix_cp)
        actor.input_qlgyros[k] = TGLFNN.InputQLGYRO()

        actor.input_qlgyros[k].N_PARALLEL = par.nky * par.cpu_per_ky
        actor.input_qlgyros[k].N_RUNS = 1
        actor.input_qlgyros[k].GAMMA_E = 0.0
        actor.input_qlgyros[k].CODE = -1
        actor.input_qlgyros[k].NKY = par.nky
        actor.input_qlgyros[k].KYGRID_MODEL = par.kygrid_model
        actor.input_qlgyros[k].KY = par.ky

        if par.sat_rule == :sat2
            actor.input_qlgyros[k].SAT_RULE = 2
        elseif par.sat_rule == :sat1
            actor.input_qlgyros[k].SAT_RULE = 1
        end

        actor.input_cgyros[k] = InputCGYRO(dd, gridpoint_cp, par.lump_ions)
        actor.input_cgyros[k].N_FIELD = par.n_field
        actor.input_cgyros[k].DELTA_T = par.delta_t
        actor.input_cgyros[k].MAX_TIME = par.max_time

	actor.input_cgyros[k].DELTA_T_METHOD = 1
	actor.input_cgyros[k].FREQ_TOL = 0.01
    end

    actor.flux_solutions = TGLFNN.run_qlgyro(actor.input_qlgyros, actor.input_cgyros)

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
    model.identifier.name = string(par.model)
    m1d = resize!(model.profiles_1d)
    m1d.grid_flux.rho_tor_norm = par.rho_transport

    GACODE.flux_gacode_to_imas((:electron_energy_flux, :ion_energy_flux,  :electron_particle_flux, :ion_particle_flux, :momentum_flux), actor.flux_solutions, m1d, eqt, cp1d)

    return actor
end

function Base.show(io::IO, ::MIME"text/plain", input::InputQLGYRO)
    for field_name in fieldnames(typeof(input))
        println(io, " $field_name = $(getfield(input,field_name))")
    end
end
