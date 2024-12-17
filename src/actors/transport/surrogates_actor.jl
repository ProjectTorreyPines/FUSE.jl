using Surrogates
using AbstractGPs #required to access different types of kernels
using SurrogatesAbstractGPs

#= =========== =#
#  ActorSurrogate  #
#= =========== =#
Base.@kwdef mutable struct FUSEparameters__ActorSurrogate{T<:Real} <: ParametersActor{T}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    _time::Float64 = NaN
    generate_surrogates::Entry{Bool} = Entry{Bool}("-", "Generate surrogates"; default=false)
end

mutable struct ActorSurrogate{D,P} <: SingleAbstractActor{D,P}
    dd::IMAS.dd{D}
    par::FUSEparameters__ActorSurrogate{P}
    all_model_inputs::Union{Vector{<:InputTGLF},Any}
    all_model_outputs::Union{Vector{<:IMAS.flux_solution},Any,Any}
    flux_surrogates::Union{Vector{<:Any},Any,Any}
    flux_solutions::Union{Vector{<:IMAS.flux_solution},Any}
end

"""
    ActorSurrogate(dd::IMAS.dd, act::ParametersAllActors; kw...)

Evaluates the QLGYRO predicted turbulence
"""
function ActorSurrogate(dd::IMAS.dd, act::ParametersAllActors; kw...)
    actor = ActorSurrogate(dd, act.ActorSurrogate; kw...)
    step(actor)
    finalize(actor)
    return actor
end


struct flux_surrogate{T<:Any}
    ENERGY_FLUX_e::T
    ENERGY_FLUX_i::T
    PARTICLE_FLUX_e::T
    PARTICLE_FLUX_i::Vector{T}
    STRESS_TOR_i::T
end


function generate_surrogates(all_input_tglfs, all_flux_solutions)

    actor.all_input_tglfs
    actor.flux_solutions
    nsims = length(flux_solutions)
    nrads = length(flux_solutions[1])
    flux_surrogates = Vector{flux_surrogate}()
    for irad in 1:nrad
        x = inputtglf_to_surrogate_input(all_inputs_tglf,irad)
        y = [flux_solutions[isim][irad].ENERGY_FLUX_e for isim in 1:nsims]
        ENERGY_FLUX_e = AbstractGPSurrogate(x,y, gp=GP(Matern32Kernel()), Σy=0.0001) 

        y = [flux_solutions[isim][irad].ENERGY_FLUX_i for isim in 1:nsims]
        ENERGY_FLUX_i = AbstractGPSurrogate(x,y, gp=GP(Matern32Kernel()), Σy=0.0001) 

        y = [flux_solutions[isim][irad].PARTICLE_FLUX_e for isim in 1:nsims]
        PARTICLE_FLUX_e = AbstractGPSurrogate(x,y, gp=GP(Matern32Kernel()), Σy=0.0001) 

        PARTICLE_FLUX_i = []

        y = [flux_solutions[isim][irad].STRESS_TOR_i for isim in 1:nsims]
        STRESS_TOR_i = AbstractGPSurrogate(x,y, gp=GP(Matern32Kernel()), Σy=0.0001) 

        push!(flux_surrogates, flux_surrogate(ENERGY_FLUX_e,ENERGY_FLUX_i,PARTICLE_FLUX_e,PARTICLE_FLUX_i,STRESS_TOR_i))
    end
    return flux_surrogates
end function

function inputtglf_to_surrogate_input(all_input_tglf,irad)
    nsims = length(all_input_tglfs)
    v = Array{Tuple{Vararg{Float64, 8}}}()

    for isim in 1:nsims
        push!(v,  (input_tglf.BETAE,
        input_tglf.P_PRIME_LOC,
        input_tglf.RLNS_1,
        input_tglf.RLTS_1,
        input_tglf.RLTS_2,
        input_tglf.VEXB_SHEAR,
        input_tglf.XNUE,
        input_tglf.ZEFF))
    end
    return v

end function

"""
    _step(actor::ActorSurrogate)

Runs QLGYRO actor to evaluate the turbulence flux on a vector of gridpoints
"""

function _step(actor::ActorSurrogate)
    par = actor.par
    dd = actor.dd

    if par.generate_surrogates
         generate_surrogates(par.all_model_inputs, par.all_model_outputs)
    end
    if par.run_surrogates

        actor.flux_solutions =  Vector{flux_solution}()
        for irad in 1:nrad
            all_input_tglf[end]
            x = inputtglf_to_surrogate_input(par.all_model_inputs,irad)[end]
            push!(actor.flux_solutions,flux_surrogates(x))
        end
    end
    return actor
end

"""
    _finalize(actor::ActorSurrogate)

Writes results to dd.core_transport
"""
function _finalize(actor::ActorSurrogate)
    dd = actor.dd
    par = actor.par
    cp1d = dd.core_profiles.profiles_1d[]
    eqt = dd.equilibrium.time_slice[]

    model = resize!(dd.core_transport.model, :anomalous; wipe=false)
    model.identifier.name = string(par.model)
    m1d = resize!(model.profiles_1d)
    m1d.grid_flux.rho_tor_norm = par.rho_transport

    IMAS.flux_gacode_to_fuse((:electron_energy_flux, :ion_energy_flux,  :electron_particle_flux, :ion_particle_flux, :momentum_flux), actor.flux_solutions, m1d, eqt, cp1d)

    return actor
end

function Base.show(io::IO, ::MIME"text/plain", input::InputQLGYRO)
    for field_name in fieldnames(typeof(input))
        println(io, " $field_name = $(getfield(input,field_name))")
    end
end


struct flux_surrogate{T<:Any}
    ENERGY_FLUX_e::T
    ENERGY_FLUX_i::T
    PARTICLE_FLUX_e::T
    PARTICLE_FLUX_i::Vector{T}
    STRESS_TOR_i::T
end
