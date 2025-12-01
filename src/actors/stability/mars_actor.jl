import CHEASE

Base.@kwdef mutable struct FUSEparameters__ActorMars{T<:Real} <: ParametersActor{T}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    _time::Float64 = NaN
    do_plot::Entry{Bool} = act_common_parameters(; do_plot=false)
    eq_type::Switch{Symbol} = Switch{Symbol}([:EFIT, :FUSE], "-", "Type of equilibrium to use: :EFIT or from FUSE"; default=:EFIT)
    EQDSK::Entry{Bool} = Entry{Bool}("-", "Enable EQDSK"; default=false)
    MHD_code::Switch{Symbol} = Switch{Symbol}([:MARS_F, :MARS_K, :M3D-C1], "-", "MHD code to use: :MARS or :MARS_F"; default=:MARS)
    tracer_type::Switch{Symbol} = Switch{Symbol}([:ORBIT, :REORBIT], "-", "Type of tracer to use: :ideal or :realistic"; default=:REORBIT)
    PEST_input::Entry{Bool} = Entry{Bool}("-", "Use PEST input files"; default=false)   
end

mutable struct ActorMars{D,P} <: SingleAbstractActor{D,P}
    dd::IMAS.dd{D}
    par::OverrideParameters{P,FUSEparameters__ActorMars{P}}
    actor_eq::ActorEquilibrium{D,P}
    wall_heat_flux::Union{Nothing,IMAS.WallHeatFlux}

    function ActorMars(dd::IMAS.dd{D}, par::FUSEparameters__ActorMars{P}; kw...) where {D<:Real,P<:Real}
        logging_actor_init(ActorMars)
        par = OverrideParameters(par; kw...)
        return new{D,P}(dd, par, nothing, nothing)
    end
end

"""
    ActorMars(dd::IMAS.dd, act::ParametersAllActors; kw...) 

"""
function ActorMars(dd::IMAS.dd, act::ParametersAllActors; kw...)
    actor = ActorMars(dd, act.ActorMars; kw...)
    step(actor)
    finalize(actor)
    return actor
end


function _step(actor::ActorMars)
    dd = actor.dd
    par = actor.par

    # Placeholder for MARS actor implementation
    # This would involve setting up the MARS simulation based on the parameters
    # and computing the wall heat flux accordingly

    
    #run_CHEASE(dd, par)
    actor.actor_eq.par.model = :CHEASE # hardcode for now
    #actor_eq = ActorEquilibrium(dd, act.ActorEquilibrium, act; ip_from=:core_profiles)
    #actor.eq_actor = ActorCHEASE(dd, act.ActorCHEASE, act)
    
    # Produce the additional inputs required for MARS
    get_additional_MARS_inputs(dd, par)
    @info "Running MARS actor with parameters: eq_type=$(par.eq_type), EQDSK=$(par.EQDSK), MHD_code=$(par.MHD_code), tracer_type=$(par.tracer_type), PEST_input=$(par.PEST_input)"
    run_MARS(dd, par)
    
    run_PARTICLE_TRACING(dd, par)

    # Here would be the implementation details for interfacing with MARS and computing the wall heat flux

    # For now, we just set wall_heat_flux to nothing
    actor.wall_heat_flux = nothing
end


function run_CHEASE(dd::IMAS.dd, par)
    # Placeholder function to run CHEASE equilibrium solver
    @info "Running CHEASE with EQDSK=$(par.EQDSK)"

    CHEASE_struct = run_chease(ϵ, z_axis, pressure_sep, Bt_center, r_center, Ip, r_bound, z_bound, mode, rho_psi, pressure, j_tor, clear_workdir=false)
end

function get_additional_MARS_inputs(dd::IMAS.dd, par)
    # Placeholder function to generate additional inputs for MARS
    # This would involve preparing files or data structures needed by MARS
    @info "Generating additional MARS inputs based on parameters."
end

function run_PARTICLE_TRACING(dd::IMAS.dd, par)
    # Placeholder function to run particle tracing simulations
    @info "Running particle tracing with tracer_type=$(par.tracer_type)."

    println("Particle tracing simulation completed.")
end