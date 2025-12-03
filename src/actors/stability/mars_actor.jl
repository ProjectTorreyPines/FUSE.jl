#import CHEASE: run_chease
#const μ_0 = 4pi * 1E-7

struct MarsEq
    OutRVAR::Vector{Float64}
    OutPVAR::Vector{Float64}
    # Placeholder for MARS equilibrium data structure
end

struct MarsInput
    # Placeholder for MARS input data structure
end


Base.@kwdef mutable struct FUSEparameters__ActorMars{T<:Real} <: ParametersActor{T}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    _time::Float64 = NaN
    do_plot::Entry{Bool} = act_common_parameters(; do_plot=false)
    eq_type::Switch{Symbol} = Switch{Symbol}([:CHEASE, :TEQUILA], "-", "Type of equilibrium to use: :CHEASE or :TEQUILA"; default=:CHEASE)
    EQDSK::Entry{Bool} = Entry{Bool}("-", "Enable EQDSK"; default=false)
    MHD_code::Switch{Symbol} = Switch{Symbol}([:MARS, :MARS_F, :MARS_K], "-", "MHD code to use: :MARS or :MARS_F"; default=:MARS)
    tracer_type::Switch{Symbol} = Switch{Symbol}([:ORBIT, :REORBIT], "-", "Type of tracer to use: :ideal or :realistic"; default=:REORBIT)
    PEST_input::Entry{Bool} = Entry{Bool}("-", "Use PEST input files"; default=false)   
end

mutable struct ActorMars{D,P} <: SingleAbstractActor{D,P}
    dd::IMAS.dd{D}
    par::OverrideParameters{P,FUSEparameters__ActorMars{P}}
    #actor_eq::ActorEquilibrium{D,P}
    wall_heat_flux::Union{Nothing,IMAS.WallHeatFlux}
    mars_equilibrium::Union{Nothing,MarsEq}
    mars_inputs::Union{Nothing,Vector{MarsInput}}
    """
        ActorMars(dd::IMAS.dd{D}, par::FUSEparameters__ActorMars{P}; kw...) where {D<:Real,P<:Real}

    """

    function ActorMars(dd::IMAS.dd{D}, par::FUSEparameters__ActorMars{P}, actor_eq, wall_heat_flux, mars_equilibrium, mars_inputs; kw...) where {D<:Real,P<:Real}
        logging_actor_init(ActorMars)
        par = OverrideParameters(par; kw...)
        return new{D,P}(dd, par, actor_eq, wall_heat_flux, mars_equilibrium, mars_inputs)
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
    if actor.eq_type == :CHEASE # hardcode for now
        run_CHEASE(dd, par)
    elseif actor.eq_type == :TEQUILA
        # run TEQUILA equilibrium solver
    end
    #actor_eq = ActorEquilibrium(dd, act.ActorEquilibrium, act; ip_from=:core_profiles)
    #actor.eq_actor = ActorCHEASE(dd, act.ActorCHEASE, act)
    
    # Produce the additional inputs required for MARS
    get_additional_MARS_inputs(dd, par)
    @info "Running MARS actor with parameters: eq_type=$(par.eq_type), EQDSK=$(par.EQDSK), MHD_code=$(par.MHD_code), tracer_type=$(par.tracer_type), PEST_input=$(par.PEST_input)"
    run_MARS(dd, par)
    
    #run_PARTICLE_TRACING(dd, par)

    
    # For now, we just set wall_heat_flux to nothing
    actor.wall_heat_flux = nothing
end


function run_CHEASE(dd::IMAS.dd, par, time_slice_index::Int=1)
    # Placeholder function to run CHEASE equilibrium solver
    @info "Running CHEASE with EQDSK=$(par.EQDSK)"
    pressure_sep = 0.0  # Placeholder value
    mode = 0            # Placeholder value
    #write_EXPEQ_file(dd.equilibrium.time_slice[time_slice_index], pressure_sep, mode)
    #CHEASE_struct = run_chease(ϵ, z_axis, pressure_sep, Bt_center, r_center, Ip, r_bound, z_bound, mode, rho_psi, pressure, j_tor, clear_workdir=false)
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


"""
  Writes a EXPEQ file for CHEASE given dd and mode number
"""
function write_EXPEQ_file(time_slice, pressure_sep::Float64, mode::Int)

    # populate the input file lines
    minor_radius = time_slice.boundary.minor_radius
    z_axis = time_slice.global_quantities.magnetic_axis.z
    Bt_center = time_slice.global_quantities.magnetic_axis.b_field_tor
    r_center = time_slice.global_quantities.magnetic_axis.r
    Ip = time_slice.global_quantities.ip
    pressure = time_slice.profiles_1d.pressure
    j_tor = time_slice.profiles_1d.j_tor

    # calculate aspect ratio
    ϵ = minor_radius / r_center

    # Normalize from SI to chease units
    pressure_sep_norm = pressure_sep / (Bt_center^2 / μ_0)
    pressure_norm = pressure / (Bt_center^2 / μ_0)
    j_tor_norm = abs.(j_tor / (Bt_center / (r_center * μ_0)))

    ip_sign = sign(Ip)
    bt_sign = sign(Bt_center)
    if ip_sign == 1 && bt_sign == 1
        j_tor_norm .*= 1
    elseif ip_sign == 1 && bt_sign == -1
        j_tor_norm .*= 1
    elseif ip_sign == -1 && bt_sign == -1
        j_tor_norm .*= 1
    else
        j_tor_norm .*= -1
    end

    r_bound_norm = r_bound / r_center
    z_bound_norm = z_bound / r_center

    write_list = [string(ϵ), string(z_axis), string(pressure_sep_norm)]
    @assert length(r_bound) == length(z_bound) "R,Z boundary arrays must have the same shape"
    write_list = vcat(write_list, string(length(r_bound)))
    for (r, z) in zip(r_bound_norm, z_bound_norm)
        write_list = vcat(write_list, "$r    $z")
    end
    @assert length(rho_pol) == length(pressure) == length(j_tor) "rho_pol, presssure and j_tor arrays must have the same shape"
    write_list = vcat(write_list, "$(length(pressure))    $(string(mode)[1])")
    write_list = vcat(write_list, "$(string(mode)[2])    0")
    write_list = vcat(write_list, map(string, rho_pol))
    write_list = vcat(write_list, map(string, pressure_norm))
    write_list = vcat(write_list, map(string, j_tor_norm))

    touch("EXPEQ")
    open("EXPEQ", "w") do file
        for line in write_list
            write(file, "$line \n")
        end
    end
end


