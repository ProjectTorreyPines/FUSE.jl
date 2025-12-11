#import CHEASE: run_chease
const μ_0 = 4pi * 1E-7

Base.@kwdef mutable struct MarsEqNamelist
    NEQDSK::Int     = 0
    NSURF::Int      = 6
    NTCASE::Int     = 0

    NBLOPT::Int     = 0
    NBSOPT::Int     = 0
    CPRESS::Float64 = 1.000
    CFBAL::Float64  = 3.0000

    NCSCAL::Int     = 2
    CSSPEC::Float64 = 0.000
    QSPEC::Float64  = 1.6185

    NTMF0::Int      = 0
    CURRT::Float64  = 0.3000

    NSTTP::Int      = 2
    NFUNC::Int      = 4
    NIPR::Int       = 1
    NISO::Int       = 100

    NPPFUN::Int     = 4
    NPP::Int        = 1
    NPPR::Int       = 30

    NSOUR::Int      = 2
    NPROPT::Int     = 2

    NS::Int         = 60
    NT::Int         = 60
    NPSI::Int       = 240
    NCHI::Int       = 200

    NV::Int         = 160
    REXT::Float64   = 6.0
    NVEXP::Int      = 8
    R0W::Float64    = 0.90
    RZ0W::Float64   = 0.0

    NMESHA::Int     = 2
    SOLPDA::Float64 = 0.60
    QWIDTH0::Float64 = 0.30
    ROTE::Float64   = 0.0000
    NTOR::Int       = 1

    NPOIDQ::Int     = 6
    QSHAVE::Float64 = 100.0

    QPLACE::Vector{Float64} =
        [2.0000, 3.0000, 4.0000, 5.0000, 6.0000, 7.0000]

    QWIDTH::Vector{Float64} =
        [0.0009, 0.0007, 0.0006, 0.0006, 0.0005, 0.0005]

    NEGP::Int       = -1
    NER::Int        = 1

    EPSLON::Float64 = 1.0e-10
    GAMMA::Float64  = 1.6666666667

    MSMAX::Int      = 40
    NINMAP::Int     = 50
    NINSCA::Int     = 50

    NOPT::Int       = 0
    NPLOT::Int      = 1
    NBAL::Int       = 0

    B0EXP::Float64  = 1.5
    R0EXP::Float64  = 3.0
end


struct MarsEq
    OutRVAR::Vector{Float64}
    OutPVAR::Vector{Float64}
    # Placeholder for MARS equilibrium data structure
end

struct MarsInput
    something::Float64
    # Placeholder for MARS input data structure
end

Base.@kwdef mutable struct FUSEparameters__ActorMars{T<:Real} <: ParametersActor{T}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    _time::Float64 = NaN
    do_plot::Entry{Bool} = act_common_parameters(; do_plot=false)
    eq_type::Switch{Symbol} = Switch{Symbol}([:CHEASE, :TEQUILA], "-", "Type of equilibrium to use: :CHEASE or :TEQUILA"; default=:CHEASE)
    EQDSK::Entry{Bool} = Entry{Bool}("-", "Enable EQDSK"; default=false)
    MHD_code::Switch{Symbol} = Switch{Symbol}([:MARS_Q, :MARS_F, :MARS_K], "-", "MHD code to use: :MARS or :MARS_F"; default=:MARS_F)
    tracer_type::Switch{Symbol} = Switch{Symbol}([:ORBIT, :REORBIT], "-", "Type of tracer to use: :ideal or :realistic"; default=:REORBIT)
    PEST_input::Entry{Bool} = Entry{Bool}("-", "Use PEST input files"; default=false)
    number_surfaces::Entry{Int} = Entry{Int}("-", "Number of surfaces to specify"; default=1)
    pressure_sep::Entry{Float64} = Entry{Float64}("-", "Pressure at separatrix in Pa"; default=0.0)
    GS_rhs::Switch{Symbol} = Switch{Symbol}([:TTpr, :Jtor, :Jpar], "-", "Specification of Grad-Shaf RHS current"; default=:TTpr)
    wall_resistivity_type::Switch{Symbol} = Switch{Symbol}([:Constant, :Variable], "-", "Wall Resistivity Model"; default=:Constant)    
end


mutable struct ActorMars{D,P} <: SingleAbstractActor{D,P}
    dd::IMAS.dd{D}
    par::OverrideParameters{P,FUSEparameters__ActorMars{P}}
    chease_inputs::Union{Nothing,MarsEqNamelist}
    mars_inputs::Union{Nothing,Vector{MarsInput}}

    function ActorMars(
        dd::IMAS.dd{D}, 
        par::FUSEparameters__ActorMars{P}; 
        namelist_overrides = NamedTuple(),
        chease_inputs = nothing,
        mars_inputs = nothing,
        kw...
    ) where {D<:Real,P<:Real}

        logging_actor_init(ActorMars)

        # Apply OverrideParameters FIRST (FUSE convention)
        par = OverrideParameters(par; kw...)

        # -------------------------
        # Create default CHEASE namelist
        # -------------------------
        nl = chease_inputs === nothing ? MarsEqNamelist() : chease_inputs

        # -------------------------
        # Apply user overrides
        # -------------------------
        for (k, v) in pairs(namelist_overrides)
            if k ∉ fieldnames(MarsEqNamelist)
                error("Unknown namelist entry '$k'")
            end
            setfield!(nl, k, v)
        end

        # -------------------------
        # Final actor construction
        # -------------------------
        return new{D,P}(dd, par, nl, mars_inputs)
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
    if par.eq_type == :CHEASE # hardcode for now
        run_CHEASE(dd, par)
    elseif par.eq_type == :TEQUILA
        # run TEQUILA equilibrium solver
    end
    
    # Produce the additional inputs required for MARS
    get_additional_MARS_inputs(dd, par)
    @info "Running MARS actor with parameters: eq_type=$(par.eq_type), EQDSK=$(par.EQDSK), MHD_code=$(par.MHD_code), tracer_type=$(par.tracer_type), PEST_input=$(par.PEST_input)"
    run_MARS(dd, par)
    
    #run_PARTICLE_TRACING(dd, par)

    # For now, we just set wall_heat_flux to nothing
    #actor.wall_heat_flux = nothing

    return actor
end


function run_CHEASE(dd::IMAS.dd, par, time_slice_index::Int=1)
    # Placeholder function to run CHEASE equilibrium solver
    @info "Running CHEASE with EQDSK=$(par.EQDSK)"
    limiter_RZ = [dd.wall.description_2d[time_slice_index].limiter.unit[1].outline.r, dd.wall.description_2d[time_slice_index].limiter.unit[1].outline.z]
    #write_EXPEQ_file(dd.equilibrium.time_slice[time_slice_index], par, limiter_RZ)
    write_EXPEQ_file(dd, par)
    
    ## Execute CHEASE
    #run(pipeline(`$(executable)`; stdout=io, stderr=io))
    
end

function get_additional_MARS_inputs(dd::IMAS.dd, par)
    # Placeholder function to generate additional inputs for MARS
    # This would involve preparing files or data structures needed by MARS
    @info "Generating additional MARS inputs based on parameters."
end

function run_MARS(dd::IMAS.dd, par)
    # Placeholder function to run MARS MHD stability code
    @info "Running MARS with MHD_code=$(par.MHD_code) and PEST_input=$(par.PEST_input)."
end

function run_PARTICLE_TRACING(dd::IMAS.dd, par)
    # Placeholder function to run particle tracing simulations
    @info "Running particle tracing with tracer_type=$(par.tracer_type)."

    println("Particle tracing simulation completed.")
end


"""
  Writes a EXPEQ file for CHEASE given dd and mode number
    This function will create the necessary input file for the CHEASE equilibrium solver.
        NWBPS: Number of wall boundary points
        NSTTP: Number of steps in pressure and current profiles
"""
function write_EXPEQ_file(dd::IMAS.dd, par, time_slice_index::Int=1)

    # initialize eqt from pulse_schedule and core_profiles
    time_slice = dd.equilibrium.time_slice[time_slice_index]
    eqt1d = time_slice.profiles_1d


    # populate the input file lines
    minor_radius = time_slice.boundary.minor_radius
    z_axis = time_slice.global_quantities.magnetic_axis.z
    Bt_center = time_slice.global_quantities.magnetic_axis.b_field_tor
    r_center = time_slice.global_quantities.magnetic_axis.r
    Ip = time_slice.global_quantities.ip
    r_bound = time_slice.boundary.outline.r
    z_bound = time_slice.boundary.outline.z
    pressure = eqt1d.pressure
    j_tor = eqt1d.j_tor
    rho_pol = eqt1d.rho_tor_norm

    wall_RZ = [dd.wall.description_2d[time_slice_index].limiter.unit[1].outline.r, dd.wall.description_2d[time_slice_index].limiter.unit[1].outline.z]

    ## get additional parameters from user
    NWBPS = par.number_surfaces
    NSTTP = if par.GS_rhs == :TTpr
        1
    elseif par.GS_rhs == :Jtor
        2
    elseif par.GS_rhs == :Jpar
        3
    else
        0
    end

    if par.wall_resistivity_type == :Constant
        NDATA = 2
        # set wall resistivity model to constant
    elseif par.wall_resistivity_type == :Variable
        NDATA = 3
        # set wall resistivity model to variable
    else
        NDATA = 1   
    end
    pressure_sep = par.pressure_sep

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

    # Remove/smooth the X-point along the boundary
    pr = r_bound
    pz = z_bound
    ab = sqrt((maximum(pr) - minimum(pr))^2 + (maximum(pz) - minimum(pz))^2) / 2.0
    pr, pz = limit_curvature(pr, pz, ab / 20.0)
    pr, pz = IMAS.resample_2d_path(pr, pz; n_points=301, method=:linear)

    plt = plot!(pr, pz; marker=:circle, aspect_ratio=:equal, title="Smoothed Boundary for EXPEQ")
    display(plt)

    r_bound_norm = pr / r_center
    z_bound_norm = pz / r_center

    
    # add the limiter/vacuum vessel outline
    # first smooth the outline
    # Remove/smooth the X-point along the boundary
    pr2 = wall_RZ[1]
    pz2 = wall_RZ[2]
    ab = sqrt((maximum(pr2) - minimum(pr2))^2 + (maximum(pz2) - minimum(pz2))^2) / 2.0
    println(ab)
    #pr2, pz2 = limit_curvature(pr2, pz2, ab/20.)
    pr2, pz2 = IMAS.resample_2d_path(pr2, pz2; n_points=301, method=:linear)
    
    
    ##----------------- Write the file -----------------##
    write_list = [string(ϵ), string(z_axis), string(pressure_sep_norm)]
    @assert length(r_bound) == length(z_bound) "R,Z boundary arrays must have the same shape"
    write_list = vcat(write_list, string(length(r_bound), " ", NWBPS, " ", NDATA))
    for (r, z) in zip(r_bound_norm, z_bound_norm)
        write_list = vcat(write_list, "$r    $z")
    end
    if NWBPS > 1. ## WHAT TO DO if > 2
        write_list = vcat(write_list, string(length(pr2)))
        for (r, z) in zip(pr2, pz2)
            write_list = vcat(write_list, "$r    $z")
        end
    end

    @assert length(rho_pol) == length(pressure) == length(j_tor) "rho_pol, presssure and j_tor arrays must have the same shape"
    write_list = vcat(write_list, "$(length(pressure))    $(string(NSTTP))")
    write_list = vcat(write_list, map(string, rho_pol))
    write_list = vcat(write_list, map(string, pressure_norm))
    write_list = vcat(write_list, map(string, j_tor_norm))

   
    # write to EXPEQ file   
    touch("EXPEQ")
    open("EXPEQ", "w") do file
        for line in write_list
            write(file, "$line \n")
        end
    end
end

"""
    write_chease_namelist(
        Bt_center::Float64,
        r_center::Float64,
        Ip::Float64,
        r_bound::Vector{Float64},
        z_bound::Vector{Float64};
        rescale_eq_to_ip::Bool=false,
        extra_box_fraction::Float64=0.33)

Writes the chease namelist to the current folder
"""
function write_chease_namelist(
    Bt_center::Float64,
    r_center::Float64,
    Ip::Float64,
    r_bound::Vector{Float64},
    z_bound::Vector{Float64};
    rescale_eq_to_ip::Bool=false,
    extra_box_fraction::Float64=0.33)

    eqdata = Dict{Symbol,Any}()
    eqdata[:R0EXP] = r_center
    eqdata[:B0EXP] = Bt_center
    eqdata[:CURRT] = abs(Ip / (r_center * Bt_center / μ_0))
    eqdata[:SIGNB0XP] = sign(Bt_center)
    eqdata[:SIGNIPXP] = sign(Ip)
    eqdata[:NT] = 80 # number of theta points
    eqdata[:COCOS_IN] = 11
    if rescale_eq_to_ip
        eqdata[:NCSCAL] = 2 
    else
        eqdata[:NCSCAL] = 4 
    end 
    eqdata[:NPROPT] = -2
    eqdata[:NPPFUN] = 8 
    eqdata[:EPSLON] = 1e-6 # convergence
    eqdata[:RELAX] = 0.9 


    touch("datain")

end


