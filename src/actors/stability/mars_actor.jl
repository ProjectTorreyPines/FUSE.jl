#import CHEASE: run_chease
using Interpolations
const μ_0 = 4pi * 1E-7

Base.@kwdef mutable struct CHEASEnamelist
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
    chease_exec::Entry{String} =
        Entry{String}("-", "Path to CHEASE executable"; default="chease.x")
    offset::Entry{Float64} = Entry{Float64}("-", "Offset for first wall (RW) in meters"; default=0.2)
    n_points::Entry{Int} = Entry{Int}("-", "Number of points for plasma boundary and surrounding walls "; default=301)
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
    chease_inputs::Union{Nothing,CHEASEnamelist}
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
        nl = chease_inputs === nothing ? CHEASEnamelist() : chease_inputs

        # -------------------------
        # Apply user overrides
        # -------------------------
        for (k, v) in pairs(namelist_overrides)
            if k ∉ fieldnames(CHEASEnamelist)
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
    nl = actor.chease_inputs
    profiles = dd.equilibrium.time_slice[].profiles_1d

    # Placeholder for MARS actor implementation
    # This would involve setting up the MARS simulation based on the parameters
    # and computing the wall heat flux accordingly

    #run_CHEASE(dd, par)
    if par.eq_type == :CHEASE # hardcode for now
        run_CHEASE(dd, par, nl)
    elseif par.eq_type == :TEQUILA
        # run TEQUILA equilibrium solver
    end
    
    # Produce the additional inputs required for MARS
    get_additional_MARS_inputs(profiles, par)
    @info "Running MARS actor with parameters: eq_type=$(par.eq_type), EQDSK=$(par.EQDSK), MHD_code=$(par.MHD_code), tracer_type=$(par.tracer_type), PEST_input=$(par.PEST_input)"
    run_MARS(dd, par)
    
    #run_PARTICLE_TRACING(dd, par)

    
    return actor
end


function run_CHEASE(dd::IMAS.dd, par, nl, time_slice_index::Int=1)
    @info "Running CHEASE with EQDSK=$(par.EQDSK)"

    chease_exec = par.chease_exec

    @assert nl !== nothing "CHEASE namelist not initialized"

    # Write EXPEQ file
    write_EXPEQ_file(dd, par)
    
    # Write CHEASE namelist file
    write_CHEASEnamelist(nl, "datain")

    # Execute CHEASE
    @assert isfile("datain") "CHEASE input file datain not found"
    @info "Executing CHEASE from $chease_exec"
    assert_executable(chease_exec)

    isfile(chease_exec) || error("CHEASE executable not found: $chease_exec")

    cmd = pipeline(
        `$(chease_exec)`,
        stdin  = "datain",
        stdout = "log_chease",
        stderr = "log_chease"
    )

    @info "Executing CHEASE from $chease_exec"
    ok = success(cmd)
    keys = ["GEXP", "Q_ZERO", "Q_EDGE"]
    extract_lines_for_keys("log_chease", keys) if ok
    ok || error("CHEASE failed — see log_chease")


    return nothing
end

function get_additional_MARS_inputs(profiles, par)
    # Placeholder function to generate additional inputs for MARS
    # This would involve preparing files or data structures needed by MARS
    @info "Generating additional MARS inputs based on parameters."
    tor_rotation = profiles.profiles_1d.rotation_frequency_tor_sonic
    pressure = profiles.profiles_1d.pressure
end

function run_MARS(dd::IMAS.dd, par)
    # Placeholder function to run MARS MHD stability code
    @info "Running MARS with MHD_code=$(par.MHD_code) and PEST_input=$(par.PEST_input)."
end

function assert_executable(path::AbstractString)
    isfile(path) || error("CHEASE executable not found: $path")
    isexecutable(path) || error("CHEASE exists but is not executable: $path")
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

    offset = par.offset  # offset for first wall (RW) in meters
    n_points = par.n_points  # number of points for first wall (RW)

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
    pprime = eqt1d.dpressure_dpsi
    rho_pol = sqrt.(eqt1d.psi_norm)

    wall_RZ = [dd.wall.description_2d[time_slice_index].limiter.unit[1].outline.r, dd.wall.description_2d[time_slice_index].limiter.unit[1].outline.z]

    if minimum(r_bound) - offset < 0
        error("Offset too large: boundary crosses R < 0 (min R = $(minimum(r_bound)))")
    end

    ## get additional parameters from user
    NWBPS = par.number_surfaces
    if par.GS_rhs == :TTpr
        NSTTP = 1
        j_tor = eqt1d.j_tor
    elseif par.GS_rhs == :Jtor
        NSTTP = 2
        j_tor = eqt1d.j_parallel
    elseif par.GS_rhs == :Jpar
        NSTTP = 3
        j_tor = eqt1d.j_parallel # NOT right!
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
    dpressure_ds = 2*pprime.*rho_pol / (Bt_center^2 / μ_0)

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
    ab = sqrt((maximum(r_bound) - minimum(r_bound))^2 + (maximum(z_bound) - minimum(z_bound))^2) / 2.0
    pr, pz = limit_curvature(r_bound, z_bound, ab / 20.0)
    rb_new, zb_new = IMAS.resample_2d_path(pr, pz; n_points=n_points, method=:linear)
    println(length(rb_new),length(pr))
    r_bound_norm = rb_new / r_center
    z_bound_norm = zb_new / r_center

    plt = plot()
    plt = plot!(r_bound_norm, z_bound_norm; linewidth=3., aspect_ratio=:equal, title="Smoothed Boundary & RW for CHEASE")
    display(plt)


    # add a smooth first wall (RW)
    r_lim, z_lim = offset_boundary(rb_new, zb_new, offset)
    r_lim, z_lim = IMAS.resample_2d_path(r_lim, z_lim; n_points=n_points, method=:linear)
    r_lim_norm = r_lim / r_center
    z_lim_norm = z_lim / r_center
    plt = plot!(r_lim_norm, z_lim_norm; linewidth=1.5, aspect_ratio=:equal, title="Smoothed Boundary for EXPEQ")
    display(plt)

    ##----------------- Write the file -----------------##
    write_list = [string(ϵ), string(z_axis), string(pressure_sep_norm)]
    @assert length(rb_new) == length(zb_new) "R,Z boundary arrays must have the same shape"
    write_list = vcat(write_list, string(length(rb_new), " ", NWBPS, " ", NDATA))
    for (r, z) in zip(r_bound_norm, z_bound_norm)
        write_list = vcat(write_list, "$r    $z")
    end
    if NWBPS > 1. ## WHAT TO DO if > 2
        for (r, z) in zip(r_lim_norm, z_lim_norm)
            write_list = vcat(write_list, "$r    $z")
        end
    end

    @assert length(rho_pol) == length(pressure) == length(j_tor) "rho_pol, presssure and j_tor arrays must have the same shape"
    write_list = vcat(write_list, "$(length(pressure))")
    write_list = vcat(write_list, "$(string(NSTTP))")
    write_list = vcat(write_list, map(string, rho_pol))
    write_list = vcat(write_list, map(string, dpressure_ds))
    write_list = vcat(write_list, map(string, j_tor_norm))

    # write to EXPEQ file   
    touch("EXPEQ")
    open("EXPEQ", "w") do file
        for line in write_list
            write(file, "$line \n")
        end
    end
end

function write_CHEASEnamelist(nl::CHEASEnamelist, filename::AbstractString="datain")
    open(filename, "w") do io
        println(io, "***")
        println(io, "***    Example Torus")
        println(io, "***")
        println(io, "***")
        println(io, "&EQDATA")

        for field in fieldnames(CHEASEnamelist)
            val = getfield(nl, field)

            if val isa Vector
                println(io, "  $(field)(1) = ", join(val, ", "), ",")
            else
                println(io, "  $(field) = ", val, ",")
            end
        end

        println(io, "&END")
    end

    return filename
end


function offset_boundary(xs, ys, d)
    ds = sqrt.(diff(xs).^2 .+ diff(ys).^2)
    t  = cumsum(vcat(0.0, ds))
    s  = t ./ maximum(t)

    itp_raw = linear_interpolation(s, xs)

    # resample onto uniform grid
    s_uniform = range(0.0, 1.0; length=length(s))
    r_uniform = itp_raw.(s_uniform)

    # now cubic spline
    itp_x = CubicSplineInterpolation(s_uniform, r_uniform)

    itp_raw = linear_interpolation(s, ys)

    # resample onto uniform grid
    s_uniform = range(0.0, 1.0; length=length(s))
    z_uniform = itp_raw.(s_uniform)

    # now cubic spline
    itp_y = CubicSplineInterpolation(s_uniform, z_uniform)

    npts = 200
    ss   = range(0,1,length=npts)

    X = Float64[]
    Y = Float64[]

    # Evaluate curve
    xv = itp_x.(ss)
    yv = itp_y.(ss)

    # Gradients (extract scalar from SVector{1})
    dxds = getindex.(Interpolations.gradient.(Ref(itp_x), ss), 1)
    dyds = getindex.(Interpolations.gradient.(Ref(itp_y), ss), 1)

    # Tangents
    normT = sqrt.(dxds.^2 .+ dyds.^2)
    Tx = dxds ./ normT
    Ty = dyds ./ normT

    # Normals
    Nx = -Ty
    Ny =  Tx

    # Offset curve
    X = xv .+ d .* Nx
    Y = yv .+ d .* Ny


    return X, Y
end
"""
    extract_lines_for_keys(logfile, keys) -> Dict

Scan `logfile` once and return a dictionary mapping each key
to the first line containing it.

Errors if any key is not found.
"""
function extract_lines_for_keys(
    logfile::AbstractString,
    keys::AbstractVector{<:AbstractString}
)
    isfile(logfile) || error("Log file not found: $logfile")

    remaining = Set(keys)
    results   = Dict{String,String}()

    for line in eachline(logfile)
        for key in remaining
            if occursin(key, line)
                results[key] = line
                delete!(remaining, key)
            end
        end
        isempty(remaining) && break
    end

    isempty(remaining) || error("Keys not found in $logfile: $(collect(remaining))")

    return results
end


