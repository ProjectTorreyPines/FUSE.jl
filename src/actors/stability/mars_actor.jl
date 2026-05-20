using Interpolations
import CHEASE

const μ_0 = 4pi * 1E-7

Base.@kwdef mutable struct CHEASEnamelist
    NEQDSK::Int     = 0
    NSURF::Int      = 6
    NTCASE::Int     = 0

    NBLOPT::Int     = 0
    NBSOPT::Int     = 0
    CPRESS::Float64 = 1.000
    CFBAL::Float64  = 1.0000 # set to 1. if NSCAL = 4

    NCSCAL::Int     = 2   # set to 4 if NOT scale q
    CSSPEC::Float64 = 0.000
    QSPEC::Float64  = 1.6185

    NTMF0::Int      = 0
    CURRT::Float64  = 0.3000

    NSTTP::Int      = 2
    NFUNC::Int      = 4
    NIPR::Int       = 1
    NISO::Int       = 100
    NIDEAL::Int     = 0

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


Base.@kwdef mutable struct MARS_BASIC
    NCASE::Int = 1
    TALPHA1::ComplexF64 = 1.e-4 + 0im
    NITMAX::Int = 100
    NSWEEP::Int = 1
    GAMMA::Float64 = 1.66667
    RNTOR::Float64 = -1.0
    M1::Int = -9
    M2::Int = 33
    NPROFN::Int = 4
    NPROFIE::Int = 0
    NPROFT::Int = 2
    ETA::Float64 = 0.0
    NV::Int = 160
    NWALL::Int = 1
    TAUW::Float64 = 1.0e4
    IWALL::Int = 22
    IWO::Int = 1
    IVISC::Int = 1
    PVISC::Float64 = 0.1
    NPROFR::Int = 0
    ROTE::Float64 = 0.0
    NPROFWE::Int = 0
    ROTWE0::Float64 = 0.0
    NPROFTTCA::Int = 0
    NPROFTTCE::Int = 0
    FRACPTH::Float64 = 1.0
    extras::Dict{Symbol,Any} = Dict{Symbol,Any}()
end

Base.@kwdef mutable struct MARS_FEEDBACK
    INCFEED::Int = 0
    NCOIL::Int = 2
    FCCHI::ComplexF64 = 0.34219 - 0.405470im
    FWCHI::ComplexF64 = 0.12187 + 0.089063im
    SCCHI::Float64 = 0.0
    SWCHI::Float64 = 0.05
    BTCHI::Float64 = 0.05
    FEEDI::Vector{ComplexF64} = ComplexF64[0. + 0.0im, 0. + 0.0im]
    IFEED::Int = 3
    ISENS::Int = 3
    NSENS::Int = 0
    KKF::Int = -1
    KREADECA::Int = 0
    extras::Dict{Symbol,Any} = Dict{Symbol,Any}()
end

Base.@kwdef mutable struct MARS_KINETIC
    INCKIN::Int = 0
    ATAU::ComplexF64 = 1.0e-4 + 0.0im
    ALTAU::Float64 = 0.5
    ALPHAP::Float64 = 0.5
    ALPHAD::Float64 = 1.0
    OMEGACI0::Float64 = 38.056
    IPERTURB::Int = 2
    KFASTRUN::Int = 0
    KENORM::Int = 2
    NSPECIES::Int = 2
    IFOWPSI0::Int = 1
    ISPECIES_F0::Vector{Int} = [0, 0, 3]
    KNBI::Int = 0
    DZETA0::Float64 = 1.0
    S0TYPE4::Float64 = 0.5
    R0TYPE4::Float64 = 1.5
    HHTYPE4C::Float64 = 0.5
    ESPECIES_M::Vector{Float64} = [2.0, 5.4463e-4, 2.0]
    ESPECIES_Z::Vector{Float64} = [1.0, -1.0, 1.0]
    PSPECIES_AP::Vector{Float64} = [0.0, 0.0, 1.0]
    PSPECIES_AT::Vector{Float64} = [0.0, 0.0, 0.0]
    PSPECIES_NP::Vector{Float64} = [0.0, 0.0, 0.0]
    PSPECIES_NTB::Vector{Float64} = [0.0, 0.0, 0.0]
    PSPECIES_NTD::Vector{Float64} = [0.0, 0.0, 0.0]
    PSPECIES_FOWP::Vector{Float64} = [0.0, 0.0, 1.0]
    PSPECIES_FOWT::Vector{Float64} = [0.0, 0.0, 0.0]
    NPROFK::Int = 0
    INUTYPE::Int = 0
    NPROFUI::Int = 0
    NUEFFIA::Float64 = 0.0
    NPROFUE::Int = 0
    NUEFFEA::Float64 = 0.0
    extras::Dict{Symbol,Any} = Dict{Symbol,Any}()
end

Base.@kwdef mutable struct MARS_QLIN
    TALPHA2::Float64 = 1.0
    TALPHA3::Float64 = 1.0
    TALPHA6::Float64 = 1.0
    TALPHA8::Float64 = 0.9
    TDELTAMIN::Float64 = 0.04
    TDELTAMAX::Float64 = 0.08
    TDELTALIM::Float64 = 100.0
    TDELTALOW::Float64 = 2.0
    CTJXB::Float64 = 1.0
    CTNTV::Float64 = 1.0
    CTREY::Float64 = 1.0
    ITSATURAT::Int = 2
    NPROFVD::Int = 3
    TCHIM0::Float64 = 8.4571e-7
    NPROFVP::Int = 0
    TVPINCH0::Float64 = 0.0
    CTEDGE::Float64 = 0.995
    KSOLSAVE::Int = 0
    KSOLREAD::Int = 0
    extras::Dict{Symbol,Any} = Dict{Symbol,Any}()
end

Base.@kwdef mutable struct MARS_NUMERIC
    PTRAPI::Float64 = 1.0
    PTRAPH::Float64 = 0.9
    NRES::Int = 7
    NFIT::Int = 3
    NTORQ::Int = 5
    IPDIVB::Int = 1
    NCONVB1::Int = 0
    NCONVCS::Int = 1
    IOMPNUM::Int = 1
    NVACJ::Int = 0
    V2XKEY::Int = 2
    extras::Dict{Symbol,Any} = Dict{Symbol,Any}()
end

Base.@kwdef mutable struct MARS_OUTOPT
    JSOUT::Int = 7
    ORMIN::Float64 = 0.0
    ORMAX::Float64 = 2.0
    OZMIN::Float64 = -2.0
    OZMAX::Float64 = 2.0
    NORR::Int = 0
    NOZZ::Int = 0
    extras::Dict{Symbol,Any} = Dict{Symbol,Any}()
end

Base.@kwdef mutable struct MARSnamelist
    BASIC::MARS_BASIC = MARS_BASIC()
    FEEDBACK::MARS_FEEDBACK = MARS_FEEDBACK()
    KINETIC::MARS_KINETIC = MARS_KINETIC()
    QLIN::MARS_QLIN = MARS_QLIN()
    NUMERIC::MARS_NUMERIC = MARS_NUMERIC()
    OUTOPT::MARS_OUTOPT = MARS_OUTOPT()
end

Base.@kwdef mutable struct MarsOverrides
    BASIC::Dict{Symbol,Any} = Dict{Symbol,Any}()
    FEEDBACK::Dict{Symbol,Any} = Dict{Symbol,Any}()
    KINETIC::Dict{Symbol,Any} = Dict{Symbol,Any}()
    QLIN::Dict{Symbol,Any} = Dict{Symbol,Any}()
    NUMERIC::Dict{Symbol,Any} = Dict{Symbol,Any}()
    OUTOPT::Dict{Symbol,Any} = Dict{Symbol,Any}()
end


struct MarsEq
    OutRVAR::Vector{Float64}
    OutPVAR::Vector{Float64}
    # Placeholder for MARS equilibrium data structure
end

Base.@kwdef mutable struct FUSEparameters__ActorMars{T<:Real} <: ParametersActor{T}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    _time::Float64 = NaN
    do_plot::Entry{Bool} = act_common_parameters(; do_plot=false)
    EQDSK::Entry{Bool} = Entry{Bool}("-", "Enable EQDSK"; default=false)
    chease_exec::Entry{String} =
        Entry{String}("-", "Path to CHEASE executable"; default="/fusion/projects/codes/mars/CHEASE/chease.x")
    mars_exec::Entry{String} =
        Entry{String}("-", "Path to MARS executable"; default="/fusion/projects/codes/mars/MARSQ/marsq.x")
    offset::Entry{Float64} = Entry{Float64}("-", "Offset for conforming first wall (RW) in units of length"; default=0.2)
    n_points::Entry{Int} = Entry{Int}("-", "Number of points for discretizing plasma boundary and surrounding walls "; default=301)
    tracer_type::Switch{Symbol} = Switch{Symbol}([:ORBIT, :REORBIT], "-", "Type of tracer to use: :ideal or :realistic"; default=:REORBIT)
    GS_rhs::Switch{Symbol} = Switch{Symbol}([:FFpr, :Jtor, :Jpar], "-", "Specification of Grad-Shaf RHS current"; default=:FFpr)
    wall_resistivity_type::Switch{Symbol} = Switch{Symbol}([:Constant, :Variable], "-", "Wall Resistivity Model"; default=:Constant)    
    wall_type::Switch{Symbol} = Switch{Symbol}([:no_wall, :conformal, :limiter], "-", "Machine wall shape to use for MARS"; default=:no_wall)
    number_surfaces::Entry{Int} = Entry{Int}("-", "Number of surfaces to specify"; default=1)
    run_equilibrium::Entry{Bool} = Entry{Bool}("-", "Whether to run equilibrium solver"; default=true)  
    restart_equilibrium::Entry{Bool} = Entry{Bool}("-", "Whether to restart from existing equilibrium"; default=false)
    run_MHD::Entry{Bool} = Entry{Bool}("-", "Whether to run MHD stability code"; default=true)  
    run_mode::Switch{Symbol} = Switch{Symbol}([:local, :batch], "-", "Whether to run MARS locally or submit to batch system"; default=:local)
    num_orbits::Entry{Int} = Entry{Int}("-", "Number of orbits to simulate in particle tracing"; default=0)  
end


mutable struct ActorMars{D,P} <: SingleAbstractActor{D,P}
    dd::IMAS.dd{D}
    par::OverrideParameters{P,FUSEparameters__ActorMars{P}}
    chease_inputs::Union{Nothing,CHEASEnamelist}
    mars_inputs::Union{Nothing,MARSnamelist}

    function ActorMars(
        dd::IMAS.dd{D}, 
        par::FUSEparameters__ActorMars{P}; 
        chease_overrides = NamedTuple(),
        mars_overrides = MarsOverrides(),
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
        for (k, v) in pairs(chease_overrides)
            if k ∉ fieldnames(CHEASEnamelist)
                error("Unknown namelist entry '$k'")
            end
            setfield!(nl, k, v)
        end

        # -------------------------
        # Create default MARS namelist and apply overrides
        # -------------------------
        mars_nl = mars_inputs === nothing ? MARSnamelist() : mars_inputs
        apply_overrides!(mars_nl, mars_overrides)

        # ------------------------- 
        # Final actor construction
        # -------------------------
        return new{D,P}(dd, par, nl, mars_nl)
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
    chease_namelist = actor.chease_inputs
    mars_namelist = actor.mars_inputs

    # plot equilibrium if requested
    if par.do_plot
        if !isempty(dd.equilibrium.time_slice)
            println("Plotting initial equilibrium...")
            plt = FUSE.Plots.plot(dd.equilibrium; label="before ActorEquilibrium")
        else
            plt = FUSE.Plots.plot()
        end
        display(plt)
    end

    #run equilibrium solver to generate initial conditions for MARS
    if par.run_equilibrium
        @info "Running CHEASE equilibrium solver with EQDSK=$(par.EQDSK)."
        run_CHEASE(dd, par, chease_namelist)
    end

    # run MARS
    if par.run_MHD
        @info "Running MARS actor with parameters: tracer_type=$(par.tracer_type)"
        run_MARS(dd, par, mars_namelist)
    end

    #run_PARTICLE_TRACING(dd, par)
    return actor
end


function run_CHEASE(dd::IMAS.dd, par, chease_namelist)
    
    chease_exec = par.chease_exec
    @assert chease_namelist !== nothing "CHEASE namelist not initialized"

    # Do the No-wall checks and get MARS wall file from the FUSE repository
    limiter = nothing
    if par.wall_type == :no_wall && par.number_surfaces > 1
        error("Invalid configuration: number_surfaces > 1 but wall_type is set to :
No_wall. Please specify a valid wall_type or set number_surfaces to 1.")
    end

    if par.restart_equilibrium
        if isfile("EXPEQ.OUT")
            @info "Restarting CHEASE from existing equilibrium files."
            # Implement logic to copy existing equilibrium files to current directory
            @info "Replacing EXPEQ with EXPEQ.OUT from previous run for the restart."
            cp("EXPEQ.OUT", "EXPEQ"; force=true)
        else
            error("Invalid CHEASE run configuration: there is NO EXPEQ.OUT to restart from. Do a CLEAN restart)")
        end
    else
        @info "Clean CHEASE run from dd."
        # extract B0 and R0 for CHEASE normalization and overwrite namelist entries
        B0, R0 = write_EXPEQ_file(dd, par)
        #CHEASE.write_EXPEQ_file(eq_chease)
        setfield!(chease_namelist, :B0EXP, B0)
        setfield!(chease_namelist, :R0EXP, R0)
    end
   
    if par.number_surfaces == 1 && chease_namelist.NVEXP == 8
        @info "Overriding NVEXP in CHEASE namelist"
        NVEXP = 1
        setfield!(chease_namelist, :NVEXP, NVEXP)
    end
    
    # Write CHEASE namelist file 
    write_CHEASEnamelist(chease_namelist, "datain")

    # Clean up any stale CHEASE output files
    for f in readdir(pwd())
        if startswith(f, "OUT") && occursin("MAR", f)
            @info "Checking deleting: $f" 
            rm(f; force=true)
        end
    end

    # Execute CHEASE
    @assert isfile("datain") "CHEASE input file datain not found"
    @info "Executing CHEASE from $chease_exec"
    isfile(chease_exec) || error("CHEASE executable not found: $chease_exec")

    cmd = pipeline(
        `$(chease_exec)`,
        stdin  = "datain",
        stdout = "log_chease",
        stderr = "log_chease"
    )

    ok = success(cmd)
    ok || error("CHEASE failed — see log_chease")

    # basic checks on CHEASE output
    keys = ["GEXP", "Q_ZERO", "Q_EDGE"]
    println(julia_grep(keys, "log_chease"))
    
    return nothing
end

function run_MARS(dd::IMAS.dd, par, mars_namelist)

    core_profiles = dd.core_profiles.profiles_1d[]
    
    # override IWALL in mars_namelist if number_surfaces > 1
    # NOTE - it's OK to overwrite IWALL here before NWALL is updated in
    # the next step because if NWALL = 0, IWALL does NOT matter
    if par.wall_type != :no_wall && par.number_surfaces > 1
        @info "Overriding IWALL in RUN.IN"
        # Implement the override logic here
        NW = julia_grep(["NW"], "log_chease"; extract_values=true)["NW"]
        setfield!(mars_namelist.BASIC, :IWALL, NW)
    end

    # Write final MARS namelist into run directory for MARS execution
    write_MARS_namelist(mars_namelist, "RUN.IN")

    # Determine which profiles to pull from dd, i.e. experiment
    write_exp_profiles(core_profiles, mars_namelist)

    # 5. Execute MARS
    mars_exec = par.mars_exec
    @assert isfile("RUN.IN") "MARS input file RUN.IN not found"
    isfile(mars_exec) || error("MARS executable not found: $mars_exec")
    
    if par.run_mode == :batch
        @info "Submitting MARS job to batch system with command: $(par.batch_submit_cmd) and script: mars_job.sh"
        run_batch_mars(mars_exec, par)
    elseif par.run_mode == :local
        @info "Running MARS interactively with command: $mars_exec"
        cmd = pipeline(
        `time $(mars_exec)`,
        stdout = "log_mars",
        stderr = "log_mars"
        )
        ok = success(cmd)
    else
        error("Unknown run mode: $(par.run_mode). Supported modes are :interactive and :batch.")
    end
    

    # Display growth rate and iteration count
    iter, growth, freq = parse_MARS_results("RESULT.OUT")
    @info "MARS iterations = $iter"
    @info "Growth rate = $growth"
    @info "Frequency = $freq"

    ok || error("MARS failed — see log_mars")

end



function run_batch_mars(mars_exec, par)

    script = """
    #!/bin/bash
    #SBATCH --job-name=mars_job
    #SBATCH --output=log_mars
    #SBATCH --error=log_mars

    cd $(pwd())

    time $mars_exec < RUN.IN
    """

    script_file = "mars_job.sh"
    open(script_file, "w") do io
        write(io, script)
    end

    run(`chmod +x $script_file`)

    submit_cmd = par.batch_submit_cmd === nothing ? "sbatch" : par.batch_submit_cmd

    run(`$submit_cmd $script_file`)
end

"""
    get_limiter_data(wall_type::Symbol) -> (r_wall::Vector{Float64}, z_wall::Vector{Float64})      
    Retrieve limiter boundary data for specified machine wall type. 
    This function reads pre-defined JSON files containing the R,Z coordinates of the limiter for 
        different machines (e.g., D3D, ITER, ASDEX, MAST, KSTAR) and returns the R and Z coordinates as vectors. 
        The wall_type argument specifies which machine's limiter data to retrieve.

    """

function get_limiter_data(machine::String)
    if machine != "D3D" && machine != "ITER"
        error("Unknown machine: $machine. Supported machines are 'D3D', 'ITER' for now.")
    end
    
    file_path = dirname(dirname(pathof(FUSE)))* "/sample/"
    file_path = joinpath(file_path, machine * "_mars_wall.json")
    dd_wall = IMAS.json2imas(file_path).wall.description_2d[].limiter.unit[1].outline
    
    return dd_wall
    
end

""" 
  apply_overrides!(MARSnamelist, MarsOverrides) -> MARSnamelist 
    Apply user-specified runtime overrides to a MARSnamelist instance. 
    The MarsOverrides struct contains Dicts for each MARS namelist block (e.g., BASIC, FEEDBACK, NUMERIC, OUTOPT, KINETIC), where the keys are the parameter names to override and the values are the new values to set. This function modifies the MARSnamelist in-place based on the provided overrides.
"""

function apply_overrides!(
    nl::MARSnamelist,
    mo::MarsOverrides
)

    for block in fieldnames(typeof(mo))

        overrides = getfield(mo, block)
        overrides === nothing && continue

        target_block = getfield(nl, block)

        for (k,v) in overrides

            if k ∈ fieldnames(typeof(target_block))

                setfield!(target_block, k, v)
            else

                target_block.extras[k] = v

            end

        end
    end
end


function write_MARS_namelist(
    nl::MARSnamelist,
    filename::AbstractString
)

    open(filename,"w") do io

        for block in fieldnames(typeof(nl))

            println(io,"&",block)

            b = getfield(nl,block)

            for k in fieldnames(typeof(b))

                k == :extras && continue

                v = getfield(b,k)

                if v isa Complex
                    sval = "($(real(v)), $(imag(v)))"

                elseif v isa AbstractVector{<:Complex}
                    sval = join(
                        ["($(real(z)), $(imag(z)))" for z in v],
                        ", "
                )

                elseif v isa AbstractVector
                    sval = "(" * join(v,", ") * ")"

                else
                    sval = v
                end

                println(io," $k = $sval,")
            end

            for (k,v) in b.extras

                if v isa Complex
                    sval = "($(real(v)), $(imag(v)))"

                elseif v isa AbstractVector{<:Complex}
                    sval = "(" * join(
                        ["($(real(z)), $(imag(z)))" for z in v],
                        ", "
                    ) * ")"

                elseif v isa AbstractVector
                    sval = "(" * join(v,", ") * ")"

                else
                    sval = string(v)
                end

                println(io," $k = $sval,")
            end

            println(io,"&END\n")
        end
    end
end



function write_exp_profiles(profiles, mars_namelist)
    @info "Generating additional PROF*.IN files from experiment or dd."
    
    s = sqrt.(profiles.grid.psi_norm)  # radial coordinate

    for (flag, outputs, msg) in (
        (
            :NPROFR,
            [("PROFROT.IN", p -> p.rotation_frequency_tor_sonic)],
            nothing
        ),

        (
            :NPROFN,
            [("PROFDEN.IN", p -> p.ion[1].density)],
            nothing
        ),

        (
            :NPROFWE,
            [("PROFNE.IN", p -> p.ion[1].rotation_frequency_tor)],
            "NO ExB rotation in dd, using bulk plasma rotation instead."
        ),

        (
            :NPROFTTCA,
            [("PROFTTCPARA.IN", p -> p.conductivity_parallel)],
            nothing
        ),

        (
            :NPROFTTCE,
            [],
            "NO χ_perpendicular in IMAS dd"
        ),

        (
            :NPROFIE,
            [
                ("PROFTI.IN", p -> p.ion[1].temperature),
                ("PROFTE.IN", p -> p.electron.temperature)
            ],
            nothing
        )
        
    )

        if getfield(mars_namelist.BASIC, flag) == 4

            msg !== nothing && @info msg

            for (file, getter) in outputs
                profile = getter(profiles)
                profile ./= maximum(profile)
                write_profile_IN(file, s, profile)
            end
        end
    end 
end

    




# write the *.IN files for each quantity
function write_profile_IN(filename::AbstractString, s::Vector{Float64}, profile::Vector{Float64})
    @assert length(s) == length(profile) "normal coordinate and profile arrays must have the same shape"
    
    KEY = "1"   # just an extra string the MARS profiles files need in the top line
    write_list = string(length(profile), "    ", KEY)
    for (r, z) in zip(s, profile)
        write_list = vcat(write_list, "$r    $z")
    end
    touch(filename)
    open(filename, "w") do file
        for line in write_list
            write(file, "$line \n")
        end
    end
end


function parse_MARS_results(filename::AbstractString)
    isfile(filename) || error("File not found: $filename")

    line = strip(readline(filename))
    vals = [strip(v, ',') for v in split(line)]

    length(vals) >= 6 || error("Unexpected format in $filename")

    iter   = round(Int, parse(Float64, vals[2]))
    growth = parse(Float64, vals[5])
    freq   = parse(Float64, vals[6])

    # Convergence warning
    if iter > 99
        @warn "MARS did not converge (iterations = $iter). A better initial guess for the eigenvalue may be required."
    end

    return iter, growth, freq
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
function write_EXPEQ_file(dd::IMAS.dd, par)

    offset = par.offset  # offset for first wall (RW) in meters
    n_points = par.n_points  # number of points for first wall (RW)
    NWBPS = par.number_surfaces
    
    # initialize eqt from pulse_schedule and core_profiles
    time_slice = dd.equilibrium.time_slice[]
    eqt1d = time_slice.profiles_1d
    
    # populate the input file lines
    minor_radius = time_slice.boundary.minor_radius
    z_axis = time_slice.global_quantities.magnetic_axis.z
    Bt_center = time_slice.global_quantities.vacuum_toroidal_field.b0
    Bt_axis = time_slice.global_quantities.magnetic_axis.b_field_tor
    r_center = time_slice.global_quantities.vacuum_toroidal_field.r0
    r0 = dd.equilibrium.vacuum_toroidal_field.r0
    Ip = time_slice.global_quantities.ip
    r_bound = time_slice.boundary.outline.r
    z_bound = time_slice.boundary.outline.z
    r_geo = time_slice.boundary.geometric_axis.r
    z_geo = time_slice.boundary.geometric_axis.z
    #Bt_geo = Bt_center * r_center / r_geo

    # choose B0 & R0 for CHEASE normalization
    B0 = abs(Bt_center)
    R0 = r0

    # inverse aspect ratio for CHEASE input
    ϵ = minor_radius / R0

    # get the normalized psi and convert d/dPsi to d/ds coordinate for CHEASE input
    psi_norm = eqt1d.psi_norm
    psi = eqt1d.psi
    s = sqrt.(psi_norm)
    #s_grid = range(0.0, 1.0; length=length(psi))
    pressure = eqt1d.pressure
    pressure_sep = pressure[end]
    pprime = eqt1d.dpressure_dpsi

    ### Currently NOT used, but may be useful later
    #wall_RZ = [dd.wall.description_2d[].limiter.unit[1].outline.r, dd.wall.description_2d[].limiter.unit[1].outline.z]

    if minimum(r_bound) - offset < 0
        error("Offset too large: boundary crosses R < 0 (min R = $(minimum(r_bound)))")
    end

    ## GS current density specification and de-dimensionalization for CHEASE input
    if par.GS_rhs == :FFpr
        NSTTP = 1
        FFpr = 2 * pi * eqt1d.f_df_dpsi
        GS_RHS_norm = FFpr / B0
    elseif par.GS_rhs == :Jtor
        NSTTP = 2
        Jtor = abs.(eqt1d.j_tor)
        GS_RHS_norm = Jtor / (B0 / R0 * μ_0)
    elseif par.GS_rhs == :Jpar
        NSTTP = 3
        Jpar = abs.(eqt1d.j_parallel) # NOT right!
        GS_RHS_norm = Jpar / (B0 / R0 * μ_0)
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
    
    # interpolate to uniform s grid for CHEASE input
    GS_RHS_final = GS_RHS_norm#IMAS.interp1d(s, abs.(j_tor_norm)).(s_grid)
    pprime_at_s = pprime #IMAS.interp1d(s, pprime).(s_grid)
    
    # Make pressure terms dimensionless for CHEASE input
    # throw in a 2pi to scale P' correctly
    pressure_sep_norm = pressure_sep / (B0^2 / μ_0)
    pprime_final = 2 * pi * pprime_at_s * R0^2 * μ_0 / B0
    
    # Remove/smooth the X-point along the boundary
    ab = sqrt((maximum(r_bound) - minimum(r_bound))^2 + (maximum(z_bound) - minimum(z_bound))^2) / 2.0
    pr, pz = limit_curvature(r_bound, z_bound, ab / 20.0)
    rb_new, zb_new = IMAS.resample_2d_path(pr, pz; n_points=n_points, method=:linear)
    r_bound_norm = rb_new / R0
    z_bound_norm = zb_new / R0

    #plt = plot()
    plt = plot!(rb_new, zb_new; linewidth=3., aspect_ratio=:equal, label="Plasma Boundary")
    display(plt)

    
    ##----------------- Write the file -----------------##
    write_list = [string(ϵ), string(z_geo/r_geo), string(pressure_sep_norm)]
    @assert length(rb_new) == length(zb_new) "R,Z boundary arrays must have the same shape"
    write_list = vcat(write_list, string(length(rb_new), " ", NWBPS, " ", NDATA))
    for (r, z) in zip(r_bound_norm, z_bound_norm)
        write_list = vcat(write_list, "$r    $z")
    end

    # if there is another surface, calclate its cooridanes given an offset and save to file
    if NWBPS > 1 ## WHAT TO DO if > 2
        if par.wall_type == :conformal 
            @info "Creating a conformal limiter offset from plasma boundary by $(par.offset)."
            r_lim, z_lim = offset_boundary(rb_new, zb_new, offset)
        elseif par.wall_type == :limiter
            @info "Using wall data .Json for CHEASE equilibrium generation."
            machine = dd.dataset_description.data_entry.machine
            limiter = get_limiter_data(machine)
            r_lim, z_lim = limiter.r, limiter.z  # put the length scale back in
        end
        # add a smooth first wall (RW)
        
        r_lim, z_lim = IMAS.resample_2d_path(r_lim, z_lim; n_points=n_points, method=:linear)
         
        r_lim_norm = r_lim / R0
        z_lim_norm = z_lim / R0
        
        plt = plot!(r_lim, z_lim; linewidth=1.5, aspect_ratio=:equal,label="MARS resistive wall", legend=:outertop)
        display(plt)
        for (r, z) in zip(r_lim_norm, z_lim_norm)
            write_list = vcat(write_list, "$r    $z")
        end
    end

    @assert length(s) == length(pprime_final) == length(GS_RHS_final) "s, presssure and GS_RHS arrays must have the same shape"
    write_list = vcat(write_list, "$(length(s))")
    write_list = vcat(write_list, "$(string(NSTTP))")
    write_list = vcat(write_list, map(string, s))
    write_list = vcat(write_list, map(string, pprime_final))
    write_list = vcat(write_list, map(string, GS_RHS_final))

    # write to EXPEQ file   
    touch("EXPEQ")
    open("EXPEQ", "w") do file
        for line in write_list
            write(file, "$line \n")
        end
    end
    return B0, R0
end

function write_EXPEQ_file2(dd::IMAS.dd, par)
    # Placeholder function to write EXPEQ file for CHEASE
    @info "Writing EXPEQ file for CHEASE equilibrium solver."

    
    offset = par.offset  # offset for first wall (RW) in meters
    n_points = par.n_points  # number of points for first wall (RW)
    NWBPS = par.number_surfaces
    
    # initialize eqt from pulse_schedule and core_profiles
    time_slice = dd.equilibrium.time_slice[]
    eqt1d = time_slice.profiles_1d
    
    # populate the input file lines
    minor_radius = time_slice.boundary.minor_radius
    z_axis = time_slice.global_quantities.magnetic_axis.z
    Bt_center = time_slice.global_quantities.vacuum_toroidal_field.b0
    Bt_axis = time_slice.global_quantities.magnetic_axis.b_field_tor
    r_center = time_slice.global_quantities.vacuum_toroidal_field.r0
    r0 = dd.equilibrium.vacuum_toroidal_field.r0
    Ip = time_slice.global_quantities.ip
    r_bound = time_slice.boundary.outline.r
    z_bound = time_slice.boundary.outline.z
    r_geo = time_slice.boundary.geometric_axis.r
    z_geo = time_slice.boundary.geometric_axis.z
    #Bt_geo = Bt_center * r_center / r_geo

    # choose B0 & R0 for CHEASE normalization
    B0 = abs(Bt_center)
    R0 = r0

    # inverse aspect ratio for CHEASE input
    ϵ = minor_radius / R0

    # get the normalized psi and convert d/dPsi to d/ds coordinate for CHEASE input
    psi_norm = eqt1d.psi_norm
    psi = eqt1d.psi
    s = sqrt.(psi_norm)
    #s_grid = range(0.0, 1.0; length=length(psi))
    pressure = eqt1d.pressure
    pressure_sep = pressure[end]
    pprime = eqt1d.dpressure_dpsi

    ### Currently NOT used, but may be useful later
    #wall_RZ = [dd.wall.description_2d[].limiter.unit[1].outline.r, dd.wall.description_2d[].limiter.unit[1].outline.z]

    if minimum(r_bound) - offset < 0
        error("Offset too large: boundary crosses R < 0 (min R = $(minimum(r_bound)))")
    end

    ## GS current density specification and de-dimensionalization for CHEASE input
    if par.GS_rhs == :FFpr
        NSTTP = 1
        FFpr = 2 * pi * eqt1d.f_df_dpsi
        GS_RHS_norm = FFpr / B0
    elseif par.GS_rhs == :Jtor
        NSTTP = 2
        Jtor = abs.(eqt1d.j_tor)
        GS_RHS_norm = Jtor / (B0 / R0 * μ_0)
    elseif par.GS_rhs == :Jpar
        NSTTP = 3
        Jpar = abs.(eqt1d.j_parallel) # NOT right!
        GS_RHS_norm = Jpar / (B0 / R0 * μ_0)
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
    
    # Make pressure terms dimensionless for CHEASE input
    # throw in a 2pi to scale P' correctly
    pressure_sep_norm = pressure_sep / (B0^2 / μ_0)
    pprime_final = 2 * pi * pprime * R0^2 * μ_0 / B0
    
    # Remove/smooth the X-point along the boundary
    ab = sqrt((maximum(r_bound) - minimum(r_bound))^2 + (maximum(z_bound) - minimum(z_bound))^2) / 2.0
    pr, pz = limit_curvature(r_bound, z_bound, ab / 20.0)
    rb_new, zb_new = IMAS.resample_2d_path(pr, pz; n_points=n_points, method=:linear)
    r_bound_norm = rb_new / R0
    z_bound_norm = zb_new / R0

    eq = CHEASE.MartianChease(
        ϵ=ϵ,
        z_axis=z_geo/r_geo,
        pressure_sep=pressure_sep,
        r_center=R0,
        Bt_center=B0,
        Ip=Ip,
        r_bound=rb_new,
        z_bound=zb_new,
        mode=NSTTP,
        rho_pol=s,
        pressure=pprime_final,
        j_tor=GS_RHS_norm,
        wall_surfaces=walls
    )

    println("EXPEQ file generation completed.")
end

function write_CHEASEnamelist(
    nl::CHEASEnamelist, 
    filename::AbstractString="datain", 
    #overrides::NamedTuple=NamedTuple()
)

    # for (k, v) in pairs(overrides)
    #     if hasproperty(nl, k)
    #         setproperty!(nl, k, v)
    #     else
    #         @warn "Unknown CHEASE namelist key: $k"
    #     end
    # end
    
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


function julia_grep(
    patterns::AbstractVector{<:AbstractString},
    filename::AbstractString;
    extract_values::Bool=false
)
    isfile(filename) || error("File not found: $filename")

    # Original behavior
    if !extract_values
        matches = String[]
        open(filename, "r") do file
            for line in eachline(file)
                any(p -> occursin(p, line), patterns) && push!(matches, line)
            end
        end
        return matches
    end

    # Extraction mode
    results = Dict{String,Any}()

    open(filename, "r") do file
        for line in eachline(file)
            for key in patterns
                # standalone KEY = value
                pattern = Regex("\\b$(key)\\b\\s*=\\s*([^,\\s]+)")
                if (m = match(pattern, line)) !== nothing
                    raw = strip(m.captures[1])

                    if (v = tryparse(Int, raw)) !== nothing
                        results[key] = v
                    elseif (v = tryparse(Float64, raw)) !== nothing
                        results[key] = v
                    else
                        results[key] = raw
                    end
                end
            end
        end
    end

    return results
end


"""
    extract_lines_for_keys(logfile, keys, DataType) -> Dict(key, value::DataType)

Scan `logfile` once and return a dictionary mapping each key
to the first line containing it.

Errors if any key is not found.
"""
function extract_lines_for_keys(
    logfile::AbstractString,
    keys::AbstractVector{<:AbstractString},
    T::Type{<:Any}
) ::Dict{String, T}
    isfile(logfile) || error("Log file not found: $logfile")

    remaining = Set(keys)
    results   = Dict{String, T}()

    for line in eachline(logfile)
        for key in remaining
            if occursin(key, line)
                # extract text after key
                val = strip(replace(line, key => ""))
                val = strip(lstrip(val, ['=', ':']))  # handle ": or =" separators
                val = strip(rstrip(val, ','))   # handle "= or :" separators
                results[key] = parse(T, val)
                delete!(remaining, key)
            end
        end
        isempty(remaining) && break
    end

    #isempty(remaining) || error("Keys not found in $logfile: $(collect(remaining))")

    return results
end


