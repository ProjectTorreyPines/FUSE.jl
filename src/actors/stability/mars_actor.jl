using Interpolations
import CHEASE

const μ_0 = 4pi * 1E-7


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
    # Resistive walls. NWALL is the number of resistive shells; IWALL/TAUW are per-shell
    # (Fortran declares them as IWALL(NWALL0), TAUW(NWALL0) with NWALL0=5, see globalm.f),
    # so they are vectors here and must both have length NWALL — enforced by
    # configure_wall!(), which also sets NWALL from par.wall_type.
    #
    # IWALL entries are indices into the VACUUM RADIAL MESH (MARS uses IV = NR + IWALL),
    # NOT pointers to a wall contour supplied via EXPEQ. That mesh spans the plasma boundary
    # (s=1) out to CHEASE's far-field REXT (s=6 in the shipped examples). MARS's only validity
    # check is a bounds test, `IWALL(J) > 1 && IWALL(J) < NV` (marsq.f:2144), which cannot tell
    # a CHEASE-computed wall position from an arbitrary vacuum grid point. A stale IWALL
    # therefore silently places a resistive shell at an unintended radius rather than erroring,
    # which is why configure_wall!() always re-derives it from log_chease.
    NWALL::Int = 1
    TAUW::Vector{Float64} = [1.0e4]
    IWALL::Vector{Int} = [22]
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
    FCCHI::Vector{Float64} = [0.34219, -0.405470]
    FWCHI::Vector{Float64} = [0.12187,  0.089063]
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
    ALTAU::ComplexF64 = 0.5 + 0.0im
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

# REORBIT particle-tracing settings. These are Fortran variables declared in the
# `OUTOPT` namelist (see MarsQ_package/MarsQ/newrun.inc), NOT a namelist of their own —
# MARS has no `&REORBIT` block. Kept as a separate Julia struct for readability; see the
# NAMELIST_GROUP mapping used by write_MARS_namelist(), which folds this struct's fields
# into the `&OUTOPT ... &END` block alongside MARS_OUTOPT's.
# Defaults below are taken from MarsQ_package/EXAMPLE/RwmFluid/Rmars_EXAMPLE_RwmFluid_EPs,
# except NREORBIT, which is driven by `par.num_orbits` (see run_MARS).
Base.@kwdef mutable struct MARS_REORBIT
    NREORBIT::Int = 0
    KRE_VAC::Int = 0
    KRE_INIT::Int = 1
    KRE1::Int = 1
    KRE2::Int = 5
    NRE1::Int = 30
    NRE2::Int = 30
    KRE_STAR::Int = 5
    KRE_STEP::Int = 10
    KRE_STEP_MAX::Int = 300000
    KRE_SIGMA0::Int = 1
    KRE_TRACE::Int = 5
    RE_TSTEP::Float64 = 1.0e-3
    RE_TIME::Float64 = 30000.0
    RE_EFIELD::Float64 = 0.0
    RE_SAC::Float64 = 0.0
    RE_SYNCH::Float64 = 0.0
    RE_BREMS::Float64 = 0.0
    RE_S0::Float64 = 0.0
    RE_CHI0::Float64 = 0.0
    RE_P0::Float64 = 2.308627075871614e-02
    RE_AN::Float64 = 2.0
    RE_CN::Float64 = 1.0
    RE_LAMBDA0::Float64 = 1.0
    RE_ATOL::Float64 = 1.0e-7
    RE_PERTURB::Float64 = 1.1995e-04
    RE_CONST::Vector{Float64} = [1.18e+20, 15.0, 10.0, 1.0]
    extras::Dict{Symbol,Any} = Dict{Symbol,Any}()
end

Base.@kwdef mutable struct MARSnamelist
    BASIC::MARS_BASIC = MARS_BASIC()
    FEEDBACK::MARS_FEEDBACK = MARS_FEEDBACK()
    KINETIC::MARS_KINETIC = MARS_KINETIC()
    QLIN::MARS_QLIN = MARS_QLIN()
    NUMERIC::MARS_NUMERIC = MARS_NUMERIC()
    OUTOPT::MARS_OUTOPT = MARS_OUTOPT()
    REORBIT::MARS_REORBIT = MARS_REORBIT()
end

# Maps each MARSnamelist field to the actual Fortran namelist group it must be written
# under in RUN.IN. Only REORBIT differs from its own field name: MARS has no `&REORBIT`
# namelist, so its fields are folded into `&OUTOPT` (see write_MARS_namelist()).
const NAMELIST_GROUP = Dict(
    :BASIC    => :BASIC,
    :FEEDBACK => :FEEDBACK,
    :KINETIC  => :KINETIC,
    :QLIN     => :QLIN,
    :NUMERIC  => :NUMERIC,
    :OUTOPT   => :OUTOPT,
    :REORBIT  => :OUTOPT,
)

Base.@kwdef mutable struct MarsOverrides
    BASIC::Dict{Symbol,Any} = Dict{Symbol,Any}()
    FEEDBACK::Dict{Symbol,Any} = Dict{Symbol,Any}()
    KINETIC::Dict{Symbol,Any} = Dict{Symbol,Any}()
    QLIN::Dict{Symbol,Any} = Dict{Symbol,Any}()
    NUMERIC::Dict{Symbol,Any} = Dict{Symbol,Any}()
    OUTOPT::Dict{Symbol,Any} = Dict{Symbol,Any}()
    REORBIT::Dict{Symbol,Any} = Dict{Symbol,Any}()
end


"""
    MarsModeStructure

Eigenfunction (mode structure) parsed from MARS `XPLASMA.OUT` (displacement ξ) and,
when available, `VPLASMA.OUT` (perturbed velocity).

- `s`         : radial coordinate `s = sqrt(ψ_norm)`, length `NRP1`
- `m_pol`     : poloidal Fourier mode numbers, length `MSMAX`
- `xi1/2/3`   : complex contravariant components of the plasma displacement ξ, size `(NRP1, MSMAX)`
                (`xi1` ≡ normal/radial component, `xi3` ≡ parallel-like component)
- `v1/2/3`    : complex contravariant components of the perturbed velocity (`nothing` if not written)
- `chi`       : geometric poloidal angle χ ∈ [0,2π], length `Nχ`
- `R`, `Z`    : real-space flux-surface geometry `R(s,χ)`, `Z(s,χ)` [m], size `(NRP1, Nχ)`,
                reconstructed from the R,Z Fourier harmonics in `RMZM_F.OUT` (for R,Z-space plotting)
"""
mutable struct MarsModeStructure
    s::Vector{Float64}
    m_pol::Vector{Float64}
    xi1::Matrix{ComplexF64}
    xi2::Matrix{ComplexF64}
    xi3::Matrix{ComplexF64}
    v1::Union{Nothing,Matrix{ComplexF64}}
    v2::Union{Nothing,Matrix{ComplexF64}}
    v3::Union{Nothing,Matrix{ComplexF64}}
    chi::Vector{Float64}
    R::Matrix{Float64}
    Z::Matrix{Float64}
end

"""
    MarsOutputs

Container for the MARS results that are stored into `dd.mhd_linear`.

- `n_tor`       : toroidal mode number (from `RNTOR`)
- `iterations`  : number of MARS eigenvalue iterations
- `growthrate`  : `Re(γ)·τ_A`, growth rate normalized to the Alfvén time
- `frequency`   : `Im(γ)·τ_A`, frequency normalized to the Alfvén time
- `ideal`       : `true` for ideal MHD (`ETA == 0`), `false` for resistive
- `mode`        : eigenfunction (`MarsModeStructure`), or `nothing` if output files are absent
"""
mutable struct MarsOutputs
    n_tor::Int
    iterations::Int
    growthrate::Float64
    frequency::Float64
    ideal::Bool
    mode::Union{Nothing,MarsModeStructure}
end

Base.@kwdef mutable struct FUSEparameters__ActorMars{T<:Real} <: ParametersActor{T}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    _time::Float64 = NaN
    do_plot::Entry{Bool} = act_common_parameters(; do_plot=false)
    EQDSK::Entry{Bool} = Entry{Bool}("-",
        "Whether CHEASE reads its equilibrium input from a g-file instead of the EXPEQ built from dd " *
        "(the current default chain). NOT YET IMPLEMENTED — reserved for future g-file input support."; default=false)
    chease_exec::Entry{String} =
        Entry{String}("-", "Path to CHEASE executable"; default="/fusion/projects/codes/mars/CHEASE/chease.x")
    mars_exec::Entry{String} =
        Entry{String}("-", "Path to MARS executable"; default="/fusion/projects/codes/mars/MARSQ/marsq.x")
    mpirun_exec::Entry{String} =
        Entry{String}("-", "mpirun/mpiexec executable used when num_orbits>0 (REORBIT tracing requires MPI even though MARS itself does not). Set to a full path if `mpirun` on PATH does not match the MPI marsq.x was built against (e.g. a conda-installed mpirun shadowing the system one)."; default="mpirun")
    offset::Entry{Float64} = Entry{Float64}("-", "Offset for conforming first wall (RW) in units of length"; default=0.2)
    n_points::Entry{Int} = Entry{Int}("-", "Number of points for discretizing plasma boundary and surrounding walls "; default=201)
    GS_rhs::Switch{Symbol} = Switch{Symbol}([:FFpr, :Jtor, :Jpar], "-", "Specification of Grad-Shaf RHS current"; default=:FFpr)
    wall_resistivity_type::Switch{Symbol} = Switch{Symbol}([:Constant, :Variable], "-", "Wall Resistivity Model"; default=:Constant)
    wall_type::Switch{Symbol} = Switch{Symbol}([:no_wall, :conformal, :limiter], "-", "Machine wall shape to use for MARS"; default=:no_wall)
    num_surfaces::Entry{Int} = Entry{Int}("-", "Number of surfaces to specify"; default=1)
    wall_time_constants::Entry{Vector{Float64}} =
        Entry{Vector{Float64}}("-", "Resistive wall time constant τ_w per wall surface, in Alfvén times (MARS TAUW). Length must equal the number of resistive walls; a single-element vector is applied to every wall. Ignored when wall_type=:no_wall"; default=[1.0e4])
    run_equilibrium::Entry{Bool} = Entry{Bool}("-", "Whether to run equilibrium solver"; default=true)
    restart_equilibrium::Entry{Bool} = Entry{Bool}("-", "Whether to restart from existing equilibrium"; default=false)
    beta_fac::Entry{Float64} = Entry{Float64}("-",
        "Sets CHEASE's CFBAL namelist parameter directly (rescales pressure gradient/edge pressure in " *
        "the EXPEQ.OUT this run writes out). NOTE: CFBAL only affects EXPEQ.OUT, not the equilibrium this " *
        "run itself solves — a beta scan needs restart_equilibrium=true on the next run to pick it up. " *
        "Cannot be combined with chease_overrides.CFBAL (use this parameter instead)."; default=1.0)
    run_MHD::Entry{Bool} = Entry{Bool}("-", "Whether to run MHD stability code"; default=true)  
    run_mode::Switch{Symbol} = Switch{Symbol}([:local, :batch], "-", "Whether to run MARS locally or submit to batch system"; default=:local)
    batch_submit_cmd::Entry{String} = Entry{String}("-", "Batch submission command used when run_mode=:batch"; default="sbatch")
    num_orbits::Entry{Int} = Entry{Int}("-", "Number of orbits to simulate in particle tracing"; default=0)
    save_dir::Entry{String} =
        Entry{String}("-", "Directory in which to run CHEASE/MARS. If empty, a temporary directory is created and recorded here so it can be reused (e.g. for restart_equilibrium or run_MHD-only follow-up runs via `save_dir=actor.par.save_dir`)"; default="")
    clear_workdir::Entry{Bool} =
        Entry{Bool}("-", "Remove the auto-created temporary run directory after the run (ignored when save_dir is provided)"; default=false)
end


mutable struct ActorMars{D,P} <: SingleAbstractActor{D,P}
    dd::IMAS.DD{D}
    par::OverrideParameters{P,FUSEparameters__ActorMars{P}}
    chease_inputs::Union{Nothing,CHEASE.CHEASEnamelist}
    mars_inputs::Union{Nothing,MARSnamelist}
    mars_outputs::Union{Nothing,MarsOutputs}

    function ActorMars(
        dd::IMAS.DD{D}, 
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
        nl = chease_inputs === nothing ? CHEASE.CHEASEnamelist() : chease_inputs

        # -------------------------
        # Apply user overrides
        # -------------------------
        if :CFBAL ∈ keys(chease_overrides)
            error("CFBAL cannot be set via chease_overrides — use the beta_fac actor parameter instead " *
                  "(par.beta_fac is transmitted to CHEASE's CFBAL; having both would silently pick one).")
        end
        for (k, v) in pairs(chease_overrides)
            if k ∉ fieldnames(CHEASE.CHEASEnamelist)
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
        return new{D,P}(dd, par, nl, mars_nl, nothing)
    end
end


"""
    ActorMars(dd::IMAS.DD, act::ParametersAllActors; kw...) 

"""
function ActorMars(dd::IMAS.DD, act::ParametersAllActors; kw...)
    actor = ActorMars(dd, act.ActorMars; kw...)
    step(actor)
    finalize(actor)
    return actor
end


"""
    _label_recipe_series!(plt, label::AbstractString)

IMAS's equilibrium recipe hard-codes `label := ""` on its own series (a forced
default, so a `label` kwarg passed to `plot(dd.equilibrium; label=...)` is silently
discarded — it never reaches the legend). Post-process the resulting plot object
instead: fill in `label` on any still-blank *primary* series (non-primary series,
e.g. the recipe's own dashed q=1 reference line, stay excluded from the legend
regardless of label, so they're left alone).
"""
function _label_recipe_series!(plt, label::AbstractString)
    for sp in plt.subplots, s in sp.series_list
        if s[:primary] && isempty(s[:label])
            s[:label] = label
        end
    end
    return plt
end

function _step(actor::ActorMars)
    dd = actor.dd
    par = actor.par
    chease_namelist = actor.chease_inputs
    mars_namelist = actor.mars_inputs

    # Fail fast on an inconsistent wall setup, before any CHEASE/MARS work. Checked here
    # rather than in run_CHEASE so it applies even when run_equilibrium=false.
    validate_wall_configuration(par)

    # plot equilibrium if requested. run_CHEASE (below) overlays the plasma boundary and
    # limiter/resistive-wall onto this same figure's flux-contour subplot, and _finalize
    # overlays the post-CHEASE profiles — all via Plots.jl's implicit current-plot state,
    # since nothing in between calls a fresh plot().
    if par.do_plot
        if !isempty(dd.equilibrium.time_slice)
            println("Plotting initial equilibrium...")
            plt = FUSE.Plots.plot(dd.equilibrium; label="before CHEASE", linewidth=3)
            _label_recipe_series!(plt, "before CHEASE")
        else
            plt = FUSE.Plots.plot()
        end
        display(plt)
    end

    # Determine the working directory for CHEASE/MARS execution.
    # If save_dir is empty we create a temporary directory and record its path in
    # par.save_dir, so the caller can reuse it (restart_equilibrium / run_MHD-only)
    # by passing `save_dir=actor.par.save_dir` to the next ActorMars call.
    created_workdir = isempty(par.save_dir)
    run_dir = created_workdir ? mktempdir() : abspath(par.save_dir)
    mkpath(run_dir)
    par.save_dir = run_dir
    @info "ActorMars running in $run_dir"

    old_dir = pwd()
    try
        cd(run_dir)

        #run equilibrium solver to generate initial conditions for MARS
        if par.run_equilibrium
            @info "Running CHEASE equilibrium solver with EQDSK=$(par.EQDSK)."
            run_CHEASE(dd, par, chease_namelist)
        end

        # run MARS
        if par.run_MHD
            actor.mars_outputs = run_MARS(dd, par, mars_namelist)
        end

    finally
        cd(old_dir)
    end

    # Only auto-clean a directory we created ourselves, and only if requested.
    # When chaining runs that depend on each other's files, set clear_workdir=false.
    if created_workdir && par.clear_workdir
        rm(run_dir; force=true, recursive=true)
    end

    return actor
end


"""
    parse_NGA(filename::AbstractString) -> (csm, p, dpdpsi, q)

Parse the `NGA` file CHEASE writes on every run (its own ASCII export of the CSM=s
mesh, pressure, dP/dψ, and q profiles — the block chease.f labels "SECTION ADDED TO
MATCH WITH GA'S EQUILIBRIUM CODE"). `p`/`dpdpsi` are in CHEASE's internal normalized
units; `q` is dimensionless already. `csm` is CHEASE's own equidistant s-mesh.
"""
function parse_NGA(filename::AbstractString)
    lines = readlines(filename)
    npsi1 = parse(Int, strip(lines[1]))

    function read_block(label::AbstractString)
        i = findfirst(l -> occursin(label, l), lines)
        i === nothing && error("NGA: block \"$label\" not found in $filename")
        vals = Float64[]
        j = i + 1
        while length(vals) < npsi1
            append!(vals, parse.(Float64, split(lines[j])))
            j += 1
        end
        return vals
    end

    csm    = read_block("CSM - MESH")
    p      = read_block("P (CSM)")
    dpdpsi = read_block("DP/DPSI(CSM)")
    q      = read_block("Q (CSM)")

    return csm, p, dpdpsi, q
end

"""
    _finalize(actor::ActorMars)

Store the MARS results into `dd.mhd_linear`:

- growth rate and frequency (`Re(γ)·τ_A`, `Im(γ)·τ_A`, normalized to the Alfvén time)
  in `time_slice[].toroidal_mode[].growthrate` / `.frequency`
- the toroidal mode number `n_tor` and the dominant poloidal mode number
- the mode structure (displacement eigenfunction, and perturbed velocity if available)
  in `time_slice[].toroidal_mode[].plasma`

Does nothing if MARS was not run (e.g. an equilibrium-only run).
"""
function _finalize(actor::ActorMars)
    dd = actor.dd
    par = actor.par

    # Overwrite dd.equilibrium with CHEASE's own recomputed pressure/q, read from NGA.
    # NGA is written unconditionally by CHEASE on every run (no namelist flag gates it),
    # so this always applies whenever CHEASE actually ran.
    if par.run_equilibrium
        nga_file = joinpath(par.save_dir, "NGA")
        if isfile(nga_file)
            csm, p_raw, _, q_new = parse_NGA(nga_file)
            psi_norm_chease = csm .^ 2   # CSM assumed = s = sqrt(psi_norm)

            B0EXP = actor.chease_inputs.B0EXP
            p_new = p_raw .* B0EXP^2 ./ μ_0

            eqt1d = dd.equilibrium.time_slice[].profiles_1d
            psi_norm_dd = eqt1d.psi_norm

            eqt1d.pressure = IMAS.interp1d(psi_norm_chease, p_new).(psi_norm_dd)
            eqt1d.q = IMAS.interp1d(psi_norm_chease, q_new).(psi_norm_dd)

            # overlay onto the "before CHEASE" figure from _step (which run_CHEASE's
            # write_EXPEQ_file already added the plasma boundary/limiter to), via
            # Plots.jl's current-plot state
            if par.do_plot
                plt = nothing
                try
                    plt = FUSE.Plots.plot!(dd.equilibrium; label="after CHEASE", linewidth=2.5)
                catch e
                    if isa(e, BoundsError)
                        plt = FUSE.Plots.plot(dd.equilibrium; label="after CHEASE", linewidth=2.5)
                    else
                        rethrow(e)
                    end
                end
                _label_recipe_series!(plt, "after CHEASE")
                display(plt)
            end
        else
            @warn "$nga_file not found (clear_workdir may have removed it before _finalize ran) — dd.equilibrium not updated from CHEASE."
        end
    end

    out = actor.mars_outputs
    out === nothing && return actor

    ml = dd.mhd_linear
    ml.code.name = "MARSQ"
    ml.model_type.index = 1            # global calculation
    ml.model_type.name = "global"
    ml.equations.index = 2             # full MHD
    ml.equations.name = "full"
    ml.ideal_flag = Int(out.ideal)
    ml.ids_properties.comment = "MARS-F/Q linear MHD stability. growthrate and frequency converted " *
                                "to SI from the MARS Alfvén-normalized eigenvalue using the on-axis Alfvén time."

    # On-axis Alfvén time τ_A(0) = R0·sqrt(μ0·ρ0)/B0, used by MARS to normalize the eigenvalue.
    # MARS reports Re(γ)·τ_A (growthrate) and Im(γ)·τ_A (angular frequency); convert back to SI.
    eqt = dd.equilibrium.time_slice[]
    B0 = abs(eqt.global_quantities.vacuum_toroidal_field.b0)
    R0 = dd.equilibrium.vacuum_toroidal_field.r0
    cp1d = dd.core_profiles.profiles_1d[]
    ρ_mass = IMAS.total_mass_density(cp1d)             # [kg/m^3] on the core_profiles grid
    s_cp = sqrt.(cp1d.grid.psi_norm)                   # MARS radial coordinate s = sqrt(ψ_norm)
    τA0 = R0 * sqrt(μ_0 * ρ_mass[1]) / B0              # on-axis Alfvén time [s]

    mhd = resize!(ml.time_slice; wipe=false)
    mode = resize!(mhd.toroidal_mode, "n_tor" => out.n_tor)
    mode.n_tor = out.n_tor
    mode.growthrate = out.growthrate / τA0             # [1/s]
    mode.frequency = out.frequency / τA0 / (2π)        # [Hz] (MARS gives angular frequency ω·τ_A)
    mode.perturbation_type.name = "MHD"
    mode.perturbation_type.description = "MARS-F/Q linear MHD eigenmode"

    # mode structure (eigenfunction)
    if out.mode !== nothing
        ms = out.mode
        pl = mode.plasma

        # radial label s = sqrt(ψ_norm) (dim1), poloidal Fourier modes m (dim2)
        pl.grid_type.index = 24
        pl.grid_type.name = "inverse_rhopolnorm_straight_field_line_fourier"
        pl.grid.dim1 = ms.s
        pl.grid.dim2 = ms.m_pol

        # Alfvén time profile τ_A(s) on the MARS radial grid
        ρ_s = IMAS.interp1d(s_cp, ρ_mass).(ms.s)
        pl.tau_alfven = R0 .* sqrt.(μ_0 .* ρ_s) ./ B0

        # Real-space flux-surface geometry R(s,χ), Z(s,χ) for plotting the mode in R,Z space.
        # The displacement harmonics ξₘ(s) live on (s, m); reconstruct ξ(s,χ)=Σₘ ξₘ(s)·e^{imχ}
        # on this same χ grid and evaluate at (R, Z) to plot.
        cs = pl.coordinate_system
        cs.grid_type.index = 2   # inverse: radial label (dim1) and poloidal angle (dim2)
        cs.grid_type.name = "inverse"
        cs.grid.dim1 = ms.s
        cs.grid.dim2 = ms.chi
        cs.r = ms.R
        cs.z = ms.Z

        # dominant poloidal harmonic = the one with the largest normal displacement
        idom = argmax(vec(maximum(abs.(ms.xi1); dims=1)))
        mode.m_pol_dominant = ms.m_pol[idom]

        # displacement ξ: normal contravariant component -> perpendicular, parallel-like -> parallel
        pl.displacement_perpendicular.real = real.(ms.xi1)
        pl.displacement_perpendicular.imaginary = imag.(ms.xi1)
        pl.displacement_parallel.real = real.(ms.xi3)
        pl.displacement_parallel.imaginary = imag.(ms.xi3)

        # perturbed velocity (3 contravariant components), if MARS wrote VPLASMA.OUT
        if ms.v1 !== nothing
            pl.velocity_perturbed.coordinate1.real = real.(ms.v1)
            pl.velocity_perturbed.coordinate1.imaginary = imag.(ms.v1)
            pl.velocity_perturbed.coordinate2.real = real.(ms.v2)
            pl.velocity_perturbed.coordinate2.imaginary = imag.(ms.v2)
            pl.velocity_perturbed.coordinate3.real = real.(ms.v3)
            pl.velocity_perturbed.coordinate3.imaginary = imag.(ms.v3)
        end
    end

    return actor
end


"""
    chease_normalization(dd::IMAS.DD) -> (B0, R0)

CHEASE/MARS normalization extracted from the equilibrium: `B0 = |B₀|` (vacuum toroidal
field at `R0`) and `R0` (vacuum toroidal field reference major radius). Used for both clean
and restart CHEASE runs so the equilibrium normalization is preserved across restarts.
"""
function chease_normalization(dd::IMAS.DD)
    time_slice = dd.equilibrium.time_slice[]
    B0 = abs(time_slice.global_quantities.vacuum_toroidal_field.b0)
    R0 = dd.equilibrium.vacuum_toroidal_field.r0
    return B0, R0
end


function run_CHEASE(dd::IMAS.DD, par, chease_namelist)

    chease_exec = par.chease_exec
    @assert chease_namelist !== nothing "CHEASE namelist not initialized"

    # wall_type/num_surfaces consistency is checked up front in validate_wall_configuration()
    limiter = nothing

    # CFBAL rescales pressure gradient/edge pressure in the EXPEQ.OUT this run writes out
    # (chease.f RPPF/PREDGE *= CFBAL at output time) — applies regardless of clean/restart.
    setfield!(chease_namelist, :CFBAL, par.beta_fac)

    if par.restart_equilibrium
        if isfile("EXPEQ.OUT")
            @info "Restarting CHEASE from existing equilibrium files."
            # Implement logic to copy existing equilibrium files to current directory
            @info "Replacing EXPEQ with EXPEQ.OUT from previous run for the restart."
            cp("EXPEQ.OUT", "EXPEQ"; force=true)
            # EXPEQ.OUT is in normalized CHEASE units; re-apply the SAME B0/R0 normalization
            # as the original equilibrium (from dd) so the restart does not fall back to the
            # CHEASE defaults (R0EXP=3.0, B0EXP=1.5) and silently rescale the equilibrium.
            B0, R0 = chease_normalization(dd)
            setfield!(chease_namelist, :B0EXP, B0)
            setfield!(chease_namelist, :R0EXP, R0)
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
   
    if par.num_surfaces == 1 && chease_namelist.NVEXP == 8
        @info "Overriding NVEXP in CHEASE namelist"
        NVEXP = 1
        setfield!(chease_namelist, :NVEXP, NVEXP)
    end
    
    # Write CHEASE namelist file
    CHEASE.write_CHEASEnamelist(chease_namelist, "datain")

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


"""
    validate_wall_configuration(par)

Check that `wall_type` and `num_surfaces` agree, in both directions, and that the implied
number of resistive walls is supported.

`num_surfaces` counts *bounding surfaces* (plasma boundary + walls), so the number of
resistive walls is `num_surfaces - 1` and MARS's `NWALL` must equal it — see
[`configure_wall!`](@ref), which derives `NWALL` from `num_surfaces` on that basis.

A resistive wall needs an actual wall contour in the equilibrium: `write_EXPEQ_file` only
builds one when `num_surfaces > 1`, so `wall_type != :no_wall` with the default
`num_surfaces == 1` would send CHEASE no wall at all and silently produce a no-wall
(far-field boundary at `REXT`) result instead of the requested resistive-wall one.
"""
function validate_wall_configuration(par)
    if par.wall_type == :no_wall && par.num_surfaces > 1
        error(
            "Invalid wall configuration: num_surfaces = $(par.num_surfaces) > 1 but " *
            "wall_type = :no_wall. Set a wall_type (:conformal or :limiter), or set " *
            "num_surfaces = 1."
        )
    end
    if par.wall_type != :no_wall && par.num_surfaces <= 1
        error(
            "Invalid wall configuration: wall_type = $(par.wall_type) requires a wall surface " *
            "in the equilibrium, but num_surfaces = $(par.num_surfaces). Set " *
            "num_surfaces > 1 so a wall contour is written to EXPEQ, or set " *
            "wall_type = :no_wall."
        )
    end
    # NWALL = num_surfaces - 1. Multi-wall (NWALL > 1) is not supported yet: CHEASE reports
    # a single wall index NW, so the positions of any further walls cannot be derived.
    if par.num_surfaces > 2
        error(
            "Unsupported wall configuration: num_surfaces = $(par.num_surfaces) implies " *
            "$(par.num_surfaces - 1) resistive walls, but only NWALL <= 1 is supported " *
            "(CHEASE reports a single wall index NW, so additional wall positions cannot be " *
            "determined). Set num_surfaces <= 2."
        )
    end
    return nothing
end


"""
    configure_wall!(mars_namelist, par)

Set the MARS resistive-wall inputs `NWALL`, `IWALL` and `TAUW` from `par`, so they can never
be left at struct defaults that do not correspond to the equilibrium actually being run.

- `wall_type == :no_wall` forces `NWALL = 0`, which is genuinely *no wall*. MARS's vacuum mesh
  runs from the plasma boundary (s=1) out to CHEASE's `REXT` (s=6 in the shipped examples), so
  the B1=0 condition BOWALL applies at the last grid point (marsq.f:3200) exists only to close
  the system and sits far enough out to approximate a wall at infinity.

  Forcing it matters because MARS validates `IWALL` only with the bounds test
  `IWALL(J) > 1 && IWALL(J) < NV` (marsq.f:2144) against that same mesh, which exists
  regardless of how many surfaces CHEASE computed. A leftover `NWALL = 1, IWALL = 22` is
  silently accepted and puts a resistive shell at s≈1.30 — just outside the plasma, nowhere
  near the far field — instead of being ignored.

- otherwise `IWALL` is re-derived from the CHEASE-reported wall index `NW` in `log_chease`
  (e.g. `NW=20` ↔ s=1.2308, matching CHEASE's own reported `RW`), and `TAUW` is taken from
  `par.wall_time_constants`.
"""
function configure_wall!(mars_namelist, par)
    basic = mars_namelist.BASIC

    if par.wall_type == :no_wall
        if basic.NWALL != 0
            @info "wall_type=:no_wall — forcing NWALL=0 (no wall; boundary condition at the far-field mesh edge REXT)"
        end
        setfield!(basic, :NWALL, 0)
        setfield!(basic, :IWALL, Int[])
        setfield!(basic, :TAUW, Float64[])
        return nothing
    end

    # Number of resistive walls is set by the equilibrium geometry, not by whatever length
    # IWALL happens to have: num_surfaces counts plasma boundary + walls. Guaranteed to be
    # exactly 1 here by validate_wall_configuration (num_surfaces == 2).
    n_walls = par.num_surfaces - 1

    # Wall position: CHEASE reports the wall grid index as `NW` in log_chease.
    iwall = copy(basic.IWALL)
    if isfile("log_chease")
        NW_dict = julia_grep(["NW"], "log_chease"; extract_values=true)
        if haskey(NW_dict, "NW")
            NW = NW_dict["NW"]
            @info "Setting IWALL from CHEASE wall index NW = $NW (log_chease)"
            isempty(iwall) ? push!(iwall, NW) : (iwall[1] = NW)
        else
            @warn "NW not found in log_chease; IWALL left at $iwall — verify it indexes the " *
                  "intended vacuum mesh position, MARS only bounds-checks it"
        end
    else
        @warn "log_chease not found (run_equilibrium=false?); IWALL left at $iwall — verify it " *
              "indexes the intended vacuum mesh position, MARS only bounds-checks it"
    end

    isempty(iwall) && error("wall_type=$(par.wall_type) but no wall position could be determined for IWALL")

    # Drop any extra user-supplied entries: CHEASE wrote only n_walls wall contour(s), and MARS
    # would read conductivity data for walls that do not exist in OUTVMAR.
    if length(iwall) > n_walls
        @warn "IWALL had $(length(iwall)) entries but num_surfaces=$(par.num_surfaces) " *
              "implies $n_walls wall(s); discarding extra entries $(iwall[n_walls+1:end])"
        iwall = iwall[1:n_walls]
    end

    # Wall time constants: a single value is broadcast to every wall.
    tauw = collect(float.(par.wall_time_constants))
    if length(tauw) == 1 && n_walls > 1
        tauw = fill(tauw[1], n_walls)
    elseif length(tauw) != n_walls
        error(
            "wall_time_constants has $(length(tauw)) entries but there are $n_walls wall " *
            "surface(s) (num_surfaces = $(par.num_surfaces)). Provide one τ_w per wall, " *
            "or a single value to apply to all."
        )
    end

    setfield!(basic, :NWALL, n_walls)
    setfield!(basic, :IWALL, iwall)
    setfield!(basic, :TAUW, tauw)

    @info "Resistive wall(s): NWALL=$(basic.NWALL), IWALL=$(basic.IWALL), TAUW=$(basic.TAUW)"
    return nothing
end


function run_MARS(dd::IMAS.DD, par, mars_namelist)

    # Resistive-wall settings (NWALL/IWALL/TAUW) are derived from par.wall_type and the
    # CHEASE-computed wall position, never left at struct defaults.
    configure_wall!(mars_namelist, par)

    # NREORBIT is driven entirely by par.num_orbits, so it can't drift out of sync with
    # the process count MARS is launched with below.
    setfield!(mars_namelist.REORBIT, :NREORBIT, par.num_orbits)

    # Write final MARS namelist into run directory for MARS execution
    write_MARS_namelist(mars_namelist, "RUN.IN")

    # Determine which profiles to pull from dd, i.e. experiment
    write_exp_profiles(dd, mars_namelist)

    # Execute MARS
    mars_exec = par.mars_exec
    @assert isfile("RUN.IN") "MARS input file RUN.IN not found"
    isfile(mars_exec) || error("MARS executable not found: $mars_exec")

    # MARS itself is serial; REORBIT particle tracing (reorbit.f) is the only part that
    # needs MPI, via a master rank plus one worker rank per orbit (see reorbit.f's
    # RE_TRACING task farm — requesting more ranks than num_orbits+1 leaves the extras
    # blocked forever, so nprocs must match exactly).
    nprocs = par.num_orbits > 0 ? par.num_orbits + 1 : 1

    if par.run_mode == :batch
        @info "Submitting MARS job to batch system with command: $(par.batch_submit_cmd) and script: mars_job.sh"
        ok = run_batch_mars(mars_exec, par, nprocs)
    elseif par.run_mode == :local
        run_cmd = nprocs > 1 ? `time $(par.mpirun_exec) -n $nprocs $(mars_exec)` : `time $(mars_exec)`
        @info "Running MARS interactively with command: $run_cmd"
        cmd = pipeline(
        run_cmd,
        stdout = "log_mars",
        stderr = "log_mars"
        )
        ok = success(cmd)
    else
        error("Unknown run mode: $(par.run_mode). Supported modes are :local and :batch.")
    end


    ok || error("MARS failed — see log_mars")

    # Display growth rate and iteration count
    iter, growth, freq = parse_MARS_results("RESULT.OUT")
    @info "MARS iterations = $iter"
    @info "Growth rate = $growth"
    @info "Frequency = $freq"

    # Collect results for storage into dd.mhd_linear (see _finalize).
    # growth/freq are Re(γ)·τ_A and Im(γ)·τ_A, normalized to the Alfvén time.
    n_tor = round(Int, mars_namelist.BASIC.RNTOR)
    ideal = iszero(mars_namelist.BASIC.ETA)

    # Mode structure (eigenfunction): needs the displacement harmonics (XPLASMA.OUT)
    # and the radial grid s (RMZM_F.OUT, written by MARS on the first parameter pass).
    mode = nothing
    if isfile("XPLASMA.OUT") && isfile("RMZM_F.OUT")
        m_pol, xi1, xi2, xi3 = read_MARS_eigenfunction("XPLASMA.OUT")
        s, R0EXP, RM, ZM, Ns1 = read_MARS_geometry("RMZM_F.OUT")
        chi, R, Z = mars_flux_surface_RZ(RM, ZM, R0EXP, Ns1)
        v1 = v2 = v3 = nothing
        if isfile("VPLASMA.OUT")
            _, v1, v2, v3 = read_MARS_eigenfunction("VPLASMA.OUT")
        end
        mode = MarsModeStructure(s, m_pol, xi1, xi2, xi3, v1, v2, v3, chi, R, Z)
    else
        @info "XPLASMA.OUT and/or RMZM_F.OUT not found; mode structure will not be stored in dd."
    end

    return MarsOutputs(n_tor, iter, growth, freq, ideal, mode)
end


"""
    read_MARS_eigenfunction(filename) -> (m_pol, c1, c2, c3)

Parse a MARS eigenfunction file (`XPLASMA.OUT` displacement, or `VPLASMA.OUT` velocity).

File layout:

    header:       MSMAX  NRP1  RNTOR  0 0 0
    MSMAX lines:  poloidal mode number m (repeated across columns)
    NRP1 lines:   dψ/ds (×3), T=R·Bφ (×3)   [equilibrium quantities, NOT the radial grid]
    MSMAX×NRP1:   Re,Im of the three complex contravariant components

Returns the poloidal Fourier mode numbers `m_pol` (length `MSMAX`) and three complex
matrices of size `(NRP1, MSMAX)`. The radial grid `s` is not in this file; read it from
`RMZM_F.OUT` with [`read_MARS_sgrid`](@ref).
"""
function read_MARS_eigenfunction(filename::AbstractString)
    isfile(filename) || error("File not found: $filename")
    lines = readlines(filename)

    hdr = split(lines[1])
    MSMAX = parse(Int, hdr[1])
    NRP1 = parse(Int, hdr[2])

    # poloidal mode numbers (first column of the MSMAX header rows)
    m_pol = [parse(Float64, split(lines[1+i])[1]) for i in 1:MSMAX]

    # complex contravariant components, ordered (MS outer, II inner) as written by MARS
    c1 = Matrix{ComplexF64}(undef, NRP1, MSMAX)
    c2 = similar(c1)
    c3 = similar(c1)
    idx = 1 + MSMAX + NRP1   # skip header + MSMAX poloidal rows + NRP1 equilibrium rows
    for ms in 1:MSMAX, ii in 1:NRP1
        idx += 1
        v = parse.(Float64, split(lines[idx]))
        c1[ii, ms] = complex(v[1], v[2])
        c2[ii, ms] = complex(v[3], v[4])
        c3[ii, ms] = complex(v[5], v[6])
    end

    return m_pol, c1, c2, c3
end

"""
    read_MARS_geometry(filename="RMZM_F.OUT") -> (s, R0EXP, RM, ZM, Ns1)

Read the MARS flux-surface geometry from the equilibrium file `RMZM_F.OUT`. Layout:

    header:        Nm0  NRP1  NVEQ1-1  R0EXP          (B0EXP in column 4 of row 2)
    Ns lines:      s (col 1), ..., q (col 3)          [Ns = NRP1 plasma + (NVEQ1-1) vacuum]
    Nm0×Ns lines:  Re,Im of R̂ₘ(s), Re,Im of Ẑₘ(s)    [R,Z Fourier harmonics, normalized by R0EXP]

Returns the plasma radial grid `s` (length `Ns1=NRP1`), the normalization `R0EXP`, the
complex R,Z Fourier harmonics `RM`/`ZM` of size `(Ns, Nm0)`, and the plasma point count `Ns1`.
"""
function read_MARS_geometry(filename::AbstractString="RMZM_F.OUT")
    isfile(filename) || error("File not found: $filename")
    lines = readlines(filename)

    h = split(lines[1])
    Nm0 = parse(Int, h[1])
    Ns1 = parse(Int, h[2])            # plasma radial points
    Ns2 = parse(Int, h[3])            # vacuum radial points
    R0EXP = parse(Float64, h[4])
    Ns = Ns1 + Ns2

    # plasma radial grid s = sqrt(ψ_pol_norm), column 1 of the first Ns1 rows
    s = [parse(Float64, split(lines[1+i])[1]) for i in 1:Ns1]

    # R,Z Fourier harmonics: Ns rows per harmonic, Nm0 harmonics, after the Ns radial rows
    RM = Matrix{ComplexF64}(undef, Ns, Nm0)
    ZM = Matrix{ComplexF64}(undef, Ns, Nm0)
    base = 1 + Ns
    k = 0
    for j in 1:Nm0, i in 1:Ns
        k += 1
        v = parse.(Float64, split(lines[base+k]))
        RM[i, j] = complex(v[1], v[2])
        ZM[i, j] = complex(v[3], v[4])
    end

    return s, R0EXP, RM, ZM, Ns1
end

"""
    mars_flux_surface_RZ(RM, ZM, R0EXP, Ns1; nchi=257) -> (chi, R, Z)

Reconstruct the real-space flux-surface geometry `R(s,χ)`, `Z(s,χ)` [m] over the plasma
region from the MARS R,Z Fourier harmonics.

`R̂ₘ`/`Ẑₘ` are the *one-sided* (m ≥ 0) Fourier spectrum of the real curves R(χ), Z(χ),
so the physical reconstruction carries the real-Fourier factor of 2 on the m ≥ 1 harmonics:

    R(s,χ) = R0EXP · [ Re(R̂₀) + 2·Σ_{m≥1} Re(R̂ₘ(s) e^{i m χ}) ],   χ ∈ [0, 2π]

(Note: this differs from the package's `MacGetRZ.m`, which omits the factor of 2; with the
factor included the reconstructed boundary matches the equilibrium boundary fed to CHEASE.)

Returns the poloidal-angle grid `chi` (length `nchi`) and `R`, `Z` of size `(Ns1, nchi)`.
"""
function mars_flux_surface_RZ(RM::Matrix{ComplexF64}, ZM::Matrix{ComplexF64}, R0EXP::Real, Ns1::Int; nchi::Int=257)
    Nm0 = size(RM, 2)
    m = 0:(Nm0-1)
    chi = collect(range(0.0, 2π; length=nchi))
    # one-sided real-Fourier basis: weight 1 for m=0, 2 for m≥1
    expmchi = [(mm == 0 ? 1.0 : 2.0) * exp(im * mm * cc) for mm in m, cc in chi]   # (Nm0 × nchi)
    R = real.(view(RM, 1:Ns1, :) * expmchi) .* R0EXP       # (Ns1 × nchi)
    Z = real.(view(ZM, 1:Ns1, :) * expmchi) .* R0EXP
    return chi, R, Z
end


function run_batch_mars(mars_exec, par, nprocs)

    run_line = nprocs > 1 ? "time $(par.mpirun_exec) -n $nprocs $mars_exec" : "time $mars_exec"

    script = """
    #!/bin/bash
    #SBATCH --job-name=mars_job
    #SBATCH --ntasks=$nprocs
    #SBATCH --output=log_mars
    #SBATCH --error=log_mars

    cd $(pwd())

    $run_line
    """

    script_file = "mars_job.sh"
    open(script_file, "w") do io
        write(io, script)
    end

    run(`chmod +x $script_file`)

    run(`$(par.batch_submit_cmd) $script_file`)

    return true
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

    format_namelist_value(v) =
        if v isa Complex
            "($(real(v)), $(imag(v)))"
        elseif v isa AbstractVector{<:Complex}
            join(["($(real(z)), $(imag(z)))" for z in v], ", ")
        elseif v isa AbstractVector
            join(v, ", ")
        else
            v
        end

    # Group MARSnamelist fields by the Fortran namelist name they must be written under
    # (NAMELIST_GROUP), in order of first appearance. Most Julia blocks map 1:1 to a
    # namelist of the same name; REORBIT is folded into OUTOPT since MARS has no
    # `&REORBIT` namelist of its own.
    group_order = Symbol[]
    for block in fieldnames(typeof(nl))
        g = NAMELIST_GROUP[block]
        g in group_order || push!(group_order, g)
    end
    blocks_in_group(group) = (getfield(nl, block) for block in fieldnames(typeof(nl)) if NAMELIST_GROUP[block] == group)

    # An empty vector has no valid Fortran namelist spelling (`IWALL = ,` is a syntax error
    # that fails the READ). Omitting the variable entirely is the correct encoding of "unset":
    # it keeps the Fortran-side default, and empty vectors only arise for inputs that are not
    # read at all in that configuration (e.g. IWALL/TAUW when NWALL=0).
    skip_value(v) = v isa AbstractVector && isempty(v)

    open(filename,"w") do io

        for group in group_order

            println(io,"&",group)

            for b in blocks_in_group(group)

                for k in fieldnames(typeof(b))
                    k == :extras && continue
                    v = getfield(b,k)
                    skip_value(v) && continue
                    println(io," $k = $(format_namelist_value(v)),")
                end

                for (k,v) in b.extras
                    skip_value(v) && continue
                    println(io," $k = $(format_namelist_value(v)),")
                end
            end

            println(io,"&END\n")
        end
    end
end



function write_exp_profiles(dd::IMAS.DD, mars_namelist)
    @info "Generating additional PROF*.IN files from experiment or dd."

    profiles = dd.core_profiles.profiles_1d[]
    s = sqrt.(profiles.grid.psi_norm)  # radial coordinate

    # core_transport has no single "chi_perp"/"chi_par" node; both are read off the
    # :combined transport model (sum of all transport models) as the closest DD proxy:
    #   chi_perp (NPROFTTCE) <- ion channel   (total_ion_energy.d, summed over ion species)
    #   chi_par  (NPROFTTCA) <- electron channel (electrons.energy.d), since parallel
    #                           thermal conduction is electron-dominated
    # Evaluated lazily (only if the corresponding flag is actually set to 4 below), so
    # this does not require a :combined model to exist in dd when these profiles are
    # left at their MARS-internal default.
    combined_1d() = dd.core_transport.model[:combined].profiles_1d[]

    from_dd = Symbol[]
    from_mars_default = Symbol[]

    for (flag, outputs, msg) in (
        (
            :NPROFR,
            [("PROFROT.IN", p -> p.ion[1].rotation_frequency_tor)],
            nothing
        ),

        (
            :NPROFN,
            [
                ("PROFDEN.IN", p -> p.ion[1].density),
                ("PROFNE.IN", p -> p.electrons.density)
            ],
            nothing
        ),

        (
            :NPROFWE,
            [("PROFWE.IN", p -> p.rotation_frequency_tor_sonic)],
            nothing
        ),

        (
            :NPROFTTCA,
            [("PROFTTCPARA.IN", _ -> combined_1d().electrons.energy.d)],
            "χ_parallel (thermal) has no dedicated dd node. Using electron energy diffusivity " *
            "from the combined transport model (core_transport.model[:combined].profiles_1d[].electrons.energy.d) " *
            "as a stand-in — not a lookup, an interpretation."
        ),

        (
            :NPROFTTCE,
            [("PROFTTCPERP.IN", _ -> combined_1d().total_ion_energy.d)],
            "χ_perpendicular has no dedicated dd node. Using combined-model total ion energy diffusivity " *
            "(core_transport.model[:combined].profiles_1d[].total_ion_energy.d) as a stand-in — not a lookup, an interpretation."
        ),

        (
            :NPROFIE,
            [
                ("PROFTI.IN", p -> p.ion[1].temperature),
                ("PROFTE.IN", p -> p.electrons.temperature)
            ],
            nothing
        )

    )

        if getfield(mars_namelist.BASIC, flag) == 4
            push!(from_dd, flag)

            msg !== nothing && @info msg

            for (file, getter) in outputs
                profile = getter(profiles)
                profile ./= maximum(profile)
                write_profile_IN(file, s, profile)
            end
        else
            push!(from_mars_default, flag)
        end
    end

    @info "Profiles written from dd: $from_dd"
    if !isempty(from_mars_default)
        @warn "Profiles NOT written from dd — MARS will use its own internal analytic profile for: $from_mars_default"
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

function run_PARTICLE_TRACING(dd::IMAS.DD, par)
    # Placeholder function to run particle tracing simulations
    @info "Running particle tracing with REORBIT."

    println("Particle tracing simulation completed.")
end


"""
  Writes a EXPEQ file for CHEASE given dd and mode number
    This function will create the necessary input file for the CHEASE equilibrium solver.
        NWBPS: Number of wall boundary points
        NSTTP: Number of steps in pressure and current profiles
"""
function write_EXPEQ_file(dd::IMAS.DD, par)

    offset = par.offset  # offset for first wall (RW) in meters
    n_points = par.n_points  # number of points for first wall (RW)
    NWBPS = par.num_surfaces
    
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

    # choose B0 & R0 for CHEASE normalization (shared with the restart path)
    B0, R0 = chease_normalization(dd)

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
        GS_RHS = FFpr / B0
    elseif par.GS_rhs == :Jtor
        NSTTP = 2
        Jtor = abs.(eqt1d.j_tor)
        GS_RHS = Jtor / (B0 / R0 * μ_0)
    elseif par.GS_rhs == :Jpar
        NSTTP = 3
        Jpar = abs.(eqt1d.j_parallel) # NOT right!
        GS_RHS = Jpar / (B0 / R0 * μ_0)
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
    
    chease_struct = CHEASE.MartianCHEASE(
        ϵ,
        z_geo/R0,
        pressure_sep,
        B0,
        R0,
        Ip,
        rb_new,
        zb_new,
        NSTTP,
        s,
        GS_RHS,
        pprime,
        NWBPS,
        NDATA,
        missing,
        missing,
    )

    r_bound_norm = rb_new / R0
    z_bound_norm = zb_new / R0

    # Overlay onto subplot 1 (the R,Z flux-contour panel) of the "before CHEASE" equilibrium
    # figure already created in _step, via Plots.jl's current-plot state — not a separate figure.
    par.do_plot && display(plot!(rb_new, zb_new; subplot=1, linewidth=3., aspect_ratio=:equal, label="Plasma Boundary"))

    @assert length(rb_new) == length(zb_new) "R,Z boundary arrays must have the same shape"

    # if there is another surface, calclate its cooridanes given an offset and save to file
    if NWBPS > 1 ## WHAT TO DO if > 2
        if par.wall_type == :conformal
            @info "Creating a conformal limiter offset from plasma boundary by $(par.offset)."
            r_lim, z_lim = offset_boundary(rb_new, zb_new, offset)
        elseif par.wall_type == :limiter
            @info "Using wall data .Json for CHEASE equilibrium generation."
            machine = dd.dataset_description.data_entry.machine
            limiter = get_limiter_data(machine)
            r_lim, z_lim = limiter.r, limiter.z
        end


        r_lim, z_lim = IMAS.resample_2d_path(r_lim, z_lim; n_points=n_points, method=:linear)

        chease_struct.r_limiter = r_lim
        chease_struct.z_limiter = z_lim

        r_lim_norm = r_lim / R0
        z_lim_norm = z_lim / R0

        par.do_plot && display(plot!(r_lim, z_lim; subplot=1, linewidth=1.5, aspect_ratio=:equal, label="MARS resistive wall"))
    end
  
    # write to EXPEQ file
    CHEASE.write_EXPEQ_file(chease_struct)   
    
    return B0, R0
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


