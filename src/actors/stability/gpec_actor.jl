import GeneralizedPerturbedEquilibrium as GPEC
import TOML

#ActorGPEC      
@actor_parameters_struct ActorGPEC{T} begin
    # Toroidal mode numbers
    nn_low::Entry{Int} = Entry{Int}("-", "Lowest toroidal mode number to analyze"; default=1)
    nn_high::Entry{Int} = Entry{Int}("-", "Highest toroidal mode number to analyze"; default=3)

    # Solver flags
    vac_flag::Entry{Bool} = Entry{Bool}("-", "Compute vacuum (free-boundary) response"; default=true)
    ode_flag::Entry{Bool} = Entry{Bool}("-", "Use ODE solver for plasma response"; default=true)
    mat_flag::Entry{Bool} = Entry{Bool}("-", "Compute FGK matrices"; default=true)
    mer_flag::Entry{Bool} = Entry{Bool}("-", "Compute Mercier criterion"; default=false)
    bal_flag::Entry{Bool} = Entry{Bool}("-", "Compute ballooning stability"; default=false)

    # Radial domain truncation
    psilow::Entry{T} = Entry{T}("-", "Lower radial boundary (normalized flux, 0=axis, 1=edge)"; default=0.01)
    psihigh::Entry{T} = Entry{T}("-", "Upper radial boundary (normalized flux, 0=axis, 1=edge)"; default=0.994)

    # Grid resolution
    mpsi::Entry{Int} = Entry{Int}("-", "Number of radial grid points"; default=64)
    mtheta::Entry{Int} = Entry{Int}("-", "Number of poloidal grid points"; default=256)

    # Coordinate system
    jac_type::Switch{String} = Switch{String}(["boozer", "sfl"], "-", "Jacobian type: boozer or straight-field-line (sfl)"; default="boozer")

    # IMAS/COCOS convention
    imas_cocos::Entry{Int} = Entry{Int}("-", "IMAS COCOS convention (11=IMAS standard, 2=GPEC internal)"; default=11)

    # Wall configuration
    wall_shape::Switch{String} = Switch{String}(["conformal", "nowall"], "-", "Wall geometry: conformal to plasma or nowall"; default="conformal")
    wall_distance::Entry{T} = Entry{T}("m", "Wall distance from plasma (for conformal wall)"; default=0.2)

    # q-surface search lower bound (q must exceed this value somewhere in domain)
    qlow::Entry{T} = Entry{T}("-", "Minimum q value for rational surface search"; default=0.0)

    # Diagnostics and output
    do_plot::Entry{Bool} = act_common_parameters(; do_plot=false)
    verbose::Entry{Bool} = act_common_parameters(; verbose=false)
end

mutable struct ActorGPEC{D,P} <: AbstractActor{D,P}
    dd::IMAS.dd{D}
    par::OverrideParameters{P,FUSEparameters__ActorGPEC{P}}
end

function ActorGPEC(dd::IMAS.dd{D}, par::FUSEparameters__ActorGPEC{P}; kw...) where {D<:Real,P<:Real}
    logging_actor_init(ActorGPEC)
    par = OverrideParameters(par; kw...)
    return ActorGPEC{D,P}(dd, par)
end

"""
    ActorGPEC(dd::IMAS.dd, act::ParametersAllActors; kw...)

Runs GPEC ideal MHD stability analysis on the equilibrium in `dd.equilibrium`.

# Physics Computed
- Vacuum energy eigenvalues for toroidal modes n = nn_low:nn_high
- Plasma response to 3D magnetic perturbations
- Optional: Mercier criterion, ballooning stability

# Key Parameters
- `nn_low`, `nn_high`: Toroidal mode number range
- `vac_flag`: Compute free-boundary (vacuum) response
- `psilow`, `psihigh`: Radial domain (edge truncation affects marginal stability)
- `mpsi`, `mtheta`: Grid resolution
- `imas_cocos`: COCOS convention (11=IMAS standard, FUSE default)

# COCOS Note
FUSE uses COCOS 11 (IMAS standard). GPEC internally uses COCOS 2. The conversion
(dividing poloidal flux by 2π) is handled automatically. For plasmas near marginal
stability, numerical sensitivity is inherent to the physics, not the conversion.

!!! note

    Stores results in `dd.mhd_linear.time_slice[].toroidal_mode[].energy_perturbed`
"""
function ActorGPEC(dd::IMAS.dd, act::ParametersAllActors; kw...)
    actor = ActorGPEC(dd, act.ActorGPEC; kw...)
    step(actor)
    finalize(actor)
    return actor
end

function _step(actor::ActorGPEC)
    dd = actor.dd
    par = actor.par

    # Verify equilibrium exists
    if length(dd.equilibrium.time_slice) == 0
        error("ActorGPEC: No equilibrium found in dd.equilibrium.time_slice")
    end

    # Create temporary directory for GPEC run
    tmpdir = mktempdir(prefix="fuse_gpec_")

    try
        # Build GPEC configuration from actor parameters
        config = Dict(
            "Equilibrium" => Dict(
                "eq_type" => "imas",
                "eq_filename" => "from_fuse_dd",  # placeholder, dd is passed directly
                "jac_type" => par.jac_type,
                "psilow" => par.psilow,
                "psihigh" => par.psihigh,
                "mpsi" => par.mpsi,
                "mtheta" => par.mtheta,
                "imas_cocos" => par.imas_cocos
            ),
            "ForceFreeStates" => Dict(
                "nn_low" => par.nn_low,
                "nn_high" => par.nn_high,
                "vac_flag" => par.vac_flag,
                "ode_flag" => par.ode_flag,
                "mat_flag" => par.mat_flag,
                "mer_flag" => par.mer_flag,
                "bal_flag" => par.bal_flag,
                "qlow" => par.qlow,
                "write_outputs_to_HDF5" => true,
                "HDF5_filename" => "gpec.h5",
                "verbose" => par.verbose
            ),
            "Wall" => Dict(
                "shape" => par.wall_shape,
                "a" => par.wall_distance
            )
        )

        # Write configuration to temporary directory
        open(joinpath(tmpdir, "gpec.toml"), "w") do io
            TOML.print(io, config)
        end

        # Run GPEC
        if par.verbose
            @info "ActorGPEC: Running GPEC in $tmpdir with modes n=$(par.nn_low):$(par.nn_high)"
        end

        # Call GPEC main - it reads equilibrium from dd, returns results
        result = GPEC.main([tmpdir]; dd=dd)

        # Write results to dd.mhd_linear
        GPEC.write_imas(dd, result)

        if par.verbose
            if result.vac_data !== nothing && length(result.vac_data.et) > 0
                @info "ActorGPEC: Eigenvalues computed for $(result.intr.npert) toroidal modes"
                @info "ActorGPEC: Dominant eigenvalue = $(result.vac_data.et[1])"
                if real(result.vac_data.et[1]) < 0
                    @warn "ActorGPEC: Free-boundary mode UNSTABLE (eigenvalue < 0)"
                else
                    @info "ActorGPEC: Free-boundary modes STABLE (all eigenvalues ≥ 0)"
                end
            end
        end

    catch e
        @error "ActorGPEC: GPEC run failed" exception=(e, catch_backtrace())
        rethrow(e)
    finally
        # Clean up temporary directory
        if isdir(tmpdir)
            rm(tmpdir, recursive=true)
        end
    end

    return actor
end

function _finalize(actor::ActorGPEC)
    dd = actor.dd

    # Validate results were written
    if length(dd.mhd_linear.time_slice) == 0
        @warn "ActorGPEC: No results written to dd.mhd_linear.time_slice"
    end

    return actor
end
