#= =========== =#
#  ActorSOLPSNN  #
#= =========== =#
import SOLPSNN

# SOLPS2imas is a declared FUSE dependency but is loaded *lazily* (only when
# grid=:solps2imas is requested). It is NOT imported at top level on purpose:
# SOLPS2imas pulls in IMASggd, whose plot recipes for IMASdd.interferometer
# collide with IMAS's, and method-overwriting is fatal during precompilation.
# A runtime load turns that collision into a harmless warning and keeps FUSE
# precompiling cleanly. See _solpsnn_load_solps2imas!.
const _SOLPS2IMAS_PKGID = Base.PkgId(Base.UUID("09becab6-0636-4c23-a92a-2b3723265c31"), "SOLPS2imas")

"""Load SOLPS2imas at runtime (activating SOLPSNN's SOLPS2imasExt) if not already loaded."""
function _solpsnn_load_solps2imas!()
    if !haskey(Base.loaded_modules, _SOLPS2IMAS_PKGID)
        try
            Base.require(_SOLPS2IMAS_PKGID)
        catch err
            error("ActorSOLPSNN: grid=:solps2imas needs SOLPS2imas.jl, which failed to load: $err")
        end
    end
    isempty(methods(SOLPSNN.build_edge_profiles_ggd_solps2imas!)) &&
        error("ActorSOLPSNN: SOLPS2imas loaded but SOLPSNN.SOLPS2imasExt did not activate")
    return nothing
end

# Map an actor-level quantity symbol onto a (SOLPS-NN item, species) pair.
const _SOLPSNN_QMAP = Dict{Symbol,Tuple{String,Union{Nothing,String}}}(
    :te => ("te", nothing),        # electron temperature field
    :ti => ("ti", nothing),        # ion temperature field
    :ne => ("na", "D1"),           # main-ion (D+) density ~ electron density (quasineutral)
    :pwmxap => ("pwmxap", nothing), # peak outer-target power flux [W/m^2]
    :fnixap => ("fnixap", nothing), # integrated D+ ion flux to outer target [1/s]
    :psol => ("psol", nothing)     # power across separatrix [W] (diagnostic)
)

# Approximate SOLPS-NN training ranges (JET-scaled B2 database) used only to warn
# when the FUSE-derived inputs fall outside the region the surrogate was trained on.
const _SOLPSNN_BOUNDS = (
    R=(1.0, 10.0),          # major radius [m]
    B=(1.0, 10.0),          # |B0| at R0 [T]
    P=(1.0e6, 3.0e8),       # power into SOL [W]
    D_puff=(1.0e18, 1.0e24),
    N_puff=(1.0e18, 1.0e23),
    D_core=(1.0e18, 1.0e24),
    D_perp=(0.05, 3.0),
    chi_perp=(0.05, 5.0)
)

@actor_parameters_struct ActorSOLPSNN{T} begin
    quantities::Entry{Vector{Symbol}} =
        Entry{Vector{Symbol}}("-", "SOLPS-NN quantities to evaluate (subset of $(sort(collect(keys(_SOLPSNN_QMAP)))))";
            default=[:te, :ti, :ne, :pwmxap, :psol])
    Psol::Entry{T} = Entry{T}("W", "Power across the separatrix (defaults to IMAS.power_sol when not set)")
    D_puff::Entry{T} = Entry{T}("particles/s", "Main-ion (D) gas puff rate"; default=1.0e22)
    N_puff::Entry{T} = Entry{T}("particles/s", "Impurity (N) seeding rate (must be > 0; model uses log10)"; default=1.0e20)
    D_core::Entry{T} = Entry{T}("particles/s", "Main-ion core fuelling rate"; default=1.0e22)
    D_perp::Entry{T} = Entry{T}("m^2/s", "Cross-field particle diffusivity"; default=0.3)
    chi_perp::Entry{T} = Entry{T}("m^2/s", "Cross-field heat diffusivity"; default=1.0)
    write_ggd::Entry{Bool} = Entry{Bool}("-", "Write 2D fields to dd.edge_profiles GGD"; default=true)
    grid::Switch{Symbol} = Switch{Symbol}([:native, :solps2imas], "-",
        "GGD grid builder: :native (lightweight SOLPSNN builder) or :solps2imas (SOLPS2imas.jl, "*
        "reproduces a real SOLPS run's mesh with the full set of physical subsets)"; default=:native)
    b2fgmtry::Entry{String} = Entry{String}("-",
        "Path to the SOLPS b2fgmtry for grid=:solps2imas (defaults to the file bundled with SOLPSNN)"; default="")
    write_divertor::Entry{Bool} = Entry{Bool}("-", "Write peak heat flux to dd.divertors outer target"; default=true)
    dir::Entry{String} = Entry{String}("-", "Override directory holding the converted ONNX artifacts"; default="")
    verify::Entry{Bool} = Entry{Bool}("-", "SHA-256 verify artifacts on load"; default=false)
    warn_out_of_bounds::Entry{Bool} = Entry{Bool}("-", "Warn when inputs fall outside the training range"; default=true)
    do_plot::Entry{Bool} = act_common_parameters(; do_plot=false)
    verbose::Entry{Bool} = act_common_parameters(; verbose=false)
end

mutable struct ActorSOLPSNN{D,P} <: SingleAbstractActor{D,P}
    dd::IMAS.dd{D}
    par::OverrideParameters{P,FUSEparameters__ActorSOLPSNN{P}}
    models::Dict{Symbol,SOLPSNN.SOLPSModel}
    geo::Union{Nothing,SOLPSNN.SOLPSGeometry}
    predictions::Dict{Symbol,Any}
    inputs::Vector{Float64}
end

"""
    ActorSOLPSNN(dd::IMAS.dd, act::ParametersAllActors; kw...)

Evaluates the SOLPS-NN edge-plasma surrogate (Dasbach & Wiesen) to predict the
divertor/SOL response for the current equilibrium and SOL power.

The surrogate maps 8 scalar inputs `[R, B, Psol, D_puff, N_puff, D_core, D_perp, chi_perp]`
onto either scalar boundary quantities (peak outer-target heat flux `pwmxap`,
integrated ion flux `fnixap`, power across the separatrix `psol`) or full 2D
fields on the fixed SOLPS-ITER B2 grid (`te`, `ti`, `ne`), rescaled to the
machine major radius. `R` and `B` are taken from `dd.equilibrium`; `Psol` defaults to `IMAS.power_sol`
when the `Psol` parameter is not set; the gas-puff/transport inputs are actor parameters.

2D fields are written to `dd.edge_profiles` as a General Grid Description (GGD)
slice; the peak heat flux is written to the outer divertor target when target
geometry is available. `psol` is kept as a diagnostic (compared against
`IMAS.power_sol`) and is *not* written back, to avoid corrupting the Psol/R
design constraint that FUSE computes from core power balance.

The GGD grid is built either with SOLPSNN's lightweight internal builder
(`grid=:native`, default) or through `SOLPS2imas.jl` (`grid=:solps2imas`), which
reproduces the exact mesh and full set of physical grid subsets (separatrix,
inner/outer target, OMP/IMP, core/SOL, …) of a real SOLPS run — making the
output a drop-in for GGDUtils / SOLPS2ctrl tooling.

!!! note

    Reads `dd.equilibrium`, `dd.core_profiles`, `dd.core_sources`.
    Stores 2D fields in `dd.edge_profiles` and peak heat flux in `dd.divertors`.
    Upstream SOLPS-NN ships TensorFlow weights, not ONNX, so on first use the
    bundled `convert/` pipeline fetches the needed quantities from SURFdrive and
    converts them to ONNX in a conda env (created on demand; needs `conda` on
    `PATH`, e.g. `module load conda` or the FUSE conda env). Artifacts cache
    under `\$PSCRATCH/solps-nn-onnx`; override with `ENV["FUSE_SOLPS_NN_DIR"]`
    or the `dir` parameter, or disable auto-conversion with
    `ENV["FUSE_SOLPS_NN_AUTOCONVERT"] = "0"`.
"""
function ActorSOLPSNN(dd::IMAS.dd, act::ParametersAllActors; kw...)
    actor = ActorSOLPSNN(dd, act.ActorSOLPSNN; kw...)
    step(actor)
    finalize(actor)
    return actor
end

function ActorSOLPSNN(dd::IMAS.dd{D}, par::FUSEparameters__ActorSOLPSNN{P}; kw...) where {D<:Real,P<:Real}
    logging_actor_init(ActorSOLPSNN)
    par = OverrideParameters(par; kw...)
    return ActorSOLPSNN{D,P}(dd, par, Dict{Symbol,SOLPSNN.SOLPSModel}(), nothing, Dict{Symbol,Any}(), Float64[])
end

"""Resolve the artifact directory honouring the `dir` parameter override."""
_solpsnn_dir(par) = isempty(par.dir) ? SOLPSNN.resolve_dir() : SOLPSNN.resolve_dir(; dir=par.dir)

"""Lazily load (and cache) the SOLPS-NN model for quantity `q`."""
function _solpsnn_model!(actor::ActorSOLPSNN, q::Symbol)
    haskey(actor.models, q) && return actor.models[q]
    haskey(_SOLPSNN_QMAP, q) || error("ActorSOLPSNN: unknown quantity :$q (valid: $(sort(collect(keys(_SOLPSNN_QMAP)))))")
    item, species = _SOLPSNN_QMAP[q]
    dir = isempty(actor.par.dir) ? nothing : actor.par.dir
    m = SOLPSNN.load_model(item; species, dir, verify=actor.par.verify)
    actor.models[q] = m
    return m
end

function _step(actor::ActorSOLPSNN)
    dd = actor.dd
    par = actor.par

    eqt = dd.equilibrium.time_slice[]

    R = eqt.boundary.geometric_axis.r
    B = abs(eqt.global_quantities.vacuum_toroidal_field.b0)
    Psol_par = getproperty(par, :Psol, missing)
    Psol = ismissing(Psol_par) ? IMAS.power_sol(dd.core_sources, dd.core_profiles.profiles_1d[]) : Float64(Psol_par)

    X = Float64[R, B, Psol, par.D_puff, par.N_puff, par.D_core, par.D_perp, par.chi_perp]
    actor.inputs = X

    # the surrogate takes log10 of the puff/fuelling rates -> require strictly positive
    for (name, val) in (("D_puff", par.D_puff), ("N_puff", par.N_puff), ("D_core", par.D_core))
        val > 0 || error("ActorSOLPSNN: $name must be > 0 (model applies log10), got $val")
    end

    if par.warn_out_of_bounds
        for (name, val) in zip(keys(_SOLPSNN_BOUNDS), X)
            lo, hi = getfield(_SOLPSNN_BOUNDS, name)
            if val < lo || val > hi
                @warn "ActorSOLPSNN: input $name=$(val) outside SOLPS-NN training range [$lo, $hi]"
            end
        end
    end

    empty!(actor.predictions)
    for q in par.quantities
        m = _solpsnn_model!(actor, q)
        actor.predictions[q] = SOLPSNN.predict(m, X)  # scalar or (nx+2, ny+2) matrix
    end

    if par.verbose
        @info "ActorSOLPSNN inputs [R,B,Psol,D_puff,N_puff,D_core,D_perp,chi_perp] = $X"
        if haskey(actor.predictions, :psol)
            @info "ActorSOLPSNN: surrogate psol=$(actor.predictions[:psol]) W vs FUSE Psol=$Psol W"
        end
    end

    return actor
end

function _finalize(actor::ActorSOLPSNN)
    dd = actor.dd
    par = actor.par
    preds = actor.predictions

    R = isempty(actor.inputs) ? dd.equilibrium.time_slice[].boundary.geometric_axis.r : actor.inputs[1]
    time0 = dd.global_time

    field_qs = [q for q in keys(preds) if preds[q] isa AbstractMatrix]

    if par.write_ggd && !isempty(field_qs)
        if par.grid == :solps2imas
            # Build the grid through SOLPS2imas so the surrogate populates the same
            # IMAS mesh (+ full physical subsets) a real SOLPS run does. The cell
            # ordering is transposed vs the native builder, hence order=:solps.
            _solpsnn_load_solps2imas!()
            b2 = isempty(par.b2fgmtry) ? SOLPSNN.bundled_b2fgmtry() : par.b2fgmtry
            # invokelatest: the ext method may have been defined in a newer world age
            grid_index = Base.invokelatest(SOLPSNN.build_edge_profiles_ggd_solps2imas!, dd, b2; R, R_JET=SOLPSNN.R_JET, time0)
            ggd_order = :solps
        else
            if actor.geo === nothing
                actor.geo = SOLPSNN.load_geometry(_solpsnn_dir(par))
            end
            grid_index = SOLPSNN.build_edge_profiles_ggd!(dd, actor.geo; R, time0)
            ggd_order = :native
        end
        ep = SOLPSNN.ggd_time_slice!(dd, time0)

        if haskey(preds, :te) && preds[:te] isa AbstractMatrix
            SOLPSNN.add_ggd_field!(ep.electrons.temperature, preds[:te]; grid_index, order=ggd_order)
        end
        if haskey(preds, :ne) && preds[:ne] isa AbstractMatrix
            SOLPSNN.add_ggd_field!(ep.electrons.density, preds[:ne]; grid_index, order=ggd_order)
        end
        if haskey(preds, :ti) && preds[:ti] isa AbstractMatrix
            resize!(ep.ion, 1)
            ep.ion[1].label = "D+"
            ep.ion[1].z_ion = 1.0
            SOLPSNN.add_ggd_field!(ep.ion[1].temperature, preds[:ti]; grid_index, order=ggd_order)
        end
    end

    if par.write_divertor && haskey(preds, :pwmxap) && preds[:pwmxap] isa Real
        _solpsnn_write_peak_heat_flux!(dd, Float64(preds[:pwmxap]); time0, verbose=par.verbose)
    end

    if par.do_plot
        _solpsnn_plot(actor)
    end

    return actor
end

"""
Write the SOLPS-NN peak heat flux to the outer divertor target's
`power_flux_peak` (best-effort; mirrors ActorDivertors' outer=first identifier
convention). Skips gracefully when no divertor target geometry is available.
"""
function _solpsnn_write_peak_heat_flux!(dd::IMAS.dd, qpeak::Float64; time0::Float64, verbose::Bool)
    isempty(dd.divertors.divertor) && return false
    try
        eqt = dd.equilibrium.time_slice[]
        sol1 = IMAS.sol(eqt, dd.wall; levels=1)[:lfs][1]
        identifiers = IMAS.identify_strike_surface(sol1, dd.divertors)
        if !isempty(identifiers) && identifiers[1] != (0, 0, 0)
            k_div, k_tgt = identifiers[1]
            target = dd.divertors.divertor[k_div].target[k_tgt]
            @ddtime(target.power_flux_peak.data = qpeak)
            return true
        end
    catch err
        verbose && @warn "ActorSOLPSNN: could not write pwmxap to a divertor target: $err"
    end
    verbose && @info "ActorSOLPSNN: pwmxap=$qpeak W/m^2 (not written to dd.divertors; no target geometry)"
    return false
end

"""Quick heatmap of the first available 2D field over the B2 index space."""
function _solpsnn_plot(actor::ActorSOLPSNN)
    for (q, v) in actor.predictions
        if v isa AbstractMatrix
            display(heatmap(permutedims(v); xlabel="poloidal index", ylabel="radial index",
                title="SOLPS-NN $(q)", color=:turbo))
            break
        end
    end
    return nothing
end
