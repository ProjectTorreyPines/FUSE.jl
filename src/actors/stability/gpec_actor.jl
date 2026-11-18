#= ========= =#
#  ActorGPEC  #
#= ========= =#
@actor_parameters_struct ActorGPEC{T} begin
    # Toroidal mode numbers
    nn_low::Entry{Int} = Entry{Int}("-", "Lowest toroidal mode number to analyze"; default=1)
    nn_high::Entry{Int} = Entry{Int}("-", "Highest toroidal mode number to analyze"; default=1)

    # Solver flags
    vac_flag::Entry{Bool} = Entry{Bool}("-", "Compute vacuum (free-boundary) response"; default=true)
    local_stability_flag::Entry{Bool} = Entry{Bool}("-", "Compute Mercier and ballooning stability across the ψ profile"; default=false)

    # Radial domain truncation
    psilow::Entry{T} = Entry{T}("-", "Lower radial boundary (normalized flux, 0=axis, 1=edge)"; default=1e-4)
    psihigh::Entry{T} = Entry{T}("-", "Upper radial boundary (normalized flux, 0=axis, 1=edge)"; default=0.994)
    psiedge::Entry{T} = Entry{T}("-", "Start of the edge δW(ψ) scan over [psiedge, psihigh]; set ≥ psihigh to disable"; default=0.99)

    # Grid resolution
    mpsi::Entry{Int} = Entry{Int}("-", "Number of radial grid points"; default=64)
    mtheta::Entry{Int} = Entry{Int}("-", "Number of poloidal grid points"; default=256)

    # Coordinate system
    jac_type::Switch{String} = Switch{String}(["hamada", "pest", "equal_arc", "boozer", "park"], "-", "Jacobian coordinate type"; default="hamada")

    # IMAS/COCOS convention
    imas_cocos::Entry{Int} = Entry{Int}("-", "IMAS COCOS convention (11=IMAS standard, 2=GPEC internal)"; default=11)

    # Wall configuration
    wall_shape::Switch{String} = Switch{String}(["conformal", "nowall"], "-", "Wall geometry: conformal to plasma or nowall"; default="nowall")
    wall_a::Entry{T} = Entry{T}("-", "Conformal wall distance from the plasma in units of the minor radius 0.5(max R - min R)")

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

Evaluates ideal MHD stability of the equilibrium in `dd.equilibrium` with GPEC.
Requires GeneralizedPerturbedEquilibrium.jl to be loaded (`using GeneralizedPerturbedEquilibrium`).

# Physics computed
- Least-stable total energy δW for each toroidal mode n = nn_low:nn_high (δW < 0 means unstable)

# Key parameters
- `nn_low`, `nn_high`: toroidal mode number range
- `psilow`, `psihigh`: radial domain in normalized poloidal flux
- `psiedge`: start of the edge δW scan
- `wall_shape`, `wall_a`: vacuum wall; `wall_a` is required for a conformal wall

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

    # Verify equilibrium exists
    if length(dd.equilibrium.time_slice) == 0
        error("ActorGPEC: No equilibrium found in dd.equilibrium.time_slice")
    end

    ext = Base.get_extension(@__MODULE__, :GeneralizedPerturbedEquilibriumExt)
    ext === nothing && throw(MissingExtensionError("ActorGPEC", "GeneralizedPerturbedEquilibrium"))
    ext.gpec_step!(actor)

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
