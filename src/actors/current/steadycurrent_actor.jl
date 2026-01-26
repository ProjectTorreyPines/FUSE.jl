#= ======================= =#
#  ActorSteadyStateCurrent  #
#= ======================= =#
@actor_parameters_struct ActorSteadyStateCurrent{T} begin
    allow_floating_plasma_current::Entry{Bool} = Entry{Bool}("-", "Zero loop voltage if non-inductive fraction exceeds 100% of the target Ip"; default=true)
    current_relaxation_radius::Entry{Float64} = Entry{Float64}(
        "-",
        "Radial position at which the artificial ohmic current profile relaxation starts to kick in";
        default=0.0,
        check=x -> @assert 1 >= x >= 0 "current_relaxation_radius must be between 0.0 and 1.0"
    )
    #== data flow parameters ==#
    ip_from::Switch{Symbol} = switch_get_from(:ip)
end

mutable struct ActorSteadyStateCurrent{D,P} <: SingleAbstractActor{D,P}
    dd::IMAS.dd{D}
    par::OverrideParameters{P,FUSEparameters__ActorSteadyStateCurrent{P}}
    function ActorSteadyStateCurrent(dd::IMAS.dd{D}, par::FUSEparameters__ActorSteadyStateCurrent{P}; kw...) where {D<:Real,P<:Real}
        logging_actor_init(ActorSteadyStateCurrent)
        par = OverrideParameters(par; kw...)
        return new{D,P}(dd, par)
    end
end

"""
    ActorSteadyStateCurrent(dd::IMAS.dd, act::ParametersAllActors; kw...)

Computes the steady-state ohmic current distribution using plasma conductivity and target plasma current.

This actor solves for the equilibrium current distribution assuming infinite current diffusion time,
where the current profile is determined by the balance between ohmic heating and resistive dissipation.
The solution uses the conductivity profile from `dd.core_profiles` and can either target a specific
plasma current or allow floating current based on non-inductive drive.

Key features:
- Computes steady-state ohmic current profile using plasma conductivity
- Supports current relaxation with adjustable radial extent
- Can float plasma current when non-inductive fraction exceeds 100%
- Updates `j_total` by combining ohmic and non-inductive current components

!!! note

    The fundamental quantitiy being solved is `j_total` in `dd.core_profiles.profiles_1d[]`
"""
function ActorSteadyStateCurrent(dd::IMAS.dd, act::ParametersAllActors; kw...)
    actor = ActorSteadyStateCurrent(dd, act.ActorSteadyStateCurrent; kw...)
    step(actor)
    finalize(actor)
    return actor
end

"""
    _step(actor::ActorSteadyStateCurrent)

Computes the steady-state ohmic current profile.

The calculation process:
1. Obtains target plasma current from specified source
2. Computes fully relaxed ohmic current profile using plasma conductivity
3. Optionally applies current relaxation blending with artificial profile
4. Checks for floating plasma current condition (non-inductive fraction > 100%)
5. Updates `j_total` by combining ohmic and non-inductive components

Current relaxation (if enabled) blends between an artificial current profile and the
physically motivated steady-state solution based on local current diffusion times.
"""
function _step(actor::ActorSteadyStateCurrent)
    dd = actor.dd
    par = actor.par

    eqt = dd.equilibrium.time_slice[]
    cp1d = dd.core_profiles.profiles_1d[]

    ip_target = IMAS.get_from(dd, Val(:ip), par.ip_from)

    # update j_ohmic
    relaxed_j_ohmic = IMAS.j_ohmic_steady_state(eqt, cp1d, ip_target, cp1d.conductivity_parallel)
    if par.current_relaxation_radius == 0.0
        j_ohmic = relaxed_j_ohmic
    else
        # blend between an initial ohmic current profile and the fully relaxed  ohmic current profile
        # the blending is proportional to the local current diffusion time and the current diffusion time
        # that is characteristic of the `current_relaxation_radius` parameter
        rho_tor_norm = cp1d.grid.rho_tor_norm
        initial_j_ohmic = IMAS.j_ohmic_steady_state(eqt, cp1d, ip_target, 1.0 ./ abs.(1.1 .- rho_tor_norm))

        j_diffusion_time = IMAS.mks.μ_0 .* eqt.boundary.minor_radius .^ 2.0 .* cp1d.conductivity_parallel

        time = IMAS.interp1d(rho_tor_norm, j_diffusion_time).(par.current_relaxation_radius) .+ 0.1 ./ par.current_relaxation_radius
        alpha = 1.0 .- exp.(-time ./ j_diffusion_time)

        interp_j = relaxed_j_ohmic .* alpha .+ initial_j_ohmic .* (1.0 .- alpha)
        j_ohmic = IMAS.j_ohmic_steady_state(eqt, cp1d, ip_target, interp_j)
    end

    # allow floating plasma current
    ip_non_inductive = IMAS.Ip_non_inductive(cp1d, eqt)
    if abs(ip_target) < abs(ip_non_inductive) && par.allow_floating_plasma_current
        j_ohmic = zeros(length(cp1d.grid.rho_tor_norm))
    end

    if ismissing(cp1d, :j_non_inductive)
        cp1d.j_total = j_ohmic
    else
        cp1d.j_total = j_ohmic .+ cp1d.j_non_inductive
    end

    IMAS.unfreeze!(cp1d, :q)

    return actor
end
