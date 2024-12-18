#= ======================= =#
#  ActorSteadyStateCurrent  #
#= ======================= =#
Base.@kwdef mutable struct FUSEparameters__ActorSteadyStateCurrent{T<:Real} <: ParametersActor{T}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    _time::Float64 = NaN
    allow_floating_plasma_current::Entry{Bool} = Entry{Bool}("-", "Zero loop voltage if non-inductive fraction exceeds 100% of the target Ip")
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
    par::FUSEparameters__ActorSteadyStateCurrent{P}
    function ActorSteadyStateCurrent(dd::IMAS.dd{D}, par::FUSEparameters__ActorSteadyStateCurrent{P}; kw...) where {D<:Real,P<:Real}
        logging_actor_init(ActorSteadyStateCurrent)
        par = par(kw...)
        return new{D,P}(dd, par)
    end
end

"""
    ActorSteadyStateCurrent(dd::IMAS.dd, act::ParametersAllActors; kw...)

Evolves the ohmic current to steady state using the conductivity from `dd.core_profiles`

!!! note

    Stores data in `dd.core_profiles.profiles_1d[].j_ohmic`
"""
function ActorSteadyStateCurrent(dd::IMAS.dd, act::ParametersAllActors; kw...)
    actor = ActorSteadyStateCurrent(dd, act.ActorSteadyStateCurrent; kw...)
    step(actor)
    finalize(actor)
    return actor
end

function _step(actor::ActorSteadyStateCurrent)
    dd = actor.dd
    par = actor.par

    eqt = dd.equilibrium.time_slice[]
    cp1d = dd.core_profiles.profiles_1d[]

    ip_target = IMAS.get_from(dd, Val{:ip}, par.ip_from)

    # update j_ohmic
    relaxed_j_ohmic = IMAS.j_ohmic_steady_state(eqt, cp1d, ip_target, cp1d.conductivity_parallel)
    if par.current_relaxation_radius == 0.0
        cp1d.j_ohmic = relaxed_j_ohmic
    else
        # blend between an initial (parabolic) ohmic current profile
        # and the fully relaxed  ohmic current profile based on the
        # current diffusion time evaluated at the 
        rho_tor_norm = cp1d.grid.rho_tor_norm
        initial_j_ohmic = IMAS.j_ohmic_steady_state(eqt, cp1d, ip_target, 1.0 ./ abs.(1.1 .- rho_tor_norm))

        j_diffusion_time = IMAS.mks.μ_0 .* eqt.boundary.minor_radius .^ 2.0 .* cp1d.conductivity_parallel

        time = IMAS.interp1d(rho_tor_norm, j_diffusion_time).(par.current_relaxation_radius)
        alpha = 1.0 .- exp.(-time ./ j_diffusion_time)

        interp_j = relaxed_j_ohmic .* alpha .+ initial_j_ohmic .* (1.0 .- alpha)
        interpo_j_ohmic = IMAS.j_ohmic_steady_state(eqt, cp1d, ip_target, interp_j)

        cp1d.j_ohmic = interpo_j_ohmic
    end

    # allow floating plasma current
    ip_non_inductive = IMAS.Ip_non_inductive(cp1d, eqt)
    if abs(ip_target) < abs(ip_non_inductive) && par.allow_floating_plasma_current
        cp1d.j_ohmic = zeros(length(cp1d.grid.rho_tor_norm))
    end

    return actor
end
