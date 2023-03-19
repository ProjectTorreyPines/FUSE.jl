#= == =#
#  LH  #
#= == =#
Base.@kwdef mutable struct FUSEparameters__ActorLHsimple{T} <: ParametersActor where {T<:Real}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    width::Entry{Union{Real,AbstractVector{<:T}}} = Entry(Union{Real,AbstractVector{<:T}}, "-", "Width of the deposition profile"; default=0.15)
    rho_0::Entry{Union{Real,AbstractVector{<:T}}} = Entry(Union{Real,AbstractVector{<:T}}, "-", "Radial location of the deposition profile"; default=0.6)
    current_efficiency::Entry{Union{Real,AbstractVector{<:T}}} = Entry(Union{Real,AbstractVector{<:T}}, "A/W", "Current drive efficiency"; default=0.4)
end

mutable struct ActorLHsimple <: HCDAbstractActor
    dd::IMAS.dd
    par::FUSEparameters__ActorLHsimple
    width::AbstractVector{<:Real}
    rho_0::AbstractVector{<:Real}
    current_efficiency::AbstractVector{<:Real}
end

"""
    ActorLHsimple(dd::IMAS.dd, act::ParametersAllActors; kw...)

Estimates the Lower-hybrid electron energy deposition and current drive as a gaussian.

!!! note 
    Stores data in `dd.lh_antennas, dd.core_sources`
"""
function ActorLHsimple(dd::IMAS.dd, act::ParametersAllActors; kw...)
    par = act.ActorLHsimple(kw...)
    actor = ActorLHsimple(dd, par)
    step(actor)
    finalize(actor)
    return actor
end

function ActorLHsimple(dd::IMAS.dd, par::FUSEparameters__ActorLHsimple; kw...)
    logging_actor_init(ActorLHsimple)
    par = par(kw...)
    n_antennas = length(dd.lh_antennas.antenna)
    _, width, rho_0, current_efficiency = same_length_vectors(1:n_antennas, par.width, par.rho_0, par.current_efficiency)
    return ActorLHsimple(dd, par, width, rho_0, current_efficiency)
end

function _step(actor::ActorLHsimple)
    for (idx, lha) in enumerate(actor.dd.lh_antennas.antenna)
        eqt = actor.dd.equilibrium.time_slice[]
        cp1d = actor.dd.core_profiles.profiles_1d[]
        cs = actor.dd.core_sources

        power_launched = @ddtime(lha.power_launched.data)

        rho_cp = cp1d.grid.rho_tor_norm
        volume_cp = IMAS.interp1d(eqt.profiles_1d.rho_tor_norm, eqt.profiles_1d.volume).(rho_cp)
        area_cp = IMAS.interp1d(eqt.profiles_1d.rho_tor_norm, eqt.profiles_1d.area).(rho_cp)

        ion_electron_fraction_cp = zeros(length(rho_cp))

        ne_vol = integrate(volume_cp, cp1d.electrons.density) / volume_cp[end]
        j_parallel = actor.current_efficiency[idx] / eqt.boundary.geometric_axis.r / (ne_vol / 1e19) * power_launched
        j_parallel *= sign(eqt.global_quantities.ip)

        source_index = IMAS.name_2_index(cs.source)[:lh]
        source = resize!(cs.source, "identifier.index" => source_index; allow_multiple_matches=true)
        gaussian_source(
            source,
            lha.name,
            source_index,
            rho_cp,
            volume_cp,
            area_cp,
            power_launched,
            ion_electron_fraction_cp,
            actor.rho_0[idx],
            actor.width[idx],
            1.0;
            j_parallel=j_parallel
        )
    end
    return actor
end