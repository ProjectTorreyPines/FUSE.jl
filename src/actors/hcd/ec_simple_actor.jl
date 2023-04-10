#= == =#
#  EC  #
#= == =#
Base.@kwdef mutable struct FUSEparameters__ActorECsimple{T} <: ParametersActor where {T<:Real}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    width::Entry{Union{Real,AbstractVector{<:T}}} = Entry(Union{Real,AbstractVector{<:T}}, "-", "Width of the deposition profile"; default=0.05)
    rho_0::Entry{Union{Real,AbstractVector{<:T}}} = Entry(Union{Real,AbstractVector{<:T}}, "-", "Radial location of the deposition profile"; default=0.5)
end

mutable struct ActorECsimple <: HCDAbstractActor
    dd::IMAS.dd
    par::FUSEparameters__ActorECsimple

    function ActorECsimple(dd::IMAS.dd, par::FUSEparameters__ActorECsimple; kw...)
        logging_actor_init(ActorECsimple)
        par = par(kw...)
        return new(dd, par)
    end
end

"""
    ActorECsimple(dd::IMAS.dd, act::ParametersAllActors; kw...)

Estimates the EC electron energy deposition and current drive as a gaussian.

NOTE: Current drive efficiency from GASC, based on "G. Tonon 'Current Drive Efficiency Requirements for an Attractive Steady-State Reactor'"

!!! note
    Reads data in `dd.ec_launchers` and stores data in `dd.core_sources`
"""
function ActorECsimple(dd::IMAS.dd, act::ParametersAllActors; kw...)
    par = act.ActorECsimple
    actor = ActorECsimple(dd, par; kw...)
    step(actor)
    finalize(actor)
    return actor
end

function _step(actor::ActorECsimple)
    dd = actor.dd
    par = actor.par

    n_launchers = length(dd.ec_launchers.beam)
    _, width, rho_0 = same_length_vectors(1:n_launchers, par.width, par.rho_0)

    for (idx, ecl) in enumerate(dd.ec_launchers.beam)
        eqt = dd.equilibrium.time_slice[]
        cp1d = dd.core_profiles.profiles_1d[]
        cs = dd.core_sources

        power_launched = @ddtime(ecl.power_launched.data)

        rho_cp = cp1d.grid.rho_tor_norm
        volume_cp = IMAS.interp1d(eqt.profiles_1d.rho_tor_norm, eqt.profiles_1d.volume).(rho_cp)
        area_cp = IMAS.interp1d(eqt.profiles_1d.rho_tor_norm, eqt.profiles_1d.area).(rho_cp)

        ion_electron_fraction_cp = zeros(length(rho_cp))

        R0 = eqt.boundary.geometric_axis.r
        ne20 = IMAS.interp1d(rho_cp, cp1d.electrons.density).(rho_0[idx]) / 1E20
        TekeV = IMAS.interp1d(rho_cp, cp1d.electrons.temperature).(rho_0[idx]) / 1E3
        zeff = IMAS.interp1d(rho_cp, cp1d.zeff).(rho_0[idx]) / 1E3

        eta = TekeV * 0.09 / (5.0 + zeff)
        j_parallel = eta / R0 / ne20 * power_launched
        j_parallel *= sign(eqt.global_quantities.ip)

        source_index = IMAS.name_2_index(cs.source)[:ec]
        source = resize!(cs.source, "identifier.index" => source_index; allow_multiple_matches=true)
        gaussian_source(
            source,
            ecl.name,
            source_index,
            rho_cp,
            volume_cp,
            area_cp,
            power_launched,
            ion_electron_fraction_cp,
            rho_0[idx],
            width[idx],
            1.0;
            j_parallel=j_parallel
        )
    end
    return actor
end