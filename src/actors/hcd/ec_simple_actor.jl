#= ========= =#
#  Simple EC  #
#= ========= =#
Base.@kwdef mutable struct _FUSEparameters__ActorSimpleECactuator{T<:Real} <: ParametersActor{T}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    _time::Float64 = NaN
    ηcd_scale::Entry{T} = Entry{T}("-", "Scaling factor for nominal current drive efficiency"; default=1.0)
    rho_0::Entry{T} = Entry{T}("-", "Desired radial location of the deposition profile"; default=0.5, check=x -> @assert x >= 0.0 "must be: rho_0 >= 0.0")
    width::Entry{T} = Entry{T}("-", "Desired width of the deposition profile"; default=0.025, check=x -> @assert x >= 0.0 "must be: width > 0.0")
end

Base.@kwdef mutable struct FUSEparameters__ActorSimpleEC{T<:Real} <: ParametersActor{T}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    _time::Float64 = NaN
    actuator::ParametersVector{_FUSEparameters__ActorSimpleECactuator{T}} = ParametersVector{_FUSEparameters__ActorSimpleECactuator{T}}()
end

mutable struct ActorSimpleEC{D,P} <: SingleAbstractActor{D,P}
    dd::IMAS.dd{D}
    par::FUSEparameters__ActorSimpleEC{P}
    function ActorSimpleEC(dd::IMAS.dd{D}, par::FUSEparameters__ActorSimpleEC{P}; kw...) where {D<:Real,P<:Real}
        logging_actor_init(ActorSimpleEC)
        par = par(kw...)
        return new{D,P}(dd, par)
    end
end

"""
    ActorSimpleEC(dd::IMAS.dd, act::ParametersAllActors; kw...)

Estimates the EC electron energy deposition and current drive as a gaussian.

NOTE: Current drive efficiency from GASC, based on "G. Tonon 'Current Drive Efficiency Requirements for an Attractive Steady-State Reactor'"

!!! note

    Reads data in `dd.ec_launchers`, `dd.pulse_schedule` and stores data in `dd.core_sources`
"""
function ActorSimpleEC(dd::IMAS.dd, act::ParametersAllActors; kw...)
    actor = ActorSimpleEC(dd, act.ActorSimpleEC; kw...)
    step(actor)
    finalize(actor)
    return actor
end

function _step(actor::ActorSimpleEC)
    dd = actor.dd
    par = actor.par

    eqt = dd.equilibrium.time_slice[]
    cp1d = dd.core_profiles.profiles_1d[]
    cs = dd.core_sources

    R0 = eqt.boundary.geometric_axis.r
    rho_cp = cp1d.grid.rho_tor_norm
    volume_cp = IMAS.interp1d(eqt.profiles_1d.rho_tor_norm, eqt.profiles_1d.volume).(rho_cp)
    area_cp = IMAS.interp1d(eqt.profiles_1d.rho_tor_norm, eqt.profiles_1d.area).(rho_cp)

    # rho interpolant
    _, _, RHO_interpolant = IMAS.ρ_interpolant(eqt)

    for (k, (ps, ecl)) in enumerate(zip(dd.pulse_schedule.ec.beam, dd.ec_launchers.beam))
        τ_th = 0.01 # what's a good time here?
        power_launched = max(0.0, IMAS.smooth_beam_power(dd.pulse_schedule.ec.time, ps.power_launched.reference, dd.global_time, τ_th))
        rho_0 = par.actuator[k].rho_0
        width = par.actuator[k].width
        ηcd_scale = par.actuator[k].ηcd_scale

        coherent_wave = resize!(dd.waves.coherent_wave, k; wipe=false)[k]

        # vacuum "ray tracing"
        if !ismissing(ecl.frequency, :data)
            position_vector = [@ddtime(ecl.launching_position.r), 0, @ddtime(ecl.launching_position.z)]
            angle_pol = @ddtime(ecl.steering_angle_pol)
            angle_tor = @ddtime(ecl.steering_angle_tor)
            resonance_layer = IMAS.ech_resonance_layer(eqt, IMAS.frequency(ecl))
            t_intersect = IMAS.toroidal_intersection(resonance_layer.r, resonance_layer.z, position_vector..., angle_pol, angle_tor)
            if t_intersect == Inf
                t_intersect = 0.0
            end
            x, y, z, r = IMAS.pencil_beam(position_vector, angle_pol, angle_tor, range(0.0, t_intersect, 100))

            # save trajectory to dd
            beam_tracing = resize!(coherent_wave.beam_tracing)
            beam = resize!(beam_tracing.beam, 1)[1]
            beam.length = cumsum(sqrt.(IMAS.gradient(x) .^ 2 .+ IMAS.gradient(y) .^ 2 .+ IMAS.gradient(z) .^ 2))
            beam.position.r = r
            beam.position.z = z
            if t_intersect != 0.0
                rho_0 = RHO_interpolant.(r[end], z[end])
            else
                power_launched = 0.0
            end
        end

        @ddtime(ecl.power_launched.data = power_launched)

        ion_electron_fraction_cp = zeros(length(rho_cp))

        ne20 = IMAS.interp1d(rho_cp, cp1d.electrons.density).(rho_0) / 1E20
        TekeV = IMAS.interp1d(rho_cp, cp1d.electrons.temperature).(rho_0) / 1E3
        zeff = IMAS.interp1d(rho_cp, cp1d.zeff).(rho_0)

        eta = ηcd_scale * TekeV * 0.09 / (5.0 + zeff)
        j_parallel = eta / R0 / ne20 * power_launched
        j_parallel *= sign(eqt.global_quantities.ip)

        source = resize!(cs.source, :ec, "identifier.name" => ecl.name; wipe=false)
        shaped_source(
            source,
            ecl.name,
            source.identifier.index,
            rho_cp,
            volume_cp,
            area_cp,
            power_launched,
            ion_electron_fraction_cp,
            ρ -> IMAS.gaus(ρ, rho_0, width, 1.0);
            j_parallel
        )
        cs1d = source.profiles_1d[]

        wv1d = resize!(coherent_wave.profiles_1d)
        wv1d.grid.rho_tor_norm = cs1d.grid.rho_tor_norm
        wv1d.grid.psi = cs1d.grid.psi
        wv1d.power_density = cs1d.electrons.energy
        wv1d.electrons.power_density_thermal = cs1d.electrons.energy
        wv1d.current_parallel_density = cs1d.j_parallel
    end
    return actor
end