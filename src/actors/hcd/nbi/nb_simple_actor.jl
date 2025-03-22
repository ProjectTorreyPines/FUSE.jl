#= ========== =#
#  Simple NBI  #
#= ========== =#
Base.@kwdef mutable struct _FUSEparameters__ActorSimpleNBactuator{T<:Real} <: ParametersActor{T}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    _time::Float64 = NaN
    dl::Entry{T} = Entry{T}("-", "Step size for pencil beam [meters]"; default=0.005, check=x -> @assert x >= 0.0 "must be: rho_0 >= 0.0")
    banana_shift_fraction::Entry{T} = Entry{T}("-", "Shift factor"; default=0.5, check=x -> @assert x >= 0.0 "must be: rho_0 >= 0.0")
    smoothing_width::Entry{T} = Entry{T}("-", "Width of the deposition profile"; default=0.12, check=x -> @assert x >= 0.0 "must be: width > 0.0")
end

Base.@kwdef mutable struct FUSEparameters__ActorSimpleNB{T<:Real} <: ParametersActor{T}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    _time::Float64 = NaN
    actuator::ParametersVector{_FUSEparameters__ActorSimpleNBactuator{T}} = ParametersVector{_FUSEparameters__ActorSimpleNBactuator{T}}()
end

mutable struct ActorSimpleNB{D,P} <: SingleAbstractActor{D,P}
    dd::IMAS.dd{D}
    par::FUSEparameters__ActorSimpleNB{P}
    function ActorSimpleNB(dd::IMAS.dd{D}, par::FUSEparameters__ActorSimpleNB{P}; kw...) where {D<:Real,P<:Real}
        logging_actor_init(ActorSimpleNB)
        par = par(kw...)
        return new{D,P}(dd, par)
    end
end

"""
    ActorSimpleNB(dd::IMAS.dd, act::ParametersAllActors; kw...)

Estimates the NBI ion/electron energy deposition, particle source, rotation and current drive source with a super-gaussian.

NOTE: Current drive efficiency from GASC, based on "G. Tonon 'Current Drive Efficiency Requirements for an Attractive Steady-State Reactor'"

!!! note

    Reads data in `dd.nbi`, `dd.pulse_schedule` and stores data in `dd.core_sources`
"""
function ActorSimpleNB(dd::IMAS.dd, act::ParametersAllActors; kw...)
    actor = ActorSimpleNB(dd, act.ActorSimpleNB; kw...)
    step(actor)
    finalize(actor)
    return actor
end


function _step(actor::ActorSimpleNB)
    dd = actor.dd
    par = actor.par

    eqt = dd.equilibrium.time_slice[]
    cp1d = dd.core_profiles.profiles_1d[]
    cs = dd.core_sources
    eqt2d = findfirst(:rectangular, eqt.profiles_2d)

    rho_eq = eqt.profiles_1d.rho_tor
    rho_cp = cp1d.grid.rho_tor_norm
    volume_cp = IMAS.interp1d(eqt.profiles_1d.rho_tor_norm, eqt.profiles_1d.volume).(rho_cp)
    area_cp = IMAS.interp1d(eqt.profiles_1d.rho_tor_norm, eqt.profiles_1d.area).(rho_cp)

    ne = cp1d.electrons.density
    Te = cp1d.electrons.temperature

    phi = eqt2d.phi
    R0 = eqt.global_quantities.vacuum_toroidal_field.r0
    B0 = eqt.global_quantities.vacuum_toroidal_field.b0
    rho2d = sqrt.(abs.((phi) ./ pi ./ B0)) ./ rho_eq[end]

    r = eqt2d.grid.dim1
    z = eqt2d.grid.dim2
    r = range(r[1], r[end], length(r))
    z = range(z[1], z[end], length(z))
    rho2d_interp = Interpolations.cubic_spline_interpolation((r, z), rho2d; extrapolation_bc=2.0)

    nenergies = 3
    ncp1d = length(rho_cp)

    q_interp = IMAS.interp1d(eqt.profiles_1d.rho_tor_norm, eqt.profiles_1d.q)
    rin = eqt.profiles_1d.r_inboard
    rout = eqt.profiles_1d.r_outboard
    eps = (rout .- rin) ./ (rout .+ rin)
    eps_interp = IMAS.interp1d(eqt.profiles_1d.rho_tor_norm, eps)

    ne_interp = IMAS.interp1d(rho_cp, ne)
    Te_interp = IMAS.interp1d(rho_cp, Te)

    nmaxgroups = maximum(length(nbu.beamlets_group) for nbu in dd.nbi.unit)
    qbeam = zeros(nmaxgroups, nenergies, ncp1d)
    sbeam = zeros(nmaxgroups, nenergies, ncp1d)
    mombeam = zeros(nmaxgroups, nenergies, ncp1d)
    qbeame = zeros(nmaxgroups, nenergies, ncp1d)
    qbeami = zeros(nmaxgroups, nenergies, ncp1d)
    curbeam = zeros(nmaxgroups, nenergies, ncp1d)

    gaus = similar(rho_cp)
    qbeamtmp = similar(rho_cp)
    IMAS.freeze!(cp1d, :zeff)
    for (ibeam, (ps, nbu)) in enumerate(zip(dd.pulse_schedule.nbi.unit, dd.nbi.unit))
        # smoothing of the instantaneous power_launched based on the NBI thermalization time, effectively turning it into a measure of the absorbed power.
        beam_mass = nbu.species.a
        beam_Z = nbu.species.z_n
        beam_energy = max(0.0, @ddtime(ps.energy.reference))
        if beam_energy == 0.0
            continue
        end
        @ddtime(nbu.energy.data = beam_energy)

        fbcur = @ddtime(nbu.beam_current_fraction.data)

        qbeam .*= 0.0
        sbeam .*= 0.0
        mombeam .*= 0.0
        qbeame .*= 0.0
        qbeami .*= 0.0
        curbeam .*= 0.0

        ngroups = length(nbu.beamlets_group)
        for igroup in 1:length(nbu.beamlets_group)
            bgroup = nbu.beamlets_group[igroup]
            if ngroups > 1
                group_power_frac = bgroup.beamlets.power_fractions
            else
                group_power_frac = 1.0
            end
            source_r = bgroup.position.r
            source_z = bgroup.position.z

            angleh = bgroup.direction * asin(bgroup.tangency_radius / (source_r))
            anglev = bgroup.angle

            # trace pencil beam
            dx = -par.actuator[ibeam].dl .* cos(angleh) .* cos(anglev)
            dy = -par.actuator[ibeam].dl .* sin(angleh) * cos(anglev)
            dZ = par.actuator[ibeam].dl .* cos(angleh) * sin(anglev)
            x = source_r
            y = 0.0
            Z = source_z
            in_box = false
            out_box = false
            istep = 0
            ngrid = 0
            Rs = []
            Zs = []
            ftors = []
            max_steps = 1e6
            while ~in_box || ~out_box && istep < max_steps
                x += dx
                y += dy
                Z += dZ
                R = sqrt(x^2 + y^2)
                phi = asin(y / R)

                ftor = abs(-dx * sin(phi) + dy * cos(phi)) / par.actuator[ibeam].dl
                istep += 1
                if ~in_box
                    if R < rout[end]
                        in_box = true
                    end
                else
                    push!(Rs, R)
                    push!(Zs, Z)
                    push!(ftors, ftor)
                    ngrid += 1
                    if R > rout[end] || R < rin[end]
                        out_box = true
                    end
                end
            end

            rho_beam = rho2d_interp.(Rs, Zs)

            dist = range(0.0, par.actuator[ibeam].dl * ngrid, ngrid)
            ne_beam = ne_interp.(rho_beam)
            ne_beam[rho_beam.>1] .= 0.0
            Te_beam = Te_interp.(rho_beam)
            Te_beam[rho_beam.>1] .= 0.0

            power_launched_allenergies = 0.0
            for (ifpow, fpow) in enumerate(fbcur)
                τ_th = IMAS.fast_ion_thermalization_time(cp1d, nbu.species, beam_energy / ifpow)
                power_launched = fpow * max(0.0, IMAS.smooth_beam_power(dd.pulse_schedule.nbi.time, ps.power.reference, dd.global_time, τ_th))
                power_launched_allenergies += power_launched

                vbeam = sqrt((IMAS.mks.e * beam_energy / ifpow) / (0.5 * beam_mass * IMAS.mks.m_p))
                times = dist ./ vbeam

                cs = zeros(ngrid)
                for (i, (R, Z)) in enumerate(zip(Rs, Zs))
                    cs1 = IMAS.imfp_electron_collisions(vbeam * 1e2, Te_beam[i], ne_beam[i] * 1e-6)
                    cs2 = IMAS.imfp_ion_collisions(beam_mass, beam_energy / beam_mass, ne_beam[i] * 1e-6, 1)
                    cs3 = IMAS.imfp_charge_exchange(beam_mass, beam_energy / beam_mass, ne_beam[i] * 1e-6)
                    cs[i] = (cs1 + cs2 + cs3)
                end
                cross_section_t = IMAS.cumtrapz(times, cs .* vbeam)

                fbeam = exp.(-cross_section_t)
                rbananas = zeros(ngrid)

                mask = rho_beam .< 1
                eps = eps_interp.(rho_beam[mask])
                rbananas[mask] .= IMAS.banana_width.(beam_energy / ifpow, B0, beam_Z, beam_mass, eps, q_interp.(rho_beam[mask]))
                for i in 1:ngrid
                    rho_beam[i] = rho2d_interp(Rs[i] - rbananas[i] * (1.0 - ftors[i]) * par.actuator[ibeam].banana_shift_fraction * bgroup.direction, Zs[i])
                end

                for itime in 1:length(times)-1
                    if rho_beam[itime] < 1
                        @. gaus .= exp.(-0.5 .* (rho_cp .- rho_beam[itime]) .^ 2 ./ par.actuator[ibeam].smoothing_width^2) ./ (par.actuator[ibeam].smoothing_width * sqrt(2 * π))
                        gaus ./= IMAS.trapz(volume_cp, gaus)
                        @. qbeamtmp .= power_launched * group_power_frac * (fbeam[itime] - fbeam[itime+1]) .* gaus
                        @. qbeam[igroup, ifpow, :] .+= qbeamtmp
                        @. sbeam[igroup, ifpow, :] .+= qbeamtmp / (beam_energy * IMAS.mks.e / ifpow)
                        @. mombeam[igroup, ifpow, :] .+= bgroup.direction .* qbeamtmp .* (beam_mass * IMAS.mks.m_p .* vbeam) .* ftors[itime] / (beam_energy * IMAS.mks.e / ifpow)
                    end
                end
            end
            @ddtime(nbu.power_launched.data = power_launched_allenergies)

            mombeam *= 2.0 # fudge factor to get momentum flux right, why is this needed?????

            eps = maximum(eps) .* cp1d.grid.rho_tor_norm
            for ifpow in eachindex(fbcur)
                frac_ie = IMAS.sivukhin_fraction(cp1d, nbu.energy.data[1] / ifpow, nbu.species.a)
                tauppff = IMAS.ion_momentum_slowingdown_time(cp1d, nbu.energy.data[1] / ifpow, nbu.species.z_n, Int(nbu.species.a))
                qbeame[igroup, ifpow, :] .= @views (1.0 .- frac_ie) .* qbeam[igroup, ifpow, :]
                qbeami[igroup, ifpow, :] .= @views frac_ie .* qbeam[igroup, ifpow, :]
                curbi = IMAS.mks.e * mombeam[igroup, ifpow, :] .* tauppff / (nbu.species.a * IMAS.mks.m_p)
                curbe = -curbi ./ cp1d.zeff
                curbet = -curbe .* ((1.55 .+ 0.85 ./ cp1d.zeff) .* sqrt.(eps) .- (0.20 .+ 1.55 ./ cp1d.zeff) .* eps)
                curbeam[igroup, ifpow, :] .= curbe .+ curbi .+ curbet
            end
        end

        electrons_energy = sum(sum(qbeame; dims=1); dims=2)[1, 1, :]
        total_ion_energy = sum(sum(qbeami; dims=1); dims=2)[1, 1, :]
        momentum_tor = sum(sum(mombeam; dims=1); dims=2)[1, 1, :]
        curbeam_tot = sum(sum(curbeam; dims=1); dims=2)[1, 1, :]
        electrons_particles = sum(sum(sbeam; dims=1); dims=2)[1, 1, :]

        j_parallel = IMAS.JtoR_2_JparB(rho_cp, curbeam_tot ./ R0, false, eqt) ./ B0

        # Convert curbeam to parallel current here
        source = resize!(dd.core_sources.source, :nbi, "identifier.name" => nbu.name; wipe=false)
        IMAS.new_source(
            source,
            source.identifier.index,
            nbu.name,
            rho_cp,
            volume_cp,
            area_cp;
            electrons_energy,
            total_ion_energy,
            electrons_particles,
            j_parallel,
            momentum_tor)

        # add nbi fast ion particles source
        source1d = source.profiles_1d[]
        ion = resize!(source1d.ion, 1)[1]
        IMAS.ion_element!(ion, 1, nbu.species.a; fast=true)
        ion.particles = source1d.electrons.particles
        ion.particles_inside = source1d.electrons.particles_inside
        ion.energy = source1d.total_ion_energy
        ion.power_inside = source1d.total_ion_power_inside
        ion.fast_particles_energy = beam_energy
    end
    IMAS.empty!(cp1d, :zeff)

    return actor
end
