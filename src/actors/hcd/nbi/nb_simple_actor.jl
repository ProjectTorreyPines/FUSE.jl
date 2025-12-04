#= ========== =#
#  Simple NBI  #
#= ========== =#
@actor_parameters_struct _ActorSimpleNBactuator{T} begin
    banana_shift_fraction::Entry{T} = Entry{T}("-", "Shift factor"; default=0.5, check=x -> @assert x >= 0.0 "must be: banana_shift_fraction >= 0.0")
    smoothing_width::Entry{T} = Entry{T}("-", "Width of the deposition profile"; default=0.12, check=x -> @assert x >= 0.0 "must be: smoothing_width > 0.0")
end

@actor_parameters_struct ActorSimpleNB{T} begin
    actuator::ParametersVector{_FUSEparameters__ActorSimpleNBactuator{T}} = ParametersVector{_FUSEparameters__ActorSimpleNBactuator{T}}()
end

mutable struct ActorSimpleNB{D,P} <: SingleAbstractActor{D,P}
    dd::IMAS.dd{D}
    par::OverrideParameters{P,FUSEparameters__ActorSimpleNB{P}}
    function ActorSimpleNB(dd::IMAS.dd{D}, par::FUSEparameters__ActorSimpleNB{P}; kw...) where {D<:Real,P<:Real}
        logging_actor_init(ActorSimpleNB)
        par = OverrideParameters(par; kw...)
        return new{D,P}(dd, par)
    end
end

"""
    ActorSimpleNB(dd::IMAS.dd, act::ParametersAllActors; kw...)

Calculates neutral beam injection (NBI) heating and current drive using a simplified pencil beam model.

This actor models the deposition of neutral beam energy and momentum through:
- Pencil beam ray tracing through the plasma
- Calculation of beam attenuation via electron collisions, ion collisions, and charge exchange
- Energy deposition to electrons and ions using Sivukhin fractions
- Toroidal momentum deposition including banana orbit effects
- Parallel current drive from beam-driven currents (beam, Ohkawa, and bootstrap contributions)

The model includes:
- Multiple beam energy components (full, half, third energy)
- Gaussian deposition profiles with configurable width
- Banana orbit shift correction for trapped particles
- Fast ion thermalization time smoothing of launched power

!!! note

    Reads data in `dd.nbi`, `dd.pulse_schedule` and stores results in `dd.core_sources` and `dd.waves.coherent_wave`
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
    eps_cp = eps_interp.(rho_cp)

    ne_interp = IMAS.interp1d(rho_cp, ne)
    Te_interp = IMAS.interp1d(rho_cp, Te)

    ngrid = 25

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

        # save ray trajectory to dd
        coherent_wave = resize!(dd.waves.coherent_wave, "identifier.antenna_name" => nbu.name; wipe=false)
        beam_tracing = resize!(coherent_wave.beam_tracing)

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

            angleh = bgroup.direction * asin(bgroup.tangency_radius / source_r)
            anglev = bgroup.angle

            # trace pencil beam

            # initial position and velocity
            px0 = source_r
            py0 = 0.0
            pz = source_z
            vx0 = -cos(angleh) .* cos(anglev)
            vy0 = -sin(angleh) * cos(anglev)
            vz = cos(angleh) * sin(anglev)

            # rotate position and velocity in phi
            px = px0 .* cos(bgroup.position.phi) .- py0 .* sin(bgroup.position.phi)
            py = px0 .* sin(bgroup.position.phi) .+ py0 .* cos(bgroup.position.phi)
            vx = vx0 .* cos(bgroup.position.phi) .- vy0 .* sin(bgroup.position.phi)
            vy = vx0 .* sin(bgroup.position.phi) .+ vy0 .* cos(bgroup.position.phi)

            t_intersects = IMAS.toroidal_intersections(eqt.boundary.outline.r, eqt.boundary.outline.z, px, py, pz, vx, vy, vz; max_intersections=2)
            if rin[end] < source_r < rout[end] && length(t_intersects) < 1
                continue
            elseif length(t_intersects) < 2
                continue
            end
            if rin[end] < source_r < rout[end]
                tt = range(0.0, t_intersects[1], ngrid)
            else
                tt = range(t_intersects[1], t_intersects[2], ngrid)
            end
            Xs, Ys, Zs, Rs = IMAS.pencil_beam([px, py, pz], [vx, vy, vz], tt)
            phi = asin.(Ys ./ Rs)
            ftors = abs.(-vx .* sin.(phi) .+ vy .* cos.(phi))
            rho_beam = rho2d_interp.(Rs, Zs)
            dist = IMAS.arc_length(Xs, Ys, Zs)

            # save ray trajectory to dd
            coherent_wave = resize!(dd.waves.coherent_wave, "identifier.antenna_name" => nbu.name; wipe=false)
            beam_tracing = resize!(coherent_wave.beam_tracing)
            beam = resize!(beam_tracing.beam, igroup)[igroup]
            pos_r = [sqrt(px^2 + py^2); sqrt.(Xs.^2 .+ Ys.^2)]
            pos_phi = [atan(py, px); atan.(Ys, Xs)]
            pos_z = [pz; Zs]
            beam.length = IMAS.arc_length_cylindrical(pos_r, pos_phi, pos_z)
            beam.position.r = pos_r
            beam.position.z = pos_z
            beam.position.phi = pos_phi

            ne_beam = ne_interp.(rho_beam)
            Te_beam = Te_interp.(rho_beam)

            power_launched_allenergies = 0.0
            for (ifpow, fpow) in enumerate(fbcur)
                # smoothing of the instantaneous power_launched based on the NBI thermalization time, effectively turning it into a measure of the absorbed power
                τ_th = IMAS.fast_ion_thermalization_time(cp1d, 1, nbu.species, beam_energy / ifpow)
                power_launched = fpow * max(0.0, IMAS.smooth_beam_power(dd.pulse_schedule.nbi.time, ps.power.reference, dd.global_time, τ_th))
                power_launched_allenergies += power_launched
                @ddtime(nbu.power_launched.data = power_launched_allenergies)
                beam.power_initial = power_launched_allenergies
                if power_launched == 0.0
                    continue
                end

                vbeam = sqrt((IMAS.mks.e * beam_energy / ifpow) / (0.5 * beam_mass * IMAS.mks.m_p))

                cs = zeros(ngrid)
                for i in 1:ngrid
                    cs1 = IMAS.imfp_electron_collisions(vbeam, Te_beam[i], ne_beam[i])
                    cs2 = IMAS.imfp_ion_collisions(beam_mass, beam_energy / beam_mass, ne_beam[i], 1)
                    cs3 = IMAS.imfp_charge_exchange(beam_mass, beam_energy / beam_mass, ne_beam[i])
                    cs[i] = cs1 + cs2 + cs3
                end
                cross_section_t = IMAS.cumtrapz(dist, cs)
                fbeam = exp.(-cross_section_t)

                rbananas = IMAS.banana_width.(beam_energy / ifpow, B0, beam_Z, beam_mass, eps_interp.(rho_beam), q_interp.(rho_beam))

                for i in 1:ngrid-1
                    rho_beam_banana = rho2d_interp(Rs[i] - rbananas[i] * (1.0 - ftors[i]) * par.actuator[ibeam].banana_shift_fraction * bgroup.direction, Zs[i])
                    @. gaus .= exp.(-0.5 .* (rho_cp .- rho_beam_banana) .^ 2 ./ par.actuator[ibeam].smoothing_width^2) ./ (par.actuator[ibeam].smoothing_width * sqrt(2 * π))
                    gaus ./= IMAS.trapz(volume_cp, gaus)
                    @. qbeamtmp .= power_launched * group_power_frac * (fbeam[i] - fbeam[i+1]) .* gaus
                    @. qbeam[igroup, ifpow, :] .+= qbeamtmp
                    @. sbeam[igroup, ifpow, :] .+= qbeamtmp / (beam_energy * IMAS.mks.e / ifpow)
                    @. mombeam[igroup, ifpow, :] .+= bgroup.direction .* qbeamtmp .* (beam_mass * IMAS.mks.m_p .* vbeam) .* Rs[i] * ftors[i] / (beam_energy * IMAS.mks.e / ifpow)
                end
            end
            for ifpow in eachindex(fbcur)
                frac_ie = IMAS.sivukhin_fraction(cp1d, beam_energy / ifpow, nbu.species.a)
                tauppff = IMAS.ion_momentum_slowingdown_time(cp1d, beam_energy / ifpow, nbu.species.a, nbu.species.z_n)
                qbeame[igroup, ifpow, :] .= @views (1.0 .- frac_ie) .* qbeam[igroup, ifpow, :]
                qbeami[igroup, ifpow, :] .= @views frac_ie .* qbeam[igroup, ifpow, :]
                # There seems to be a factor of 0.1 pull from freya, Maybe going from momentum of g*cm to kg*m? 
                # from freya: charge/(2.99792458e9*atwb*xmassp)[A/cm^2] = 47894.15 * 1e4 =  0.1*IMAS.mks.e/(nbu.species.a * IMAS.mks.m_p)[A/m^2]
                curbi = @views 0.1 * IMAS.mks.e * mombeam[igroup, ifpow, :] .* tauppff / (nbu.species.a * IMAS.mks.m_p)
                curbe = -curbi ./ cp1d.zeff
                curbet = -curbe .* ((1.55 .+ 0.85 ./ cp1d.zeff) .* sqrt.(eps_cp) .- (0.20 .+ 1.55 ./ cp1d.zeff) .* eps_cp)
                curbeam[igroup, ifpow, :] .= curbi .+ curbe .+ curbet
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
        resize!(source1d.ion, 2)
        ion_fast = source1d.ion[1]
        IMAS.ion_element!(ion_fast, 1, nbu.species.a; fast=true)
        ion_fast.particles = source1d.electrons.particles
        ion_fast.particles_inside = source1d.electrons.particles_inside
        ion_fast.energy = source1d.total_ion_energy
        ion_fast.power_inside = source1d.total_ion_power_inside
        ion_fast.fast_particles_energy = beam_energy

        ion_thermal = source1d.ion[2]
        IMAS.ion_element!(ion_thermal, 1, nbu.species.a; fast=false)
        ion_thermal.particles        = source1d.electrons.particles
        ion_thermal.particles_inside = source1d.electrons.particles_inside

    end
    IMAS.unfreeze!(cp1d, :zeff)

    return actor
end
