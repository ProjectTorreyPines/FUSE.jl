import RABBIT
using Plots

#= =========== =#
#  ActorRABBIT  #
#= =========== =#
@actor_parameters_struct ActorRABBIT{T} begin
    remove_inputs::Entry{Bool} = Entry{Bool}("-", "Delete directory containing RABBIT input files after run"; default=true)
    Δt_history::Entry{Float64} = Entry{Float64}("s", "Amount of history to include such that simulation proceeds from (dd.global_time - Δt_history) to dd.global_time")
end

mutable struct ActorRABBIT{D,P} <: SingleAbstractActor{D,P}
    dd::IMAS.dd{D}
    par::OverrideParameters{P,FUSEparameters__ActorRABBIT{P}}
    outputs::Union{RABBIT.RABBIToutput,Vector{<:RABBIT.RABBIToutput}}
end

function ActorRABBIT(dd::IMAS.dd, par::FUSEparameters__ActorRABBIT; kw...)
    logging_actor_init(ActorRABBIT)
    par = OverrideParameters(par; kw...)
    return ActorRABBIT(dd, par, RABBIT.RABBIToutput[])
end

"""
    ActorRABBIT(dd::IMAS.dd, act::ParametersAllActors; kw...)

Calculates neutral beam injection (NBI) heating, current drive, and fast ion physics
using the RABBIT Monte Carlo code. RABBIT provides detailed modeling of beam-plasma
interactions including collisional processes, beam thermalization, and current drive.

The actor interfaces with the external RABBIT code to perform:
- Detailed beam ionization and thermalization calculations
- Power deposition to electrons and ions
- Current drive from beam-driven currents
- Toroidal momentum input from NBI
- Fast ion particle source generation

The calculation includes temporal history to properly account for beam slowing-down
physics and requires equilibrium data over the specified time window.

!!! note

    Requires RABBIT external code. Reads data from `dd.nbi`, `dd.pulse_schedule` 
    and equilibrium data, stores results in `dd.core_sources`
"""
function ActorRABBIT(dd::IMAS.dd, act::ParametersAllActors; kw...)
    actor = ActorRABBIT(dd, act.ActorRABBIT; kw...)
    step(actor)
    finalize(actor)
    return actor
end

function _step(actor::ActorRABBIT)
    dd = actor.dd
    par = actor.par

    vessel_hfs = IMAS.get_build_layer(dd.build.layer; type=IMAS._vessel_, fs=IMAS._hfs_).start_radius
    vessel_lfs = IMAS.get_build_layer(dd.build.layer; type=IMAS._vessel_, fs=IMAS._lfs_).start_radius

    all_inputs = FUSEtoRABBITinput(dd, par.Δt_history)
    actor.outputs = RABBIT.run_RABBIT(all_inputs, vessel_hfs, vessel_lfs; par.remove_inputs)
    return actor
end

"""
    _finalize(actor::ActorRABBIT)

Processes RABBIT simulation outputs and populates IMAS core_sources with calculated
heating, current drive, and particle sources. Creates fast ion source terms for
each NBI unit with appropriate energy and particle fluxes.
"""
function _finalize(actor::ActorRABBIT)
    dd = actor.dd
    cs = dd.core_sources
    outputs = actor.outputs

    eqt = dd.equilibrium.time_slice[]

    rho = [0.0; outputs[1].rho_data; 1.0]
    volume = IMAS.interp1d(eqt.profiles_1d.rho_tor_norm, eqt.profiles_1d.volume).(rho)
    area = IMAS.interp1d(eqt.profiles_1d.rho_tor_norm, eqt.profiles_1d.area).(rho)

    for (idx, (ps, nbu)) in enumerate(zip(dd.pulse_schedule.nbi.unit, dd.nbi.unit))
        power_launched = max(0.0, @ddtime(ps.power.reference))
        beam_energy = max(0.0, @ddtime(nbu.energy.data))

        @ddtime(nbu.power_launched.data = power_launched)

        # evaluate various source channels
        # here we also extend the RABBIT grid to the edges
        # so that the linear interpolation does not extrapolate in weird ways
        electrons_energy = outputs[idx].powe_data[:, end]
        electrons_energy = [electrons_energy[1]; electrons_energy; electrons_energy[end]]
        total_ion_energy = outputs[idx].powi_data[:, end]
        total_ion_energy = [total_ion_energy[1]; total_ion_energy; total_ion_energy[end]]
        electrons_particles = vec(sum(outputs[idx].bdep_data[:, end, :]; dims=2))
        electrons_particles = [electrons_particles[1]; electrons_particles; electrons_particles[end]]
        j_parallel = outputs[idx].jnbcd_data[:, end]
        j_parallel = [j_parallel[1]; j_parallel; j_parallel[end]]
        momentum_tor = vec(sum(outputs[idx].torqdepo_data[:, end, :]; dims=2))
        momentum_tor = [momentum_tor[1]; momentum_tor; momentum_tor[end]]

        source = resize!(cs.source, :nbi, "identifier.name" => nbu.name; wipe=false)
        IMAS.new_source(
            source,
            source.identifier.index,
            nbu.name,
            rho,
            volume,
            area;
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

    return actor
end

function FUSEtoRABBITinput(dd::IMAS.dd, Δt_history::Float64)
    eV_to_keV = 1e-3
    cm3_to_m3 = 1e-6

    all_inputs = RABBIT.RABBITinput[]

    index_start = IMAS.nearest_causal_time(dd.equilibrium.time, (dd.global_time - Δt_history); bounds_error=false).index
    index_end = IMAS.nearest_causal_time(dd.equilibrium.time, dd.global_time).index
    eqts = [dd.equilibrium.time_slice[index] for index in index_start:index_end]
    selected_equilibrium_times = [eqt.time for eqt in eqts]

    for eqt in eqts
        time = eqt.time

        eqt2d = findfirst(:rectangular, eqt.profiles_2d)
        if eqt2d === nothing
            continue
        end

        inp = RABBIT.RABBITinput()
        inp.time = time

        inp.nw = length(eqt2d.grid.dim1)
        inp.nh = length(eqt2d.grid.dim2)

        inp.psirz = eqt2d.psi ./ 2pi
        inp.npsi1d = length(eqt.profiles_1d.psi)
        inp.qpsi = abs.(eqt.profiles_1d.q) # RABBIT assumes positive q 
        inp.fpol = eqt.profiles_1d.f
        inp.sibry = eqt.global_quantities.psi_boundary / 2pi
        inp.simag = eqt.global_quantities.psi_axis / 2pi
        inp.signip = sign(eqt.global_quantities.ip)
        inp.rmaxis = eqt.global_quantities.magnetic_axis.r
        inp.zmaxis = eqt.global_quantities.magnetic_axis.z

        inp.r = eqt2d.grid.dim1
        inp.z = eqt2d.grid.dim2

        phi = eqt.profiles_1d.phi
        rhorz_n = eqt2d.phi ./ last(phi)
        rhorz_sq = map(x -> x < 0 ? 0 : x, rhorz_n)
        rhorz = sqrt.(rhorz_sq)

        if inp.simag == inp.sibry
            psirz_norm = abs.(inp.psirz .- inp.simag)
        else
            psirz_norm = abs.(inp.psirz .- inp.simag) ./ (inp.sibry - inp.simag)
        end

        rhoprz = sqrt.(psirz_norm)

        # as is done in OMFITrabbitEq class, use rho_tor inside lcfs, rho_pol outside (omfit_rabbit.py line 193)
        rhorz[findall(rhorz .> 1)] .= rhoprz[findall(rhorz .> 1)]
        inp.rhorz = rhorz

        inp.psi = eqt.profiles_1d.psi ./ 2pi
        inp.vol = eqt.profiles_1d.volume
        inp.area = eqt.profiles_1d.area

        cp1d = try
            dd.core_profiles.profiles_1d[time]
        catch
            dd.core_profiles.profiles_1d[1]
        end

        inp.rho = cp1d.grid.rho_tor_norm
        inp.eq_rho = eqt.profiles_1d.rho_tor_norm
        inp.n_rho = length(inp.rho)

        inp.te = IMAS.interp1d(cp1d.grid.rho_tor_norm, cp1d.electrons.temperature).(inp.rho) .* eV_to_keV
        inp.dene = IMAS.interp1d(cp1d.grid.rho_tor_norm, cp1d.electrons.density).(inp.rho) .* cm3_to_m3
        inp.rot_freq_tor = inp.rho .* 0.0
        inp.zeff = IMAS.interp1d(cp1d.grid.rho_tor_norm, cp1d.zeff).(inp.rho)
        inp.ti = IMAS.interp1d(cp1d.grid.rho_tor_norm, cp1d.t_i_average).(inp.rho) .* eV_to_keV

        push!(all_inputs, inp)
    end

    function gather_beams(dd::IMAS.dd)
        nbeams = length(dd.nbi.unit)
        nv = 3

        xyz_src = zeros(3, nbeams)
        xtan = zeros(3, nbeams)
        xyz_vec = zeros(3, nbeams)
        beamwidthpoly = zeros(3, nbeams)
        part_frac = zeros(3, nbeams)

        Einj = Vector{Float64}(undef, nbeams)
        abeam = Vector{Float64}(undef, nbeams)

        for n in 1:nbeams
            unit = dd.nbi.unit[n]
            bgrp = unit.beamlets_group[1]
            pos = bgrp.position
            dir = bgrp.direction
            angle = bgrp.angle
            Rt = bgrp.tangency_radius
            R = pos.r
            z = pos.z
            phi = 2π - pos.phi

            x = R * cos(phi)
            y = R * sin(phi)

            xyz_src[:, n] .= (x, y, z)

            l2d = sqrt(R^2 - Rt^2)
            delta = atan(l2d, Rt)
            phit = phi + delta * dir
            zt = z + tan(angle) * l2d

            xtan[:, n] .= (Rt * cos(phit), Rt * sin(phit), zt)

            for i in 1:3
                xyz_vec[i, n] = xtan[i, n] - xyz_src[i, n]
            end
            norm = sqrt(sum(xyz_vec[:, n] .^ 2))
            xyz_vec[:, n] ./= norm

            Einj[n] = maximum(unit.energy.data)
            abeam[n] = unit.species.a

            for i in 1:3
                part_frac[i, n] = maximum(unit.beam_current_fraction.data[i, :])
            end
            s = sum(part_frac[:, n])
            part_frac[:, n] ./= s

            beamwidthpoly[2, n] = bgrp.divergence_component[1].vertical
        end

        return nbeams, nv, xyz_src, xyz_vec, beamwidthpoly, Einj, part_frac, abeam
    end

    function get_pnbi(dd::IMAS.dd, selected_equilibrium_times::Vector{Float64})
        pnbis = Vector{Float64}[]
        if length(selected_equilibrium_times) == 1
            for ps in dd.pulse_schedule.nbi.unit
                power = @ddtime ps.power.reference
                push!(pnbis, [power])
            end
        else
            for ps in dd.pulse_schedule.nbi.unit
                power_downsampled = IMAS.moving_average(dd.pulse_schedule.nbi.time, ps.power.reference, selected_equilibrium_times)
                push!(pnbis, power_downsampled)
            end
        end
        return pnbis
    end

    pnbis = get_pnbi(dd, selected_equilibrium_times)
    all_inputs[1].pnbi = pnbis

    # beam info isn't time dependent so store it in the first timeslice
    all_inputs[1].n_sources,
    all_inputs[1].nv,
    all_inputs[1].start_pos,
    all_inputs[1].beam_unit_vector,
    all_inputs[1].beam_width_polynomial_coefficients,
    all_inputs[1].injection_energy,
    all_inputs[1].particle_fraction,
    all_inputs[1].a_beam = gather_beams(dd)

    if length(all_inputs) == 1
        inp = deepcopy(all_inputs[1])
        inp.time = -1e6
        push!(all_inputs, inp)
        reverse!(all_inputs)
        i = 1
        while i <= length(all_inputs[1].pnbi) # pnbi is written per beam, per timeslice e.g. beam1@time1, beam1@time2, beam2@time1, beam2@time2, etc.
            copy_elem = deepcopy(all_inputs[1].pnbi[i])
            insert!(all_inputs[1].pnbi, i + 1, copy_elem)
            i += 2
        end
    end

    return all_inputs

end