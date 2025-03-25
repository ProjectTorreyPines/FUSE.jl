import RABBIT
using Plots

#= =========== =#
#  ActorRABBIT  #
#= =========== =#
Base.@kwdef mutable struct FUSEparameters__ActorRABBIT{T<:Real} <: ParametersActor{T}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    _time::Float64 = NaN
    remove_inputs::Entry{Bool} = Entry{Bool}("-", "Delete directory containing RABBIT input files after run"; default=true)
end

mutable struct ActorRABBIT{D,P} <: SingleAbstractActor{D,P}
    dd::IMAS.dd{D}
    par::FUSEparameters__ActorRABBIT{P}
    outputs::Union{RABBIT.RABBIToutput,Vector{<:RABBIT.RABBIToutput}}
end

function ActorRABBIT(dd::IMAS.dd, par::FUSEparameters__ActorRABBIT; kw...)
    logging_actor_init(ActorRABBIT)
    par = par(kw...)
    return ActorRABBIT(dd, par, RABBIT.RABBIToutput[])
end

"""
    ActorRABBIT(dd::IMAS.dd, act::ParametersAllActors; kw...)
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

    all_inputs = FUSEtoRABBITinput(dd)
    actor.outputs = RABBIT.run_RABBIT(all_inputs; par.remove_inputs)
    return actor
end

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

function FUSEtoRABBITinput(dd::IMAS.dd)
    eV_to_keV = 1e-3
    cm3_to_m3 = 1e-6

    eq = dd.equilibrium

    all_inputs = RABBIT.RABBITinput[]

    for eqt in eq.time_slice
        time = eqt.time

        eqt2d = findfirst(:rectangular, eqt.profiles_2d)
        if eqt2d === nothing
            continue
        end

        inp = RABBIT.RABBITinput()
        inp.time = time * 1e3

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
    
        Rs = [unit.beamlets_group[1].position.r for unit in dd.nbi.unit]
        zs = [unit.beamlets_group[1].position.z for unit in dd.nbi.unit]
        phis = [unit.beamlets_group[1].position.phi for unit in dd.nbi.unit]
        Rts = [unit.beamlets_group[1].tangency_radius for unit in dd.nbi.unit]
        angles = [unit.beamlets_group[1].angle for unit in dd.nbi.unit]
        direcs = [unit.beamlets_group[1].direction for unit in dd.nbi.unit]
    
        phis .= 2 * π .- phis
    
        xyz_src[1, :] .= Rs .* cos.(phis)
        xyz_src[2, :] .= Rs .* sin.(phis)
        xyz_src[3, :] .= zs
    
        l2d = sqrt.(Rs .^ 2 .- Rts .^ 2)
        delta = atan.(Rts, l2d)
        phit = phis .+ delta .* direcs
        zt = xyz_src[3, :] .+ tan.(angles) .* l2d
    
        xtan[1, :] .= Rts .* cos.(phit)
        xtan[2, :] .= Rts .* sin.(phit)
        xtan[3, :] .= zt
    
        xyz_vec .= xtan .- xyz_src
    
        for n in 1:nbeams
            xyz_vec[:, n] .= xyz_vec[:, n] ./ sqrt(sum(xyz_vec[:, n] .^ 2))
        end
    
        Einj = [maximum(unit.energy.data) for unit in dd.nbi.unit]
        abeam = [unit.species.a for unit in dd.nbi.unit]
    
        for (n,unit) in enumerate(dd.nbi.unit)
            for i in 1:3
                part_frac[i,n] = maximum(unit.beam_current_fraction.data[i,:])
            end
        end
    
        for (n,unit) in enumerate(dd.nbi.unit)
            part_frac[:,n] = part_frac[:,n] ./ sum(part_frac[:,n])
        end
    
        for (n,unit) in enumerate(dd.nbi.unit)
            beamwidthpoly[2,n] = unit.beamlets_group[1].divergence_component[1].vertical
        end
    
        return nbeams, nv, xyz_src, xyz_vec, beamwidthpoly, Einj, part_frac, abeam
        
    end

    function get_pnbi(dd::IMAS.dd)
        pnbis = Vector{Float64}[]
        for ps in dd.pulse_schedule.nbi.unit
            power_downsampled = IMAS.interp1d(dd.pulse_schedule.nbi.time, ps.power.reference).(dd.equilibrium.time)
            push!(pnbis, power_downsampled)
        end
        return pnbis
    end

    pnbis = get_pnbi(dd)
    all_inputs[1].pnbi = pnbis

    # beam info isn't time dependent so store it in the first timeslice
    all_inputs[1].n_sources, all_inputs[1].nv, all_inputs[1].start_pos, all_inputs[1].beam_unit_vector, all_inputs[1].beam_width_polynomial_coefficients, all_inputs[1].injection_energy, all_inputs[1].particle_fraction, all_inputs[1].a_beam = gather_beams(dd)

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