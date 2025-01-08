import RABBIT
using Plots

#= =========== =#
#  ActorRABBIT  #
#= =========== =#
Base.@kwdef mutable struct FUSEparameters__ActorRABBIT{T<:Real} <: ParametersActor{T}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    _time::Float64 = NaN
end

mutable struct ActorRABBIT{D,P} <: SingleAbstractActor{D,P}
    dd::IMAS.dd{D}
    par::FUSEparameters__ActorRABBIT{P}
    outputs::Union{RABBIT.RABBIToutput,Vector{<:RABBIT.RABBIToutput}}
end

function ActorRABBIT(dd::IMAS.dd, par::FUSEparameters__ActorRABBIT; kw...)
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

    @assert length(dd.nbi.unit) == 1 "For now only one NBI unit is supported"

    all_inputs = FUSEtoRABBITinput(dd)
    actor.outputs = RABBIT.run_RABBIT(all_inputs; remove_inputs=true)
    return actor
end

function _finalize(actor::ActorRABBIT)
    dd = actor.dd
    cs = dd.core_sources
    output = actor.outputs

    eqt = dd.equilibrium.time_slice[]

    rho = [0.0; output.rho_data; 1.0]
    volume = IMAS.interp1d(eqt.profiles_1d.rho_tor_norm, eqt.profiles_1d.volume).(rho)
    area = IMAS.interp1d(eqt.profiles_1d.rho_tor_norm, eqt.profiles_1d.area).(rho)

    for (ps, nbu) in zip(dd.pulse_schedule.nbi.unit, dd.nbi.unit)
        power_launched = @ddtime(ps.power.reference)
        beam_energy = @ddtime(nbu.energy.data)

        @ddtime(nbu.power_launched.data = power_launched)

        # evaluate various source channels
        # here we also extend the RABBIT grid to the edges
        # so that the linear interpolation does not extrapolate in weird ways
        electrons_energy = output.powe_data[:, end]
        electrons_energy = [electrons_energy[1]; electrons_energy; electrons_energy[end]]
        total_ion_energy = output.powi_data[:, end]
        total_ion_energy = [total_ion_energy[1]; total_ion_energy; total_ion_energy[end]]
        electrons_particles = vec(sum(output.bdep_data[:, end, :]; dims=2))
        electrons_particles = [electrons_particles[1]; electrons_particles; electrons_particles[end]]
        j_parallel = output.jnbcd_data[:, end]
        j_parallel = [j_parallel[1]; j_parallel; j_parallel[end]]
        momentum_tor = vec(sum(output.torqdepo_data[:, end, :]; dims=2))
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
            momentum_tor
        )

        # add nbi fast ion particles source
        ion = resize!(source.profiles_1d[].ion, 1)[1]
        IMAS.ion_element!(ion, 1, nbu.species.a; fast=true)
        ion.particles = source.profiles_1d[].electrons.particles
        ion.particles_inside = source.profiles_1d[].electrons.particles_inside
        ion.fast_particles_energy = beam_energy

        # add neutron rate profile 
        neutral = resize!(source.profiles_1d[].neutral, 1)[1]
        neutral_particles = output.nrate_data[:, end]
        neutral_particles = [neutral_particles[1]; neutral_particles; neutral_particles[end]]
        neutral.particles = neutral_particles
        neutral.label = "neutrons"

        element = resize!(neutral.element, 1)[1]
        element.z_n = 0.0
        
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

        cp1d = dd.core_profiles.profiles_1d[time]

        inp.rho = cp1d.grid.rho_tor_norm
        inp.eq_rho = eqt.profiles_1d.rho_tor_norm
        inp.n_rho = length(inp.rho)

        inp.te = IMAS.interp1d(cp1d.grid.rho_tor_norm, cp1d.electrons.temperature).(inp.rho) .* eV_to_keV
        inp.dene = IMAS.interp1d(cp1d.grid.rho_tor_norm, cp1d.electrons.density).(inp.rho) .* cm3_to_m3
        inp.rot_freq_tor = inp.rho .* 0.0
        inp.zeff = IMAS.interp1d(cp1d.grid.rho_tor_norm, cp1d.zeff).(inp.rho)
        inp.ti = IMAS.interp1d(cp1d.grid.rho_tor_norm, cp1d.t_i_average).(inp.rho) .* eV_to_keV

        pnbis = Float64[]
        for ps in dd.pulse_schedule.nbi.unit
            push!(pnbis, IMAS.get_time_array(ps.power, :reference, time))
        end
        inp.pnbi = pnbis

        inp.n_sources = length(dd.nbi.unit)
        inp.injection_energy = dd.nbi.unit[1].energy.data
        inp.a_beam = [dd.nbi.unit[1].species.a]

        # the settings below reflect the default beams.dat input file for DIII-D from OMFIT
        inp.nv = 3
        inp.start_pos = [5.804921, 5.6625959, 0.0000000]
        inp.beam_unit_vector = [-0.80732277, -0.59011012, 0.0000000]
        inp.beam_width_polynomial_coefficients = [0.0000000, 0.023835, 0.0000000]
        inp.particle_fraction = [0.52422392, 0.3088602, 0.16691588]

        push!(all_inputs, inp)
    end

    if length(all_inputs) == 1
        inp = deepcopy(all_inputs[1])
        inp.time = -1e6
        push!(all_inputs, inp)
        reverse!(all_inputs)
    end

    return all_inputs

end