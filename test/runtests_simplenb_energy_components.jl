using FUSE
using Test

@testset "ActorSimpleNB energy components" begin
    ini, act = FUSE.case_parameters(:ITER; init_from=:scalars)
    dd = IMAS.dd()
    FUSE.init(dd, ini, act)
    for nbu in dd.nbi.unit
        nbu.beam_current_fraction.data = repeat([0.6, 0.3, 0.1], 1, length(nbu.beam_current_fraction.time))
    end
    FUSE.ActorSimpleNB(dd, act)
    s = findfirst(:nbi, dd.core_sources.source).profiles_1d[]
    # one fast ion per component at E, E/2, E/3, together carrying every beam particle and all the ion heating
    E = s.ion[1].fast_particles_energy
    @test length(s.ion) == 4 && [ion.fast_particles_energy for ion in s.ion[1:3]] ≈ [E, E / 2, E / 3]
    @test sum(ion.particles for ion in s.ion[1:3]) ≈ s.electrons.particles
    @test sum(ion.energy for ion in s.ion[1:3]) ≈ s.total_ion_energy
end
