using FUSE
using Test

@testset "ActorSimpleNB energy components" begin
    ini, act = FUSE.case_parameters(:ITER; init_from=:scalars)
    dd = IMAS.dd()
    FUSE.init(dd, ini, act)
    ρ = dd.core_profiles.profiles_1d[].grid.rho_tor_norm
    for nbu in dd.nbi.unit
        nbu.beam_current_fraction.data = repeat([0.6, 0.3, 0.1], 1, length(nbu.beam_current_fraction.time))
    end
    function nb(; fast_ion_components=false, component_cross_sections=false)
        act.ActorSimpleNB.fast_ion_components = fast_ion_components
        act.ActorSimpleNB.component_cross_sections = component_cross_sections
        FUSE.ActorSimpleNB(dd, act)
        return findfirst(:nbi, dd.core_sources.source).profiles_1d[]
    end
    # one fast ion per component at E, E/2, E/3, together carrying every beam particle
    s = nb(; fast_ion_components=true)
    E = s.ion[1].fast_particles_energy
    @test length(s.ion) == 4 && [ion.fast_particles_energy for ion in s.ion[1:3]] ≈ [E, E / 2, E / 3]
    @test sum(ion.particles for ion in s.ion[1:3]) ≈ s.electrons.particles
    # ion impact and charge exchange at E/3 stop the E/3 component further out
    centroid(p) = IMAS.trapz(ρ, ρ .* p) / IMAS.trapz(ρ, p)
    c3 = centroid(s.ion[3].particles) # nb() refills the same source
    s_xs = nb(; fast_ion_components=true, component_cross_sections=true)
    @test centroid(s_xs.ion[3].particles) > c3
end
