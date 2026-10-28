using FUSE
using Test

@testset "ActorSawteethSource.flatten_bootstrap" begin
    ini, act = FUSE.case_parameters(:ITER; init_from=:scalars)
    dd = IMAS.dd()
    FUSE.init(dd, ini, act)
    cp1d = dd.core_profiles.profiles_1d[]
    ρ = cp1d.grid.rho_tor_norm
    act.ActorSawteethSource.flat_factor = 1.0
    # spread of the non-inductive current inside the inversion radius, relative to its peak
    function spread(flatten_bootstrap)
        act.ActorSawteethSource.flatten_bootstrap = flatten_bootstrap
        FUSE.ActorSawteethSource(dd, act)
        j = IMAS.total_sources(dd.core_sources, cp1d; time0=dd.global_time, exclude_indexes=[7, 409], fields=[:j_parallel]).j_parallel
        ρ0 = ρ[findlast(abs.(cp1d.q) .< act.ActorSawteethSource.qmin_desired)]
        inner = ρ .< ρ0 - min(ρ0 / 4, 0.05)
        return (maximum(j[inner]) - minimum(j[inner])) / maximum(abs, j)
    end
    # the bootstrap current is left out of the flattening by default, and vanishes on axis
    @test spread(false) > 0.1
    @test spread(true) < 0.01
end
