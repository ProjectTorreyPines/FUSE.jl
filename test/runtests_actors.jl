using FUSE
using Test

@testset "fluxmatcher" begin
    dd, ini, act = FUSE.init(:ITER, init_from=:scalars)
    act.ActorFluxMatcher.max_iterations = 3
    act.ActorFluxMatcher.evolve_pedestal = true
    act.ActorFluxMatcher.evolve_densities = :flux_match
    FUSE.ActorFluxMatcher(dd, act)
end
