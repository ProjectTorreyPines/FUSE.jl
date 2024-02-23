using FUSE
using Test

@testset "fluxmatcher" begin
    dd, ini, act = FUSE.init(:ITER, init_from=:scalars)
    act.ActorFluxMatcher.max_iterations = 2
    act.ActorFluxMatcher.evolve_pedestal = true
    act.ActorFluxMatcher.evolve_densities = :flux_match
    act.ActorTGLF.model = :TJLF
    FUSE.ActorFluxMatcher(dd, act)
    act.ActorTGLF.model = :TGLFNN
    FUSE.ActorFluxMatcher(dd, act)
end
