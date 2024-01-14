using FUSE
using Test

@testset "yaml_workflow" begin
    flow="""
    - init
    - ActorEquilibrium:
        symmetrize: true
        ip_from: :pulse_schedule
    """

    ini,act = FUSE.case_parameters(:KDEMO);
    actors = FUSE.yaml_workflow(flow, IMAS.dd(), ini, act)

    @test length(actors) == 1
end
