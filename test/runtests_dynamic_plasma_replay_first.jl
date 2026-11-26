using FUSE
using Test

@testset "ActorDynamicPlasma.replay_before_sources" begin
    ini, act = FUSE.case_parameters(:ITER; init_from=:scalars)
    dd = IMAS.dd()
    FUSE.init(dd, ini, act)
    # replayed profiles: 20% hotter one step later
    δt = 0.1
    replay_dd = deepcopy(dd)
    IMAS.new_timeslice!(replay_dd.core_profiles, dd.global_time + δt)
    cp1d = replay_dd.core_profiles.profiles_1d[end]
    cp1d.electrons.temperature = 1.2 .* cp1d.electrons.temperature
    act.ActorReplay.replay_dd = replay_dd
    act.ActorCoreTransport.model = :replay
    # bootstrap source at the end of the step, relative to the bootstrap current of the profiles there
    function jbs_lag(replay_before_sources)
        dd1 = deepcopy(dd)
        FUSE.ActorDynamicPlasma(dd1, act; Δt=δt, Nt=1, replay_before_sources,
            evolve_current=false, evolve_equilibrium=false, evolve_pedestal=false, evolve_sawteeth=false)
        jbs = findfirst(:bootstrap_current, dd1.core_sources.source).profiles_1d[].j_parallel
        jbs_t = IMAS.Sauter_neo2021_bootstrap(dd1.equilibrium.time_slice[], dd1.core_profiles.profiles_1d[])
        return maximum(abs, jbs .- jbs_t) / maximum(abs, jbs_t)
    end
    @test jbs_lag(false) > 0.01
    @test jbs_lag(true) < 1e-6
end
