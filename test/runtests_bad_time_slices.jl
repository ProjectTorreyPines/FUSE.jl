using FUSE
using Test

# an equilibrium + core_profiles dd with slices at 0.1:0.1:0.8 s; the `bad` equilibrium slices get a
# rectangular grid with no ψ, so flux_surfaces throws on them (the others have no profiles_2d and are skipped)
function dd_with_bad_eq(bad=Int[])
    dd = IMAS.dd()
    eq = dd.equilibrium
    eq.time = collect(0.1:0.1:0.8)
    eq.vacuum_toroidal_field.b0 = collect(1.0:8.0)
    resize!(eq.time_slice, length(eq.time))
    for (k, eqt) in enumerate(eq.time_slice)
        eqt.time = eq.time[k]
        if k in bad
            resize!(eqt.profiles_2d, 1)
            eqt.profiles_2d[1].grid_type.index = 1
        end
    end
    for t in eq.time
        resize!(dd.core_profiles.profiles_1d, t)
    end
    dd.core_profiles.time = eq.time
    return dd
end

@testset "droppable_time_slices" begin
    # runs at the ends and single ones in between are dropped
    @test FUSE.droppable_time_slices([1, 2, 5, 8], 8) == ([1, 2, 5, 8], nothing)
    # two in a row in between end the shot there
    @test FUSE.droppable_time_slices([1, 4, 5], 8) == ([1, 4, 5, 6, 7, 8], 4)
    @test_throws ErrorException FUSE.droppable_time_slices([2, 3], 8)
    @test_throws ErrorException FUSE.droppable_time_slices(1:8, 8)
end

@testset "flux_surfaces!" begin
    dd = dd_with_bad_eq([1, 2, 5, 8])
    # a failed slice's magnetic axis is cleared, so the next slice does not start its search there
    eqt5 = dd.equilibrium.time_slice[5]
    eqt5.global_quantities.magnetic_axis.r, eqt5.global_quantities.magnetic_axis.z = -1e70, 0.0
    bad_time_slices = FUSE.BadTimeSlices()
    FUSE.flux_surfaces!(dd, bad_time_slices)
    @test ismissing(eqt5.global_quantities.magnetic_axis, :r)
    @test dd.equilibrium.time ≈ [0.3, 0.4, 0.6, 0.7] && dd.core_profiles.time ≈ [0.3, 0.4, 0.6, 0.7]
    @test [x.time for x in bad_time_slices.records] ≈ [0.1, 0.2, 0.5, 0.8]
    @test all(x -> x.stage == "equilibrium flux surfaces", bad_time_slices.records)
    # off: IMAS.flux_surfaces as before, its threaded error (TaskFailedException) thrown
    @test_throws Exception FUSE.flux_surfaces!(dd_with_bad_eq([4, 5]), nothing)
end

@testset "init_dropping_bad_start!" begin
    ini = FUSE.ParametersInits()
    ini.general.dd = dd_with_bad_eq()
    ini.time.simulation_start = 0.2
    dd = IMAS.dd()
    setup() = (ini.time.simulation_start < 0.35 ? error("bad start") : (dd.global_time = ini.time.simulation_start))
    bad_time_slices = FUSE.BadTimeSlices()
    FUSE.init_dropping_bad_start!(setup, dd, ini, bad_time_slices)
    @test ini.time.simulation_start ≈ 0.4 && dd.global_time ≈ 0.4
    @test [x.time for x in bad_time_slices.records] ≈ [0.2, 0.3]
    @test ini.general.dd.equilibrium.time ≈ [0.1, 0.4, 0.5, 0.6, 0.7, 0.8]
    @test_throws ErrorException FUSE.init_dropping_bad_start!(() -> error("bad start"), IMAS.dd(), ini, nothing)
    # a prepare() failure is not the start slice's: thrown at once, nothing dropped
    bad_time_slices = FUSE.BadTimeSlices()
    @test_throws ErrorException FUSE.init_dropping_bad_start!(() -> nothing, IMAS.dd(), ini, bad_time_slices; prepare=() -> error("no data"))
    @test isempty(bad_time_slices.records)
end
