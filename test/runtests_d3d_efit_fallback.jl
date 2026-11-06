using FUSE
using Test

@testset "d3d_efit_interior_gaps" begin
    t = collect(0.1:0.02:3.0)
    gaps, dt = FUSE.d3d_efit_interior_gaps(t, 2.0, Inf)
    @test isempty(gaps) && dt ≈ 0.02
    t2 = vcat(t[t.<1.5], t[t.>2.0])
    gaps, _ = FUSE.d3d_efit_interior_gaps(t2, 2.0, Inf)
    @test length(gaps) == 1 && gaps[1][1] < 1.5 && gaps[1][2] > 2.0
    # the first interval is never a gap (the simulation starts at the second slice)
    t3 = vcat(0.0, t)
    @test isempty(FUSE.d3d_efit_interior_gaps(t3, 2.0, Inf)[1])
    # hard gap: an interval >= hard_gap_s even when the median spacing is large
    @test length(FUSE.d3d_efit_interior_gaps(t2, 100.0, 0.035)[1]) == 1
end

@testset "d3d_efit_fallback!" begin
    # an equilibrium IDS with slices at `times`, toroidal field b0
    function eqdd(times, b0)
        dd = IMAS.dd()
        eq = dd.equilibrium
        eq.time = times
        eq.vacuum_toroidal_field.r0 = 1.6955
        eq.vacuum_toroidal_field.b0 = fill(b0, length(times))
        resize!(eq.time_slice, length(times))
        for (ts, t) in zip(eq.time_slice, times)
            ts.time = t
            ts.global_quantities.ip = 1e6
            resize!(ts.profiles_2d, 1)
            ts.profiles_2d[1].grid.dim1 = collect(range(0.8, 2.5, 33))
            ts.profiles_2d[1].grid.dim2 = collect(range(-1.6, 1.6, 65))
            ts.profiles_1d.psi = collect(range(0.0, 1.0, 21))
        end
        return dd
    end
    primary = vcat(0.1:0.02:1.0, 1.5:0.02:2.0)
    fallback = collect(0.1:0.025:2.0)

    dd1 = eqdd(primary, -2.0)
    prov = FUSE.d3d_efit_fallback!(dd1, () -> eqdd(fallback, -2.1), "EFIT02", "EFIT01"; gap_factor=2.0, hard_gap_s=0.035)
    eq = dd1.equilibrium
    filled = [t for t in fallback if 1.0 < t < 1.5]
    # the gap is filled with the fallback slices strictly inside it, each with its own tree's b0
    @test eq.time == sort(vcat(primary, filled)) && issorted(eq.time)
    @test eq.vacuum_toroidal_field.b0 == [t in filled ? -2.1 : -2.0 for t in eq.time]
    @test eq.code.name == "EFIT02+EFIT01"
    @test prov["slice_trees"] == [t in filled ? "EFIT01" : "EFIT02" for t in eq.time]
    @test prov["gaps"][1]["n_filled"] == length(filled)

    # trees that cannot be spliced (different r0): the primary is kept as fetched
    dd1 = eqdd(primary, -2.0)
    dd_fb = eqdd(fallback, -2.1)
    dd_fb.equilibrium.vacuum_toroidal_field.r0 = 1.7
    FUSE.d3d_efit_fallback!(dd1, () -> dd_fb, "EFIT02", "EFIT01"; gap_factor=2.0, hard_gap_s=0.035)
    @test dd1.equilibrium.time == primary

    # too few primary slices: the fallback replaces it
    dd1 = eqdd([0.5, 0.6], -2.0)
    FUSE.d3d_efit_fallback!(dd1, () -> eqdd(fallback, -2.1), "EFIT02", "EFIT01"; gap_factor=2.0, hard_gap_s=0.035)
    @test dd1.equilibrium.time == fallback
end
