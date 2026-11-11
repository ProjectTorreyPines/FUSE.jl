using FUSE
using Test
import FUSE: IMAS

@testset "init_expressions" begin
    all_expr_fields = intersect(FUSE.init_expressions()...)
    all_expr_fields_from_file = FUSE.load_init_expressions()
    @assert isempty(symdiff(all_expr_fields, all_expr_fields_from_file)) "expressions after init() have changed! run `make init_expressions` or `fusebot init_expressions` or `FUSE.init_expressions(;save=true)`"
end

@testset "restore_init_expressions keeps the stored density" begin
    dd = IMAS.dd()
    resize!(dd.core_profiles.profiles_1d, 2)
    rho = collect(range(0.0, 1.0, 11))
    n = 1E19 .* (1.0 .- 0.5 .* rho .^ 2)
    for (k, cp1d) in enumerate(dd.core_profiles.profiles_1d)
        cp1d.time = 0.1 * k
        cp1d.grid.rho_tor_norm = rho
        cp1d.electrons.density_thermal = 2.0 .* n
        resize!(cp1d.ion, 2)
        # main ion in the state left by quasi-neutrality without fast ions: total density stored, thermal density derived
        cp1d.ion[1].density = 1.5 .* n
        cp1d.ion[1].density_fast = zero(n)
        cp1d.ion[2].density_thermal = 0.1 .* n
    end
    dd.core_profiles.time = [0.1, 0.2]
    dd.global_time = 0.2
    FUSE.restore_init_expressions!(dd; verbose=false)
    for cp1d in dd.core_profiles.profiles_1d
        @test IMAS.hasdata(cp1d.ion[1], :density_thermal)
        @test cp1d.ion[1].density_thermal ≈ 1.5 .* n
        @test cp1d.ion[1].density ≈ 1.5 .* n
        @test cp1d.ion[2].density_thermal ≈ 0.1 .* n
    end
    # the first ion of the first slice is where the stored density used to be discarded
    @test !IMAS.hasdata(dd.core_profiles.profiles_1d[1].ion[1], :density)
end
