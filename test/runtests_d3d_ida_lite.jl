using FUSE
using Test

@testset "IDA-lite" begin
    # Phi_N from q: q = 1 gives Phi_N = psi_N
    ψ = collect(range(0.0, 1.1, 56))
    @test FUSE.ida_lite_phi_n(ψ, ones(length(ψ)))[1:51] ≈ ψ[1:51]
    @test FUSE.ida_lite_phi_n(ψ[2:end], ones(length(ψ) - 1)) === nothing

    # Zeff/ni route ladder
    @test FUSE.ida_lite_routes(:all, true, true, 6.0) == (true, true, "VB+CER")
    @test FUSE.ida_lite_routes(:all, false, false, 6.0) == (true, false, "VB (ne-fraction sigma)")
    @test FUSE.ida_lite_routes(:CER, true, false, 6.0) == (false, false, "none")

    # a synthetic direct-layout file
    nr, nt = 56, 7
    path = joinpath(mktempdir(), "IDA_123456_.cdf")
    FUSE.HDF5.h5open(path, "w") do f
        f["time"] = collect(range(1000.0, 2200.0, nt))
        f["psi_n"] = ψ
        f["n_e"] = [5e19 * (1 - 0.8x) for x in ψ, _ in 1:nt]
        f["T_e"] = [3e3 * (1 - 0.9x) + 10 for x in ψ, _ in 1:nt]
        f["T_12C6"] = [2.5e3 * (1 - 0.9x) + 10 for x in ψ, _ in 1:nt]
        f["Zeff"] = fill(2.0, nr, nt)
        f["Zeff_err"] = fill(0.2, nr, nt)
        f["q"] = [1.0 + 2x^2 for x in ψ, _ in 1:nt]
    end
    @test FUSE.ida_lite_file(dirname(path), 123456) == path
    @test FUSE.ida_lite_time(path) ≈ collect(range(1.0, 2.2, nt))
    @test FUSE.ida_lite_rho_grid(path) == count(<=(1.0), ψ)
    ida = FUSE.read_ida_lite(path; need_rotation=false)
    @test ida.layout === :direct && size(ida.n_e) == (nr, nt) && all(ida.vb_live) && !any(ida.cer_live)
    ni, tier = FUSE.ida_lite_ni(ida; ni_source=:all, impurity_Z=6.0)
    @test ni ≈ ida.n_e .* (6 - 2) / 5 && all(==("VB"), tier)
end

@testset "postdiction with external profiles" begin
    eq = collect(1.0:0.05:3.0)
    src = collect(1.2:0.02:2.8)
    sim = FUSE.postdictive_instants(1.05, 40 * 0.04, 40)
    @test sim[1] == 1.05 && sim[end] ≈ 1.05 + 40 * 0.04 && length(sim) == 81

    base = FUSE.core_profiles_time_base(eq, src, sim)
    @test issorted(base) && allunique(base)
    @test first(base) >= first(eq) && last(base) <= last(eq)
    # the zero-order hold at each instant lands on the tick made for it
    for t in sim
        t - 1e-9 >= first(eq) && @test base[searchsortedlast(base, t)] == t - 1e-9
    end

    a = FUSE.assign_profile_slices(base, src)
    # nearest slice; distances within 1 μs are a tie
    @test all(abs(base[i] - src[a.k[i]]) <= 0.01 + 1e-6 for i in eachindex(base) if src[1] <= base[i] <= src[end])
    # an exact tie goes to the slice not yet used
    b = FUSE.assign_profile_slices([1.01, 1.03], [1.0, 1.02, 1.04])
    @test b.k == [1, 2] && b.n_deadheat == 2

    k_of = Dict(zip(base, a.k))
    replaced = [t for t in base if src[1] <= t <= src[end]]
    p = FUSE.replay_pairing(base, sim[1:2:end], k_of, src, replaced; t_lo=src[1], t_hi=src[end])
    @test all(r.outcome === :source for r in p if src[1] + 1e-6 <= r.t_sim <= src[end])
    @test all(abs(r.source_time - r.t_sim) <= 0.011 for r in p if r.outcome === :source)
end
