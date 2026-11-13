using FUSE
using Test
import FUSE.ForwardDiff
import FUSE.ForwardDiff: value, partials
import FUSE.LinearAlgebra: norm

# Actors on a dd whose element type is a ForwardDiff.Dual: every derivative is compared with central
# finite differences of the same actor on a Float64 dd. Inputs: scale factors of the electron
# temperature (partial 1) and of the densities (partial 2).

const D2 = IMAS.InputsDual{2}

function scale_inputs!(dd, s)
    for cp1d in dd.core_profiles.profiles_1d
        cp1d.electrons.temperature = cp1d.electrons.temperature .* s[1]
        cp1d.electrons.density_thermal = cp1d.electrons.density_thermal .* s[2]
        for ion in cp1d.ion
            ion.density_thermal = ion.density_thermal .* s[2]
        end
    end
    return dd
end
seeded(dd0) = scale_inputs!(IMAS.retype(D2, dd0), [IMAS.seed(D2, 1), IMAS.seed(D2, 2)])
perturbed(dd0, k, h) = scale_inputs!(deepcopy(dd0), [k == 1 ? 1.0 + h : 1.0, k == 2 ? 1.0 + h : 1.0])

profiles(dd) = (cp1d = dd.core_profiles.profiles_1d[]; vcat(cp1d.electrons.temperature, cp1d.t_i_average, cp1d.electrons.density_thermal))

"""
AD against central finite differences for `run!(dd)` followed by `out(dd)`: returns the largest difference
over the outputs, relative to the largest derivative of each input
"""
function ad_vs_fd(dd0, run!, out; h=1e-4)
    dd = seeded(dd0)
    run!(dd)
    y = out(dd)
    ddf = deepcopy(dd0)
    run!(ddf)
    @test value.(y) ≈ out(ddf) rtol = 1e-8
    @test all(v -> all(isfinite, partials(v)), y)
    worst = 0.0
    for k in 1:2
        fd = (out(run!(perturbed(dd0, k, +h))) .- out(run!(perturbed(dd0, k, -h)))) ./ (2h)
        ad = partials.(y, k)
        worst = max(worst, maximum(abs, ad .- fd) / max(maximum(abs, fd), floatmin()))
    end
    return worst
end

ini, act = FUSE.case_parameters(:D3D, :default)
const dd0 = IMAS.dd()
FUSE.init(dd0, ini, act)
act.ActorTGLF.model = :TGLFNN
act.ActorTGLF.warn_nn_train_bounds = false
act.ActorFluxMatcher.verbose = false

@testset "ad helpers" begin
    x = [IMAS.seed(D2, 1) * 2.0, IMAS.seed(D2, 2) * 3.0]
    @test FUSE.dual_values(x) == [2.0, 3.0]
    @test FUSE.partials_matrix(x) == [2.0 0.0; 0.0 3.0]
    @test FUSE.with_partials(D2, [2.0, 3.0], [2.0 0.0; 0.0 3.0]) == x
    @test all(iszero, partials(FUSE.strip_partials(x[1])))
    # with duals a type error is never swallowed; other element types are unaffected
    @test_throws MethodError try
        throw(MethodError(sin, ("a",)))
    catch e
        FUSE.rethrow_if_dual_type_error(D2, e)
    end
    @test (try
        throw(MethodError(sin, ("a",)))
    catch e
        FUSE.rethrow_if_dual_type_error(Float64, e)
    end) === nothing
end

@testset "sources with duals" begin
    run!(dd) = (IMAS.intrinsic_sources!(dd); dd)
    out(dd) = (ts = IMAS.total_sources(dd.core_sources, dd.core_profiles.profiles_1d[]; time0=dd.global_time); vcat(ts.electrons.energy, ts.total_ion_energy))
    @test ad_vs_fd(dd0, run!, out) < 1e-4

    run2!(dd) = (FUSE.ActorNeutralFueling(dd, act); dd)
    out2(dd) = (i = findfirst(s -> s.identifier.index == 108, dd.core_sources.source); dd.core_sources.source[i].profiles_1d[].electrons.particles)
    @test ad_vs_fd(dd0, run2!, out2) < 1e-3
end

@testset "pedestal with duals" begin
    tb, ts = IMAS.pedestal_search_tolerance.brent[], IMAS.pedestal_search_tolerance.simplex[]
    IMAS.pedestal_search_tolerance.brent[] = 1e-12
    IMAS.pedestal_search_tolerance.simplex[] = 1e-14
    try
        for model in (:WPED, :EPED)
            a = deepcopy(act)
            a.ActorPedestal.model = model
            run!(dd) = (FUSE.ActorPedestal(dd, a; ip_from=:core_profiles, βn_from=:core_profiles, ne_from=:core_profiles, zeff_from=:core_profiles, rho_nml=0.775, rho_ped=0.8); dd)
            @test ad_vs_fd(dd0, run!, profiles) < 2e-2
        end
    finally
        IMAS.pedestal_search_tolerance.brent[] = tb
        IMAS.pedestal_search_tolerance.simplex[] = ts
    end
end

@testset "current diffusion with duals" begin
    a = deepcopy(act)
    a.ActorCurrent.model = :QED
    a.ActorQED.Δt = 0.05
    run!(dd) = (FUSE.ActorCurrent(dd, a; ip_from=:core_profiles, vloop_from=:equilibrium); dd)
    out(dd) = dd.core_profiles.profiles_1d[].j_tor
    @test ad_vs_fd(dd0, run!, out; h=1e-3) < 1e-3
end

@testset "flux matching with duals" begin
    a = deepcopy(act)
    a.ActorPedestal.model = :none
    a.ActorFluxMatcher.evolve_pedestal = false
    a.ActorFluxMatcher.consistent_fast_pressure = true

    # residual and its derivatives at fixed unknowns
    function context(dd; norms=nothing)
        actor = FUSE.ActorFluxMatcher(dd, a.ActorFluxMatcher, a)
        cp1d = dd.core_profiles.profiles_1d[]
        IMAS.intrinsic_sources!(dd)
        IMAS.refreeze!(cp1d, :j_non_inductive)
        IMAS.refreeze!(cp1d, :j_ohmic)
        IMAS.freeze!(cp1d, :rotation_frequency_tor_sonic)
        initial = FUSE.cp1d_copy_primary_quantities(cp1d)
        z = Vector{Float64}(FUSE.dual_values(FUSE.scale_z_profiles(FUSE.pack_z_profiles(cp1d, actor.par).z_profiles)))
        nch = length(z) ÷ length(actor.par.rho_transport)
        actor.norms = norms === nothing ? fill(eltype(dd)(NaN), nch) : eltype(dd).(norms)
        return actor, initial, z
    end
    c0 = context(deepcopy(dd0))
    z = c0[3]
    F0 = FUSE.flux_match_errors(c0[1], copy(z), c0[2]).errors
    norms = copy(c0[1].norms)
    residual(dd) = (c = context(dd; norms); FUSE.flux_match_errors(c[1], eltype(dd).(z), c[2]).errors)

    Fd = residual(seeded(dd0))
    @test value.(Fd) ≈ F0 rtol = 1e-9
    h = 1e-4
    for k in 1:2
        fd = (residual(perturbed(dd0, k, +h)) .- residual(perturbed(dd0, k, -h))) ./ (2h)
        @test norm(partials.(Fd, k) .- fd) < 1e-3 * norm(fd)
    end

    # Jacobian with respect to the unknowns
    J, FJ = FUSE.flux_match_jacobian(c0[1], z, c0[2])
    @test FJ ≈ F0 rtol = 1e-8
    for col in (1, length(z) ÷ 2, length(z))
        δ = 1e-4 * max(abs(z[col]), 1.0)
        zp = copy(z); zp[col] += δ
        zm = copy(z); zm[col] -= δ
        fd = (FUSE.flux_match_errors(c0[1], zp, c0[2]).errors .- FUSE.flux_match_errors(c0[1], zm, c0[2]).errors) ./ (2δ)
        @test norm(J[:, col] .- fd) < 1e-3 * norm(fd)
    end

    # the implicit derivative of the solution moves the unknowns so that the residual stays where it is
    Fθ = FUSE.partials_matrix(Fd)
    dz = -(FUSE.LinearAlgebra.pinv(J; rtol=1e-10) * Fθ)
    hh = 1e-3
    for k in 1:2
        c = context(perturbed(dd0, k, hh); norms)
        fixed = norm(FUSE.flux_match_errors(c[1], copy(z), c[2]).errors .- F0)
        shifted = norm(FUSE.flux_match_errors(c[1], z .+ hh .* dz[:, k], c[2]).errors .- F0)
        @test shifted < 0.05 * fixed
    end

    # the actor: same values as a Float64 run when the solver follows the same iterates, finite derivatives
    a.ActorFluxMatcher.algorithm = :simple_dfsane
    a.ActorFluxMatcher.max_iterations = 30
    dd = seeded(dd0)
    actor = FUSE.ActorFluxMatcher(dd, a)
    ddf = deepcopy(dd0)
    FUSE.ActorFluxMatcher(ddf, a)
    @test value.(profiles(dd)) ≈ profiles(ddf) rtol = 1e-6
    @test all(v -> all(isfinite, partials(v)), profiles(dd))
    @test actor.ad_info[:jacobian_copy_mismatch] < 1e-8
end

@testset "equilibrium with duals" begin
    a = deepcopy(act)
    a.ActorEquilibrium.model = :FRESCO
    a.ActorFRESCO.nR = a.ActorFRESCO.nZ = 33
    run!(dd) = (FUSE.ActorEquilibrium(dd, a; ip_from=:core_profiles); dd)
    out(dd) = (g = dd.equilibrium.time_slice[].global_quantities; [g.beta_normal, g.beta_pol, g.magnetic_axis.r])
    dd = seeded(dd0)
    run!(dd)
    y = out(dd)
    @test value.(y) ≈ out(run!(deepcopy(dd0))) rtol = 1e-4
    h = 0.05 # the equilibrium solve converges to 1e-4: a large step keeps the finite difference above its noise
    for k in 1:2
        fd = (out(run!(perturbed(dd0, k, +h))) .- out(run!(perturbed(dd0, k, -h)))) ./ (2h)
        @test partials.(y, k)[1:2] ≈ fd[1:2] rtol = 0.1
    end
end

@testset "TEQUILA refuses duals" begin
    a = deepcopy(act)
    a.ActorEquilibrium.model = :TEQUILA
    @test_throws ArgumentError FUSE.ActorEquilibrium(seeded(dd0), a; ip_from=:core_profiles)
end

@testset "ActorStationaryPlasma with duals" begin
    a = deepcopy(act)
    a.ActorPedestal.model = :none
    a.ActorFluxMatcher.evolve_pedestal = false
    a.ActorFluxMatcher.algorithm = :simple_dfsane
    a.ActorFluxMatcher.max_iterations = 30
    a.ActorEquilibrium.model = :FRESCO
    a.ActorFRESCO.nR = a.ActorFRESCO.nZ = 33
    a.ActorCurrent.model = :none
    a.ActorStationaryPlasma.max_iterations = 2

    dd = seeded(dd0)
    actor = FUSE.ActorStationaryPlasma(dd, a)
    ddf = deepcopy(dd0)
    FUSE.ActorStationaryPlasma(ddf, a)

    @test eltype(dd) === D2
    @test value.(profiles(dd)) ≈ profiles(ddf) rtol = 5e-2 # the loosely converged solves amplify rounding differences
    @test all(v -> all(isfinite, partials(v)), profiles(dd))
    @test any(v -> !iszero(partials(v, 1)), profiles(dd)) && any(v -> !iszero(partials(v, 2)), profiles(dd))
    eqt = dd.equilibrium.time_slice[]
    @test all(v -> all(isfinite, partials(v)), eqt.profiles_1d.q)
    @test any(v -> !iszero(partials(v, 2)), eqt.profiles_1d.q) # the equilibrium responds to the density
    @test haskey(actor.ad_info, :sensitivity_error) && all(isfinite, actor.ad_info[:sensitivity_error])
    @test actor.actor_tr.tr_actor.ad_info[:jacobian_copy_mismatch] < 1e-8
end

@testset "ActorDynamicPlasma with duals" begin
    ini_t, act_t = FUSE.case_parameters(:ITER; init_from=:scalars, time_dependent=true)
    ini_t.time.pulse_shedule_time_basis = range(0, 300; step=1.0)
    ini_t.time.simulation_start = 50.0
    ddt = IMAS.dd()
    FUSE.init(ddt, ini_t, act_t)

    act_t.ActorTGLF.model = :TGLFNN
    act_t.ActorTGLF.warn_nn_train_bounds = false
    act_t.ActorFluxMatcher.verbose = false
    act_t.ActorFluxMatcher.max_iterations = -20
    act_t.ActorFluxMatcher.relax = 1.0
    act_t.ActorSawteethSource.flat_factor = 1.0
    for m in (:ic_model, :lh_model, :pellet_model)
        setproperty!(act_t.ActorSources, m, :none)
    end
    act_t.ActorDynamicPlasma.Nt = 2
    act_t.ActorDynamicPlasma.Δt = 2.0
    for f in (:evolve_current, :evolve_equilibrium, :evolve_transport, :evolve_sources, :evolve_pedestal, :evolve_sawteeth)
        setproperty!(act_t.ActorDynamicPlasma, f, true)
    end
    act_t.ActorDynamicPlasma.evolve_pf_active = false

    dd = seeded(ddt)
    FUSE.ActorDynamicPlasma(dd, act_t; verbose=false)
    ddf = deepcopy(ddt)
    FUSE.ActorDynamicPlasma(ddf, act_t; verbose=false)

    @test dd.global_time ≈ ddf.global_time
    @test value.(profiles(dd)) ≈ profiles(ddf) rtol = 5e-2
    @test all(v -> all(isfinite, partials(v)), profiles(dd))
    @test any(v -> !iszero(partials(v, 2)), profiles(dd))
    @test all(v -> all(isfinite, partials(v)), dd.core_profiles.profiles_1d[].q)
end
