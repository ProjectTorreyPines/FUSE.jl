using FUSE
using Test

#include("test_zmq_actor.jl")  # requires GSLite connection — disabled until GSLite integration is merged

@testset "fluxmatcher" begin
    ini, act = FUSE.case_parameters(:ITER; init_from=:scalars)
    dd = IMAS.dd()
    FUSE.init(dd, ini, act)
    act.ActorFluxMatcher.max_iterations = 2
    act.ActorFluxMatcher.evolve_pedestal = true

    act.ActorFluxMatcher.algorithm = :simple
    act.ActorTGLF.model = :TJLF
    FUSE.ActorFluxMatcher(dd, act)
    act.ActorFluxMatcher.algorithm = :anderson
    act.ActorTGLF.model = :TGLFNN
    FUSE.ActorFluxMatcher(dd, act)
end

@testset "Mars Actor CHEASE run" begin
    ini, act = FUSE.case_parameters(:D3D, :default)
    dd = IMAS.dd()
    FUSE.init(dd, ini, act)

    # control execution of CHEASE & MARS for testing purposes
    act.ActorMars.run_equilibrium = true
    act.ActorMars.run_MHD = false
    act.ActorMars.wall_type = :no_wall # Begin with NO wall
    #act.ActorMars.chease_exec = "/path/to/MarsQ_package/CheaseMerge/chease.x"
    #act.ActorMars.mars_exec = "/path/to/MarsQ_package/MarsQ/marsq.x"

    # Configure CHEASE parameters for testing
    chease_overrides = (NPSI=64, NVEXP=2, NCSCAL=4)

    # configure MARS parameters for testing
    mars_overrides = FUSE.MarsOverrides()
    mars_overrides.BASIC[:M1] = -10
    mars_overrides.BASIC[:NV] = 120 # moves the IW in from the original CHEASE location

    if !isfile(act.ActorMars.chease_exec)
        @test_skip "CHEASE executable not found at $(act.ActorMars.chease_exec). Skipping CHEASE & MARS test."
    else
        # The actor creates and manages its own run directory; keep it across the
        # chained runs below (restart / MHD-only depend on each other's files).
        act.ActorMars.clear_workdir = false

        @info "Test 1: Clean CHEASE equilibrium run"
        actor = FUSE.ActorMars(dd, act; chease_overrides)
        run_dir = actor.par.save_dir
        @testset "CHEASE output files" begin
            @test isfile(joinpath(run_dir, "EXPEQ"))
            @test isfile(joinpath(run_dir, "datain"))
        end

        @info "Test 2: CHEASE with resistive wall & higher beta"
        act.ActorMars.number_surfaces = 2
        act.ActorMars.wall_type = :limiter
        FUSE.ActorMars(dd, act; chease_overrides=(CFBAL=1.1,), save_dir=run_dir)

        @info "Test 3: CHEASE restart from existing equilibrium"
        act.ActorMars.restart_equilibrium = true
        FUSE.ActorMars(dd, act; chease_overrides=(NCSCAL=2,), save_dir=run_dir) # keep Ip fixed

        @info "Test 4: MARS MHD stability run"
        act.ActorMars.run_equilibrium = false # do NOT rerun CHEASE, just run MARS on the existing equilibrium
        act.ActorMars.run_MHD = true
        FUSE.ActorMars(dd, act; mars_overrides, save_dir=run_dir)

        @testset "MARS output files" begin
            @test isfile(joinpath(run_dir, "log_mars"))
            @test isfile(joinpath(run_dir, "RESULT.OUT"))
        end

        @testset "MARS outputs in dd.mhd_linear" begin
            mode = dd.mhd_linear.time_slice[].toroidal_mode[1]

            # growth rate & frequency (converted to SI in _finalize)
            @test mode.n_tor == -1
            @test isfinite(mode.growthrate)
            @test isfinite(mode.frequency)

            # displacement eigenfunction on the (s, m) grid
            ns = length(mode.plasma.grid.dim1)
            nm = length(mode.plasma.grid.dim2)
            @test size(mode.plasma.displacement_perpendicular.real) == (ns, nm)
            @test size(mode.plasma.displacement_perpendicular.imaginary) == (ns, nm)

            # Alfvén-time profile on the radial grid
            @test length(mode.plasma.tau_alfven) == ns
            @test all(mode.plasma.tau_alfven .> 0)

            # real-space R,Z geometry on the (s, χ) grid
            cs = mode.plasma.coordinate_system
            @test size(cs.r) == size(cs.z)
            @test size(cs.r, 1) == length(cs.grid.dim1)
            @test size(cs.r, 2) == length(cs.grid.dim2)

            # plot recipe for the mode structure works
            @test (FUSE.plot(dd.mhd_linear); true)
        end
    end
end
        

@testset "ActorTJLFEP" begin
    # End-to-end TGLF-EP -> ALPHA EP transport on an ITER dd:
    #   TJLFEP.runTHD finds the critical EP density/pressure gradients (the scan's last
    #   point is the separatrix ir=NR (rho~1), where a singular TGLF Hermite matrix now
    #   degrades gracefully to "stable" instead of erroring), ALPHA.run_alpha integrates
    #   them into the EP profiles/flux, and _finalize writes the fast-ion population to
    #   core_profiles and the EP flux to core_transport.
    # Kept small (ngrid=51, SCAN_N=2, N_BASIS=2, :marginal solver) for a fast, robust
    # integration test. The TGLF-EP eigenvalue solve itself is also covered by TJLFEP's
    # own nb6 regression.
    ini, act = FUSE.case_parameters(:ITER; init_from=:ods)
    ini.core_profiles.ngrid = 51
    dd = IMAS.dd()
    FUSE.init(dd, ini, act)

    is_ep = 3
    act.ActorTJLFEP.is_ep = is_ep
    act.ActorTJLFEP.rho_scan = [0.41, 0.61]   # SCAN_N=2; INPUT_PROFILE_METHOD=2 makes the last point ir=NR (separatrix)
    act.ActorTJLFEP.n_basis = 2
    act.ActorTJLFEP.solver = :grid            # fast Fortran-equivalent path for this smoke test (actor default is now :ad)
    act.ActorTJLFEP.alpha_solver = :marginal

    actor = FUSE.ActorTJLFEP(dd, act)         # full pipeline: runTHD -> run_alpha -> finalize

    @testset "TGLF-EP stability metrics" begin
        @test length(actor.SFmin) == length(act.ActorTJLFEP.rho_scan)
        @test all(isfinite, actor.SFmin)
        @test all(actor.SFmin .> 0)
        @test actor.alpha !== nothing
        @test length(actor.rho_grid) == length(dd.core_profiles.profiles_1d[].grid.rho_tor_norm)
    end

    @testset "EP profiles in dd.core_profiles" begin
        cp1d = dd.core_profiles.profiles_1d[]
        ep_ion = cp1d.ion[is_ep]
        nfast = ep_ion.density_fast
        @test length(nfast) == length(cp1d.grid.rho_tor_norm)
        @test all(isfinite, nfast)
        @test all(nfast .>= 0.0)
        @test maximum(nfast) > 0.0
        pfast = ep_ion.pressure_fast_parallel .+ 2.0 .* ep_ion.pressure_fast_perpendicular
        @test all(isfinite, pfast)
        @test all(pfast .>= 0.0)
        @test maximum(pfast) > 0.0
    end

    @testset "EP flux in dd.core_transport" begin
        model = dd.core_transport.model[:anomalous]
        m1d = model.profiles_1d[]
        @test occursin("TJLFEP-ALPHA", model.identifier.name)
        @test length(m1d.ion) >= 1
        @test all(isfinite, m1d.ion[1].particles.flux)
        @test all(isfinite, m1d.ion[1].energy.flux)
    end
end

@testset "SOL-pedestal separatrix coupling" begin
    ini, act = FUSE.case_parameters(:D3D, :default)
    dd = FUSE.init(ini, act)
    rho_nml, rho_ped = act.ActorFluxMatcher.rho_transport[end-1:end]
    actor_ped = FUSE.ActorPedestal(dd, act.ActorPedestal, act; ip_from=:core_profiles, βn_from=:core_profiles, ne_from=:pulse_schedule, zeff_from=:pulse_schedule, rho_nml, rho_ped)

    # SOL model output: separatrix point of edge_profiles
    ep1d = resize!(dd.edge_profiles.profiles_1d)
    ep1d.grid.rho_pol_norm = [1.0]
    ep1d.electrons.density = [1e19]
    ep1d.electrons.temperature = [120.0]
    ep1d.t_i_average = [200.0]

    # full update, then under-relaxed update
    actor_ped.par.Te_sep = 80.0
    sep = FUSE.couple_separatrix!(actor_ped, dd)
    @test sep.Te_sep ≈ 120.0
    @test ismissing(sep.Ti_sep) # Ti_sep is not coupled by default
    ep1d.electrons.temperature[1] = 100.0
    @test FUSE.couple_separatrix!(actor_ped, dd; relaxation=0.5).Te_sep ≈ 110.0
    @test actor_ped.par.Te_sep ≈ 110.0
    # Ti_sep: first coupling takes the SOL value, then it is under-relaxed (SOL Te held at the current Te_sep)
    ep1d.electrons.temperature[1] = 110.0
    sep = FUSE.couple_separatrix!(actor_ped, dd; relaxation=0.5, couple_ti=true)
    @test sep.Ti_sep ≈ 200.0
    @test sep.Te_sep ≈ 110.0
    ep1d.t_i_average[1] = 160.0
    @test FUSE.couple_separatrix!(actor_ped, dd; relaxation=0.5, couple_ti=true).Ti_sep ≈ 180.0
    @test_throws AssertionError FUSE.couple_separatrix!(actor_ped, dd; relaxation=0.0)

    # Te_sep changed after construction must reach the sub-pedestal actors at the next step
    FUSE.sync_separatrix!(actor_ped)
    for sub_actor in (actor_ped.eped_actor, actor_ped.wped_actor, actor_ped.analytic_actor)
        @test sub_actor.par.Te_sep ≈ 110.0
        @test sub_actor.par.Ti_sep ≈ 180.0
    end
    FUSE.finalize(FUSE.step(actor_ped))
    cp1d = dd.core_profiles.profiles_1d[]
    @test cp1d.electrons.temperature[end] ≈ 110.0 rtol = 0.2
    @test cp1d.t_i_average[end] ≈ 180.0 rtol = 0.2
    @test FUSE.separatrix_ti(actor_ped.par, 1.0) ≈ 180.0
    actor_ped.par.Ti_sep = missing
    @test FUSE.separatrix_ti(actor_ped.par, 1.5) ≈ 1.5 * 110.0

    # SOL models that do not populate the separatrix are rejected
    act.ActorStationaryPlasma.couple_sol = true
    act.ActorSOL.model = :none
    @test_throws AssertionError FUSE.ActorStationaryPlasma(dd, act.ActorStationaryPlasma, act)
end

@testset "pedestal Te/Ti split" begin
    # pure closure functions
    ne, ni, ae, ai = 4e19, 3.5e19, 300.0, 100.0
    Pstar = 4e19 * 500.0
    for (chi_i_neo, f) in ((0.0, 1.0), (0.2, 1.0), (0.2, 0.1), (0.2, 0.0))
        χe = FUSE.pedestal_chi_e(Pstar, ne, ni, ae, ai, chi_i_neo, f)
        @test χe > 0.0
        @test ne * ae / χe + ni * ai / (chi_i_neo + f * χe) ≈ Pstar
    end
    @test isnan(FUSE.pedestal_chi_e(-1.0, ne, ni, ae, ai, 0.2, 1.0))
    @test isnan(FUSE.pedestal_chi_e(Pstar, ne, ni, ae, 1e6, 0.2, 0.0)) # neoclassical ions alone exceed Pstar

    split0 = (model=:ratio, Te_sep=80.0, Ti_sep=200.0, Ti_over_Te=1.2, rho_top=0.93, f=NaN, chi_i_neo=0.0, ae, ai, Qe=NaN, Qi=NaN, De=NaN, De_over_chi_e_expected=NaN)
    p_ped = 1.5e4 # Pa
    for (model, f, chi_i_neo) in ((:ratio, NaN, 0.0), (:increment, NaN, 0.0), (:heat_flux, 2.0, 0.0), (:fingerprint, 0.1, 0.2))
        split = merge(split0, (; model, f, chi_i_neo))
        T = FUSE.pedestal_temperatures(split, p_ped, ne, ni)
        @test FUSE.IMAS.mks.e * (ne * T.Te_ped + ni * T.Ti_ped) ≈ p_ped
        if model == :ratio
            @test T.Ti_ped / T.Te_ped ≈ 1.2
        elseif model == :increment
            @test (T.Ti_ped - 200.0) / (T.Te_ped - 80.0) ≈ 1.2
        elseif model == :heat_flux
            @test (T.Ti_ped - 200.0) / (T.Te_ped - 80.0) ≈ (ai / ae) / f
        end
        # closure on Te_ped (WPED) is consistent with closure on pressure (EPED)
        @test FUSE.pedestal_ti(split, T.Te_ped).Ti_ped ≈ T.Ti_ped
    end
    # bracketed bisection
    sol = FUSE.solve_pedestal_ti(T -> 500.0 - 0.5 * T - T, 100.0, 1000.0)
    @test sol.converged
    @test sol.T ≈ 500.0 / 1.5 atol = 0.5
    @test !FUSE.solve_pedestal_ti(T -> 1.0, 100.0, 1000.0).converged

    # pressure below separatrix temperatures at pedestal density: separatrix Ti/Te is kept
    T = FUSE.pedestal_temperatures(merge(split0, (model=:increment,)), 100.0, ne, ni)
    @test T.Ti_ped / T.Te_ped ≈ 200.0 / 80.0

    # through ActorPedestal
    ini, act = FUSE.case_parameters(:D3D, :default)
    dd0 = FUSE.init(ini, act)
    rho_nml, rho_ped = act.ActorFluxMatcher.rho_transport[end-1:end]
    kw = (ip_from=:core_profiles, βn_from=:core_profiles, ne_from=:pulse_schedule, zeff_from=:pulse_schedule, rho_nml, rho_ped, Te_sep=80.0, Ti_sep=200.0)
    ratio = Dict{Tuple{Symbol,Symbol},Float64}()
    for model in (:EPED, :WPED), (Ti_ped_model, fingerprint_mode) in ((:ratio, :MHD), (:increment, :MHD), (:heat_flux, :MHD), (:fingerprint, :MHD), (:fingerprint, :MTM))
        dd = deepcopy(dd0)
        actor = FUSE.ActorPedestal(dd, act.ActorPedestal, act; model, Ti_ped_model, fingerprint_mode, kw...)
        FUSE.finalize(FUSE.step(actor))
        s = actor.ped_actor.ti_split
        @test s.model == Ti_ped_model
        @test s.Te_ped > 80.0 && s.Ti_ped > 0.0
        @test all(isfinite, dd.core_profiles.profiles_1d[].t_i_average)
        if Ti_ped_model == :ratio
            @test s.Ti_ped ≈ s.Te_ped
        elseif Ti_ped_model == :fingerprint
            @test s.chi_e > 0.0 && s.chi_i_neo > 0.0
            ratio[(model, fingerprint_mode)] = s.Ti_ped / s.Te_ped
        end
    end
    # ions closer to neoclassical (MTM) hold a higher Ti_ped than MHD-like transport, for the same pedestal pressure
    @test ratio[(:EPED, :MTM)] > ratio[(:EPED, :MHD)]
    @test ratio[(:WPED, :MTM)] > ratio[(:WPED, :MHD)]

    # Ti_ped consistent with collisional exchange and χi_neo recomputed on the resulting profiles
    for fingerprint_mode in (:MHD, :MTM)
        dd = deepcopy(dd0)
        actor = FUSE.ActorPedestal(dd, act.ActorPedestal, act; model=:EPED, Ti_ped_model=:fingerprint, fingerprint_mode, Ti_ped_self_consistent=true, kw...)
        FUSE.finalize(FUSE.step(actor))
        s = actor.ped_actor.ti_split
        @test length(s.Ti_ped_history) > 2
        @test s.Ti_ped ≈ s.Ti_ped_history[end] atol = 5.0
        # exchange pulls the pedestal towards Ti = Te
        @test abs(log(s.Ti_ped / s.Te_ped)) < abs(log(ratio[(:EPED, fingerprint_mode)]))
        @test dd.core_profiles.profiles_1d[].t_i_average[end] ≈ 200.0 rtol = 0.2
    end
end
