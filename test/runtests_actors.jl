using FUSE
using Test

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
         
