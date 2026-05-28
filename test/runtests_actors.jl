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

    #control execution of CHEASE & MARS for testing purposes
    act.ActorMars.run_equilibrium = true
    act.ActorMars.run_MHD = false
    act.ActorMars.wall_type = :no_wall # Begin with NO wall
    #act.ActorMars.chease_exec = "/Users/akcay/Codes/MarsQ_package/CheaseMerge/chease.x"
    #act.ActorMars.mars_exec = "/Users/akcay/Codes/MarsQ_package/MarsQ/marsq.x"
    
    # Configure CHEASE parameters for testing
    chease_overrides = (NPSI=64, NVEXP=2, NCSCAL=4)

    # configure MARS parameters for testing
    mars_overrides = FUSE.MarsOverrides()
    mars_overrides.BASIC[:M1]=-10
    mars_overrides.BASIC[:NV]=120 # moves the IW in from the original CHEASE lcation
    

    mktempdir() do tempdir
        @info "Running CHEASE & MARS test in temporary directory: $tempdir"

        cd(tempdir) do  
            @info "=============================================================="
            @info "       Test 1: Clean CHEASE equilibrium run"
            @info "=============================================================="
            if !isfile(act.ActorMars.chease_exec)
                @test_skip "CHEASE executable not found at $(act.ActorMars.chease_exec). Skipping CHEASE & MARS test."
            else
                FUSE.ActorMars(dd, act; chease_overrides=chease_overrides)
                @testset "CHEASE output files" begin
                    @test isfile("EXPEQ")|| error("CHEASE equilibrium file EXPEQ not found.")
                    #@test filesize("EXPEQ") > 0 "CHEASE equilibrium file EXPEQ is empty."
                    @test isfile("datain") || error("CHEASE datain file not found.")
                end

                # next test the restart capability and RW set-up
                @info "=============================================================="
                @info "    Test 2: CHEASE capability with resistive wall & higher beta"
                @info "=============================================================="
                act.ActorMars.number_surfaces = 2
                act.ActorMars.wall_type = :limiter
                chease_overrides_RW = (CFBAL=1.1,)
                FUSE.ActorMars(dd, act; chease_overrides=chease_overrides_RW)

                @info "=============================================================="
                @info "    Test 3: CHEASE capability with resistive wall & higher beta"
                @info "=============================================================="
                act.ActorMars.restart_equilibrium = true
                chease_overrides_restart = (NCSCAL=2,) # keep Ip fixed
                FUSE.ActorMars(dd, act; chease_overrides=chease_overrides_restart)

                # test MARS MHD stability run
                @info "=============================================================="
                @info "       Test 4: Testing MARS MHD stability run"
                @info "=============================================================="
                act.ActorMars.run_equilibrium = false # do NOT rerun CHEASE, just run MARS on the existing equilibrium
                act.ActorMars.run_MHD = true

                FUSE.ActorMars(dd, act; mars_overrides=mars_overrides)
                @testset "MARS output files" begin
                    @test isfile("log_mars") || error("MARS output file MARS_OUTPUT not found.")
                    #@test filesize("MARS_OUTPUT") > 0 "MARS output file MARS_OUTPUT is empty."
                end
                

            end
        end
    end
end

         
