using FUSE
using Test

# @testset "fluxmatcher" begin
#     ini, act = FUSE.case_parameters(:ITER; init_from=:scalars)
#     dd = IMAS.dd()
#     FUSE.init(dd, ini, act)
#     act.ActorFluxMatcher.max_iterations = 2
#     act.ActorFluxMatcher.evolve_pedestal = true

#     act.ActorFluxMatcher.algorithm = :simple
#     act.ActorTGLF.model = :TJLF
#     FUSE.ActorFluxMatcher(dd, act)
#     act.ActorFluxMatcher.algorithm = :anderson
#     act.ActorTGLF.model = :TGLFNN
#     FUSE.ActorFluxMatcher(dd, act)
# end

@testset "Mars Actor CHEASE run" begin
    ini, act = FUSE.case_parameters(:D3D, :default)
    dd = IMAS.dd()
    FUSE.init(dd, ini, act)

    #control execution of CHEASE & MARS for testing purposes
    act.ActorMars.run_equilibrium = true
    act.ActorMars.restart_equilibrium = false
    act.ActorMars.run_MHD = false
    act.ActorMars.run_mode = :local
    act.ActorMars.chease_exec = "/Users/akcay/Codes/MarsQ_package/CheaseMerge/chease.x"  # <-- Update this path to your CHEASE executable

    # Configure CHEASE parameters for testing
    act.ActorMars.number_surfaces = 1
    chease_overrides = (NPSI=64, NVEXP=2, NCSCAL=4)

    mktempdir() do tempdir
        @info "Running CHEASE & MARS test in temporary directory: $tempdir"

        cd(tempdir) do  
            if !isfile(act.ActorMars.chease_exec)
                @test_skip "CHEASE executable not found at $(act.ActorMars.chease_exec). Skipping CHEASE & MARS test."
            else
                FUSE.ActorMars(dd, act; chease_overrides=chease_overrides)
                @testset "CHEASE output files" begin
                    @test isfile("EXPEQ")|| error("CHEASE equilibrium file EXPEQ not found.")
                    #@test filesize("EXPEQ") > 0 "CHEASE equilibrium file EXPEQ is empty."
                    @test isfile("datain") || error("CHEASE datain file not found.")
                end

                # next test the restart capability
                @info "Testing CHEASE restart capability"
                act.ActorMars.restart_equilibrium = true
                FUSE.ActorMars(dd, act; chease_overrides=chease_overrides)
            end
        end
    end
end

@testset "MARS actor full CHEASE-MARS MHD stability run" begin
    ini, act = FUSE.case_parameters(:D3D, :default)
    dd = IMAS.dd()
    FUSE.init(dd, ini, act)

    #control execution of CHEASE & MARS for testing purposes
    act.ActorMars.run_equilibrium = true
    act.ActorMars.restart_equilibrium = false
    act.ActorMars.run_MHD = true
    act.ActorMars.run_mode = :local
    act.ActorMars.chease_exec = "/Users/akcay/Codes/MarsQ_package/CheaseMerge/chease.x"  # <-- Update this path to your CHEASE executable

    # Configure CHEASE & MARS parameters for testing
    chease_overrides = (NPSI=64, NVEXP=2, NCSCAL=4)
    mars_overrides=(BASIC=(NWALL=0, M1=-10, NPROFN=0, 
                    NV=160, TALPHA1=(0.05, 0.)),)

    mktempdir() do tempdir
        @info "Running MARS test in temporary directory: $tempdir"

        cd(tempdir) do  
            if !isfile(act.ActorMars.chease_exec)
                @test_skip "CHEASE executable not found at $(act.ActorMars.chease_exec). Skipping MARS test."
            else
                if !isfile("RUN.IN")
                    @test_skip "MARS input file MARS_INPUT not found after CHEASE run. Skipping MARS test."
                else
                    FUSE.ActorMars(dd, act; chease_overrides=chease_overrides)
                    @testset "MARS output files" begin
                        @test isfile("log_mars") || error("MARS output file MARS_OUTPUT not found.")
                        #@test filesize("MARS_OUTPUT") > 0 "MARS output file MARS_OUTPUT is empty."
                    end
                end
            end
        end
    end
end
