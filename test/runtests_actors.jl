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

@testset "Martian CHEASE" begin
    ini, act = FUSE.case_parameters(:D3D; :default)
    dd = IMAS.dd()
    FUSE.init(dd, ini, act)

    # control execution of CHEASE & MARS for testing purposes
    act.ActorMARS.run_equilibrium = true
    act.ActorMARS.restart_equilibrium = false
    act.ActorMARS.run_MHD = false
    act.ActorMARS.run_tracer = false
    act.ActorMARS.run_mode = :local

    # Configure CHEASE parameters for testing
    act.ActorMARS.number_surfaces = 1
    act.ActorMARS.chease_overrides = (NPSI=64, NVEXP=2, NCSCAL=4,)

    mktempdir() do tempdir
        @info "Running CHEASE & MARS test in temporary directory: $tempdir"

        cd(tempdir) do  
            @test isfile("datain") "datain file not found in current directory"
            content = read("datain", String)
            @test occursin("number_of_surfaces = 1", content) "datain file does not contain expected case"  
            
            if !isfile(act.ActorMARS.chease_exec)
                @warn "CHEASE executable not found at $(act.ActorMARS.chease_exec). Skipping CHEASE & MARS test."
            else
                FUSE.ActorMARS(dd, act)
            end
        end
    end
end
