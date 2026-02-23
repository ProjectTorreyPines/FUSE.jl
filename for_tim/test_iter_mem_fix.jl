# Test script for the generator-based DatabaseGenerator (memory fix)

using FUSE

# ── Step 1: Hardware init from ODS ───────────────────────────────────────────
ini, act = FUSE.case_parameters(:ITER; init_from=:ods)
dd = IMAS.dd()
@time FUSE.init(dd, ini, act)
@checkin :hw_init dd ini act

# ── Step 2: Base ini/act from scalars ────────────────────────────────────────
ini_base, _ = FUSE.case_parameters(:ITER; init_from=:scalars, boundary_from=:scalars, time_dependent=false)

ini_base.core_profiles.ne_setting = :greenwald_fraction

act_base = deepcopy(act)
act_base.ActorPedestal.density_match = :ne_line
act_base.ActorStationaryPlasma.max_iterations = 2

# ── Step 3: Define parameter grid ───────────────────────────────────────────
ne_values      = [0.2, 0.4, 0.6, 0.8, 1.0]
impurity_species = [:Kr, :Ne, :Xe]
plasma_modes   = [:H_mode, :L_mode]
zeff_values    = [1.1, 2.0, 3.0, 4.0]

grid = collect(Iterators.product(ne_values, impurity_species, plasma_modes, zeff_values))
@info "Grid has $(length(grid)) parameter combinations"

# ── Step 4: Generator function ───────────────────────────────────────────────
function make_ini_act(item::Int)
    (ne, imp, mode, zeff) = grid[item]

    ini_i = deepcopy(ini_base)
    act_i = deepcopy(act_base)

    ini_i.core_profiles.ne_value = ne
    ini_i.core_profiles.impurity = imp
    ini_i.core_profiles.plasma_mode = mode
    ini_i.core_profiles.zeff = zeff

    return ini_i, act_i
end

# ── Step 5: Configure and create study ───────────────────────────────────────
sty = FUSE.study_parameters(:DatabaseGenerator)
sty.server = "localhost"
sty.n_workers = 2
sty.file_save_mode = :safe_write
sty.save_folder = "/Users/tims/source_codes/run_FUSEDesign/test_ITER_after_mem_fix"
sty.n_simulations = length(grid)

study = FUSE.StudyDatabaseGenerator(sty, make_ini_act)

# ── Step 6: Define workflow ──────────────────────────────────────────────────
using Distributed
@everywhere import FUSE
@everywhere ProgressMeter = FUSE.ProgressMeter

@everywhere function workflow_DatabaseGenerator(dd::FUSE.IMAS.dd, ini::FUSE.ParametersAllInits, act::FUSE.ParametersAllActors)
    if ini.core_profiles.plasma_mode == :L_mode
        act.ActorPedestal.model = :WPED
        act.ActorWPED.ped_to_core_fraction = 0.2
    else
        act.ActorPedestal.model = :EPED
        act.ActorEPED.ped_factor = 1.0
    end
    FUSE.init(dd, ini, act; initialize_hardware=false)
    #FUSE.ActorStationaryPlasma(dd, act)
    return nothing
end

study.workflow = workflow_DatabaseGenerator

# ── Step 7: Run ──────────────────────────────────────────────────────────────
@info "Starting test run..."
@info "  Runs: $(sty.n_simulations)"
@info "  Workers: $(sty.n_workers)"
@info "  Save folder: $(sty.save_folder)"

FUSE.run(study)

@info "Test complete!"
