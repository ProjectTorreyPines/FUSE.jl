# ============================================================================
# FUSE ITER Digital Twin - Database Generator
# ============================================================================
# Parameter sweep for the ITER digital twin database.
# Total nominal combinations: 311,040 runs
# Expected viable (passes quasineutrality): ~180k-200k
# Expected converged: ~100k-150k
#
# Uses a lazy generator function so that ini/act objects are created on-the-fly
# per worker, instead of pre-allocating 311k deepcopy'd objects in memory.
#
# Run a small local test first (LOCAL_TEST = true), then switch to NERSC.
# ============================================================================

using Pkg
Pkg.activate("/pscratch/sd/s/slendebt/.julia/dev/FUSE_Design")

using FUSE # This is FUSE_Design

LOCAL_TEST = false  # Set to false for full NERSC submission
FOLDER_NAME = "ITER_digital_twin_scan_v1.1"

# ── Step 1: Initialize hardware from ODS (done once, reused for all runs) ────
ini, act = FUSE.case_parameters(:ITER; init_from=:ods)
dd = IMAS.dd()
@time FUSE.init(dd, ini, act)
@checkin :hw_init dd ini act

# ── Step 2: Base ini/act from scalars ────────────────────────────────────────
ini_base, _ = FUSE.case_parameters(:ITER;
    init_from=:scalars,
    boundary_from=:scalars,
    time_dependent=false
)

# Fixed settings for all runs
ini_base.core_profiles.ne_setting = :greenwald_fraction

# Base act settings
act_base = deepcopy(act)
act_base.ActorPedestal.density_match = :ne_line
act_base.ActorStationaryPlasma.max_iterations = 2

# ── Step 3: Define parameter grid ───────────────────────────────────────────
# Continuous parameters
ip_values       = [10e6, 12e6, 15e6, 18e6]         # [A] Plasma current
fgw_values      = [0.4, 0.6, 0.8, 1.0]              # Greenwald fraction
zeff_values     = [1.5, 2.0, 3.0]                    # Effective charge
ne_sep_ped      = [0.15, 0.25, 0.40]                 # ne separatrix-to-pedestal ratio

# Heating power levels (mapped to ini launcher units)
# NBI: 2 units, each 16.7 MW max => total 0, 16.7, 25.05, 33.4 MW
nbi_configs = [
    (0.0,    0.0),      # 0 MW total
    (16.7e6, 0.0),      # 16.7 MW total (1 beam)
    (16.7e6, 8.35e6),   # 25 MW total
    (16.7e6, 16.7e6),   # 33.4 MW total (2 beams)
]
ec_values  = [0.0, 10e6, 20e6]   # [W] EC power (1 launcher)
ic_values  = [0.0, 12e6, 24e6]   # [W] IC power (1 antenna)

# Discrete parameters
impurity_species = [:Ne, :Ar, :Xe]
he_fractions     = [0.0, 0.01, 0.02, 0.05, 0.10]
w_densities      = [0.0, 5e15, 1.5e16, 5e16]         # [m^-3]
bulk_species     = [:DT, :D]
plasma_modes     = [:H_mode, :L_mode]

# ── Step 4: Build grid and generator function ────────────────────────────────
# Use Iterators.product to define the full parameter grid as a compact array of
# tuples. This is only ~311k tuples of scalars/symbols (~few MB), NOT 311k
# deepcopy'd ini/act objects (~300GB).
grid = collect(Iterators.product(
    bulk_species, plasma_modes, ip_values, fgw_values,
    nbi_configs, ec_values, ic_values, zeff_values,
    impurity_species, he_fractions, w_densities, ne_sep_ped
))

@info "Grid has $(length(grid)) parameter combinations"

# The generator function creates one (ini, act) pair at a time.
# It captures only ini_base, act_base, and grid (all small).
# Each pmap worker receives this closure once and calls it per case.
function make_ini_act(item::Int)
    (bulk, mode, ip, fgw, (nbi1, nbi2), ec, ic, zeff, imp, he, w, sep_ped) = grid[item]

    ini_i = deepcopy(ini_base)
    act_i = deepcopy(act_base)

    # Plasma parameters
    ini_i.equilibrium.ip = ip
    ini_i.core_profiles.ne_value = fgw
    ini_i.core_profiles.zeff = zeff
    ini_i.core_profiles.ne_sep_to_ped_ratio = sep_ped

    # Heating
    ini_i.nb_unit[1].power_launched = nbi1
    ini_i.nb_unit[2].power_launched = nbi2
    ini_i.ec_launcher[1].power_launched = ec
    ini_i.ic_antenna[1].power_launched = ic

    # Species
    ini_i.core_profiles.impurity = imp
    ini_i.core_profiles.helium_fraction = he
    ini_i.core_profiles.tungsten_density = w
    ini_i.core_profiles.bulk = bulk
    ini_i.core_profiles.plasma_mode = mode

    return ini_i, act_i
end

# ── Step 5: Configure the study ──────────────────────────────────────────────
sty = FUSE.study_parameters(:DatabaseGenerator)

sty.server = "localhost"
sty.n_workers = 110
sty.save_folder = "/pscratch/sd/s/slendebt/FUSE/ITERSimulations/$(FOLDER_NAME)"

sty.file_save_mode = :safe_write
sty.n_simulations = length(grid)

# ── Step 6: Create the study with a lazy generator ───────────────────────────
# This passes only the make_ini_act function, NOT 311k pre-allocated objects.
# Memory on the manager: ~few MB (grid + base ini/act)
# Memory per worker: one ini/act pair at a time
study = FUSE.StudyDatabaseGenerator(sty, make_ini_act)

using Distributed
@everywhere import FUSE
@everywhere ProgressMeter = FUSE.ProgressMeter

@everywhere function workflow_DatabaseGenerator(
    dd::FUSE.IMAS.dd,
    ini::FUSE.ParametersAllInits,
    act::FUSE.ParametersAllActors
)
    # Set pedestal model based on plasma mode
    if ini.core_profiles.plasma_mode == :L_mode
        act.ActorPedestal.model = :WPED
        act.ActorWPED.ped_to_core_fraction = 0.2
    else
        act.ActorPedestal.model = :EPED
        act.ActorEPED.ped_factor = 1.0
    end

    # Initialize profiles (hardware already loaded from @checkin)
    FUSE.init(dd, ini, act; initialize_hardware=false)

    # Run the stationary plasma solver for self-consistent solution
    FUSE.ActorStationaryPlasma(dd, act)

    return nothing
end

study.workflow = workflow_DatabaseGenerator

# ── Step 7: Run ──────────────────────────────────────────────────────────────
@info "Starting ITER Digital Twin database generation..."
@info "  Runs: $(sty.n_simulations)"
@info "  Workers: $(sty.n_workers)"
@info "  Save folder: $(sty.save_folder)"

FUSE.run(study)

@info "Database generation complete!"
