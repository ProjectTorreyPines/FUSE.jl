# CLAUDE.md - Development Notes

## Changes Made (2026-05-20)

### 1. Ohmic Power Clamp (IMAS_Design)
**File:** `IMAS_Design/src/physics/sources.jl`

Clamped `powerDensityOhm` to non-negative values in `ohmic_source!()`. At very low densities,
the Spitzer conductivity breaks down and bootstrap current can overshoot `ip`, causing `j_ohmic`
to flip sign. The product `j_tor * j_ohmic` then goes negative, producing unphysical negative
ohmic heating. The clamp prevents this while preserving valid ohmic contributions.

### 2. Allow `density_match = :fixed` in Pedestal Initialization
**File:** `src/ddinit/init_pulse_schedule.jl`

Relaxed the assertion that enforced strict pairing between `ini.core_profiles.ne_setting` and
`act.ActorPedestal.density_match`. The `:fixed` option now bypasses the check, since the pedestal
actor doesn't read density from the pulse schedule when using `:fixed`. The mismatch protection
between `:ne_ped` and `:ne_line` remains active.

### 3. EPED-NN Extrapolation Distance (EPEDNN + FUSE_Design)
**Files:**
- `EPEDNN/src/EPEDNN.jl` - Added `extrapolation_distance()` function
- `src/actors/pedestal/EPED_actor.jl` - Store and warn on extrapolation
- `src/actors/pedestal/pedestal_actor.jl` - Write to `dd.summary`

Added a normalized extrapolation distance metric to detect when EPED-NN inputs are outside
the training range. For each input, distance is measured as a fraction of the training range
width beyond the bounds (0 = within bounds, 1.0 = one full range outside).

The `ActorEPED` now:
- Computes extrapolation distance after each EPED call
- Warns with per-input details when any input is out of bounds
- Stores `extrapolation` on the actor struct

The `ActorPedestal._finalize` writes `Te_ped * extrapolation_distance` to
`dd.summary.local.pedestal.t_e.value_σ`, allowing scan post-processing to filter
on `value_σ / value` as a reliability metric.

### 4. Explicit Pedestal Rotation Parameter (`rot_ped`)
**Files:**
- `src/parameters/parameters_inits.jl` - Added `rot_ped` entry
- `src/ddinit/init_core_profiles.jl` - Use `rot_ped` in profile construction
- All case files in `src/cases/` - Set `rot_ped` for each case

Previously the pedestal rotation was implicitly derived as `rot_core / ne_core_to_ped_ratio`,
borrowing a density parameter for rotation shaping. Now `rot_ped` is an explicit required
parameter. If not set, the init errors with a message suggesting the old default value.

Case file settings:
- Zero rotation cases: `rot_ped = 0.0`
- MASTU, HDB5, D3D: `rot_ped = rot_core / 1.4` (preserves old behavior)
- ITER: `rot_ped = 4e3` rad/s (intrinsic rotation from dimensionless empirical scaling,
  Chrystal et al., Nucl. Fusion 60 (2020) 036003, https://doi.org/10.1088/1741-4326/ab5c86).
  Core rotation scales linearly with NBI power fraction.

## Architecture Notes

### Rotation Profile
- `rotation_frequency_tor_sonic` is the primary rotation variable set at init and used by transport
- `ion.rotation_frequency_tor` is a derived expression via `sonic2ωtor()` that subtracts
  the diamagnetic correction `dp_i/dψ / (n_i * Z_i * e)` — this can produce visual artifacts
  near rho=0 due to numerical gradient issues but does not affect physics (transport starts at rho>=0.2)
- With `evolve_rotation = :fixed` (default), the sonic rotation profile is never modified by the flux matcher

### EPED-NN Behavior at Low Density
- The power law component scales as `10^(c0 + c1*log10(neped) + ...) ^ 2 * neped`
- The NN correction is additive on top of the power law (before squaring)
- Both extrapolate badly outside training range — power law diverges, NN produces arbitrary corrections
- The `only_powerlaw = true` default avoids the NN but does not fix the power law extrapolation
- Use `extrapolation_distance` to detect and filter these cases in scans
