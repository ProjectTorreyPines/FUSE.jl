# locking_actor.jl commit message

EF phase, error_field Gauss, MDS+ robustness, frequency-unit C1, plot labeling

## EF phase propagation

- deg2rad(par.EF_phase) passed to simulate_one_case and solve_and_classify.
- par.EF_phase = -90.0 when overwrite_params = true (PoP2024 convention).

## error_field Gauss conversion

- error_field in ODEparams is now in Gauss (default 10 G).
- set_control_parameters! converts once: error_field *= 1e-4/b0 * rc/m_pol.
- rc and m_pol extracted before the if/else block.
- Removed stale EpsUp inverse-conversion line.

## MDS+ robustness

- Negative μ guard: abs(muSI) with warning when torque/rotation < 0.
- Mass density fallback: n_e/Z_i via quasi-neutrality when ion density_thermal
  is missing or zero (common in MDS+ loads).

## Frequency-unit C1 convention

- Control1_min/max are in kHz; grid conversion: C1_dim = f_kHz * 1e3 * t0.
- Control1_max auto-set: max(default, rot_at_rs_kHz) from dd.
- Removed spurious 2π from old Omega0_calc.
- C1_op in :evaluate_probability uses rot/(2π) * t0 (frequency, not angular).

## Physical parameter printout

- set_phys_params! logs all parameters: torque, rotation (kHz at core and q=2),
  μ_SI, inertia, mass density, R0, ψ₀, U0, μ, I, Δ', ΔW, l12, l21, l32,
  error_field (Gauss), EF_phase (degrees).

## Operating-point evaluation (:evaluate_probability)

- Runs set_up_ode_params! on dd to extract C1_op (rotation) and C2_op
  (error_field / stability_index / saturation_param).
- Prints P(locked) at operating point.
- plot_probability shows yellow star at (C2_op, C1_op).

## Plot labeling

- Shot number + time (ms) displayed on all plots via _shot_label(dd).
- Appends "(OW)" when overwrite_params = true.
- save_locking_plots: filenames include shot/time and method prefix
  (NN_probability, Conv5x5_probability, KDE5x5_probability).

## Conv/KDE probability (actor wrappers)

- conv_locking_probability(actor) and kde_locking_probability(actor) wrappers.
- save_prob_model/load_prob_model dispatch on model type.
- :calc_conv_prob, :calc_kde_prob tasks added.
- :evaluate_probability tries KDE → conv → NN fallback.
- conv_window_C1/C2 par entries (shared by conv and KDE).

## _finalize() remains commented out

- finalize(actor) call (line 113) still commented out due to unresolved
  errors writing to dd.mhd_linear and dd.limits.model fields.
