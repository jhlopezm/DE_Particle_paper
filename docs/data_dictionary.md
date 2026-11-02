# Data dictionary

The authoritative, immutable inputs are `data/raw/runs.csv`, `data/raw/targets.csv`, and `data/raw/production_config.json`.
Each row of `runs.csv` is one optimization, identified by the unique pair `(target_index, restart)`; each target has three restarts, with base seeds 11, 23 and 37 and per-target seed `base_seed + 100003 * target_index`.

- `target_index`: zero-based target index, 0–999.
- `restart`: zero-based restart, 0–2; `seed`: effective optimizer seed.
- `lambda_1`, `lambda_2`, `lambda_3`, `lambda_4`, `lambda_sigma`, `lambda_1sigma`, `lambda_2sigma`: dimensionless potential couplings.
- `v_1`, `v_sigma`, `mu`: dimensionful parameters in GeV.
- `target_H1`, `target_H2`, `target_H3`: sampled target masses in GeV.
- `mass_H1`, `mass_H2`, `mass_H3`: predicted CP-even masses in GeV (ordered increasingly).
- `mc`: charged-scalar mass in GeV; `kappa2V`: reduced vector coupling (nonnegative eigenvector sign convention).
- `z_H1`, `z_H2`, `z_H3`: *signed* mass residuals `(prediction - target)/sigma_i`.
- `z_charged`: magnitude of charged mass residual divided by 1 GeV.
- `z_kappa`: magnitude of coupling residual divided by 0.05.
- `loss`: normalized squared mass reconstruction loss minimized by DE.
- `optimizer_success`: convergence status from SciPy (distinct from the selection criteria).
- `nit`, `nfev`: iteration and function-evaluation counts.
- `mass_pass`, `aux_pass`, `all_pass`: original run-time diagnostics (historical flags). **Do not use them for the final 5×1σ selection**; recalculate all five residuals strictly with absolute values below one.
- `minimum_eigenvalue`: smallest eigenvalue of the CP-even mass-squared matrix, in GeV².

`data/processed/accepted_2971.csv` and `data/processed/best_per_target_1000.csv` are *derived* from the original 3000-row CSV via final 5×1σ cuts. The representative solution is the **lowest-loss accepted** restart for each target; it is not a posterior sample.

These numerical σ widths are chosen reconstruction/selection scales; they are **not claimed to be independent experimental standard deviations**. All masses, VEVs and μ are in GeV unless indicated otherwise.
