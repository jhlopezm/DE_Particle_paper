# Reproducible numerical protocol

## Forward problem

The ten cubic-potential parameters are

`x = (lambda_1, lambda_2, lambda_3, lambda_4, lambda_sigma, lambda_1sigma, lambda_2sigma, v_1, v_sigma, mu)`.

Set `v = 246.22 GeV`, `v_2 = sqrt(v^2-v_1^2)`. The real symmetric 3×3 CP-even mass-squared matrix is implemented explicitly in `src/model.py`. For ordered nonnegative eigenvalues `m_i²` and eigenvector matrix `R`, `R.T @ M2 @ R = diag(m_i²)`, `phi = R @ H`, and `H = R.T @ phi`. The H2 eigenvector is oriented such that `kappa_2V >= 0`.

## Sampling and optimization

1000 independent synthetic target triplets were generated with NumPy seed `20261020` from normals centered at (95.40, 125.09, 151.50) GeV with **numerical scales** `sigma=(0.05, 0.09, 0.09)` GeV. Three DE restarts per target: base seeds (11,23,37), per-target seed `base_seed + 100003*target_index`.

`scipy.optimize.differential_evolution`: `best1bin`, `popsize=50`, `maxiter=300`, `tol=0.01`, `mutation=(0.5,1.0)`, `recombination=0.7`, `polish=True`. Original optimizer constraints: `|m_charged² - (130 GeV)²| <= 260 GeV²` and `0.98 <= kappa_2V <= 1.08`. These are **not identical** to the tighter final cuts.

## Final selection (post-processing)

For each run, define

- `z_i = |m_i(x)-m_i_star|/sigma_i` for i=1,2,3,
- `z_plus = |m_charged(x)-130 GeV|/(1 GeV)`,
- `z_kappa = |kappa_2V(x)-1.03|/0.05`.

Accept if **all five** satisfy `z < 1` (strict inequality). Numerically, 2977 satisfy all three mass criteria; 2994 the charged scalar criterion; 3000 the kappa criterion; 2971 satisfy all five. Every target has at least two accepted restarts (971 targets have three, 29 have two).

## Interpretation

A restart is a reconstruction of a sampled target, not a new experimental observation. Empirical reconstruction distributions are optimizer- and bounds-dependent; they must not be represented as Bayesian posterior densities. The DE runs reached the configured `maxiter`; SciPy's formal `success=False` remains in the raw data even when reconstruction residuals satisfy all final acceptance tests. Correlations are measured separately from analytic approximations; singlet-dominance approximations are checked against eigenvector composition.

The frozen original run is a record of the procedure executed with the **original optimizer constraints**; re-running DE with final 5×1σ conditions imposed internally would define a *different* numerical experiment.
