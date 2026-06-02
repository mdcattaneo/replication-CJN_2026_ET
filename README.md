# Pairwise Estimator Simulation

Simulation code for bootstrap confidence intervals for pairwise difference estimators
Both simulations study four CI variants that combine two bandwidth choices (B = 1 and B = 3^(1/d)) with two bias correction levels (L = 0 no bias correction, L = 1 jackknife).

---

## Partially linear regression (PLR)

### Core functions

| File | Description |
|------|-------------|
| `plr.cpp` | Rcpp source compiled by the R scripts. Exports: `generate_data_plr` (simulate from DGP models 1-3), `estimate_plr` (closed-form pairwise-difference estimator at a single bandwidth), and `estimate_plr_grid` (efficient multi-bandwidth version using a single O(n²) pair scan). |
| `function_plr_boot.R` | Helper functions sourced by `boot_plr.R`. `bootstrap_once` runs a full bootstrap for one dataset: it evaluates the estimator on the union bandwidth grid, stores quantiles for both CI variants (B = 1 and B = 3^(1/d)), and returns point estimates and bootstrap quantiles. |

### Simulation script

| File | Description |
|------|-------------|
| `boot_plr.R` | Bootstrap CI simulation. Same distributed design as `boot_logit.R`: each SLURM node writes a partial RDS file. Arguments: `n`, `B_boot`, `R`, `model`, `node_id`, `n_nodes`, `outdir`, `jk_c`. |

### Post-processing scripts

| File | Description |
|------|-------------|
| `combine_plr.R` | Reads all partial RDS files from `boot_plr.R`, validates metadata consistency, combines arrays, computes empirical coverage and mean CI length, and saves a compiled RDS. |
| `mse_plr.R` | Standalone MSE simulation (runs all three models sequentially, no distributed setup). Reports bias, SD, RMSE, and MSE for the plain and jackknife estimators across a bandwidth grid. Saves one RDS per model. |
| `make_table_plr.R` | Reads the compiled RDS from `combine_plr.R` and writes a LaTeX `tabular` block (plain text) with coverage and CI-length columns for all four CI variants. |
| `plot_coverage_plr.R` | Reads the compiled RDS from `combine_plr.R` and plots empirical coverage probability vs. bandwidth for all four CI variants. Saves a PDF. |

---

## Pairwise logit

### Core functions

| File | Description |
|------|-------------|
| `logit.cpp` | Rcpp source compiled by the R scripts. Exports: `generate_data` (simulate from DGP models 1-3), `precompute_active_pairs` (collect kernel-weighted discordant pairs), `fast_wlogit` (Newton-Raphson weighted logistic regression, no intercept), and `gradient_fnc` / `objective_fnc` (gradient and value of the pairwise logistic objective). |
| `function_logit_boot.R` | Helper functions sourced by `boot_logit.R`. `thetahat` computes all four estimators; `boot_result` assembles the bootstrap CIs. |
| `function_logit_mse.R` | Reference `glm()`-based implementations of the L=0 (`estL0`) and jackknife-debiased (`estL1`) estimators. Used by the MSE simulation scripts. |

### Simulation script

| File | Description |
|------|-------------|
| `boot_logit.R` | Bootstrap CI simulation. Each SLURM node handles a disjoint subset of replications and writes a partial RDS file. Arguments: `n`, `R`, `B_boot`, `c1`, `model`, `node_id`, `n_nodes`, `outdir`. |

### Post-processing scripts

| File | Description |
|------|-------------|
| `combine_logit.R` | Reads all partial RDS files from `boot_logit.R`, stacks them, computes empirical coverage and mean CI length for each bandwidth and CI variant, and saves a compiled RDS. Prints a summary table to stdout. |
| `mse_combine_logit.R` | Aggregator for MSE simulation output. Reports bias, SD, RMSE, and MSE per bandwidth. |
| `make_table_logit.R` | Reads the compiled RDS from `combine_logit.R` and writes a LaTeX `tabular` block (plain text) with coverage and CI-length columns for all four CI variants. |
| `plot_coverage_logit.R` | Reads the compiled RDS from `combine_logit.R` and plots empirical coverage probability vs. bandwidth for all four CI variants. Saves a PDF. |
