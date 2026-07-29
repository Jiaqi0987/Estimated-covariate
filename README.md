# Simulation Code for Regression With Estimated Covariates

The code studies nonparametric regression when the covariate used in the second-stage regression is itself estimated in a first stage. 

## Methods Implemented

The simulation compares the following estimator families:

- `OR`: oracle benchmark that regresses on the true propensity score `r(X)`.
- `PI`: plug-in estimator that regresses `Y` on the estimated propensity score `rhat(X)`.
- `CPI`: corrected plug-in estimator that subtracts an estimated first-stage bias correction before the second-stage regression.
- `BC`: influence-function-based bias-corrected estimator.
- `PI_cali`: plug-in estimator after isotonic calibration of `rhat(X)`.
- `BC_cali`: bias-corrected estimator after isotonic calibration of `rhat(X)`.
- `CPI_or`, `BC_or`, and `BC_or_cali`: oracle variants used for benchmarking pieces of the correction.

Each estimator is evaluated with multiple second-stage regression bases:

- Fourier/cosine basis
- Polynomial basis
- B-spline basis
- Local-linear smoothing via `nprobust::lprobust`

## Requirements

The code requires R and the following R packages:

```r
install.packages(c(
  "nprobust",
  "splines2",
  "gbm",
  "caret",
  "glmnet",
  "mvtnorm",
  "mgcv",
  "future.apply",
  "ranger",
  "ggplot2",
  "dplyr",
  "tidyr",
  "patchwork"
))
```

## CV-Only Simulation Workflow

The CV-only workflow runs each replication with cross-validated tuning. For local-polynomial estimators, `--bandwidth` controls how the bandwidth is chosen:

- `cv`: each estimator uses its own default/CV bandwidth.
- `oracle`: non-oracle local-polynomial estimators use the bandwidth selected by the oracle estimator in that replication.
- `fixed`: all local-polynomial estimators use user-specified fixed bandwidths.
- `all`: runs `cv`, `oracle`, and each fixed bandwidth in one invocation.

Example:

```bash
Rscript run_simulation.R \
  --mode cv_only \
  --setting 1 \
  --g_type null \
  --n 2500 \
  --n_tr 2500 \
  --nsim 500 \
  --bandwidth all \
  --h_fixed_values 0.10,0.15,0.20 \
  --outdir results_cv_ntr2500
```

## Pilot-Fixed Simulation Workflow

The pilot-fixed workflow is used for coverage and confidence-interval comparisons with fixed tuning parameters:

1. Run `nsim_pilot` pilot replications with cross-validated `df` or `h`.
2. Take the median selected tuning value separately for each estimator and basis.
3. Run `nsim` final replications with those estimator-specific fixed tuning choices.

Example:

```bash
Rscript run_simulation.R \
  --mode pilot_fixed \
  --setting 1 \
  --g_type null \
  --n 2500 \
  --n_tr 2500 \
  --nsim_pilot 50 \
  --nsim 500 \
  --outdir results
```

## Running on Slurm

The included Slurm script runs simulation arrays on a cluster.

Pilot-fixed runs:

```bash
sbatch submit_simulation.sbatch
```

CV-only runs:

```bash
MODE=cv_only BANDWIDTH=all H_FIXED_VALUES=0.10,0.15,0.20 \
  sbatch --array=1-6 submit_simulation.sbatch
```

The script uses `SLURM_CPUS_PER_TASK` to set the number of parallel workers for `parallel::mclapply`.

## Plotting

After result files exist in run:

```bash
Rscript generate_plots.R
```

For pilot-fixed result files, the coverage plot uses the fixed tuning table saved in the result file and targets the corresponding smoothed estimand. For CV-only results, no fixed tuning table is available, so the coverage plot falls back to realized CV-selected `df/h` values from each replication.
