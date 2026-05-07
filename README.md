# Estimated Covariates Simulation Code

This folder contains a small, self-contained version of the simulation code for the estimated-covariates project. 

The simulation studies estimation of

```tex
m(t) = E\{Y \mid r(X)=t\},
```

where the index `r(X)` is unknown and must be estimated from data.

## File Structure

```text
R/
  estimators.R
  plot_functions.R
  run_simulation_and_plot.R
```

## What Each File Does

### `R/estimators.R`

This file contains the main estimator functions:

- plug-in estimator,
- corrected plug-in estimator,
- influence-function bias-corrected estimator,
- local-linear implementation,
- sieve/basis implementation,
- nuisance estimation helpers for `rhat` and `muhat`,
- standard-error calculations.

This is the file to inspect if you want to understand or modify the estimators.

### `R/plot_functions.R`

This file contains the plotting and result-summary functions:

- pointwise semioracle coverage,
- pointwise CI length,
- pointwise MSE relative to the true curve,
- CSV summaries,
- PDF plot generation.

This is the file to modify if you want to change colors, estimator labels, facets, or output format.

### `R/run_simulation_and_plot.R`

This is the only file users need to run directly. It does three things:

- defines the simulation design,
- runs one simulation configuration,
- optionally generates plots from saved result files.

It has three modes:

```text
--mode run   : run one simulation configuration
--mode plot  : plot all existing result files in R/results
--mode both  : run one configuration, then regenerate plots
```

## R Dependencies

Install the required R packages with:

```r
install.packages(c(
  "nprobust", "splines2", "gbm", "caret", "glmnet", "mvtnorm",
  "mgcv", "future.apply", "ranger", "ggplot2", "dplyr", "tidyr"
))
```

## Run A Small Example

From the repository root:

```bash
Rscript R/run_simulation_and_plot.R \
  --mode run \
  --setting 3 \
  --g_type null \
  --n 2500 \
  --n_tr 2500 \
  --nsim 5 \
  --outdir results
```

This saves an `.rds` result file under:

```text
R/results/
```

For a real simulation, increase `--nsim`, for example:

```bash
Rscript R/run_simulation_and_plot.R \
  --mode run \
  --setting 3 \
  --g_type null \
  --n 2500 \
  --n_tr 10000 \
  --nsim 500 \
  --outdir results
```

To use multiple local cores, set `SIM_CORES`:

```bash
SIM_CORES=4 Rscript R/run_simulation_and_plot.R \
  --mode run \
  --setting 1 \
  --g_type linear \
  --n 2500 \
  --n_tr 5000 \
  --nsim 500 \
  --outdir results
```

## Generate Plots

After one or more result files have been created, run:

```bash
Rscript R/run_simulation_and_plot.R --mode plot --outdir results
```

This creates:

```text
R/results/plots/coverage_pointwise.csv
R/results/plots/coverage_summary.csv
R/results/plots/ci_length_pointwise.csv
R/results/plots/ci_length_summary.csv
R/results/plots/mse_pointwise.csv
R/results/plots/mse_summary.csv
R/results/plots/*_coverage.pdf
R/results/plots/*_ci_length.pdf
R/results/plots/*_mse_truth.pdf
```

You can also run one configuration and immediately regenerate plots:

```bash
Rscript R/run_simulation_and_plot.R \
  --mode both \
  --setting 2 \
  --g_type linear \
  --n 2500 \
  --n_tr 2500 \
  --nsim 20 \
  --outdir results
```

