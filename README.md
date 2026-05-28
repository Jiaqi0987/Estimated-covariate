This github repo contains the code for the pilot-fixed and CV-only simulation runs, plus plotting scripts for the generated results.

The core estimator/function files are copied from the working simulation folder and should be treated as source files:

- `estimator_functions.R`: estimator, basis, nuisance, and helper functions.
- `simulation_functions.R`: simulation data generation and per-replication estimator execution.
- `plot_helpers.R`: plotting helper functions used by the main-plot scripts.

The run script sources those files directly. Generated results and plots are intentionally ignored by Git.

## Requirements

R packages used by the scripts include:

```r
install.packages(c(
  "nprobust", "splines2", "gbm", "caret", "glmnet",
  "mvtnorm", "mgcv", "future.apply", "ranger",
  "ggplot2", "dplyr", "tidyr", "patchwork"
))
```

Some clusters manage R packages through modules or a shared library; use the local cluster setup if that is preferred.

## Pilot-Fixed Simulation

This is the main two-stage workflow:

1. Run `nsim_pilot` pilot replications with CV-selected `df/h`.
2. Take the median selected `df/h` separately for each estimator and basis.
3. Run `nsim` final replications with those fixed choices.

For each setting, the pilot result and final fixed-tuning result are saved together in one `fixed_*.rds` file.

Run one combination manually:

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

Run all pilot-fixed jobs on Slurm:

```bash
sbatch submit_simulation.sbatch
```

The Slurm script runs:

- `setting = 1, 2, 3`
- `g_type = null, linear`
- `n.tr = 2500, 10000`
- `n = 2500`

## CV-Only Simulation

This runs 500 simulation replications without fixing `df/h`; each replication uses CV tuning.

Run one combination manually:

```bash
Rscript run_simulation.R \
  --mode cv_only \
  --setting 1 \
  --g_type null \
  --n 2500 \
  --n_tr 2500 \
  --nsim 500 \
  --outdir results_cv_ntr2500
```

Run all CV-only jobs for `n.tr = 2500` on Slurm:

```bash
MODE=cv_only sbatch --array=1-6 submit_simulation.sbatch
```

## Plotting

After CV-only `.rds` files exist in `results_cv_ntr2500/` and pilot-fixed `.rds` files exist in `results/`, run:

```bash
Rscript generate_plots.R
```

Outputs:

- `plots/main_plots/`: main plots from CV-only results.
- `plots/smoothing_target_coverage/`: smoothing-target coverage from pilot-fixed results.
- `plots/ci_length_ntr2500/`: CI length from pilot-fixed results with `n.tr = 2500`.

Smoothing-target coverage uses `fixed_tuning$fixed_value` from each final result file, so the smoothing target is fixed across replications for each estimator/basis.

For the CV-only main plots, the upper-right pointwise-error panels use a 99th-percentile y-axis cap for readability. The lower MSE boxplot panel is not clipped and uses a shared y-axis across bases.

## Expected Output Directories

The scripts create these directories as needed:

- `results/`
- `results_cv_ntr2500/`
- `logs/`
- `plots/`

These generated directories are ignored by Git.
