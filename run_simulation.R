#!/usr/bin/env Rscript

# ============================================================
# run_simulation.R
#
# Modes:
#   pilot_fixed:
#     1. Run a pilot simulation with CV tuning.
#     2. Take the median selected df/h separately for each estimator and basis.
#     3. Run the final simulation using estimator-specific fixed choices.
#   cv_only:
#     Run the final simulation with CV tuning in every replication.
#     Local-polynomial bandwidth can be chosen by --bandwidth:
#       cv     = each estimator uses its own default/CV bandwidth.
#       oracle = all non-oracle local-polynomial estimators use the OR bandwidth.
#       fixed  = all local-polynomial estimators use --h_fixed.
#       all    = run cv, oracle, and fixed in one invocation.
#
# Example:
#   Rscript run_simulation.R --mode pilot_fixed \
#     --setting 2 --g_type linear --n 2500 --n_tr 2500 \
#     --nsim_pilot 50 --nsim 500 --outdir results
#   Rscript run_simulation.R --mode cv_only \
#     --setting 2 --g_type linear --n 2500 --n_tr 2500 \
#     --nsim 500 --outdir results_cv_ntr2500
# ============================================================

# -------------------------------
# Helper: parse command-line args
# -------------------------------
parse_args <- function(x) {
  out <- list()
  i <- 1
  while (i <= length(x)) {
    key <- sub("^--", "", x[i])
    val <- x[i + 1]
    out[[key]] <- val
    i <- i + 2
  }
  out
}

`%||%` <- function(x, y) {
  if (!is.null(x) && !is.na(x) && nzchar(x)) x else y
}

args <- parse_args(commandArgs(trailingOnly = TRUE))

mode       <- args$mode %||% "pilot_fixed"
setting    <- as.integer(args$setting %||% stop("Need --setting"))
g_type     <- args$g_type %||% stop("Need --g_type")
n          <- as.integer(args$n %||% stop("Need --n"))
n.tr       <- as.integer(args$n_tr %||% args$ntr %||% as.character(n))
d          <- as.integer(args$d %||% "200")
nsim_pilot <- as.integer(args$nsim_pilot %||% "50")
nsim       <- as.integer(args$nsim %||% "500")
outdir     <- args$outdir %||% if (mode == "cv_only") "results_cv_ntr2500" else "results"
rate       <- as.numeric(args$rate %||% "0.3")  # unused when est_r = TRUE
bandwidth  <- args$bandwidth %||% if (mode == "pilot_fixed") "pilot_fixed" else "cv"
h_fixed    <- as.numeric(args$h_fixed %||% args$h %||% NA_real_)
h_fixed_values_arg <- args$h_fixed_values %||% args$h_values %||% ""
h_fixed_values <- if (nzchar(h_fixed_values_arg)) {
  as.numeric(strsplit(h_fixed_values_arg, ",", fixed = TRUE)[[1]])
} else {
  h_fixed
}

if (!(mode %in% c("pilot_fixed", "cv_only"))) {
  stop("--mode must be one of: pilot_fixed, cv_only")
}
if (!(bandwidth %in% c("cv", "fixed", "oracle", "all", "pilot_fixed"))) {
  stop("--bandwidth must be one of: cv, fixed, oracle, all, pilot_fixed")
}
if (mode == "cv_only" && bandwidth == "pilot_fixed") {
  stop("--bandwidth pilot_fixed is only valid with --mode pilot_fixed")
}
if (bandwidth %in% c("fixed", "all") && any(!is.finite(h_fixed_values) | h_fixed_values <= 0)) {
  stop("--bandwidth fixed/all requires positive --h_fixed or comma-separated --h_fixed_values")
}
if (!(setting %in% c(1, 2, 3))) {
  stop("--setting must be one of 1, 2, 3")
}
if (!(g_type %in% c("null", "linear", "interact"))) {
  stop("--g_type must be one of: null, linear, interact")
}
if (is.na(n) || is.na(n.tr)) {
  stop("--n and --n_tr must be integers")
}
if (is.na(nsim_pilot) || nsim_pilot < 1) {
  stop("--nsim_pilot must be a positive integer")
}
if (is.na(nsim) || nsim < 1) {
  stop("--nsim must be a positive integer")
}

# -------------------------------
# Set working directory = script dir
# -------------------------------
cmdArgs <- commandArgs(trailingOnly = FALSE)
fileArg <- cmdArgs[grep("^--file=", cmdArgs)]
if (length(fileArg) > 0) {
  script_path <- sub("^--file=", "", fileArg[1])
  script_dir  <- normalizePath(dirname(script_path))
  setwd(script_dir)
}

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# -------------------------------
# Cores
# -------------------------------
cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "1"))
if (is.na(cores) || cores < 1) cores <- 1

message("SLURM_CPUS_PER_TASK = ", Sys.getenv("SLURM_CPUS_PER_TASK"),
        " -> using mc.cores = ", cores)

# -------------------------------
# Sources & packages
# -------------------------------
# Keep these file names in the same folder as this script.
source("./estimator_functions.R")
source("./simulation_functions.R")

suppressPackageStartupMessages({
  library(parallel)
  library(MASS)
  library(mvtnorm)
  library(splines)
  library(splines2)
})

# -------------------------------
# Model setup
# -------------------------------
set.seed(15)
Sigma <- 0.7^abs(row(diag(d)) - col(diag(d)))

poly_deg <- 1
beta_g <- rnorm(d, 0.5 / (1:d), 0.5 / (1:d))
beta_r <- 0.5 / (1:d)
eval.pts <- seq(0.15, 0.85, length.out = 50)
k1 <- 10
k2 <- 2

get_g <- function(g_type, beta_g = NULL) {
  if (g_type == "null") {
    return(function(beta_g, x) 0)
  }
  if (g_type == "linear") {
    return(function(beta_g, x) as.vector(x %*% beta_g))
  }
  if (g_type == "interact") {
    return(function(beta_g, x) beta_g[1] * x[, 1] * x[, 2])
  }
  stop("Unknown g_type: ", g_type)
}

r.x <- function(beta_r, x) as.vector(expit(x %*% beta_r))

make_f_setting <- function(setting, ...) {
  dots <- list(...)

  k1      <- dots$k1
  k2      <- dots$k2
  beta_r  <- dots$beta_r
  beta_g  <- dots$beta_g
  Sigma   <- dots$Sigma
  g_type  <- dots$g_type
  r.x     <- dots$r.x
  g.x     <- get_g(dots$g_type, beta_g = dots$beta_g)

  out <- switch(
    as.character(setting),

    "1" = list(
      m.t = function(t) {
        0.15 * sin(k1 * pi * t) + cos(k2 * pi * t) +
          truth_fun(t, g_type, beta_r, beta_g, Sigma)
      },
      der.m.t = function(t) {
        0.15 * k1 * pi * cos(k1 * pi * t) - k2 * pi * sin(k2 * pi * t) +
          deriv_fun(t, g_type, beta_r, beta_g, Sigma)
      },
      mu.x = function(x) {
        0.15 * sin(k1 * pi * r.x(beta_r, x)) +
          cos(k2 * pi * r.x(beta_r, x)) + g.x(beta_g, x)
      }
    ),

    "2" = list(
      m.t = function(t) {
        4 + 3 * t + 0.15 * sin(4 * pi * t) +
          truth_fun(t, g_type, beta_r, beta_g, Sigma)
      },
      der.m.t = function(t) {
        3 + 0.6 * pi * cos(4 * pi * t) +
          deriv_fun(t, g_type, beta_r, beta_g, Sigma)
      },

      mu.x = function(x) {
        4 + 3 * r.x(beta_r, x) +
          0.15 * sin(4 * pi * r.x(beta_r, x)) + g.x(beta_g, x)
      }
    ),

    "3" = list(
      m.t = function(t) {
        1 + 0.1 * sin(2 * pi * t) +
          sin(2 * pi * t) * exp(-40 * (t - 0.3)^2) +
          truth_fun(t, g_type, beta_r, beta_g, Sigma)
      },
      mu.x = function(x) {
        rr <- r.x(beta_r, x)
        1 + 0.1 * sin(2 * pi * rr) +
          sin(2 * pi * rr) * exp(-40 * (rr - 0.3)^2) +
          g.x(beta_g, x)
      },
      der.m.t = function(t) {
        0.2 * pi * cos(2 * pi * t) +
          2 * pi * cos(2 * pi * t) * exp(-40 * (t - 0.3)^2) +
          sin(2 * pi * t) * (-80) * (t - 0.3) * exp(-40 * (t - 0.3)^2) +
          deriv_fun(t, g_type, beta_r, beta_g, Sigma)
      }
    )
  )

  if (is.null(out)) stop("Unknown setting: ", setting)
  out
}

make_config <- function(fix.df.or.h = FALSE,
                        df.fourier = 1:10,
                        df.poly = 1:10,
                        df.bs = 1:10,
                        h.seq = seq(0.05, 0.35, length.out = 15),
                        fixed_tuning = NULL,
                        bandwidth_mode = "cv",
                        fixed_h = NA_real_) {
  config <- list(
    n = n,
    n.tr = n.tr,
    rate = rate,
    est_r = TRUE,
    eval.pts = eval.pts,
    poly_deg = poly_deg,

    d = d,
    beta_r = beta_r,
    beta_g = beta_g,
    Sigma = Sigma,
    r.x = r.x,
    g.x = NULL,
    mu.x = NULL,
    m.t = NULL,
    der.m.t = NULL,
    der2.m.t = NULL,
    g_type = g_type,
    truth_fun = truth_fun,

    df.seq.fourier = df.fourier,
    df.seq.poly = df.poly,
    df.seq.bs = df.bs,
    h.seq = h.seq,
    fix.df.or.h = fix.df.or.h,
    fixed_tuning = fixed_tuning,
    bandwidth_mode = bandwidth_mode,
    fixed_h = fixed_h
  )

  config$g.x <- get_g(beta_g = beta_g, g_type = g_type)

  fs <- make_f_setting(
    setting = setting,
    k1 = k1,
    k2 = k2,
    g_type = g_type,
    r.x = r.x,
    beta_r = beta_r,
    beta_g = beta_g,
    Sigma = Sigma
  )

  config$m.t      <- fs$m.t
  config$mu.x     <- fs$mu.x
  config$der.m.t  <- fs$der.m.t
  config$der2.m.t <- fs$der2.m.t

  config
}

run_sims <- function(config, nsim, seed_offset, label) {
  message("Running ", label, ": setting=", setting,
          " g_type=", g_type,
          " n=", config$n,
          " n.tr=", config$n.tr,
          " nsim=", nsim,
          " fixed_tuning=", isTRUE(config$fix.df.or.h),
          " bandwidth=", config$bandwidth_mode,
          " cores=", cores)

  mclapply(
    X = seq_len(nsim),
    FUN = function(s) simu_rslt(seed = seed_offset + s, config = config),
    mc.cores = cores,
    mc.preschedule = FALSE
  )
}

get_dfh <- function(one_result, estimator, basis) {
  x <- try(one_result[[estimator]][[basis]]$df.or.h, silent = TRUE)
  if (inherits(x, "try-error") || is.null(x)) return(NA_real_)
  as.numeric(x)
}

median_discrete <- function(x) {
  x <- sort(na.omit(as.numeric(x)))
  if (length(x) == 0) return(NA_real_)
  # For even nsim, this chooses an observed value rather than a half-integer.
  x[ceiling(length(x) / 2)]
}

choose_fixed_tuning <- function(pilot_result) {
  estimators <- c(
    "OR", "PI", "CPI_or", "CPI", "PI_cali",
    "BC_or_cali", "BC_cali", "BC", "BC_or"
  )
  bases <- c("fourier", "poly", "bs", "lpoly")

  rows <- list()
  k <- 1
  for (estimator in estimators) {
    for (basis in bases) {
      vals <- vapply(
        pilot_result,
        get_dfh,
        numeric(1),
        estimator = estimator,
        basis = basis
      )

      raw_median <- median(na.omit(vals))
      if (basis == "lpoly") {
        fixed_value <- raw_median
      } else {
        fixed_value <- as.numeric(median_discrete(vals))
      }

      rows[[k]] <- data.frame(
        estimator = estimator,
        basis = basis,
        median_raw = as.numeric(raw_median),
        fixed_value = as.numeric(fixed_value),
        n_nonmissing = sum(is.finite(vals)),
        stringsAsFactors = FALSE
      )
      k <- k + 1
    }
  }

  fixed_tbl <- do.call(rbind, rows)

  bad <- !is.finite(fixed_tbl$fixed_value)
  if (any(bad)) {
    print(fixed_tbl[bad, ])
    stop("Failed to extract finite pilot df/h choices for some estimator/basis rows.")
  }

  fixed_tbl
}

stamp <- format(Sys.time(), "%m%d_%H%M%S")
base_name <- sprintf("setting%d_%s_n%d_ntr%d", setting, g_type, n, n.tr)
h_grid <- seq(0.05, 0.35, length.out = 15)

format_h_tag <- function(h) gsub("[.]", "p", sprintf("%.3f", h))

run_cv_bandwidth <- function(bw_mode, fixed_h_value = h_fixed) {
  file_bw_tag <- bw_mode
  if (bw_mode == "fixed") {
    file_bw_tag <- paste0("fixedh", format_h_tag(fixed_h_value))
  }

  cv_config <- make_config(
    fix.df.or.h = FALSE,
    df.fourier = 1:10,
    df.poly = 1:10,
    df.bs = 1:10,
    h.seq = h_grid,
    fixed_tuning = NULL,
    bandwidth_mode = bw_mode,
    fixed_h = fixed_h_value
  )

  cv_result <- run_sims(
    config = cv_config,
    nsim = nsim,
    seed_offset = 100000,
    label = paste0("cv-only final, bandwidth=", bw_mode)
  )

  cv_file <- file.path(
    outdir,
    sprintf("cv_%s_nsim%d_bw%s_%s.rds", base_name, nsim, file_bw_tag, stamp)
  )

  saveRDS(
    list(
      result = cv_result,
      config = cv_config,
      fixed_tuning = NULL,
      bandwidth_mode = bw_mode,
      bandwidth_label = file_bw_tag,
      fixed_h = fixed_h_value,
      note = "CV-only run; local-polynomial bandwidth mode recorded in bandwidth_mode."
    ),
    cv_file
  )

  message("CV-only file: ", cv_file)
  invisible(cv_file)
}

if (mode == "cv_only") {
  if (bandwidth == "all") {
    out_files <- c(
      cv = run_cv_bandwidth("cv"),
      oracle = run_cv_bandwidth("oracle"),
      setNames(
        vapply(h_fixed_values, function(h) run_cv_bandwidth("fixed", fixed_h_value = h), character(1)),
        paste0("fixed_h_", format_h_tag(h_fixed_values))
      )
    )
  } else if (bandwidth == "fixed") {
    out_files <- setNames(
      vapply(h_fixed_values, function(h) run_cv_bandwidth("fixed", fixed_h_value = h), character(1)),
      paste0("fixed_h_", format_h_tag(h_fixed_values))
    )
  } else {
    out_files <- setNames(run_cv_bandwidth(bandwidth), bandwidth)
  }

  message("Finished.")
  message("CV-only file(s):")
  print(out_files)
  quit(save = "no", status = 0)
}

# -------------------------------
# 1. Pilot run: CV tuning
# -------------------------------
pilot_config <- make_config(
  fix.df.or.h = FALSE,
  df.fourier = 1:10,
  df.poly = 1:10,
  df.bs = 1:10,
  h.seq = h_grid
)

pilot_result <- run_sims(
  config = pilot_config,
  nsim = nsim_pilot,
  seed_offset = 0,
  label = "pilot"
)

fixed_tuning <- choose_fixed_tuning(pilot_result)

message("Pilot fixed tuning choices by estimator and basis:")
print(fixed_tuning)

# -------------------------------
# 2. Final run: fixed median df/h
# -------------------------------
final_config <- make_config(
  fix.df.or.h = TRUE,
  # Keep default grids as fallback only. The actual fixed choices are read
  # from fixed_tuning inside simu_rslt(), estimator by estimator.
  df.fourier = 1:10,
  df.poly = 1:10,
  df.bs = 1:10,
  h.seq = h_grid,
  fixed_tuning = fixed_tuning,
  bandwidth_mode = "pilot_fixed",
  fixed_h = NA_real_
)

final_result <- run_sims(
  config = final_config,
  nsim = nsim,
  seed_offset = 100000,
  label = "final fixed-tuning"
)

final_file <- file.path(
  outdir,
  sprintf("fixed_%s_pilot%d_nsim%d_bwpilot_fixed_%s.rds", base_name, nsim_pilot, nsim, stamp)
)
saveRDS(
  list(
    result = final_result,
    config = final_config,
    pilot_result = pilot_result,
    fixed_tuning = fixed_tuning,
    pilot_config = pilot_config,
    bandwidth_mode = "pilot_fixed",
    fixed_h = NA_real_,
    note = "Pilot and fixed-tuning final results are stored together."
  ),
  final_file
)
message("Saved final: ", final_file)

message("Finished.")
message("Final file: ", final_file)
