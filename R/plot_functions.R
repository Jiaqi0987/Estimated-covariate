#!/usr/bin/env Rscript

# Plot and summary functions for the estimated-covariates simulation.
# These functions assume result files are saved as .rds objects produced by
# run_simulation_and_plot.R.

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

`%||%` <- function(x, y) if (!is.null(x)) x else y

BASIS <- c("fourier", "bs", "poly", "lpoly")
BASIS_LABELS <- c(
  fourier = "Fourier",
  bs = "B-spline",
  poly = "Polynomial",
  lpoly = "Local linear"
)

ESTIMATORS <- c(
  "PI", "OR", "CPI", "CPI_or", "BC", "BC_or",
  "PI_cali", "BC_cali", "BC_or_cali"
)

COVERAGE_ESTIMATORS <- setdiff(ESTIMATORS, "OR")

PALETTE <- c(
  PI = "#FC8D62",
  OR = "#999999",
  CPI = "#1F78B4",
  CPI_or = "#A6CEE3",
  BC = "#33A02C",
  BC_or = "#B2DF8A",
  PI_cali = "#E7298A",
  BC_cali = "#FCC5C0",
  BC_or_cali = "#C994C7"
)

parse_result_meta <- function(path) {
  b <- basename(path)
  data.frame(
    setting = as.integer(sub(".*setting([0-9]+)_.*", "\\1", b)),
    g_type = sub(".*setting[0-9]+_([^_]+)_n.*", "\\1", b),
    n = as.integer(sub(".*_n([0-9]+)_ntr.*", "\\1", b)),
    n_tr = as.integer(sub(".*_ntr([0-9]+)_nsim.*", "\\1", b)),
    file = b,
    path = path,
    stringsAsFactors = FALSE
  )
}

list_result_files <- function(result_dir) {
  files <- list.files(
    result_dir,
    pattern = "^sim_semior_setting[0-9]+_[^_]+_n[0-9]+_ntr[0-9]+_nsim[0-9]+_.*[.]rds$",
    full.names = TRUE
  )
  if (length(files) == 0) {
    stop("No simulation .rds files found in: ", result_dir)
  }
  files
}

get_est_mat <- function(results_list, est, basis) {
  vals <- lapply(results_list, function(res) {
    obj <- res[[est]][[basis]]
    if (is.null(obj) || is.null(obj$est)) return(rep(NA_real_, 0))
    as.numeric(obj$est)
  })
  do.call(rbind, vals)
}

get_se_mat <- function(results_list, est, basis) {
  vals <- lapply(results_list, function(res) {
    obj <- res[[est]][[basis]]
    if (is.null(obj)) return(rep(NA_real_, 0))
    se <- if (basis == "lpoly") obj$se %||% obj$se.us else obj$se.adj %||% obj$se
    if (is.null(se)) return(rep(NA_real_, 0))
    as.numeric(se)
  })
  do.call(rbind, vals)
}

get_semior_mat <- function(results_list, est, basis) {
  vals <- lapply(results_list, function(res) {
    obj <- res$semi_or[[est]][[basis]]
    if (is.null(obj) || is.null(obj$est)) return(rep(NA_real_, 0))
    as.numeric(obj$est)
  })
  do.call(rbind, vals)
}

read_result <- function(path) {
  out <- readRDS(path)
  out$result <- Filter(Negate(is.null), out$result)
  out
}

coverage_rows_one_file <- function(path, basis = BASIS, estimators = COVERAGE_ESTIMATORS) {
  out <- read_result(path)
  meta <- parse_result_meta(path)
  eval_pts <- out$config$eval.pts

  bind_rows(lapply(estimators, function(est) {
    bind_rows(lapply(basis, function(b) {
      tryCatch({
        est_mat <- get_est_mat(out$result, est, b)
        se_mat <- get_se_mat(out$result, est, b)
        so_mat <- get_semior_mat(out$result, est, b)
        if (nrow(est_mat) == 0 || nrow(se_mat) == 0 || nrow(so_mat) == 0) {
          return(NULL)
        }
        cov <- colMeans(
          est_mat - 1.96 * se_mat <= so_mat &
            so_mat <= est_mat + 1.96 * se_mat,
          na.rm = TRUE
        )
        data.frame(
          setting = meta$setting,
          g_type = meta$g_type,
          n = meta$n,
          n_tr = meta$n_tr,
          basis = b,
          estimator = est,
          t = eval_pts,
          coverage = cov,
          file = meta$file
        )
      }, error = function(e) NULL)
    }))
  }))
}

ci_length_rows_one_file <- function(path, basis = BASIS, estimators = ESTIMATORS) {
  out <- read_result(path)
  meta <- parse_result_meta(path)
  eval_pts <- out$config$eval.pts

  bind_rows(lapply(estimators, function(est) {
    bind_rows(lapply(basis, function(b) {
      tryCatch({
        se_mat <- get_se_mat(out$result, est, b)
        if (nrow(se_mat) == 0) return(NULL)
        data.frame(
          setting = meta$setting,
          g_type = meta$g_type,
          n = meta$n,
          n_tr = meta$n_tr,
          basis = b,
          estimator = est,
          t = eval_pts,
          ci_length = colMeans(2 * 1.96 * se_mat, na.rm = TRUE),
          file = meta$file
        )
      }, error = function(e) NULL)
    }))
  }))
}

mse_rows_one_file <- function(path, basis = BASIS, estimators = ESTIMATORS) {
  out <- read_result(path)
  meta <- parse_result_meta(path)
  eval_pts <- out$config$eval.pts
  truth <- as.numeric(out$config$m.t(eval_pts))

  bind_rows(lapply(estimators, function(est) {
    bind_rows(lapply(basis, function(b) {
      tryCatch({
        est_mat <- get_est_mat(out$result, est, b)
        if (nrow(est_mat) == 0) return(NULL)
        mse <- colMeans((est_mat - matrix(truth, nrow(est_mat), length(truth), byrow = TRUE))^2,
                        na.rm = TRUE)
        data.frame(
          setting = meta$setting,
          g_type = meta$g_type,
          n = meta$n,
          n_tr = meta$n_tr,
          basis = b,
          estimator = est,
          t = eval_pts,
          mse = mse,
          file = meta$file
        )
      }, error = function(e) NULL)
    }))
  }))
}

add_plot_labels <- function(df) {
  df %>%
    mutate(
      tag = paste0("setting", setting, "_", g_type),
      basis_f = factor(BASIS_LABELS[basis], levels = unname(BASIS_LABELS[BASIS])),
      estimator = factor(estimator, levels = ESTIMATORS),
      ntr_f = factor(paste0("n.tr = ", n_tr), levels = paste0("n.tr = ", sort(unique(n_tr))))
    )
}

plot_pointwise <- function(df, yvar, ylab, outfile, hline = NULL) {
  p <- ggplot(df, aes(x = t, y = .data[[yvar]], color = estimator)) +
    geom_line(linewidth = 0.55, alpha = 0.95) +
    facet_grid(ntr_f ~ basis_f, drop = TRUE) +
    scale_color_manual(values = PALETTE, drop = FALSE, name = "") +
    labs(x = "t", y = ylab) +
    theme_bw(base_size = 12) +
    theme(
      legend.position = "bottom",
      legend.text = element_text(size = 10),
      strip.text = element_text(face = "bold"),
      plot.title = element_blank()
    ) +
    guides(color = guide_legend(nrow = 2))

  if (!is.null(hline)) {
    p <- p + geom_hline(yintercept = hline, linetype = "dashed", color = "gray35")
  }
  if (yvar == "coverage") {
    p <- p + scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 0.95, 1))
  }
  ggsave(outfile, p, width = 15, height = 10, bg = "white")
  invisible(p)
}

summarise_pointwise <- function(df, value_col) {
  df %>%
    group_by(setting, g_type, n, n_tr, basis, estimator) %>%
    summarise(
      mean_value = mean(.data[[value_col]], na.rm = TRUE),
      min_value = min(.data[[value_col]], na.rm = TRUE),
      max_value = max(.data[[value_col]], na.rm = TRUE),
      .groups = "drop"
    )
}

generate_plots <- function(result_dir = "results",
                           out_dir = file.path(result_dir, "plots"),
                           basis = BASIS) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  files <- list_result_files(result_dir)

  coverage <- bind_rows(lapply(files, coverage_rows_one_file, basis = basis)) %>% add_plot_labels()
  ci_length <- bind_rows(lapply(files, ci_length_rows_one_file, basis = basis)) %>% add_plot_labels()
  mse <- bind_rows(lapply(files, mse_rows_one_file, basis = basis)) %>% add_plot_labels()

  write.csv(coverage, file.path(out_dir, "coverage_pointwise.csv"), row.names = FALSE)
  write.csv(ci_length, file.path(out_dir, "ci_length_pointwise.csv"), row.names = FALSE)
  write.csv(mse, file.path(out_dir, "mse_pointwise.csv"), row.names = FALSE)
  write.csv(summarise_pointwise(coverage, "coverage"), file.path(out_dir, "coverage_summary.csv"), row.names = FALSE)
  write.csv(summarise_pointwise(ci_length, "ci_length"), file.path(out_dir, "ci_length_summary.csv"), row.names = FALSE)
  write.csv(summarise_pointwise(mse, "mse"), file.path(out_dir, "mse_summary.csv"), row.names = FALSE)

  for (tag_i in sort(unique(coverage$tag))) {
    plot_pointwise(
      filter(coverage, tag == tag_i),
      "coverage",
      "Coverage of matched semi-oracle target",
      file.path(out_dir, paste0(tag_i, "_coverage.pdf")),
      hline = 0.95
    )
    plot_pointwise(
      filter(ci_length, tag == tag_i),
      "ci_length",
      "Mean CI length",
      file.path(out_dir, paste0(tag_i, "_ci_length.pdf"))
    )
    plot_pointwise(
      filter(mse, tag == tag_i),
      "mse",
      "Pointwise MSE relative to truth",
      file.path(out_dir, paste0(tag_i, "_mse_truth.pdf"))
    )
  }

  invisible(list(coverage = coverage, ci_length = ci_length, mse = mse))
}

