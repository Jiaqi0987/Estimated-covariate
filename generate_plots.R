#!/usr/bin/env Rscript

# Generate PDF plots from simulation result files:
#   1. Combined main plots from CV-only and/or pilot-fixed results.
#   2. CI-length plots.
#   3. Smoothing-target coverage plots.
#
# This script writes PDFs only. Smoothing targets are cached as .rds files.

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(splines)
  library(splines2)
})

THIS_DIR <- normalizePath(getwd(), mustWork = FALSE)
CV_RESULT_DIR <- Sys.getenv("CV_RESULT_DIR", unset = file.path(THIS_DIR, "results_cv_ntr2500"))
FIXED_RESULT_DIR <- Sys.getenv("FIXED_RESULT_DIR", unset = file.path(THIS_DIR, "results"))
OUT_ROOT <- Sys.getenv("OUT_DIR", unset = file.path(THIS_DIR, "plots"))
MAIN_OUT <- file.path(OUT_ROOT, "main_plots")
CI_OUT <- file.path(OUT_ROOT, "ci_length")
COV_OUT <- file.path(OUT_ROOT, "smoothing_target_coverage")
CACHE_DIR <- file.path(COV_OUT, "target_cache")

dir.create(MAIN_OUT, showWarnings = FALSE, recursive = TRUE)
dir.create(CI_OUT, showWarnings = FALSE, recursive = TRUE)
dir.create(COV_OUT, showWarnings = FALSE, recursive = TRUE)
dir.create(CACHE_DIR, showWarnings = FALSE, recursive = TRUE)

source(file.path(THIS_DIR, "estimator_functions.R"))
source(file.path(THIS_DIR, "simulation_functions.R"))
source(file.path(THIS_DIR, "plot_helpers.R"))

`%||%` <- function(x, y) if (!is.null(x)) x else y

EST_ALL <- c("PI", "OR", "CPI", "CPI_or",
             "BC", "BC_or", "PI_cali", "BC_cali", "BC_or_cali")
EST_ERR <- c("PI", "OR", "CPI", "BC", "PI_cali", "BC_cali")
BASIS <- c("fourier", "bs", "poly", "lpoly")
BASIS_LABELS <- c(fourier = "Cosine", bs = "B-spline",
                  poly = "Polynomial", lpoly = "Local linear")
PALETTE <- c(
  PI = "#FC8D62", OR = "#999999", CPI = "#1F78B4", CPI_or = "#A6CEE3",
  BC = "#33A02C", BC_or = "#B2DF8A", PI_cali = "#E7298A",
  BC_cali = "#FCC5C0", BC_or_cali = "#C994C7"
)
EST_LABELS <- c(
  PI = expression(bold(widehat(m)["pi"])),
  CPI = expression(bold(widehat(m)[cpi])),
  BC = expression(bold(widehat(m)[bc])),
  OR = expression(bold(widehat(m)[or])),
  PI_cali = expression(bold(widehat(m)["pi"~"."~cali])),
  BC_cali = expression(bold(widehat(m)[bc~"."~cali])),
  CPI_or = expression(bold(widehat(m)[cpi~"."~or])),
  BC_or = expression(bold(widehat(m)[bc~"."~or])),
  BC_or_cali = expression(bold(widehat(m)[bc~"."~or~"."~cali]))
)
CURVE_LABELS <- c(Truth = "Truth", PI = EST_LABELS[["PI"]],
                  CPI = EST_LABELS[["CPI"]], BC = EST_LABELS[["BC"]])

parse_meta <- function(path) {
  b <- basename(path)
  info <- file.info(path)
  bandwidth <- if (startsWith(b, "fixed_")) {
    "pilot_fixed"
  } else if (grepl("_bw[^_]+_", b)) {
    sub(".*_bw([^_]+)_.*", "\\1", b)
  } else {
    "cv"
  }
  data.frame(
    setting = as.integer(sub(".*setting([0-9]+)_.*", "\\1", b)),
    g_type = sub(".*setting[0-9]+_([^_]+)_n.*", "\\1", b),
    n = as.integer(sub(".*_n([0-9]+)_ntr.*", "\\1", b)),
    n_tr = as.integer(sub(".*_ntr([0-9]+).*", "\\1", b)),
    nsim = suppressWarnings(as.integer(sub(".*_nsim([0-9]+).*", "\\1", b))),
    bandwidth = bandwidth,
    file = b,
    path = path,
    mtime = info$mtime,
    stringsAsFactors = FALSE
  )
}

get_truth <- function(config) as.numeric(config$m.t(config$eval.pts))

result_files <- function() {
  cv_files <- list.files(
    CV_RESULT_DIR,
    pattern = "^cv_setting[0-9]+_.*[.]rds$",
    full.names = TRUE
  )
  fixed_files <- if (dir.exists(FIXED_RESULT_DIR)) {
    list.files(FIXED_RESULT_DIR, pattern = "^fixed_setting[0-9]+_.*[.]rds$", full.names = TRUE)
  } else {
    character(0)
  }
  c(cv_files, fixed_files)
}

latest_by_group <- function(files) {
  if (length(files) == 0) return(data.frame())
  bind_rows(lapply(files, parse_meta)) %>%
    group_by(setting, g_type, n, n_tr, nsim, bandwidth) %>%
    arrange(desc(mtime), .by_group = TRUE) %>%
    slice(1) %>%
    ungroup() %>%
    arrange(g_type, setting, n_tr, bandwidth)
}

tag_for <- function(meta) {
  sprintf("setting%d_%s_ntr%d_nsim%d_bw%s",
          meta$setting, meta$g_type, meta$n_tr, meta$nsim, meta$bandwidth)
}

make_curve_df <- function(item) {
  config <- item$res$config
  results_list <- Filter(Negate(is.null), item$res$result)
  eval_pts <- config$eval.pts
  truth <- get_truth(config)
  res1 <- results_list[[1]]
  rbind(
    data.frame(t = eval_pts, value = truth, label = "Truth"),
    data.frame(t = eval_pts, value = as.numeric(res1$PI$lpoly$est), label = "PI"),
    data.frame(t = eval_pts, value = as.numeric(res1$CPI$lpoly$est), label = "CPI"),
    data.frame(t = eval_pts, value = as.numeric(res1$BC$lpoly$est), label = "BC")
  )
}

make_scatter_plot <- function(item) {
  config <- item$res$config
  curve_df <- make_curve_df(item)
  curve_df$label <- factor(curve_df$label, levels = c("Truth", "PI", "CPI", "BC"))

  set.seed(42)
  samp <- gen_data(beta_r = config$beta_r, beta_g = config$beta_g,
                   n = config$n, d = config$d, Sigma = config$Sigma,
                   r.x = config$r.x, mu.x = config$mu.x)
  scatter_df <- data.frame(r = samp$dat$r, y = samp$dat$y)
  yr <- quantile(scatter_df$y, c(0.01, 0.99), na.rm = TRUE)

  ggplot() +
    geom_point(data = scatter_df, aes(r, y), color = "gray55",
               alpha = 0.18, size = 0.45) +
    geom_line(data = curve_df, aes(t, value, color = label, linewidth = label, alpha = label)) +
    scale_color_manual(values = c(Truth = "#B2182B", PI = "#FC8D62",
                                  CPI = "#1F78B4", BC = "#1B7837"),
                       labels = CURVE_LABELS, name = NULL) +
    scale_linewidth_manual(values = c(Truth = 1.05, PI = 0.75, CPI = 0.75, BC = 0.75),
                           guide = "none") +
    scale_alpha_manual(values = c(Truth = 1, PI = 0.95, CPI = 0.85, BC = 0.85),
                       guide = "none") +
    coord_cartesian(xlim = range(config$eval.pts), ylim = yr) +
    labs(x = expression(bold(t)), y = expression(bold(Y))) +
    theme_bw(base_size = 18) +
    theme(
      legend.position = "top",
      legend.text = element_text(size = 18, face = "bold"),
      legend.key.width = unit(1.7, "cm"),
      legend.key.height = unit(0.8, "cm"),
      axis.title = element_text(size = 18, face = "bold"),
      axis.text = element_text(size = 15, face = "bold")
    )
}

make_ptwise_err_plot <- function(item) {
  results_list <- Filter(Negate(is.null), item$res$result)
  config <- item$res$config
  eval_pts <- config$eval.pts
  truth <- get_truth(config)
  all_err <- unlist(lapply(EST_ERR, function(est) {
    as.numeric(get_pointwise_mse(results_list, est, "lpoly", truth))
  }), use.names = FALSE)
  y_cap <- as.numeric(quantile(all_err[is.finite(all_err)], 0.99, na.rm = TRUE))
  if (!is.finite(y_cap) || y_cap <= 0) y_cap <- 1

  make_plot_df <- function(err_mat) {
    nsim <- nrow(err_mat)
    as.data.frame(err_mat) %>%
      mutate(sim = seq_len(nsim)) %>%
      pivot_longer(cols = -sim, names_to = "t_index", values_to = "err") %>%
      mutate(t_index = as.integer(gsub("V", "", t_index)),
             t = eval_pts[t_index])
  }
  plot_panel <- function(plot_df, title_label) {
    mean_df <- plot_df %>% group_by(t) %>%
      summarise(mean_err = mean(err, na.rm = TRUE), .groups = "drop")
    ggplot() +
      geom_line(data = plot_df, aes(t, err, group = sim),
                color = "gray80", alpha = 0.18, linewidth = 0.2) +
      geom_line(data = mean_df, aes(t, mean_err), color = "red", linewidth = 0.95) +
      labs(title = title_label, x = "t", y = NULL) +
      coord_cartesian(xlim = range(eval_pts), ylim = c(0, y_cap)) +
      theme_bw(base_size = 16) +
      theme(plot.title = element_text(size = 19, face = "bold", hjust = 0.5),
            axis.title.x = element_text(size = 15),
            axis.text = element_text(size = 13))
  }
  plot_list <- lapply(EST_ERR, function(est) {
    plot_panel(make_plot_df(get_pointwise_mse(results_list, est, "lpoly", truth)),
               EST_LABELS[[est]])
  })
  (plot_list[[1]] | plot_list[[2]]) /
    (plot_list[[3]] | plot_list[[4]]) /
    (plot_list[[5]] | plot_list[[6]])
}

make_mse_plot <- function(item) {
  results_list <- Filter(Negate(is.null), item$res$result)
  truth <- get_truth(item$res$config)
  rows <- list()
  for (est in EST_ALL) {
    for (basis in BASIS) {
      vals <- tryCatch(get_average_mse(results_list, est, basis, truth),
                       error = function(e) NULL)
      if (!is.null(vals)) rows[[length(rows) + 1]] <- data.frame(estimator = est, basis = basis, mse = vals)
    }
  }
  mse_df <- bind_rows(rows) %>%
    mutate(estimator = factor(estimator, levels = EST_ALL),
           basis = factor(basis, levels = BASIS))
  y_cap <- as.numeric(quantile(mse_df$mse, 0.99, na.rm = TRUE))
  if (!is.finite(y_cap) || y_cap <= 0) y_cap <- max(mse_df$mse, na.rm = TRUE)
  ggplot(mse_df, aes(estimator, mse, fill = estimator)) +
    geom_boxplot(outlier.size = 0.35, outlier.alpha = 0.35, linewidth = 0.35) +
    facet_wrap(~ basis, nrow = 1, scales = "fixed",
               labeller = labeller(basis = BASIS_LABELS)) +
    scale_fill_manual(values = PALETTE, drop = FALSE) +
    coord_cartesian(ylim = c(0, y_cap)) +
    labs(x = NULL, y = "Average MSE") +
    theme_bw(base_size = 16) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14, face = "bold"),
          axis.text.y = element_text(size = 12),
          axis.title.y = element_text(size = 15, face = "bold"),
          legend.position = "none",
          strip.text = element_text(face = "bold", size = 15))
}

make_main_plot <- function(item) {
  top_row <- (make_scatter_plot(item) | wrap_elements(full = make_ptwise_err_plot(item))) +
    plot_layout(widths = c(1, 2.1))
  (top_row / wrap_elements(make_mse_plot(item))) +
    plot_layout(heights = c(1.05, 1))
}

get_se_vec <- function(res, est, basis) {
  obj <- res[[est]][[basis]]
  if (is.null(obj)) return(NULL)
  se <- if (basis == "lpoly") obj$se %||% obj$se.us else obj$se.adj %||% obj$se
  if (is.null(se)) return(NULL)
  as.numeric(se)
}

get_se_mat <- function(results_list, est, basis) {
  vals <- lapply(results_list, get_se_vec, est = est, basis = basis)
  if (any(vapply(vals, is.null, logical(1)))) return(NULL)
  do.call(rbind, vals)
}

get_realized_dfh <- function(results_list, est, basis) {
  vals <- vapply(results_list, function(res) {
    obj <- res[[est]][[basis]]
    val <- obj$df.or.h %||% obj$df %||% obj$h
    if (is.null(val)) return(NA_real_)
    as.numeric(val)[1]
  }, numeric(1))
  if (basis %in% c("fourier", "bs", "poly")) vals <- as.integer(round(vals))
  vals
}

get_fixed_dfh <- function(fixed_tuning, est, basis) {
  if (is.null(fixed_tuning)) return(NA_real_)
  idx <- which(fixed_tuning$estimator == est & fixed_tuning$basis == basis)
  if (length(idx) == 0) return(NA_real_)
  val <- as.numeric(fixed_tuning$fixed_value[idx[1]])
  if (!is.finite(val)) return(NA_real_)
  if (basis %in% c("fourier", "bs", "poly")) val <- as.integer(round(val))
  val
}

make_ci_length_plot <- function(item) {
  results_list <- Filter(Negate(is.null), item$res$result)
  eval_pts <- item$res$config$eval.pts
  rows <- bind_rows(lapply(seq_along(results_list), function(rep_id) {
    res <- results_list[[rep_id]]
    bind_rows(lapply(BASIS, function(basis) {
      bind_rows(lapply(EST_ALL, function(est) {
        se <- tryCatch(get_se_vec(res, est, basis), error = function(e) NULL)
        if (is.null(se)) return(NULL)
        data.frame(rep = rep_id, basis = basis, estimator = est,
                   t = eval_pts, ci_length = 2 * 1.96 * se)
      }))
    }))
  }))
  pointwise <- rows %>%
    group_by(basis, estimator, t) %>%
    summarise(mean_ci_length = mean(ci_length, na.rm = TRUE), .groups = "drop") %>%
    mutate(basis_f = factor(BASIS_LABELS[basis], levels = unname(BASIS_LABELS[BASIS])),
           estimator = factor(estimator, levels = EST_ALL))
  ggplot(pointwise, aes(t, mean_ci_length, color = estimator)) +
    geom_line(linewidth = 0.6, alpha = 0.95) +
    facet_wrap(~ basis_f, nrow = 1, ncol = 4, scales = "free_y") +
    scale_color_manual(values = PALETTE, labels = EST_LABELS, drop = FALSE, name = NULL) +
    labs(x = "t", y = "Mean 95% CI length") +
    theme_bw(base_size = 15) +
    theme(legend.position = "bottom",
          legend.text = element_text(size = 15, face = "bold"),
          legend.key.width = grid::unit(2.0, "lines"),
          legend.key.height = grid::unit(1.1, "lines"),
          strip.text = element_text(face = "bold", size = 14),
          axis.title = element_text(size = 15, face = "bold"),
          axis.text = element_text(size = 13)) +
    guides(color = guide_legend(nrow = 2))
}

make_sieve_targets_one <- function(r_pop, mr_pop, eval_pts, df_seq, basis) {
  basis_str <- switch(basis, fourier = "fourier", bs = "bspline", poly = "polynomial")
  targets <- lapply(df_seq, function(k) {
    Phi <- cbind(1, basis_fun(r_pop, k, basis = basis_str))
    Phi_eval <- cbind(1, basis_fun(eval_pts, k, basis = basis_str))
    XtX <- crossprod(Phi)
    beta_k <- tryCatch(
      solve(XtX, crossprod(Phi, mr_pop)),
      error = function(e) solve(XtX + 1e-8 * diag(ncol(XtX)), crossprod(Phi, mr_pop))
    )
    as.numeric(Phi_eval %*% beta_k)
  })
  names(targets) <- as.character(df_seq)
  targets
}

get_sieve_targets <- function(config, tag, n_pop = as.integer(Sys.getenv("N_POP_SIEVE", "50000"))) {
  cache_file <- file.path(CACHE_DIR, sprintf("sieve_targets_%s_npop%d.rds", tag, n_pop))
  if (file.exists(cache_file)) return(readRDS(cache_file))
  message("  computing sieve smoothing targets for ", tag, " with n_pop = ", n_pop)
  set.seed(0)
  pop <- gen_data(beta_r = config$beta_r, beta_g = config$beta_g,
                  n = n_pop, d = config$d, Sigma = config$Sigma,
                  r.x = config$r.x, mu.x = config$mu.x)
  r_pop <- pop$dat$r
  mr_pop <- as.numeric(config$m.t(r_pop))
  targets <- list(
    fourier = make_sieve_targets_one(r_pop, mr_pop, config$eval.pts,
                                     config$df.seq.fourier, "fourier"),
    bs = make_sieve_targets_one(r_pop, mr_pop, config$eval.pts,
                                config$df.seq.bs, "bs"),
    poly = make_sieve_targets_one(r_pop, mr_pop, config$eval.pts,
                                  config$df.seq.poly, "poly")
  )
  saveRDS(targets, cache_file)
  targets
}

get_lpoly_population <- function(config, tag, n_pop = as.integer(Sys.getenv("N_POP_LPOLY", "50000"))) {
  cache_file <- file.path(CACHE_DIR, sprintf("lpoly_population_%s_npop%d.rds", tag, n_pop))
  if (file.exists(cache_file)) return(readRDS(cache_file))
  message("  computing local-linear population for ", tag, " with n_pop = ", n_pop)
  set.seed(1)
  pop <- gen_data(beta_r = config$beta_r, beta_g = config$beta_g,
                  n = n_pop, d = config$d, Sigma = config$Sigma,
                  r.x = config$r.x, mu.x = config$mu.x)
  out <- list(r = as.numeric(pop$dat$r), m = as.numeric(config$m.t(pop$dat$r)))
  saveRDS(out, cache_file)
  out
}

make_lpoly_targets <- function(config, tag, h_values) {
  pop <- get_lpoly_population(config, tag)
  h_key <- function(h) sprintf("%.3f", as.numeric(h))
  local_linear_gaussian <- function(r, y, eval_pts, h) {
    vapply(eval_pts, function(t0) {
      x <- r - t0
      w <- dnorm(x / h)
      s0 <- sum(w); s1 <- sum(w * x); s2 <- sum(w * x * x)
      y0 <- sum(w * y); y1 <- sum(w * x * y)
      denom <- s0 * s2 - s1 * s1
      if (!is.finite(denom) || abs(denom) < 1e-12) return(NA_real_)
      (s2 * y0 - s1 * y1) / denom
    }, numeric(1))
  }
  do.call(rbind, lapply(h_values, function(h) {
    h_round <- round(as.numeric(h), 3)
    cache_file <- file.path(CACHE_DIR, sprintf("lpoly_target_%s_h%s.rds",
                                               tag, gsub("[.]", "p", h_key(h_round))))
    if (file.exists(cache_file)) return(readRDS(cache_file))
    target <- local_linear_gaussian(pop$r, pop$m, config$eval.pts, h_round)
    saveRDS(target, cache_file)
    target
  }))
}

make_coverage_plot <- function(item, meta) {
  results_list <- Filter(Negate(is.null), item$res$result)
  config <- item$res$config
  fixed_tuning <- item$res$fixed_tuning
  eval_pts <- config$eval.pts
  cov_tag <- sprintf("setting%d_%s_ntr%d_%s", meta$setting, meta$g_type, meta$n_tr, meta$bandwidth)
  sieve_targets <- get_sieve_targets(config, cov_tag)

  cov_df <- bind_rows(lapply(EST_ALL, function(est) {
    bind_rows(lapply(BASIS, function(basis) {
      tryCatch({
        est_mat <- get_preds_mat(results_list, est, basis)
        se_mat <- get_se_mat(results_list, est, basis)
        if (is.null(se_mat) || !identical(dim(est_mat), dim(se_mat))) return(NULL)
        fixed_dfh <- get_fixed_dfh(fixed_tuning, est, basis)
        df_h <- if (is.finite(fixed_dfh)) rep(fixed_dfh, nrow(est_mat)) else get_realized_dfh(results_list, est, basis)
        target_mat <- if (basis == "lpoly") {
          make_lpoly_targets(config, cov_tag, df_h)
        } else {
          do.call(rbind, lapply(as.character(df_h), function(k) sieve_targets[[basis]][[k]]))
        }
        data.frame(
          basis = basis,
          estimator = est,
          t = eval_pts,
          coverage = colMeans(est_mat - 1.96 * se_mat <= target_mat &
                                target_mat <= est_mat + 1.96 * se_mat,
                              na.rm = TRUE)
        )
      }, error = function(e) {
        message("  coverage skip ", est, " / ", basis, ": ", conditionMessage(e))
        NULL
      })
    }))
  })) %>%
    mutate(ntr_f = factor(paste0("n.tr = ", meta$n_tr)),
           basis_f = factor(BASIS_LABELS[basis], levels = unname(BASIS_LABELS[BASIS])),
           estimator = factor(estimator, levels = EST_ALL))

  ggplot(cov_df, aes(t, coverage, color = estimator)) +
    geom_hline(yintercept = 0.95, linetype = "dashed", color = "gray35", linewidth = 0.45) +
    geom_line(linewidth = 0.6, alpha = 0.95) +
    facet_grid(ntr_f ~ basis_f) +
    scale_color_manual(values = PALETTE, labels = EST_LABELS, drop = FALSE, name = NULL) +
    scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 0.95, 1)) +
    labs(x = "t", y = "Coverage") +
    theme_bw(base_size = 15) +
    theme(legend.position = "bottom",
          legend.text = element_text(size = 15, face = "bold"),
          legend.key.width = grid::unit(2.0, "lines"),
          legend.key.height = grid::unit(1.1, "lines"),
          strip.text = element_text(face = "bold", size = 14),
          axis.title = element_text(size = 15, face = "bold"),
          axis.text = element_text(size = 13)) +
    guides(color = guide_legend(nrow = 2))
}

all_meta <- latest_by_group(result_files())
if (nrow(all_meta) == 0) {
  stop("No result files found in ", CV_RESULT_DIR, " or ", FIXED_RESULT_DIR)
}
message("Using ", nrow(all_meta), " latest result file(s)")

for (i in seq_len(nrow(all_meta))) {
  meta <- all_meta[i, ]
  tag <- tag_for(meta)
  message("Building plots: ", tag)
  res_obj <- readRDS(meta$path)
  item <- list(setting = meta$setting, g_type = meta$g_type, n_tr = meta$n_tr, res = res_obj)

  ggsave(file.path(MAIN_OUT, sprintf("main_%s.pdf", tag)),
         make_main_plot(item), width = 18.5, height = 12, bg = "white")
  ggsave(file.path(CI_OUT, sprintf("ci_length_%s.pdf", tag)),
         make_ci_length_plot(item), width = 18, height = 5, bg = "white")
  ggsave(file.path(COV_OUT, sprintf("coverage_smoothing_target_%s.pdf", tag)),
         make_coverage_plot(item, meta), width = 15, height = 8, bg = "white")
}

message("All PDF plots saved to: ", OUT_ROOT)
