#!/usr/bin/env Rscript

# Generate all publication plots from saved simulation results:
#   - main plots from CV-only results in results_cv_ntr2500/
#   - smoothing-target coverage from pilot-fixed results in results/
#   - CI length for n.tr = 2500 from pilot-fixed results in results/

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(splines)
  library(splines2)
  library(nprobust)
  library(glmnet)
  library(mvtnorm)
  library(MASS)
})

THIS_DIR <- normalizePath(getwd(), mustWork = FALSE)
CV_RESULT_DIR <- Sys.getenv("CV_RESULT_DIR", unset = file.path(THIS_DIR, "results_cv_ntr2500"))
FIXED_RESULT_DIR <- Sys.getenv("FIXED_RESULT_DIR", unset = file.path(THIS_DIR, "results"))
OUT_ROOT <- Sys.getenv("OUT_DIR", unset = file.path(THIS_DIR, "plots"))
MAIN_OUT <- file.path(OUT_ROOT, "main_plots")
COV_OUT <- file.path(OUT_ROOT, "smoothing_target_coverage")
CI_OUT <- file.path(OUT_ROOT, "ci_length_ntr2500")
CACHE_DIR <- file.path(COV_OUT, "target_cache")
dir.create(MAIN_OUT, showWarnings = FALSE, recursive = TRUE)
dir.create(COV_OUT, showWarnings = FALSE, recursive = TRUE)
dir.create(CI_OUT, showWarnings = FALSE, recursive = TRUE)
dir.create(CACHE_DIR, showWarnings = FALSE, recursive = TRUE)

source(file.path(THIS_DIR, "estimator_functions.R"))
source(file.path(THIS_DIR, "simulation_functions.R"))
source(file.path(THIS_DIR, "plot_helpers.R"))

`%||%` <- function(x, y) if (!is.null(x)) x else y

EST_ALL <- c("PI", "OR", "CPI", "CPI_or",
             "BC", "BC_or", "PI_cali", "BC_cali", "BC_or_cali")
EST_ERR <- c("PI", "OR", "CPI", "BC", "PI_cali", "BC_cali")
EST_COV <- EST_ALL
BASIS <- c("fourier", "bs", "poly", "lpoly")
BASIS_LABELS <- c(fourier = "Cosine", bs = "B-spline",
                  poly = "Polynomial", lpoly = "Local linear")
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
EST_LABELS <- c(
  PI      = expression(bold(widehat(m)["pi"])),
  CPI     = expression(bold(widehat(m)[cpi])),
  BC      = expression(bold(widehat(m)[bc])),
  OR      = expression(bold(widehat(m)[or])),
  PI_cali = expression(bold(widehat(m)["pi"~"."~cali])),
  BC_cali = expression(bold(widehat(m)[bc~"."~cali])),
  CPI_or  = expression(bold(widehat(m)[cpi~"."~or])),
  BC_or   = expression(bold(widehat(m)[bc~"."~or])),
  BC_or_cali = expression(bold(widehat(m)[bc~"."~or~"."~cali]))
)
SCATTER_LABELS <- c(
  Truth = "Truth",
  PI = EST_LABELS[["PI"]],
  CPI = EST_LABELS[["CPI"]],
  BC = EST_LABELS[["BC"]]
)

parse_meta <- function(path) {
  b <- basename(path)
  info <- file.info(path)
  data.frame(
    setting = as.integer(sub(".*setting([0-9]+)_.*", "\\1", b)),
    g_type = sub(".*setting[0-9]+_([^_]+)_n.*", "\\1", b),
    n = as.integer(sub(".*_n([0-9]+)_ntr.*", "\\1", b)),
    n_tr = as.integer(sub(".*_ntr([0-9]+).*", "\\1", b)),
    file = b,
    path = path,
    mtime = info$mtime,
    stringsAsFactors = FALSE
  )
}

get_truth <- function(config) as.numeric(config$m.t(config$eval.pts))

get_se_mat_local <- function(results_list, est, basis) {
  vals <- lapply(results_list, function(res) {
    obj <- res[[est]][[basis]]
    if (is.null(obj)) return(rep(NA_real_, 0))
    se <- if (basis == "lpoly") obj$se %||% obj$se.us else obj$se.adj %||% obj$se
    if (is.null(se)) return(rep(NA_real_, 0))
    as.numeric(se)
  })
  do.call(rbind, vals)
}

get_main_plot_y_max <- function(item) {
  results_list <- Filter(Negate(is.null), item$res$result)
  truth <- get_truth(item$res$config)

  pointwise_vals <- unlist(lapply(EST_ERR, function(est) {
    as.numeric(get_pointwise_mse(results_list, est, "lpoly", truth))
  }), use.names = FALSE)

  avg_mse_vals <- unlist(lapply(EST_ALL, function(est) {
    unlist(lapply(BASIS, function(basis) {
      tryCatch(get_average_mse(results_list, est, basis, truth),
               error = function(e) numeric(0))
    }), use.names = FALSE)
  }), use.names = FALSE)

  vals <- c(pointwise_vals, avg_mse_vals)
  vals <- vals[is.finite(vals)]
  max(vals, na.rm = TRUE)
}

make_scatter_plot <- function(item) {
  config <- item$res$config
  results_list <- Filter(Negate(is.null), item$res$result)
  eval.pts <- config$eval.pts
  truth <- get_truth(config)
  res1 <- results_list[[1]]

  curve_df <- rbind(
    data.frame(t = eval.pts, value = truth, label = "Truth"),
    data.frame(t = eval.pts, value = as.numeric(res1$PI$lpoly$est), label = "PI"),
    data.frame(t = eval.pts, value = as.numeric(res1$CPI$lpoly$est), label = "CPI"),
    data.frame(t = eval.pts, value = as.numeric(res1$BC$lpoly$est), label = "BC")
  )
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
    geom_line(data = curve_df,
              aes(t, value, color = label, linewidth = label, alpha = label)) +
    scale_color_manual(values = c(Truth = "#B2182B", PI = "#FC8D62",
                                  CPI = "#1F78B4", BC = "#1B7837"),
                       labels = SCATTER_LABELS,
                       name = NULL) +
    scale_linewidth_manual(values = c(Truth = 1.05, PI = 0.75, CPI = 0.75, BC = 0.75),
                           guide = "none") +
    scale_alpha_manual(values = c(Truth = 1, PI = 0.95, CPI = 0.85, BC = 0.85),
                       guide = "none") +
    coord_cartesian(xlim = range(eval.pts), ylim = yr) +
    labs(x = expression(bold(t)), y = expression(bold(Y))) +
    theme_bw(base_size = 18) +
    theme(
      legend.position = "top",
      legend.text = element_text(size = 18, face = "bold"),
      legend.title = element_text(face = "bold"),
      legend.key.width = unit(1.7, "cm"),
      legend.key.height = unit(0.8, "cm"),
      axis.title = element_text(size = 18, face = "bold"),
      axis.text = element_text(size = 15, face = "bold")
    )
}

make_ptwise_err_plot <- function(item) {
  results_list <- Filter(Negate(is.null), item$res$result)
  config <- item$res$config
  eval.pts <- config$eval.pts
  truth <- get_truth(config)
  ptwise_vals <- unlist(lapply(EST_ERR, function(est) {
    as.numeric(get_pointwise_mse(results_list, est, "lpoly", truth))
  }), use.names = FALSE)
  ptwise_vals <- ptwise_vals[is.finite(ptwise_vals)]
  y_cap <- as.numeric(quantile(ptwise_vals, 0.99, na.rm = TRUE))
  if (!is.finite(y_cap) || y_cap <= 0) y_cap <- max(ptwise_vals, na.rm = TRUE)
  if (!is.finite(y_cap) || y_cap <= 0) y_cap <- 1

  make_plot_df <- function(err_mat) {
    nsim <- nrow(err_mat)
    as.data.frame(err_mat) %>%
      mutate(sim = seq_len(nsim)) %>%
      pivot_longer(cols = -sim, names_to = "t_index", values_to = "err") %>%
      mutate(t_index = as.integer(gsub("V", "", t_index)),
             t = eval.pts[t_index])
  }

  plot_panel <- function(plot_df, title_label) {
    mean_df <- plot_df %>%
      group_by(t) %>%
      summarise(mean_err = mean(err, na.rm = TRUE), .groups = "drop")
    ggplot() +
      geom_line(data = plot_df, aes(t, err, group = sim),
                color = "gray80", alpha = 0.18, linewidth = 0.2) +
      geom_line(data = mean_df, aes(t, mean_err), color = "red", linewidth = 0.95) +
      labs(title = title_label, x = "t", y = NULL) +
      coord_cartesian(xlim = range(eval.pts), ylim = c(0, y_cap)) +
      scale_y_continuous(breaks = pretty(c(0, y_cap), n = 3)) +
      theme_bw(base_size = 16) +
      theme(
        plot.title = element_text(size = 19, face = "bold", hjust = 0.5),
        axis.title.x = element_text(size = 15),
        axis.text = element_text(size = 13)
      )
  }

  plot_list <- lapply(EST_ERR, function(est) {
    err_mat <- get_pointwise_mse(results_list, est, "lpoly", truth)
    plot_panel(make_plot_df(err_mat), EST_LABELS[[est]])
  })
  (plot_list[[1]] | plot_list[[2]]) /
    (plot_list[[3]] | plot_list[[4]]) /
    (plot_list[[5]] | plot_list[[6]])
}

make_mse_plot <- function(item) {
  results_list <- Filter(Negate(is.null), item$res$result)
  truth <- get_truth(item$res$config)

  mse_rows <- list()
  for (est in EST_ALL) {
    for (basis in BASIS) {
      vals <- tryCatch(get_average_mse(results_list, est, basis, truth),
                       error = function(e) NULL)
      if (!is.null(vals)) {
        mse_rows[[length(mse_rows) + 1]] <- data.frame(
          estimator = est, basis = basis, mse = vals
        )
      }
    }
  }

  mse_df <- bind_rows(mse_rows) %>%
    mutate(estimator = factor(estimator, levels = EST_ALL),
           basis = factor(basis, levels = BASIS))

  lower_max <- max(mse_df$mse, na.rm = TRUE)
  if (!is.finite(lower_max) || lower_max <= 0) lower_max <- 1

  ggplot(mse_df, aes(estimator, mse, fill = estimator)) +
    geom_boxplot(outlier.size = 0.35, outlier.alpha = 0.35, linewidth = 0.35) +
    facet_wrap(~ basis, nrow = 1, scales = "fixed",
               labeller = labeller(basis = BASIS_LABELS)) +
    scale_fill_manual(values = PALETTE, drop = FALSE) +
    coord_cartesian(ylim = c(0, lower_max)) +
    scale_y_continuous(breaks = pretty(c(0, lower_max), n = 3)) +
    labs(x = NULL, y = "Average MSE") +
    theme_bw(base_size = 16) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 14, face = "bold"),
      axis.text.y = element_text(size = 12),
      axis.title.y = element_text(size = 15, face = "bold"),
      legend.position = "none",
      strip.text = element_text(face = "bold", size = 15)
    )
}

make_sieve_targets_one <- function(r_pop, mr_pop, eval.pts, df_seq, basis) {
  basis_str <- switch(basis, fourier = "fourier", bs = "bspline", poly = "polynomial")
  targets <- lapply(df_seq, function(k) {
    Phi <- cbind(1, basis_fun(r_pop, k, basis = basis_str))
    Phi_eval <- cbind(1, basis_fun(eval.pts, k, basis = basis_str))
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

get_sieve_targets <- function(config, tag, n_pop = 200000) {
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

get_lpoly_population <- function(config, tag, n_pop = 50000) {
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

make_lpoly_targets <- function(config, tag, h_values, n_pop = 50000) {
  pop <- get_lpoly_population(config, tag, n_pop = n_pop)
  h_key <- function(h) sprintf("%.3f", as.numeric(h))
  unique_h <- sort(unique(round(as.numeric(h_values), 3)))
  local_linear_gaussian <- function(r, y, eval_pts, h) {
    vapply(eval_pts, function(t0) {
      x <- r - t0
      w <- dnorm(x / h)
      s0 <- sum(w)
      s1 <- sum(w * x)
      s2 <- sum(w * x * x)
      y0 <- sum(w * y)
      y1 <- sum(w * x * y)
      denom <- s0 * s2 - s1 * s1
      if (!is.finite(denom) || abs(denom) < 1e-12) return(NA_real_)
      (s2 * y0 - s1 * y1) / denom
    }, numeric(1))
  }
  targets <- lapply(unique_h, function(h) {
    cache_file <- file.path(
      CACHE_DIR,
      sprintf("lpoly_target_%s_npop%d_h%s.rds",
              tag, n_pop, gsub("[.]", "p", h_key(h)))
    )
    if (file.exists(cache_file)) return(readRDS(cache_file))
    target <- local_linear_gaussian(pop$r, pop$m, config$eval.pts, h)
    saveRDS(target, cache_file)
    target
  })
  names(targets) <- h_key(unique_h)
  targets[h_key(round(as.numeric(h_values), 3))]
}

coverage_rows_one_file <- function(path) {
  out <- readRDS(path)
  meta <- parse_meta(path)
  results_list <- Filter(Negate(is.null), out$result)
  config <- out$config
  fixed_tuning <- out$fixed_tuning
  eval_pts <- config$eval.pts
  tag <- sprintf("setting%d_%s", meta$setting, meta$g_type)
  sieve_targets <- get_sieve_targets(config, tag)

  get_fixed_dfh <- function(estimator, basis) {
    if (is.null(fixed_tuning)) return(NA_real_)
    idx <- which(fixed_tuning$estimator == estimator & fixed_tuning$basis == basis)
    if (length(idx) == 0) return(NA_real_)
    val <- as.numeric(fixed_tuning$fixed_value[idx[1]])
    if (!is.finite(val)) return(NA_real_)
    if (basis %in% c("fourier", "poly", "bs")) val <- as.integer(round(val))
    val
  }

  bind_rows(lapply(EST_COV, function(est) {
    bind_rows(lapply(BASIS, function(basis) {
      tryCatch({
        est_mat <- get_preds_mat(results_list, est, basis)
        se_mat <- get_se_mat_local(results_list, est, basis)
        if (nrow(est_mat) == 0 || nrow(se_mat) == 0) return(NULL)
        fixed_dfh <- get_fixed_dfh(est, basis)
        if (!is.finite(fixed_dfh)) {
          message("Missing fixed df/h for ", basename(path), " / ", est, " / ", basis,
                  "; falling back to realized df.or.h.")
          df_h <- get_df_or_h(results_list, est, basis)
        } else {
          df_h <- rep(fixed_dfh, nrow(est_mat))
        }

        if (basis == "lpoly") {
          target_mat <- do.call(rbind, make_lpoly_targets(config, tag, df_h))
        } else {
          target_mat <- do.call(rbind, lapply(as.character(df_h), function(k) {
            sieve_targets[[basis]][[k]]
          }))
        }

        cov_vec <- colMeans(est_mat - 1.96 * se_mat <= target_mat &
                              target_mat <= est_mat + 1.96 * se_mat,
                            na.rm = TRUE)
        data.frame(
          setting = meta$setting,
          g_type = meta$g_type,
          n = meta$n,
          n_tr = meta$n_tr,
          basis = basis,
          estimator = est,
          t = eval_pts,
          coverage = cov_vec,
          file = meta$file,
          stringsAsFactors = FALSE
        )
      }, error = function(e) {
        message("Coverage skip ", basename(path), " / ", est, " / ", basis, ": ",
                conditionMessage(e))
        NULL
      })
    }))
  }))
}

get_se_vec <- function(res, est, basis) {
  obj <- res[[est]][[basis]]
  if (is.null(obj)) return(NULL)
  se <- if (basis == "lpoly") obj$se %||% obj$se.us else obj$se.adj %||% obj$se
  if (is.null(se)) return(NULL)
  as.numeric(se)
}

ci_rows_one_file <- function(path) {
  out <- readRDS(path)
  meta <- parse_meta(path)
  results_list <- Filter(Negate(is.null), out$result)
  eval_pts <- out$config$eval.pts

  bind_rows(lapply(seq_along(results_list), function(rep_id) {
    res <- results_list[[rep_id]]
    bind_rows(lapply(BASIS, function(basis) {
      bind_rows(lapply(EST_ALL, function(est) {
        se <- tryCatch(get_se_vec(res, est, basis), error = function(e) NULL)
        if (is.null(se)) return(NULL)
        data.frame(
          setting = meta$setting,
          g_type = meta$g_type,
          n = meta$n,
          n_tr = meta$n_tr,
          rep = rep_id,
          basis = basis,
          estimator = est,
          t = eval_pts,
          ci_length = 2 * 1.96 * se,
          file = meta$file,
          stringsAsFactors = FALSE
        )
      }))
    }))
  }))
}

cv_files <- list.files(
  CV_RESULT_DIR,
  pattern = "^cv_setting[0-9]+_(null|linear)_n2500_ntr2500_nsim500_.*[.]rds$",
  full.names = TRUE
)
if (length(cv_files) == 0) stop("No CV n.tr=2500 result files found in ", CV_RESULT_DIR)

cv_meta <- bind_rows(lapply(cv_files, parse_meta)) %>%
  group_by(setting, g_type, n, n_tr) %>%
  arrange(desc(mtime), .by_group = TRUE) %>%
  slice(1) %>%
  ungroup() %>%
  arrange(g_type, setting, n_tr)

write.csv(cv_meta, file.path(MAIN_OUT, "selected_cv_main_plot_files.csv"), row.names = FALSE)

message("Using ", nrow(cv_meta), " CV n.tr=2500 result file(s) for main plots")

for (i in seq_len(nrow(cv_meta))) {
  f <- cv_meta$path[i]
  setting <- cv_meta$setting[i]
  g_type <- cv_meta$g_type[i]
  n_tr <- cv_meta$n_tr[i]
  tag <- sprintf("setting%d_%s_ntr%d", setting, g_type, n_tr)
  message("Building main plot: ", tag)

  res_obj <- readRDS(f)
  item <- list(setting = setting, g_type = g_type, n_tr = n_tr, res = res_obj)

  p_scatter <- make_scatter_plot(item)
  p_err <- make_ptwise_err_plot(item)
  p_mse <- make_mse_plot(item)

  top_row <- (p_scatter | wrap_elements(full = p_err)) +
    plot_layout(widths = c(1, 2.1))
  combined <- (top_row / wrap_elements(p_mse)) +
    plot_layout(heights = c(1.05, 1))

  pdf_file <- file.path(MAIN_OUT, sprintf("main_%s.pdf", tag))
  png_file <- file.path(MAIN_OUT, sprintf("main_%s.png", tag))
  ggsave(pdf_file, combined, width = 18.5, height = 12, bg = "white")
  ggsave(png_file, combined, width = 18.5, height = 12, dpi = 300, bg = "white")
  message("  saved: ", pdf_file)
}

message("All CV main plots saved to: ", MAIN_OUT)

fixed_files <- list.files(
  FIXED_RESULT_DIR,
  pattern = "^fixed_setting[0-9]+_(null|linear)_n2500_ntr[0-9]+_pilot50_nsim500_.*[.]rds$",
  full.names = TRUE
)
if (length(fixed_files) == 0) {
  stop("No fixed pilot result files found in ", FIXED_RESULT_DIR)
}

fixed_meta <- bind_rows(lapply(fixed_files, parse_meta)) %>%
  group_by(setting, g_type, n, n_tr) %>%
  arrange(desc(mtime), .by_group = TRUE) %>%
  slice(1) %>%
  ungroup() %>%
  arrange(g_type, setting, n_tr)

write.csv(fixed_meta, file.path(COV_OUT, "selected_fixed_coverage_files.csv"), row.names = FALSE)
message("Using ", nrow(fixed_meta), " pilot-fixed result file(s) for coverage")

cov_df <- bind_rows(lapply(fixed_meta$path, coverage_rows_one_file)) %>%
  mutate(
    tag = paste0("setting", setting, "_", g_type),
    ntr_f = factor(paste0("n.tr = ", n_tr),
                   levels = paste0("n.tr = ", sort(unique(fixed_meta$n_tr)))),
    basis_f = factor(BASIS_LABELS[basis], levels = unname(BASIS_LABELS[BASIS])),
    estimator = factor(estimator, levels = EST_COV)
  )

write.csv(cov_df, file.path(COV_OUT, "smoothing_target_coverage_pointwise.csv"),
          row.names = FALSE)

summary_df <- cov_df %>%
  group_by(setting, g_type, n, n_tr, basis, estimator) %>%
  summarise(
    mean_coverage = mean(coverage, na.rm = TRUE),
    min_coverage = min(coverage, na.rm = TRUE),
    max_coverage = max(coverage, na.rm = TRUE),
    .groups = "drop"
  )
write.csv(summary_df, file.path(COV_OUT, "smoothing_target_coverage_summary.csv"),
          row.names = FALSE)

for (tag_i in sort(unique(cov_df$tag))) {
  plot_df <- filter(cov_df, tag == tag_i)
  p <- ggplot(plot_df, aes(t, coverage, color = estimator)) +
    geom_hline(yintercept = 0.95, linetype = "dashed",
               color = "gray35", linewidth = 0.45) +
    geom_line(linewidth = 0.6, alpha = 0.95) +
    facet_grid(ntr_f ~ basis_f) +
    scale_color_manual(values = PALETTE, labels = EST_LABELS,
                       drop = FALSE, name = NULL) +
    scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 0.95, 1)) +
    labs(x = "t", y = "Coverage") +
    theme_bw(base_size = 15) +
    theme(
      legend.position = "bottom",
      legend.text = element_text(size = 15, face = "bold"),
      legend.title = element_text(face = "bold"),
      legend.key.width = grid::unit(2.0, "lines"),
      legend.key.height = grid::unit(1.1, "lines"),
      strip.text = element_text(face = "bold", size = 14),
      axis.title = element_text(size = 15, face = "bold"),
      axis.text = element_text(size = 13)
    ) +
    guides(color = guide_legend(nrow = 2))

  pdf_file <- file.path(COV_OUT, paste0("coverage_smoothing_target_", tag_i, ".pdf"))
  png_file <- file.path(COV_OUT, paste0("coverage_smoothing_target_", tag_i, ".png"))
  ggsave(pdf_file, p, width = 15, height = 8, bg = "white")
  ggsave(png_file, p, width = 15, height = 8, dpi = 300, bg = "white")
  message("Saved: ", pdf_file)
}

ci_meta <- fixed_meta %>%
  filter(n_tr == 2500)
if (nrow(ci_meta) == 0) {
  stop("No n.tr=2500 fixed result files found in ", FIXED_RESULT_DIR)
}

write.csv(ci_meta, file.path(CI_OUT, "selected_ci_length_files.csv"), row.names = FALSE)
message("Using ", nrow(ci_meta), " pilot-fixed n.tr=2500 result file(s) for CI length")

ci_df <- bind_rows(lapply(ci_meta$path, ci_rows_one_file)) %>%
  mutate(
    tag = paste0("setting", setting, "_", g_type),
    basis_f = factor(BASIS_LABELS[basis], levels = unname(BASIS_LABELS[BASIS])),
    estimator = factor(estimator, levels = EST_ALL)
  )

ci_pointwise <- ci_df %>%
  group_by(setting, g_type, n, n_tr, basis, estimator, t, file) %>%
  summarise(
    mean_ci_length = mean(ci_length, na.rm = TRUE),
    median_ci_length = median(ci_length, na.rm = TRUE),
    q90_ci_length = quantile(ci_length, 0.90, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    tag = paste0("setting", setting, "_", g_type),
    basis_f = factor(BASIS_LABELS[basis], levels = unname(BASIS_LABELS[BASIS])),
    estimator = factor(estimator, levels = EST_ALL)
  )

ci_summary <- ci_df %>%
  group_by(setting, g_type, n, n_tr, basis, estimator) %>%
  summarise(
    mean_ci_length = mean(ci_length, na.rm = TRUE),
    median_ci_length = median(ci_length, na.rm = TRUE),
    q90_ci_length = quantile(ci_length, 0.90, na.rm = TRUE),
    max_ci_length = max(ci_length, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(ci_pointwise, file.path(CI_OUT, "ci_length_pointwise_ntr2500.csv"), row.names = FALSE)
write.csv(ci_summary, file.path(CI_OUT, "ci_length_summary_ntr2500.csv"), row.names = FALSE)

for (tag_i in sort(unique(ci_pointwise$tag))) {
  plot_df <- filter(ci_pointwise, tag == tag_i)
  p <- ggplot(plot_df, aes(x = t, y = mean_ci_length, color = estimator)) +
    geom_line(linewidth = 0.6, alpha = 0.95) +
    facet_wrap(~ basis_f, nrow = 1, ncol = 4, scales = "free_y") +
    scale_color_manual(values = PALETTE, labels = EST_LABELS,
                       drop = FALSE, name = NULL) +
    labs(x = "t", y = "Mean 95% CI length") +
    theme_bw(base_size = 15) +
    theme(
      legend.position = "bottom",
      legend.text = element_text(size = 15, face = "bold"),
      legend.title = element_text(face = "bold"),
      legend.key.width = grid::unit(2.0, "lines"),
      legend.key.height = grid::unit(1.1, "lines"),
      strip.text = element_text(face = "bold", size = 14),
      axis.title = element_text(size = 15, face = "bold"),
      axis.text = element_text(size = 13)
    ) +
    guides(color = guide_legend(nrow = 2))

  pdf_file <- file.path(CI_OUT, paste0("ci_length_", tag_i, "_ntr2500.pdf"))
  png_file <- file.path(CI_OUT, paste0("ci_length_", tag_i, "_ntr2500.png"))
  ggsave(pdf_file, p, width = 18, height = 5, bg = "white")
  ggsave(png_file, p, width = 18, height = 5, dpi = 300, bg = "white")
  message("Saved: ", pdf_file)
}

message("All outputs saved to: ", OUT_ROOT)
