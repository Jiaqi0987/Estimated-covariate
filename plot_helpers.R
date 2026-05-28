
get_est_vector <- function(res, est, basis, sim = NULL) {
  if (!is.null(sim)) {
    res <- res[[sim]]
  }
  
  if (!est %in% names(res)) {
    stop(sprintf("Estimator '%s' not found.", est))
  }
  if (!basis %in% names(res[[est]])) {
    stop(sprintf("Basis '%s' not found under estimator '%s'.", basis, est))
  }
  
  obj <- res[[est]][[basis]]
  
  if (is.null(obj$est)) {
    stop(sprintf("No $est found for estimator '%s', basis '%s'.", est, basis))
  }
  
  x <- obj$est
  
  if (is.data.frame(x)) {
    x <- as.matrix(x)
  }
  
  if (is.matrix(x)) {
    if (ncol(x) == 1) {
      x <- x[, 1]
    } else if (nrow(x) == 1) {
      x <- x[1, ]
    } else {
      stop(sprintf(
        "est for '%s' / '%s' is a %d x %d matrix, not a vector.",
        est, basis, nrow(x), ncol(x)
      ))
    }
  }
  
  x <- as.numeric(x)
  
  if (anyNA(x)) {
    stop(sprintf("est for '%s' / '%s' contains NA.", est, basis))
  }
  
  x
}

# Extract the tuning parameter (df or h) for a given estimator and basis
# from each data.frame in the result list.
get_df_or_h <- function(results_list, est, basis) {
  
  vals <- sapply(results_list, function(res) {
    
    if (!est %in% names(res)) {
      stop(sprintf("Estimator '%s' not found.", est))
    }
    
    if (!basis %in% names(res[[est]])) {
      stop(sprintf("Basis '%s' not found under estimator '%s'.", basis, est))
    }
    
    obj <- res[[est]][[basis]]
    
    if (!is.null(obj$df.or.h)) {
      return(obj$df.or.h)
    } else {
      stop(sprintf(
        "Neither 'df' nor 'h' found for estimator '%s' and basis '%s'.",
        est, basis
      ))
    }
  })
  
  return(vals)
}

get_pointwise_mse <- function(results_list, est, basis, truth) {
  truth <- as.numeric(truth)
  
  preds_list <- lapply(seq_along(results_list), function(i) {
    x <- get_est_vector(results_list, est, basis, sim = i)
    
    if (length(x) != length(truth)) {
      stop(sprintf(
        "Length mismatch at simulation %d for %s / %s: est has length %d, truth has length %d.",
        i, est, basis, length(x), length(truth)
      ))
    }
    
    x
  })
  
  preds_mat <- do.call(rbind, preds_list)
  sweep(preds_mat, 2, truth, FUN = "-")^2
}

get_average_mse <- function(results_list, est, basis, truth) {
  mse_mat <- get_pointwise_mse(results_list, est, basis, truth)
  rowMeans(mse_mat)
}
get_df_or_h <- function(results_list, est, basis) {
  sapply(results_list, function(res) {
    obj <- res[[est]][[basis]]
    if (!is.null(obj$df)) return(obj$df)
    if (!is.null(obj$h)) return(obj$h)
    stop(sprintf("Neither df nor h found for estimator '%s', basis '%s'.", est, basis))
  })
}
make_plot_mse <- function(results_list, truth,
                          ylab = "",
                          palette = c("#1b9e77", "#d95f02", "#7570b3",
                                      "#e7298a", "#66a61e", "#e6ab02",
                                      "pink", "brown", "grey"),
                          ymax = NULL,
                          coef = 1.5) {
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  
  estimators <- c("PI", "OR", "CPI", "CPI_or",
                  "BC", "BC_or", "PI_cali", "BC_cali", "BC_or_cali")
  bases <- c("fourier", "bs", "poly", "lpoly")
  
  basis_labels <- c(
    fourier = "Fourier",
    bs      = "B-spline",
    poly    = "Polynomial",
    lpoly   = "Local linear"
  )
  
  rows <- list()
  k <- 1
  
  for (est in estimators) {
    for (basis in bases) {
      tmp <- tryCatch({
        mse_vec <- get_average_mse(results_list, est, basis, truth)
        data.frame(
          estimator = est,
          basis = basis,
          sim = seq_along(mse_vec),
          mse = mse_vec
        )
      }, error = function(e) {
        message(sprintf("Skipping %s / %s: %s", est, basis, e$message))
        NULL
      })
      
      if (!is.null(tmp)) {
        rows[[k]] <- tmp
        k <- k + 1
      }
    }
  }
  
  mse_df <- bind_rows(rows) %>%
    mutate(
      basis = factor(basis, levels = bases),
      estimator = factor(estimator, levels = estimators)
    )
  
  mse_report <- mse_df %>%
    group_by(basis, estimator) %>%
    summarize(
      mse_mean = mean(mse, na.rm = TRUE),
      mse_sd   = sd(mse, na.rm = TRUE),
      mse_max  = max(mse, na.rm = TRUE),
      .groups = "drop"
    )
  
  print(mse_report, n = Inf)
  
  # helper: compute boxplot stats manually
  get_box_stats <- function(x, coef = 1.5) {
    x <- x[is.finite(x)]
    if (length(x) == 0) {
      return(data.frame(
        ymin = NA, lower = NA, middle = NA,
        upper = NA, ymax = NA
      ))
    }
    
    qs <- quantile(x, probs = c(0.25, 0.5, 0.75), na.rm = TRUE, names = FALSE)
    iqr <- qs[3] - qs[1]
    lower_whisker <- min(x[x >= qs[1] - coef * iqr], na.rm = TRUE)
    upper_whisker <- max(x[x <= qs[3] + coef * iqr], na.rm = TRUE)
    
    data.frame(
      ymin   = lower_whisker,
      lower  = qs[1],
      middle = qs[2],
      upper  = qs[3],
      ymax   = upper_whisker
    )
  }
  
  box_df <- mse_df %>%
    group_by(basis, estimator) %>%
    group_modify(~ get_box_stats(.x$mse, coef = coef)) %>%
    ungroup() %>%
    mutate(
      basis = factor(basis, levels = bases),
      estimator = factor(estimator, levels = estimators)
    )
  
  # per-basis ymax (whisker max × 1.15) so each panel uses its own scale
  ymax_per_basis <- sapply(bases, function(b) {
    vals <- box_df$ymax[box_df$basis == b]
    max(vals, na.rm = TRUE) * 1.15
  })

  draw_panel <- function(b) {
    ymax_b    <- ymax_per_basis[[b]]
    raw_b     <- filter(mse_df, basis == b)
    box_b     <- filter(box_df, basis == b)

    ggplot() +
      geom_jitter(data = raw_b,
                  aes(x = estimator, y = mse, color = estimator),
                  width = 0.25, alpha = 0.25, size = 0.6, show.legend = FALSE) +
      geom_boxplot(data = box_b,
                   aes(x = estimator,
                       ymin = ymin, lower = lower, middle = middle,
                       upper = upper, ymax = ymax,
                       fill = estimator),
                   stat = "identity", width = 0.55, alpha = 0.7,
                   outlier.shape = NA) +
      scale_y_continuous(limits = c(0, ymax_b)) +
      scale_fill_manual(values  = palette) +
      scale_color_manual(values = palette) +
      labs(title = basis_labels[[b]], x = NULL, y = ylab) +
      theme_bw() +
      theme(
        plot.title    = element_text(hjust = 0.5, face = "bold"),
        axis.text.x   = element_text(size = 14, angle = 45, hjust = 1),
        axis.text.y   = element_text(size = 16),
        axis.title    = element_text(size = 20),
        legend.position = "none"
      )
  }
  
  wrap_plots(lapply(bases, draw_panel), nrow = 1)
}
get_preds_mat <- function(results_list, est, basis) {
  preds_list <- lapply(results_list, function(res) {
    x <- res[[est]][[basis]]$est
    
    if (is.matrix(x)) {
      if (ncol(x) == 1) x <- x[,1]
      else if (nrow(x) == 1) x <- x[1,]
      else stop("est not vector")
    }
    
    as.numeric(x)
  })
  
  do.call(rbind, preds_list)
}
get_bias_var <- function(results_list, est, basis, truth) {
  preds_mat <- get_preds_mat(results_list, est, basis)
  
  mean_est <- colMeans(preds_mat)
  
  bias <- mean_est - truth
  var  <- apply(preds_mat, 2, var)
  mse  <- bias^2 + var
  
  return(list(
    bias = bias,
    var  = var,
    mse  = mse,
    mean = mean_est
  ))
}

if (FALSE) {
  out_PI  <- get_bias_var(res$result, "CPI", "lpoly", truth = truth)
  out_CPI <- get_bias_var(res$result, "CPI_or", "lpoly", truth)

  plot(eval.pts, out_PI$var, type = "l", col = "blue", ylim = range(c(out_PI$var, out_CPI$var)))
  lines(eval.pts, out_CPI$var, col = "red")
  legend("topright", legend = c("PI", "CPI"), col = c("blue", "red"), lty = 1)

  plot(eval.pts, out_PI$bias^2, type = "l", col = "blue", ylim = range(c(out_PI$bias^2, out_CPI$bias^2)))
  lines(eval.pts, out_CPI$bias^2, col = "red")
  legend("topright", legend = c("PI", "CPI"), col = c("blue", "red"), lty = 1)

  plot(eval.pts, out_PI$mse, type = "l", col = "blue", ylim = range(c(out_PI$mse, out_CPI$mse)))
  lines(eval.pts, out_CPI$mse, col = "red")
  legend("topright", legend = c("PI", "CPI"), col = c("blue", "red"), lty = 1)
}



get_preds_at_t <- function(results_list, est, basis, eval.pts, t0) {
  j <- which.min(abs(eval.pts - t0))
  
  vals <- sapply(results_list, function(res) {
    x <- res[[est]][[basis]]$est
    if (is.matrix(x)) {
      if (ncol(x) == 1) x <- x[, 1]
      else if (nrow(x) == 1) x <- x[1, ]
      else stop("est is not a vector")
    }
    as.numeric(x)[j]
  })
  
  return(list(values = vals, index = j, t = eval.pts[j]))
}

plot_pointwise_density <- function(results_list, basis, eval.pts, truth = NULL,
                                   t_grid = c(0.25, 0.5, 0.75),
                                   estimators = c("PI", "CPI", "BC")) {
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  
  par(mfrow = c(1, length(t_grid)))
  
  cols <- c("blue", "red", "darkgreen")
  names(cols) <- estimators
  
  for (t0 in t_grid) {
    vals_list <- lapply(estimators, function(est) {
      get_preds_at_t(results_list, est, basis, eval.pts, t0)$values
    })
    names(vals_list) <- estimators
    
    dens_list <- lapply(vals_list, density)
    xr <- range(sapply(dens_list, function(d) range(d$x)))
    yr <- range(sapply(dens_list, function(d) range(d$y)))
    
    j <- which.min(abs(eval.pts - t0))
    plot(NA, xlim = xr, ylim = yr,
         xlab = "estimate", ylab = "density",
         main = paste0("t = ", round(eval.pts[j], 3)))
    
    for (est in estimators) {
      lines(dens_list[[est]], col = cols[est], lwd = 2)
    }
    
    if (!is.null(truth)) {
      abline(v = truth[j], col = 2, lwd = 2, lty = 2)
    }
    
    legend("topright", legend = estimators, col = cols[estimators],
           lwd = 2, bty = "n")
  }
}



plot_pointwise_err = function(results_list, basis = "lpoly", r ,ymax) {
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  
  rslt    <- results_list$result
  config  <- results_list$config
  eval.pts <- config$eval.pts
  truth   <- config$m.t(eval.pts)
  
  
  estimators <- c("PI", "OR", "CPI", "BC", "BC_cali","PI_cali")
  
  # -------------------------------------------------
  # helper: compute pointwise squared error 500×50
  # -------------------------------------------------
  get_err_mat <- function(est) {
    get_pointwise_mse(rslt, est = est, basis = basis, truth = truth)
  }
  
  # -------------------------------------------------
  # helper: convert to long format
  # -------------------------------------------------
  make_plot_df <- function(err_mat) {
    nsim <- nrow(err_mat)
    err_mat %>%
      as.data.frame() %>%
      mutate(sim = 1:nsim) %>%
      pivot_longer(cols = -sim, names_to = "t_index", values_to = "err") %>%
      mutate(
        t_index = as.integer(gsub("V", "", t_index)),
        t = eval.pts[t_index]
      )
  }
  
  # -------------------------------------------------
  # global y-range across all 4 estimators
  # -------------------------------------------------
  all_err <- do.call(rbind, lapply(estimators, get_err_mat))
  global_max <- max(all_err, na.rm = TRUE)
  
  # -------------------------------------------------
  # panel plot with background density of r
  # -------------------------------------------------
  plot_panel <- function(plot_df, title) {
    mean_df <- plot_df %>% group_by(t) %>% summarise(mean_err = mean(err), .groups = "drop")
    
    # 1. density of r
    hist_r <- hist(r, breaks = 30, plot = FALSE)
    hist_df <- data.frame(
      t = hist_r$mids,          # midpoints of bins
      d = hist_r$density
    )
    hist_df[hist_df$t<0.85& hist_df$t>0.15,]
    scale_factor <- (global_max) / max(hist_df$d)
    
    ggplot() +
      
      # spaghetti
      geom_line(data = plot_df,
                aes(x = t, y = err, group = sim),
                color = "gray80", alpha = 0.20) +
      # mean
      geom_line(data = mean_df,
                aes(x = t, y = mean_err),
                color = "red", linewidth = 0.8,alpha = 0.8) +
      # background density
      #geom_col(
      #  data = hist_df,
      #  aes(x = t, y = d * scale_factor),
      #  width = diff(hist_r$breaks)[1],   # bin width
      #  fill = "skyblue", alpha = 0.25, color = NA
      #) +
      labs(
        title = title,
        x = "t",
        y = ""
      ) +
      coord_cartesian(
        xlim = range(eval.pts),   # force x-range matches eval.pts
        ylim = c(0, 1)
      )+
      scale_y_continuous(breaks = c(0, 0.5, 1)) +
      theme_bw(base_size = 20) +
      theme(
        plot.title      = element_text(size = 35, face = "bold", hjust = 0.5),
        axis.title.x    = element_text(size = 35, face = "bold", hjust = 0.5),
        axis.title.y    = element_text(size = 20),
        axis.text.x     = element_text(size = 20),
        axis.text.y     = element_text(size = 20),
      )
  }
  
  
  # -------------------------------------------------
  # build all 6 panels
  # -------------------------------------------------
  plot_list <- lapply(estimators, function(est) {
    err_mat  <- get_err_mat(est)
    plot_df  <- make_plot_df(err_mat)
    plot_panel(plot_df, est)
  })
  
  # -------------------------------------------------
  # assemble into 2×2
  # -------------------------------------------------
  (plot_list[[1]] | plot_list[[2]]) /
    (plot_list[[3]] | plot_list[[4]])/
    (plot_list[[6]] | plot_list[[5]])
}
