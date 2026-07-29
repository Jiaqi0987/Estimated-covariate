# simulation_fns_parallel_Feb4_pilot_fixed.R
# Simulation functions used by run_simulation.R.
# Supports two-stage tuning: pilot CV choices first, then fixed median df/h in final run.

library(nprobust)
`%||%` <- function(x, y) if (!is.null(x)) x else y
library(splines)
library(splines2)
library(gbm)
library(caret)
library(glmnet)
library(mvtnorm)
library(mgcv)
library(future.apply)

expit <- function(u) 1/(exp(-u) + 1)
logit <- function(u) log(u/(1-u))

# functions to compute the true E(g(x)|r(X) = t) and its derivative when \rho!= 0
truth_fun <- function(u, g_type, beta_r, beta_g, Sigma) {

  if (g_type == "interact"){
    den <- as.numeric(t(beta_r) %*% Sigma %*% beta_r)
    v   <- as.vector(Sigma %*% beta_r)                  # Σβ
    v12 <- v[1] * v[2]
    c12 <- Sigma[1, 2]
    s   <- logit(u)
   beta_g[1]*(c12 - v12/den + (v12/den^2) * s^2)
  }else if (g_type =="linear"){
    den <- as.numeric(t(beta_r) %*% Sigma %*% beta_r)  
    lambda <- as.numeric(t(beta_g) %*% Sigma %*% beta_r / den)
    lambda * logit(u)
  }else if (g_type =="null"){
    0
  }
}

deriv_fun <- function(t, g_type, beta_r, beta_g, Sigma){
  if(g_type =="interact"){
    den <- as.numeric(t(beta_r) %*% Sigma %*% beta_r)
    v   <- as.vector(Sigma %*% beta_r)
    v12 <- v[1] * v[2]
    s   <- logit(t)
    2 * beta_g[1]*(v12/den^2) * s * dlogit(t)
  }else if (g_type =="linear"){
    den <- as.numeric(t(beta_r) %*% Sigma %*% beta_r)
    lambda <- as.numeric(t(beta_g) %*% Sigma %*% beta_r / den)
    lambda * dlogit(t)
  }else if(g_type =="null"){
    0
  }
}


logit <- function(x) log(x / (1 - x))
dlogit <- function(p) {
  if (any(p <= 0 | p >= 1)) {
    stop("p must be strictly between 0 and 1")
  }
  return(1 / (p * (1 - p)))
}


gen_data <- function(beta_r, beta_g, n, d, Sigma, r.x, mu.x) {
  x <- rmvnorm(n, mean = rep(0, d), sigma = Sigma)
  r <- r.x(beta_r,x)
  a <- rbinom(n, p = r, size = 1)
  y <- rnorm(n, mean = mu.x(x), sd = 1)
  dat <- data.frame(y = y, a = a, r = r)
  return(list(dat = dat, x = x))
}

simu_rslt = function(seed, config){
  
  set.seed(seed)
  # extract everything from config
  n = config$n
  n.tr = config$n.tr
  df.seq.bs = config$df.seq.bs
  df.seq.poly = config$df.seq.poly
  df.seq.fourier = config$df.seq.fourier
  g_type = config$g_type
  der.m.t = config$der.m.t
  m.t = config$m.t
  r.x = config$r.x
  g.x = config$g.x
  mu.x = config$mu.x
  beta_g = config$beta_g
  beta_r = config$beta_r
  Sigma = config$Sigma
  d = config$d
  poly_deg = config$poly_deg
  eval.pts = config$eval.pts
 # rate = config$rate
 # est_r = config$est_r
  h.seq = config$h.seq
  fix.df.or.h = isTRUE(config$fix.df.or.h)

  # When fix.df.or.h = TRUE, allow each estimator/basis to use its own
  # pilot-selected median df/h. config$fixed_tuning should be a data.frame
  # with columns: estimator, basis, fixed_value.
  fixed_tuning <- config$fixed_tuning

  get_tuning_value <- function(estimator, basis, default_value) {
    if (!fix.df.or.h) return(default_value)
    if (is.null(fixed_tuning)) return(default_value)

    idx <- which(fixed_tuning$estimator == estimator & fixed_tuning$basis == basis)
    if (length(idx) == 0) return(default_value)

    val <- fixed_tuning$fixed_value[idx[1]]
    if (is.na(val) || !is.finite(as.numeric(val))) return(default_value)

    if (basis %in% c("fourier", "poly", "bs")) {
      return(as.integer(round(as.numeric(val))))
    }
    as.numeric(val)
  }

  df_fourier_for <- function(estimator) get_tuning_value(estimator, "fourier", df.seq.fourier)
  df_poly_for    <- function(estimator) get_tuning_value(estimator, "poly",    df.seq.poly)
  df_bs_for      <- function(estimator) get_tuning_value(estimator, "bs",      df.seq.bs)
  h_for          <- function(estimator) get_tuning_value(estimator, "lpoly",   h.seq)
  bandwidth_mode <- config$bandwidth_mode %||% if (fix.df.or.h) "pilot_fixed" else "cv"
  fixed_h <- suppressWarnings(as.numeric(config$fixed_h %||% NA_real_))

  lpoly_h_arg <- function(estimator, oracle_h = NULL) {
    if (bandwidth_mode == "oracle" && !is.null(oracle_h) && is.finite(oracle_h)) {
      return(oracle_h)
    }
    if (bandwidth_mode == "fixed") {
      if (!is.finite(fixed_h) || fixed_h <= 0) {
        stop("bandwidth_mode='fixed' requires config$fixed_h > 0")
      }
      return(fixed_h)
    }
    if (fix.df.or.h) {
      return(h_for(estimator))
    }
    NULL
  }

  fit_lprobust <- function(y, x, estimator, oracle_h = NULL) {
    h_use <- lpoly_h_arg(estimator, oracle_h = oracle_h)
    if (is.null(h_use)) {
      return(lprobust(y = y, x = x, p = poly_deg, eval = eval.pts,
                      kernel = "gau"))
    }
    lprobust(y = y, x = x, p = poly_deg, eval = eval.pts,
             kernel = "gau", h = h_use)
  }

  bc_hseq_for <- function(estimator, oracle_h = NULL) {
    h_use <- lpoly_h_arg(estimator, oracle_h = oracle_h)
    if (is.null(h_use)) return(h.seq)
    h_use
  }

  # generate sample for inference
  inference_data <- gen_data(beta_r = beta_r, beta_g = beta_g, n = n, d = d, Sigma = Sigma,
                             r.x = r.x, mu.x = mu.x)
  
  # generate auxiliary sample to estimate r.x of size 
  train_data <- gen_data(beta_r = beta_r, beta_g = beta_g, n = n.tr, d = d, Sigma = Sigma,
                          r.x = r.x, mu.x = mu.x)
  

  dat <- inference_data$dat
  y <- dat$y; r <- dat$r; a <- dat$a; x <- inference_data$x
  y.tr <- train_data$dat$y; a.tr <- train_data$dat$a; x.tr <- train_data$x
  # estimate rhat or generate rhat given the convergence rate
  #if(est_r == TRUE){
    rhat_fit <- train_rhat(a=a.tr, x=x.tr, new.x=x)
    rhat <- rhat_fit$rhat
    rhat.tr <- rhat_fit$rhat_train
  #}else{
  #  rhat <- expit(logit(r) + rnorm(n, mean =1 / (n)^rate, sd = 1 / (n)^rate)*sample(c(-1,1),size = n, replace = T))
  #  rhat.tr= expit(logit(train_data$dat$r) + rnorm(n.tr, mean = 1 / (n.tr)^rate, sd = 1 / (n.tr)^rate)*sample(c(-1,1),size = n.tr, replace = T))
  #}
  
  #calibrated rhat
  rhat_star_fit = isoreg(y = a, x = rhat)
  rhat_star  = numeric(length(a))
  rhat_star[ order(rhat)] = rhat_star_fit$yf
  eps <- 1e-6
  rhat_cali <- pmin(pmax(rhat_star, eps), 1 - eps)

  # estimate derivative
  eval.pts2 <- seq(0, 1, length.out=150)
  tmp<- lprobust(y = train_data$dat$y, x = rhat.tr, eval = eval.pts2,
                      deriv = 1, kernel = "gau")$Estimate[, "tau.us"]
  #estimated derivative for rhat
  der.est <- approx(x = eval.pts2, y = tmp, xout = rhat,      rule = 2)$y
  #estimated derivative for calibrated rhat
  der.est_cali = approx(x = eval.pts2, y = tmp, xout = rhat_cali, rule = 2)$y
  # estimtate muhat
  muhat_fit <- train_muhat(y=y.tr, x=x.tr, new.x=x)
  muhat <- muhat_fit$muhat
  mu <- as.vector(mu.x(x))

  ## oracle estimator knowing r.x(x)
  or.fourier <-  cpi.spline(y = y, a = a, rhat = r, eval = eval.pts, 
                            der.est = 0, df.seq = df_fourier_for("OR"), basis = "fourier")
  or.poly <-  cpi.spline(y = y, a = a, rhat = r, eval = eval.pts, 
                         der.est = 0, df.seq = df_poly_for("OR"), basis = "polynomial")
  or.bs <-   cpi.spline(y = y, a = a, rhat = r, eval = eval.pts, 
                        der.est = 0, df.seq = df_bs_for("OR"), basis = "bspline")
  or.lpoly <- fit_lprobust(y = y, x = r, estimator = "OR")
  oracle_lpoly_h <- mean(or.lpoly$Estimate[, "h"], na.rm = TRUE)

  ##PI
  pi.fourier <- cpi.spline(y = y, rhat = rhat, a = a, eval=eval.pts, der.est = 0,
                           df.seq = df_fourier_for("PI") , basis = "fourier")
  pi.poly <- cpi.spline(y = y, rhat = rhat, a=a ,eval=eval.pts, der.est= 0,
                        df.seq = df_poly_for("PI") , basis ="polynomial")
  pi.bs <-cpi.spline(y = y, rhat = rhat, a=a, eval=eval.pts, der.est = 0,
                     df.seq = df_bs_for("PI") , basis = "bspline")
  pi.lpoly <- fit_lprobust(y = y, x = rhat, estimator = "PI",
                           oracle_h = oracle_lpoly_h)

  #if(fix.df.or.h == TRUE){
  #  df.seq.fourier = pi.fourier$df.or.h
  #  df.seq.bs = pi.bs$df.or.h
  #  df.seq.poly = pi.poly$df.or.h
  #  h.seq = mean(pi.lpoly$Estimate[,"h"])
  #}

  ### Corrected plug-ins oracle (knowledge of derivative) ###
  f.der_p = der.m.t(rhat)
  cpi.or.fourier <- cpi.spline(y = y-f.der_p*(a-rhat), rhat = rhat, a = a, eval=eval.pts, der.est = f.der_p,
                               df.seq = df_fourier_for("CPI_or") , basis = "fourier")
  cpi.or.poly <- cpi.spline(y = y-f.der_p*(a-rhat), rhat = rhat, a = a, eval=eval.pts, der.est = f.der_p,
                            df.seq = df_poly_for("CPI_or") , basis = "polynomial")
  cpi.or.bs <- cpi.spline(y = y-f.der_p*(a-rhat), rhat = rhat, a = a, eval=eval.pts, der.est = f.der_p,
                          df.seq = df_bs_for("CPI_or") , basis = "bspline")
  cpi.or.lpoly <- fit_lprobust(y = y - f.der_p * (a - rhat), x = rhat,
                               estimator = "CPI_or", oracle_h = oracle_lpoly_h)
  
  ### Corrected plug-ins ###

  cpi.fourier <- cpi.spline(y = y-der.est*(a-rhat), rhat = rhat, a = a, eval=eval.pts, der.est = der.est,
                            df.seq = df_fourier_for("CPI") , basis = "fourier")
  cpi.poly <- cpi.spline(y = y-der.est*(a-rhat), rhat = rhat, a = a, eval=eval.pts, der.est = der.est,
                         df.seq = df_poly_for("CPI") , basis = "polynomial")
  cpi.bs <- cpi.spline(y = y-der.est*(a-rhat), rhat = rhat, a = a, eval=eval.pts, der.est = der.est,
                       df.seq = df_bs_for("CPI") , basis = "bspline")
  cpi.lpoly <- fit_lprobust(y = y - der.est * (a - rhat), x = rhat,
                            estimator = "CPI", oracle_h = oracle_lpoly_h)
  
  ##PI with calibration (corrected outcome: same as CPI but using rhat_cali)
  y_cali <- y - der.est_cali * (a - rhat_cali)
  pi.cali.fourier <- cpi.spline(y = y_cali, rhat = rhat_cali, a = a, eval=eval.pts, der.est = der.est_cali,
                            df.seq = df_fourier_for("PI_cali"), basis = "fourier")
  pi.cali.poly <- cpi.spline(y = y_cali, rhat = rhat_cali, a = a, eval=eval.pts, der.est = der.est_cali,
                         df.seq = df_poly_for("PI_cali"), basis = "polynomial")
  pi.cali.bs <- cpi.spline(y = y_cali, rhat = rhat_cali, a = a, eval=eval.pts, der.est = der.est_cali,
                       df.seq = df_bs_for("PI_cali"), basis = "bspline")
  # first try automatic bandwidth
  pi.cali.lpoly <- try(
    fit_lprobust(y = y_cali, x = rhat_cali, estimator = "PI_cali",
                 oracle_h = oracle_lpoly_h),
    silent = TRUE
  )
  # if failed, retry with chosen bandwidth
  if (inherits(pi.cali.lpoly, "try-error")) {
    message("lprobust failed for pi.cali.lpoly; retrying with CPI bandwidth.")
    
    h_retry <- lpoly_h_arg("CPI", oracle_h = oracle_lpoly_h)
    if (is.null(h_retry)) h_retry <- mean(cpi.lpoly$Estimate[, "h"], na.rm = TRUE)
    
    pi.cali.lpoly <- lprobust(
      y = y_cali,
      x = rhat_cali,
      p = poly_deg,
      eval = eval.pts,
      kernel = "gau",
      h = h_retry
    )
  }
  
  ### IF-based bias correction oracle (knowledge of mu) ### Calibration+debias
  bc.fourier.or.cali = bc.spline(y = y, rhat = rhat_cali, a = a, muhat = mu, eval = eval.pts,
                  df.seq = df_fourier_for("BC_or_cali") , basis = "fourier", g_type = g_type)
  bc.poly.or.cali = bc.spline(y = y, rhat = rhat_cali, a = a, muhat = mu, eval = eval.pts,
                  df.seq = df_poly_for("BC_or_cali"), basis ="polynomial", g_type = g_type)
  bc.bs.or.cali = bc.spline(y = y, rhat = rhat_cali, a = a, muhat = mu, eval = eval.pts,
                  df.seq = df_bs_for("BC_or_cali"), basis = "bspline", g_type = g_type)
  bc.lpoly.or.cali <- bc.estimator(y = y, rhat = rhat_cali, a = a, muhat = mu,
                         eval = eval.pts, hseq = bc_hseq_for("BC_or_cali", oracle_lpoly_h),
                         deg = poly_deg, g_type = g_type)

  ### IF-based bias correction ### Calibration+debias
  bc.fourier.cali = bc.spline(y = y, rhat = rhat_cali, a = a, muhat = muhat, eval = eval.pts,
                  df.seq = df_fourier_for("BC_cali") , basis = "fourier", g_type = g_type)
  bc.poly.cali = bc.spline(y = y, rhat = rhat_cali, a = a, muhat = muhat, eval = eval.pts,
                  df.seq = df_poly_for("BC_cali"), basis ="polynomial", g_type = g_type)
  bc.bs.cali = bc.spline(y = y, rhat = rhat_cali, a = a, muhat = muhat, eval = eval.pts,
                  df.seq = df_bs_for("BC_cali"), basis = "bspline", g_type = g_type)

  bc.lpoly.cali <- bc.estimator(y = y, rhat = rhat_cali, a = a, muhat = muhat,
                                   eval = eval.pts, hseq = bc_hseq_for("BC_cali", oracle_lpoly_h),
                                   deg = poly_deg, g_type = g_type)


  ## BC without cailbration
  bc.fourier = bc.spline(y = y, rhat = rhat, a = a, muhat = muhat, eval = eval.pts,
                  df.seq = df_fourier_for("BC") , basis = "fourier", g_type = g_type)
  bc.poly = bc.spline(y = y, rhat = rhat, a = a, muhat = muhat, eval = eval.pts,
                  df.seq = df_poly_for("BC"), basis ="polynomial", g_type = g_type)
  bc.bs = bc.spline(y = y, rhat = rhat, a = a, muhat = muhat, eval = eval.pts,
                  df.seq = df_bs_for("BC"), basis = "bspline", g_type = g_type)
  bc.lpoly <- bc.estimator(y = y, rhat = rhat, a = a, muhat = muhat,
                                eval = eval.pts, hseq = bc_hseq_for("BC", oracle_lpoly_h),
                                deg = poly_deg, g_type = g_type)

  ## BC.or without cailbration (with knowledga of mu)
  bc.fourier.or = bc.spline(y = y, rhat = rhat, a = a, muhat = mu, eval = eval.pts,
                      df.seq = df_fourier_for("BC_or") , basis = "fourier", g_type = g_type)
  bc.poly.or = bc.spline(y = y, rhat = rhat, a = a, muhat = mu, eval = eval.pts,
                      df.seq = df_poly_for("BC_or"), basis ="polynomial", g_type = g_type)
  bc.bs.or = bc.spline(y = y, rhat = rhat, a = a, muhat = mu, eval = eval.pts,
                      df.seq = df_bs_for("BC_or"), basis = "bspline", g_type = g_type)
  bc.lpoly.or <- bc.estimator(y = y, rhat = rhat, a = a, muhat = mu,
                              eval = eval.pts, hseq = bc_hseq_for("BC_or", oracle_lpoly_h),
                              deg = poly_deg, g_type = g_type)
  

  rslt <- list(
    OR = list(
      poly = or.poly,
      bs = or.bs,
      fourier = or.fourier,
      lpoly = list(est = or.lpoly$Estimate[,"tau.us"], se = or.lpoly$Estimate[,"se.us"],
                   df.or.h = mean(or.lpoly$Estimate[,"h"]))

    ),
    
    PI = list(
      fourier = pi.fourier,
      poly = pi.poly,
      bs = pi.bs,
      lpoly  = list(est = pi.lpoly$Estimate[,"tau.us"], se  = pi.lpoly$Estimate[,"se.us"],
                    df.or.h = mean(pi.lpoly$Estimate[,"h"]))
    ),
    
    CPI_or = list(
      fourier =cpi.or.fourier,
      poly = cpi.or.poly,
      bs = cpi.or.bs,
      lpoly = list(est = cpi.or.lpoly$Estimate[,"tau.us"], se = cpi.or.lpoly$Estimate[,"se.us"],
                   df.or.h = mean(cpi.or.lpoly$Estimate[,"h"]))

    ),
    
    CPI = list(
      fourier = cpi.fourier,
      poly = cpi.poly,
      bs = cpi.bs,
      lpoly = list(est = cpi.lpoly$Estimate[,"tau.us"], se  = cpi.lpoly$Estimate[,"se.us"],
                   df.or.h = mean(cpi.lpoly$Estimate[,"h"]))

    ),
    
    PI_cali = list(
      fourier = pi.cali.fourier,
      bs = pi.cali.bs,
      poly = pi.cali.poly,
      lpoly = list(est = pi.cali.lpoly$Estimate[,"tau.us"], se  = pi.cali.lpoly$Estimate[,"se.us"],
                   df.or.h = mean(pi.cali.lpoly$Estimate[,"h"]))

    ),
    
    BC_or_cali = list(
      poly = bc.poly.or.cali,
      fourier = bc.fourier.or.cali,
      bs = bc.bs.or.cali,
      lpoly = bc.lpoly.or.cali
    ),
    
    BC_cali = list(
      poly = bc.poly.cali,
      fourier = bc.fourier.cali,
      bs = bc.bs.cali,
      lpoly = bc.lpoly.cali   
    ),
    
    BC = list(
      poly = bc.poly,
      fourier = bc.fourier,
      bs = bc.bs,
      lpoly = bc.lpoly
    ),
    
    BC_or = list(
      poly = bc.poly.or,
      fourier = bc.fourier.or,
      bs = bc.bs.or,
      lpoly = bc.lpoly.or
    )
  )
  
  
  return(rslt)
} 
