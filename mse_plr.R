rm(list = ls())

library(Rcpp)
sourceCpp("plr.cpp")


# -----------------------------------------------------------
# Command-line arguments
# Usage:
#   Rscript mse_plr.R n=2000 R=2000 outdir=results_plr seed=123
#
# Arguments:
#   n      : sample size              [default 2000]
#   R      : Monte Carlo replications [default 10000]
#   outdir : directory for output RDS [default "results_plr"]
#   seed   : base RNG seed            [default 123]
#
# Runs all three models sequentially and saves one RDS per model.
# Reports MSE for both the plain estimator theta_hat(h) and the
# jackknife L=2 debiased estimator theta_tilde(h) with fixed
# scaling constant c=2:
#
#   theta_tilde(h) = (4/3)*theta_hat(h) - (1/3)*theta_hat(2h)
#
# Weights: lambda0 = c^2/(c^2-1) = 4/3,  lambda1 = -1/(c^2-1) = -1/3.
# -----------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(key, default) {
  hit <- grep(paste0("^", key, "="), args, value = TRUE)
  if (length(hit) == 0) return(default)
  sub(paste0("^", key, "="), "", hit[1])
}

n      <- as.integer(get_arg("n",      "2000"))
R      <- as.integer(get_arg("R",      "10000"))
outdir <- get_arg("outdir", "results_plr")
seed   <- as.integer(get_arg("seed",   "123"))

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
if (!dir.exists(outdir))
  stop(sprintf("Cannot create output directory '%s'.", outdir))

theta0 <- 1.0   # true parameter

# Fixed jackknife scaling constant c = 2
JK_C    <- 2
JK_LAM0 <-  JK_C^2 / (JK_C^2 - 1)   #  4/3
JK_LAM1 <- -1      / (JK_C^2 - 1)    # -1/3


# -----------------------------------------------------------
# Model metadata
# -----------------------------------------------------------
model_info <- list(
  "1" = list(
    d     = 1L,
    hvals = seq(0.1, 1.5, by = 0.05),
    desc  = "Model (6.1): y = x + (1+w^2)         + eps  [d=1]"
  ),
  "2" = list(
    d     = 2L,
    hvals = seq(0.10, 1.50, by = 0.05),
    desc  = "Model (new): y = x + (1+sum_j w_j^2) + eps  [d=2, q=2, equicorrelated N(1,3)]"
  ),
  "3" = list(
    d     = 3L,
    hvals = seq(0.10, 1.50, by = 0.05),
    desc  = "Model (new): y = x + (1+sum_j w_j^2) + eps  [d=3, q=3, equicorrelated N(1,3)]"
  )
)


# -----------------------------------------------------------
# Helper: run MC for one model.
#
# For each replication r and each bandwidth h_i:
#   est_hat[r, i] = theta_hat(h_i)
#   est_jk[r, i]  = (4/3)*theta_hat(h_i) - (1/3)*theta_hat(2*h_i)
#
# Both arrays have the same dimensions (hlen rows).
# The partner 2*h is evaluated in the same loop at no extra overhead.
# -----------------------------------------------------------
run_mse <- function(model, hvals, R, n, theta0, seed) {

  hlen    <- length(hvals)
  hvalsJK <- JK_C * hvals        # partner bandwidths: 2*h

  # Union grid (sorted) passed to estimate_plr_grid in one C++ call per rep.
  # idx_h and idx_2h map back to the two sub-grids.
  h_all  <- sort(unique(c(hvals, hvalsJK)))
  key    <- function(x) sprintf("%.17g", x)
  idx_h  <- match(key(hvals),   key(h_all))
  idx_2h <- match(key(hvalsJK), key(h_all))

  est_hat <- matrix(NA_real_, nrow = R, ncol = hlen)
  est_jk  <- matrix(NA_real_, nrow = R, ncol = hlen)

  set.seed(seed)
  t_start <- proc.time()

  for (r in seq_len(R)) {
    dat       <- generate_data_plr(n, theta0 = theta0, model = model)
    theta_all <- estimate_plr_grid(dat$y, dat$x, dat$w, h_all)
    th_h      <- theta_all[idx_h]
    th_2h     <- theta_all[idx_2h]
    est_hat[r, ] <- th_h
    est_jk[r, ]  <- JK_LAM0 * th_h + JK_LAM1 * th_2h
    if (r %% 500 == 0) {
      elapsed <- (proc.time() - t_start)["elapsed"]
      cat(sprintf("    rep %4d / %4d   elapsed: %.1f s\n", r, R, elapsed))
    }
  }

  elapsed_total <- (proc.time() - t_start)["elapsed"]
  cat(sprintf("  Done.  Total: %.1f s (%.1f min)\n\n",
              elapsed_total, elapsed_total / 60))

  summarise <- function(est) {
    bias <- colMeans(est, na.rm = TRUE) - theta0
    vrn  <- apply(est, 2L, var, na.rm = TRUE)
    mse  <- bias^2 + vrn
    data.frame(h = hvals, bias = bias, sd = sqrt(vrn), rmse = sqrt(mse), mse = mse)
  }

  list(
    est_hat = est_hat,
    est_jk  = est_jk,
    res_hat = summarise(est_hat),
    res_jk  = summarise(est_jk)
  )
}


# -----------------------------------------------------------
# Helper: print one model's result table (both estimators)
# -----------------------------------------------------------
print_table <- function(out, desc, n, R) {
  header <- function(label) {
    cat(sprintf("\n  --- %s ---\n", label))
    cat(sprintf("%-8s  %9s  %9s  %9s  %9s\n", "h", "Bias", "SD", "RMSE", "MSE"))
    cat(strrep("-", 52), "\n")
  }
  print_rows <- function(res) {
    for (i in seq_len(nrow(res))) {
      r <- res[i, ]
      cat(sprintf("%-8.2f  %9.4f  %9.4f  %9.4f  %9.4f\n",
                  r$h, r$bias, r$sd, r$rmse, r$mse))
    }
    best <- which.min(res$mse)
    cat(sprintf("  => MSE-minimising h: %.2f  (MSE = %.5f)\n",
                res$h[best], res$mse[best]))
  }

  cat("=================================================================\n")
  cat("  MSE of CJN PLR Pairwise-Difference Estimator\n")
  cat(sprintf("  %s\n", desc))
  cat(sprintf("  n = %d,  R = %d replications\n", n, R))
  cat("=================================================================\n")

  header("Plain: theta_hat(h)")
  print_rows(out$res_hat)

  header(sprintf("Jackknife L=2: theta_tilde(h) = (%.4f)*theta_hat(h) + (%.4f)*theta_hat(%.1f*h)",
                 JK_LAM0, JK_LAM1, JK_C))
  print_rows(out$res_jk)
  cat("\n")
}


# -----------------------------------------------------------
# Main loop: run all three models
# -----------------------------------------------------------
cat(sprintf(
  "PLR MSE simulation | n=%d, R=%d, seed=%d\n",
  n, R, seed))
cat(sprintf(
  "Jackknife: c=%.1f, lambda0=%.4f, lambda1=%.4f\n\n",
  JK_C, JK_LAM0, JK_LAM1))

summary_rows <- vector("list", 3)

for (model in 1:3) {
  info <- model_info[[as.character(model)]]
  cat(sprintf("--- Model %d (d=%d) ---\n", model, info$d))
  cat(sprintf("  %s\n", info$desc))
  cat(sprintf("  h grid: [%.2f, %.2f], step %.2f, %d points\n\n",
              min(info$hvals), max(info$hvals),
              info$hvals[2] - info$hvals[1], length(info$hvals)))

  out <- run_mse(model, info$hvals, R, n, theta0, seed)
  print_table(out, info$desc, n, R)

  outfile <- file.path(outdir, sprintf("mse_plr_n%04d_m%d.rds", n, model))
  saveRDS(
    list(
      n       = as.integer(n),
      R       = as.integer(R),
      model   = as.integer(model),
      d       = info$d,
      theta0  = theta0,
      jk_c    = JK_C,
      jk_lam  = c(JK_LAM0, JK_LAM1),
      hvals   = info$hvals,
      est_hat = out$est_hat,
      est_jk  = out$est_jk,
      res_hat = out$res_hat,
      res_jk  = out$res_jk
    ),
    file = outfile
  )
  cat("Saved:", outfile, "\n\n")

  bh <- which.min(out$res_hat$mse)
  bj <- which.min(out$res_jk$mse)
  summary_rows[[model]] <- data.frame(
    model       = model,
    d           = info$d,
    h_opt_hat   = out$res_hat$h[bh],
    mse_opt_hat = out$res_hat$mse[bh],
    h_opt_jk    = out$res_jk$h[bj],
    mse_opt_jk  = out$res_jk$mse[bj]
  )
}


# -----------------------------------------------------------
# Cross-model summary
# -----------------------------------------------------------
cat("=================================================================\n")
cat("  Summary: MSE-minimising bandwidths\n")
cat(sprintf("  n=%d,  R=%d,  JK c=%g (lam0=%.4f, lam1=%.4f)\n",
            n, R, JK_C, JK_LAM0, JK_LAM1))
cat("=================================================================\n")
cat(sprintf("%-7s  %-4s  %8s  %10s  %8s  %10s\n",
            "Model", "d", "h*(hat)", "MSE(hat)", "h*(jk)", "MSE(jk)"))
cat(strrep("-", 58), "\n")
for (model in 1:3) {
  s <- summary_rows[[model]]
  cat(sprintf("%-7d  %-4d  %8.2f  %10.5f  %8.2f  %10.5f\n",
              s$model, s$d,
              s$h_opt_hat, s$mse_opt_hat,
              s$h_opt_jk,  s$mse_opt_jk))
}
cat("\n")
