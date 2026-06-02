library(Rcpp)
sourceCpp("plr.cpp")





# -----------------------------------------------------------
# Point estimate at bandwidth h.
# -----------------------------------------------------------
estimate_beta <- function(obs, h) {
  estimate_plr(obs$y, obs$x, obs$w, h)
}


# -----------------------------------------------------------
# One bootstrap draw.
#
# Evaluates the bootstrap estimator at every point in h_all
# (union of h, 2h, r3*h, 2*r3*h) then maps to the four grids.
#
# Returns a hlen x 4 matrix with columns:
#   1: theta_hat(h)          B=1, no bias correction
#   2: theta_tilde(h)        B=1, jackknife c=2
#   3: theta_hat(r3*h)       B=3^(1/d), no bias correction
#   4: theta_tilde(r3*h)     B=3^(1/d), jackknife c=2
# -----------------------------------------------------------
boot_PLR_union <- function(dat, n, h_all,
                           idx_h, idx_2h, idx_3h, idx_23h) {
  idx      <- sample.int(n, size = n, replace = TRUE)
  boot.dat <- list(y = dat$y[idx], x = dat$x[idx],
                   w = dat$w[idx, , drop = FALSE])

  theta_all <- estimate_plr_grid(boot.dat$y, boot.dat$x, boot.dat$w, h_all)

  th_h   <- theta_all[idx_h]
  th_2h  <- theta_all[idx_2h]
  th_3h  <- theta_all[idx_3h]
  th_23h <- theta_all[idx_23h]

  hlen       <- length(idx_h)
  boot.naive <- matrix(NA_real_, nrow = hlen, ncol = 2)
  boot.3h    <- matrix(NA_real_, nrow = hlen, ncol = 2)

  for (i in seq_len(hlen)) {
    boot.naive[i, 1] <- th_h[i]
    boot.naive[i, 2] <- JK_LAM0 * th_h[i]  + JK_LAM1 * th_2h[i]
    boot.3h[i,    1] <- th_3h[i]
    boot.3h[i,    2] <- JK_LAM0 * th_3h[i] + JK_LAM1 * th_23h[i]
  }

  cbind(boot.naive, boot.3h)
}


# -----------------------------------------------------------
# Full bootstrap for one dataset.
#
# Union grid: sort(unique(c(h, 2h, r3*h, 2*r3*h)))
#   where r3 = 3^(1/d), the CJN bandwidth rescaling factor.
#
# Full-sample estimates:
#   t0[i, 1]    = theta_hat(h_i)            no bc
#   t0[i, 2]    = theta_tilde(h_i)          JK c=2
#   t0.3h[i, 1] = theta_hat(r3*h_i)         no bc  (used for B=3^(1/d) centering)
#   t0.3h[i, 2] = theta_tilde(r3*h_i)       JK c=2 (used for B=3^(1/d) centering)
#
# Bootstrap quantiles (centred):
#   q_naive[i, , j] = quantile of boot_col_j(h)   - t0[i, j]
#   q[i, , j]       = quantile of boot_col_j(r3*h) - t0.3h[i, j]
#
# Returns list(est, boot.quant.naive, boot.quant)
#   est             : hlen x 2  (t0, point estimates at h)
#   boot.quant.naive: hlen x 2 x 2  [bandwidth, quantile side, estimator]
#   boot.quant      : hlen x 2 x 2
# -----------------------------------------------------------
bootstrap_once <- function(dat, B_boot, hvecs, alpha = 0.05) {

  n    <- length(dat$y)
  hlen <- length(hvecs)
  d    <- ncol(dat$w)
  r3   <- 3^(1 / d)

  # ---- union grid ----
  h_all <- sort(unique(c(
    hvecs,
    JK_C        * hvecs,
    r3          * hvecs,
    JK_C * r3   * hvecs
  )))
  key    <- function(x) sprintf("%.17g", x)
  idx_h   <- match(key(hvecs),              key(h_all))
  idx_2h  <- match(key(JK_C * hvecs),       key(h_all))
  idx_3h  <- match(key(r3 * hvecs),         key(h_all))
  idx_23h <- match(key(JK_C * r3 * hvecs),  key(h_all))
  if (anyNA(idx_h) || anyNA(idx_2h) || anyNA(idx_3h) || anyNA(idx_23h))
    stop("Bandwidth matching failed.")

  # ---- full-sample estimates ----
  theta_all <- estimate_plr_grid(dat$y, dat$x, dat$w, h_all)

  th_h   <- theta_all[idx_h]
  th_2h  <- theta_all[idx_2h]
  th_3h  <- theta_all[idx_3h]
  th_23h <- theta_all[idx_23h]

  t0    <- matrix(NA_real_, nrow = hlen, ncol = 2)
  t0.3h <- matrix(NA_real_, nrow = hlen, ncol = 2)

  for (i in seq_len(hlen)) {
    t0[i, 1]    <- th_h[i]
    t0[i, 2]    <- JK_LAM0 * th_h[i]  + JK_LAM1 * th_2h[i]
    t0.3h[i, 1] <- th_3h[i]
    t0.3h[i, 2] <- JK_LAM0 * th_3h[i] + JK_LAM1 * th_23h[i]
  }

  # ---- B_boot bootstrap draws: hlen x 4 x B_boot ----
  boot.out <- replicate(
    B_boot,
    boot_PLR_union(dat, n, h_all, idx_h, idx_2h, idx_3h, idx_23h),
    simplify = "array"
  )

  # ---- bootstrap quantiles ----
  q_naive <- array(NA_real_, dim = c(hlen, 2, 2))
  q       <- array(NA_real_, dim = c(hlen, 2, 2))
  probs   <- c(alpha / 2, 1 - alpha / 2)

  for (i in seq_len(hlen)) {
    for (j in 1:2) {
      q_naive[i, , j] <- quantile(boot.out[i, j,      ] - t0[i, j],    probs, na.rm = TRUE)
      q[i,       , j] <- quantile(boot.out[i, j + 2L, ] - t0.3h[i, j], probs, na.rm = TRUE)
    }
  }

  list(est = t0, boot.quant.naive = q_naive, boot.quant = q)
}

