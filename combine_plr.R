rm(list = ls())


# -----------------------------------------------------------
# combine_plr.R
#
# Reads all partial RDS files produced by boot_plr.R, combines
# them, and computes coverage probabilities and average CI
# lengths.
#
# Usage:
#   Rscript combine_plr.R model=1 n=2000 jk_c=2 outdir=partials_plr \
#                         savedir=results_plr
#
# Arguments:
#   model   : 1, 2, or 3, must match the files  [default 1]
#   n       : sample size                        [default 2000]
#   jk_c    : jackknife scaling constant used    [default 2]
#   outdir  : directory containing partial RDS   [default "partials_plr"]
#   savedir : directory for combined output RDS  [default "results_plr"]
#   alpha   : nominal level                      [default 0.05]
# -----------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(key, default) {
  hit <- grep(paste0("^", key, "="), args, value = TRUE)
  if (length(hit) == 0) return(default)
  sub(paste0("^", key, "="), "", hit[1])
}

model   <- as.integer(get_arg("model",   "1"))
n       <- as.integer(get_arg("n",       "2000"))
jk_c    <- as.numeric(get_arg("jk_c",   "2"))
outdir  <- get_arg("outdir",  "partials_plr")
savedir <- get_arg("savedir", "results_plr")
alpha   <- as.numeric(get_arg("alpha",   "0.05"))

stopifnot(model %in% c(1L, 2L, 3L))
dir.create(savedir, showWarnings = FALSE, recursive = TRUE)


# -----------------------------------------------------------
# Discover partial files
# Filename format (from boot_plr.R):
#   plr_part_node{NNN}_of_{KKK}_n{NNNN}_m{M}_c{jk_c}.rds
# -----------------------------------------------------------
pattern <- sprintf("plr_part_node[0-9]+_of_[0-9]+_n%04d_m%d_c%g\\.rds",
                   n, model,jk_c)
files   <- list.files(outdir, pattern = pattern, full.names = TRUE)

if (length(files) == 0)
  stop(sprintf("No partial files found in '%s' for model=%d, n=%d.\n  Pattern: %s",
               outdir, model, n, pattern))

cat(sprintf("Found %d partial file(s) for model=%d, n=%d:\n", length(files), model, n))
for (f in files) cat("  ", f, "\n")
cat("\n")


# -----------------------------------------------------------
# Load and validate all partial files
# -----------------------------------------------------------
parts <- lapply(files, readRDS)

ref <- parts[[1]]
for (k in seq_along(parts)) {
  p <- parts[[k]]
  if (!isTRUE(all.equal(p$hvals, ref$hvals))) stop(sprintf("File %d: hvals mismatch.", k))
  if (p$jk_c  != ref$jk_c)                   stop(sprintf("File %d: jk_c mismatch.", k))
  if (p$model != ref$model)                   stop(sprintf("File %d: model mismatch.", k))
  if (p$n     != ref$n)                       stop(sprintf("File %d: n mismatch.", k))
  if (p$jk_c  != jk_c)
    stop(sprintf("File %d: jk_c=%g but requested jk_c=%g.", k, p$jk_c, jk_c))
}

hvals   <- ref$hvals
hlen    <- length(hvals)
jk_c    <- ref$jk_c
jk_lam  <- ref$jk_lam
theta0  <- ref$theta0
R_total <- ref$R_total


# -----------------------------------------------------------
# Check rep_ids coverage
# -----------------------------------------------------------
all_rep_ids <- unlist(lapply(parts, `[[`, "rep_ids"))
if (anyDuplicated(all_rep_ids))
  warning("Duplicate replication IDs — some replications may be counted twice.")
missing <- setdiff(seq_len(R_total), all_rep_ids)
if (length(missing) > 0)
  warning(sprintf("%d replication(s) missing: %s%s",
                  length(missing),
                  paste(head(missing, 10), collapse = ", "),
                  if (length(missing) > 10) " ..." else ""))

R_found <- length(all_rep_ids)
cat(sprintf("Replications: %d found / %d expected  (%d missing, %d duplicate)\n\n",
            R_found, R_total, length(missing), sum(duplicated(all_rep_ids))))


# -----------------------------------------------------------
# Combine arrays across nodes.
#
# Each partial file has:
#   est     [R_i, hlen, 2]      point estimates (t0 at h)
#   ci      [R_i, hlen, 2, 2]  B=3^(1/d) bootstrap quantiles
#   ci.naive[R_i, hlen, 2, 2]  B=1 bootstrap quantiles
#
# Dimensions:
#   dim 1 : local replication index
#   dim 2 : bandwidth index  (1..hlen)
#   dim 3 : estimator type   1=no bc, 2=jackknife c=2
#   dim 4 : quantile side    1=lower alpha/2, 2=upper 1-alpha/2
# -----------------------------------------------------------
est_all      <- do.call(rbind, lapply(parts, function(p)
  matrix(p$est,      nrow = nrow(p$est),      ncol = hlen * 2)))
ci_all       <- do.call(rbind, lapply(parts, function(p)
  matrix(p$ci,       nrow = nrow(p$ci),       ncol = hlen * 2 * 2)))
ci_naive_all <- do.call(rbind, lapply(parts, function(p)
  matrix(p$ci.naive, nrow = nrow(p$ci.naive), ncol = hlen * 2 * 2)))

est_all      <- array(est_all,      dim = c(R_found, hlen, 2))
ci_all       <- array(ci_all,       dim = c(R_found, hlen, 2, 2))
ci_naive_all <- array(ci_naive_all, dim = c(R_found, hlen, 2, 2))


# -----------------------------------------------------------
# Compute coverage and average CI length.
#
# CI construction:
#   lower = est[r,i,j] - quants[r,i,2,j]   (subtract upper quantile)
#   upper = est[r,i,j] - quants[r,i,1,j]   (subtract lower quantile)
#   covered iff quants[r,i,1,j] <= est[r,i,j] - theta0 <= quants[r,i,2,j]
# -----------------------------------------------------------
compute_coverage_length <- function(est, quants, theta0) {
  nh <- dim(est)[2]

  cov_mat <- matrix(NA_real_, nrow = nh, ncol = 2)
  len_mat <- matrix(NA_real_, nrow = nh, ncol = 2)

  for (i in seq_len(nh)) {
    for (j in 1:2) {
      bias    <- est[, i, j] - theta0
      q_lo    <- quants[, i, 1, j]
      q_hi    <- quants[, i, 2, j]
      cov_mat[i, j] <- mean((q_lo <= bias) & (bias <= q_hi), na.rm = TRUE)
      len_mat[i, j] <- mean(q_hi - q_lo, na.rm = TRUE)
    }
  }
  list(coverage = cov_mat, length = len_mat)
}

res_naive <- compute_coverage_length(est_all, ci_naive_all, theta0)
res_main  <- compute_coverage_length(est_all, ci_all,       theta0)


# -----------------------------------------------------------
# Print table
# -----------------------------------------------------------
model_str <- if (model == 1L) {
  "Model (6.1): y = x + (1+w^2)         + eps  [d=1]"
} else if (model == 2L) {
  "Model (new): y = x + (1+sum_j w_j^2) + eps  [d=2, q=2, equicorrelated N(1,3)]"
} else {
  "Model (new): y = x + (1+sum_j w_j^2) + eps  [d=3, q=3, equicorrelated N(1,3)]"
}

d <- c("1" = 1L, "2" = 2L, "3" = 3L)[[as.character(model)]]

cat(strrep("=", 90), "\n")
cat(sprintf("  Bootstrap %.0f%% Confidence Intervals for theta_0\n", (1 - alpha) * 100))
cat(sprintf("  %s\n", model_str))
cat(sprintf("  n = %d,  R = %d,  Jackknife c=%g (lam0=%.4f, lam1=%.4f)\n",
            n, R_found, jk_c, jk_lam[1], jk_lam[2]))
cat(strrep("=", 90), "\n\n")

cat(sprintf("%-7s  %19s  %19s  %19s  %19s\n",
            "", "   Naive CI: no-bc  ", "  Naive CI: jk c=2  ",
            "   Main CI:  no-bc  ", "   Main CI:  jk c=2 "))
cat(sprintf("%-7s  %9s %9s  %9s %9s  %9s %9s  %9s %9s\n",
            "h", "Coverage", "Length", "Coverage", "Length",
            "Coverage", "Length", "Coverage", "Length"))
cat(strrep("-", 90), "\n")

for (i in seq_len(hlen)) {
  cat(sprintf(
    "%-7.3f  %9.3f %9.3f  %9.3f %9.3f  %9.3f %9.3f  %9.3f %9.3f\n",
    hvals[i],
    res_naive$coverage[i, 1], res_naive$length[i, 1],
    res_naive$coverage[i, 2], res_naive$length[i, 2],
    res_main$coverage[i, 1],  res_main$length[i, 1],
    res_main$coverage[i, 2],  res_main$length[i, 2]
  ))
}

cat("\nNotes:\n")
cat(sprintf("  'no-bc'  = theta_hat (no bias correction)\n"))
cat(sprintf("  'jk c=%g' = theta_tilde = (%.4f)*theta_hat(h) + (%.4f)*theta_hat(%g*h)\n",
            jk_c, jk_lam[1], jk_lam[2], jk_c))
cat(sprintf("  Naive CI : bootstrap resamples at same bandwidth h\n"))
cat(sprintf("  Main CI  : bootstrap resamples at rescaled bandwidth 3^{1/d}*h\n"))
cat(sprintf("             (d=%d, factor = %.4f)\n", d, 3^(1/d)))
cat("\n")


# -----------------------------------------------------------
# Save combined results
# -----------------------------------------------------------
outfile <- file.path(savedir, sprintf("combined_plr_n%04d_m%d_c%g.rds", n, model, jk_c))

saveRDS(
  list(
    n           = as.integer(n),
    R           = R_found,
    R_total     = R_total,
    model       = as.integer(model),
    alpha       = alpha,
    theta0      = theta0,
    hvals       = hvals,
    jk_c        = jk_c,
    jk_lam      = jk_lam,
    est         = est_all,        # R x hlen x 2
    ci          = ci_all,         # R x hlen x 2 x 2
    ci.naive    = ci_naive_all,   # R x hlen x 2 x 2
    res_naive   = res_naive,      # list(coverage, length): hlen x 2
    res_main    = res_main        # list(coverage, length): hlen x 2
  ),
  file = outfile
)

cat("Combined results saved to:", outfile, "\n")
