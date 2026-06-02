################################################################################
# mse_combine_logit.R
#
# Aggregates partial MSE simulation results into a bias / variance / MSE table.
# Reads partial files produced by mse_logit.R (single jackknife constant).
#
#   files : mse_part_model<M>_node<NNN>_of_<MMM>_n<N>_c<c1>.rds
#   hatL0 : R_i x K0          plain estimator (L=0)
#   hatL1 : R_i x K1          debiased estimator (L=1)
#   hvals : list(hL0, hL1)
#   c1    : jackknife scaling constant
#
# Usage:
#   Rscript mse_combine_logit.R model=1 n=2000 c1=2 indir=mse_partials outdir=mse_partials
#
# Arguments:
#   model  : DGP model index (1, 2, or 3)        [default 1]
#   n      : sample size                         [default 2000]
#   c1     : jackknife scaling constant          [default "2"]
#   indir  : directory with partial RDS files    [default "mse_partials"]
#   outdir : directory for compiled RDS          [default same as indir]
################################################################################

BETA_TRUE <- 1.0     # true value of the first coefficient

# ── Argument parsing ───────────────────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(key, default = NULL) {
  hit <- grep(paste0("^", key, "="), args, value = TRUE)
  if (length(hit) == 0) return(default)
  sub(paste0("^", key, "="), "", hit[1])
}

model  <- as.integer(get_arg("model", "1"))
N      <- as.integer(get_arg("n", "2000"))
indir  <- get_arg("indir",  "mse_partials")
outdir <- get_arg("outdir", indir)
c1_arg <- as.numeric(get_arg("c1", "2"))

if (is.na(N)) stop("Invalid n= argument.")
if (is.na(model) || !model %in% 1:3) stop("model= must be 1, 2, or 3.")
if (is.na(c1_arg)) stop("Invalid c1= argument.")

c1_tag_re <- gsub(".", "\\.", sprintf("%g", c1_arg), fixed = TRUE)
pattern   <- sprintf("mse_part_model%d_node.*_n%d_c%s\\.rds", model, N, c1_tag_re)

# ── Discover partial files ─────────────────────────────────────────────────────
files <- sort(list.files(indir, pattern = pattern, full.names = TRUE))

if (length(files) == 0)
  stop(sprintf("No partial files found matching '%s' in '%s'.", pattern, indir))

cat(sprintf("Model: %d\n", model))
cat(sprintf("Found %d partial file(s) for n=%d.\n\n", length(files), N))

# ── Load and stack partial results ─────────────────────────────────────────────
hL0_grid         <- NULL
hL1_grid         <- NULL
c1_ref           <- NULL
R_total_declared <- NA_integer_

hatL0_list <- list()
hatL1_list <- list()

for (i in seq_along(files)) {
  obj <- tryCatch(readRDS(files[i]), error = function(e) {
    warning(sprintf("Cannot read %s: %s", files[i], e$message))
    NULL
  })
  if (is.null(obj)) next

  # On first valid file, record reference metadata
  if (is.null(hL1_grid)) {
    R_total_declared <- obj$R_total
    hL0_grid <- obj$hvals[[1]]
    hL1_grid <- obj$hvals[[2]]
    c1_ref   <- obj$c1
  } else {
    # Validate consistency across files
    if (!is.null(obj$model) && obj$model != model) {
      warning(sprintf("model mismatch in %s — skipping.", files[i])); next
    }
    if (!isTRUE(all.equal(hL0_grid, obj$hvals[[1]])) ||
        !isTRUE(all.equal(hL1_grid, obj$hvals[[2]]))) {
      warning(sprintf("h-grid mismatch in %s — skipping.", files[i])); next
    }
    if (!isTRUE(all.equal(c1_ref, obj$c1))) {
      warning(sprintf("c1 mismatch in %s — skipping.", files[i])); next
    }
  }

  if (nrow(obj$hatL1) > 0) {
    hatL0_list[[i]] <- obj$hatL0   # R_i x K0
    hatL1_list[[i]] <- obj$hatL1   # R_i x K1
  }
}

# Drop NULLs
hatL0_list <- Filter(Negate(is.null), hatL0_list)
hatL1_list <- Filter(Negate(is.null), hatL1_list)

if (length(hatL1_list) == 0)
  stop("No valid replications found across all partial files.")

# Stack along replications
hatL0_all <- do.call(rbind, hatL0_list)   # R_tot x K0
hatL1_all <- do.call(rbind, hatL1_list)   # R_tot x K1

total_reps <- nrow(hatL1_all)
K0 <- length(hL0_grid)
K1 <- length(hL1_grid)

if (!is.na(R_total_declared) && total_reps < R_total_declared)
  warning(sprintf(
    "%d of %d declared replications collected — some task files missing or empty.",
    total_reps, R_total_declared
  ))

# ── Print loading summary ──────────────────────────────────────────────────────
cat(sprintf("Total replications stacked : %d\n", total_reps))
cat(sprintf("L=0 bandwidths (K0=%d)     : %s\n", K0,
            paste(sprintf("%.3f", hL0_grid), collapse = ", ")))
cat(sprintf("L=1 bandwidths (K1=%d)     : %s\n", K1,
            paste(sprintf("%.3f", hL1_grid), collapse = ", ")))
cat(sprintf("c1                         : %g\n\n", c1_ref))

# ── Summary statistics ─────────────────────────────────────────────────────────
mse_stats <- function(estimates, true_val = BETA_TRUE) {
  ok   <- !is.na(estimates)
  n_ok <- sum(ok)
  if (n_ok == 0) return(c(bias=NA, sd=NA, rmse=NA, mse=NA, n_ok=0L))
  e        <- estimates[ok]
  bias     <- mean(e) - true_val
  variance <- var(e)
  mse      <- bias^2 + variance
  c(bias = bias, sd = sqrt(variance), rmse = sqrt(mse), mse = mse, n_ok = n_ok)
}

# ── Print helper ───────────────────────────────────────────────────────────────
print_mse_block <- function(stats_mat, h_grid, label) {
  cat(sprintf("\n  --- %s ---\n", label))
  cat(sprintf("  %-6s  %9s  %9s  %9s  %9s  %6s\n",
              "h", "Bias", "SD", "RMSE", "MSE", "n_ok"))
  cat(strrep("-", 60), "\n")
  for (hh in seq_len(nrow(stats_mat))) {
    cat(sprintf("  %-6.3f  %9.5f  %9.5f  %9.5f  %9.5f  %6d\n",
                h_grid[hh],
                stats_mat[hh, "bias"],
                stats_mat[hh, "sd"],
                stats_mat[hh, "rmse"],
                stats_mat[hh, "mse"],
                as.integer(stats_mat[hh, "n_ok"])))
  }
  best <- which.min(stats_mat[, "mse"])
  cat(sprintf("  => min-MSE h: %.3f  (MSE=%.5f, Bias=%.5f, RMSE=%.5f)\n",
              h_grid[best], stats_mat[best, "mse"],
              stats_mat[best, "bias"], stats_mat[best, "rmse"]))
}

# ── Compute and print ──────────────────────────────────────────────────────────
cat("=================================================================\n")
cat("  MSE Simulation Results\n")
cat(sprintf("  model = %d,  n = %d,  R = %d replications,  beta0[1] = %.1f\n",
            model, N, total_reps, BETA_TRUE))
cat("=================================================================\n")

statsL0 <- t(apply(hatL0_all, 2L, mse_stats))   # K0 x 5
statsL1 <- t(apply(hatL1_all, 2L, mse_stats))   # K1 x 5

print_mse_block(statsL0, hL0_grid,
                "Plain estimator, L=0  [hat_theta(h)]")
print_mse_block(statsL1, hL1_grid,
                sprintf("Debiased estimator, L=1  [theta_tilde(h)],  c = %g", c1_ref))

cat("\n")

# ── Save compiled results ──────────────────────────────────────────────────────
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

outfile <- file.path(outdir,
                     sprintf("compiled_mse_model%d_n%d_c%g.rds", model, N, c1_ref))
saveRDS(
  list(
    model      = model,
    n          = N,
    R          = total_reps,
    beta_true  = BETA_TRUE,
    c1         = c1_ref,
    hL0        = hL0_grid,
    hL1        = hL1_grid,
    statsL0    = statsL0,     # K0 x 5
    statsL1    = statsL1,     # K1 x 5
    hatL0_all  = hatL0_all,   # R x K0
    hatL1_all  = hatL1_all    # R x K1
  ),
  file = outfile
)

cat("Compiled results saved to:", outfile, "\n")
