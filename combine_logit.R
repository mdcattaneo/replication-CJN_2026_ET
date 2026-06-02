################################################################################
# combine_logit.R
#
# Aggregates partial bootstrap simulation results produced by boot_logit.R
# (SLURM job-array tasks) into a single coverage / CI-length table.
#
# Each partial file (boot_part_model<M>_node<NNN>_of_<MMM>_n<N>_c<c1>.rds):
#   cis    : R_i x hlen x 2 x 4 array  (reps x bandwidths x lo/hi x CI type)
#   hvals  : hlen-vector of bandwidths
#   CI types: 1: B=1 L=0, 2: B=1 L=1, 3: B=3^(1/d) L=0, 4: B=3^(1/d) L=1
#
# Usage:
#   Rscript combine_logit.R model=1 n=2000 indir=partials outdir=partials c1=2
#
# Arguments:
#   model  : DGP model index (1, 2, or 3)        [default 1]
#   n      : sample size                          [required]
#   indir  : directory with partial RDS files     [default "partials"]
#   outdir : directory for compiled RDS           [default same as indir]
#   c1     : jackknife scaling constant           [default 2]
################################################################################

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(key, default = NULL) {
  hit <- grep(paste0("^", key, "="), args, value = TRUE)
  if (length(hit) == 0) return(default)
  sub(paste0("^", key, "="), "", hit[1])
}

model  <- as.integer(get_arg("model", "1"))
N      <- as.integer(get_arg("n"))
indir  <- get_arg("indir", "partials")
outdir <- get_arg("outdir", indir)
c1     <- as.numeric(get_arg("c1", "2"))

a.theta0 <- 1

if (is.na(model) || !model %in% 1:3) stop("model= must be 1, 2, or 3.")
if (is.na(N)) stop("Required argument missing: n=<integer>")

# ── Discover partial files ─────────────────────────────────────────────────────
pattern <- sprintf("boot_part_model%d_node.*_n%d_c%g\\.rds", model, N, c1)
files   <- sort(list.files(indir, pattern = pattern, full.names = TRUE))

if (length(files) == 0)
  stop(sprintf("No partial files found matching '%s' in '%s'.", pattern, indir))

cat(sprintf("Found %d partial file(s).\n\n", length(files)))

# ── Load and stack all partial arrays ─────────────────────────────────────────
ci_list          <- vector("list", length(files))
H_GRID           <- NULL
meta             <- NULL
n_tasks_expected <- NA_integer_

for (i in seq_along(files)) {
  obj <- tryCatch(readRDS(files[i]), error = function(e) {
    warning(sprintf("Cannot read %s: %s", files[i], e$message))
    NULL
  })
  if (is.null(obj)) next

  if (is.null(H_GRID)) {
    H_GRID           <- obj$hvals
    meta             <- obj[c("n", "R_total")]
    n_tasks_expected <- obj$n_nodes
  } else {
    if (!is.null(obj$model) && obj$model != model) {
      warning(sprintf("model mismatch in %s — skipping.", files[i])); next
    }
    if (!isTRUE(all.equal(H_GRID, obj$hvals))) {
      warning(sprintf("H_GRID mismatch in %s — skipping.", files[i])); next
    }
  }

  ci_list[[i]] <- obj$cis
}

# Drop NULLs and stack along dim 1 (replications)
ci_list <- Filter(Negate(is.null), ci_list)

if (length(ci_list) == 0)
  stop("No valid replications found across all partial files.")

hlen     <- dim(ci_list[[1]])[2]
total_reps <- sum(sapply(ci_list, function(a) dim(a)[1]))
full_arr <- array(NA_real_, dim = c(total_reps, hlen, 2L, 4L))
row <- 1L
for (a in ci_list) {
  r <- dim(a)[1L]
  full_arr[row:(row + r - 1L), , , ] <- a
  row <- row + r
}

n_tasks_found <- length(ci_list)
if (!is.na(n_tasks_expected) && n_tasks_found < n_tasks_expected)
  warning(sprintf("%d of %d task file(s) missing — results based on %d completed tasks.",
                  n_tasks_expected - n_tasks_found, n_tasks_expected, n_tasks_found))

cat(sprintf("Total replications: %d\n", total_reps))
cat(sprintf("Bandwidths (hlen=%d): %s\n\n",
            hlen, paste(sprintf("%.3f", H_GRID), collapse = ", ")))

# ── Compute summary statistics ─────────────────────────────────────────────────
# full_arr dims: (total_reps, hlen, 2, 4) — dim3 = lo/hi, dim4 = CI type
# apply over dims 2 and 4: function receives a (total_reps x 2) matrix per (h, CI-type)
coverage_mat <- apply(full_arr, c(2, 4), function(a)
  sum(a[, 1] <= a.theta0 & a[, 2] >= a.theta0, na.rm = TRUE))
length_mat   <- apply(full_arr, c(2, 4), function(a)
  mean(a[, 2] - a[, 1], na.rm = TRUE))

# Interleave coverage (odd cols) and length (even cols) into hlen x 8
summary_mat                       <- matrix(NA, nrow = hlen, ncol = 8)
summary_mat[, seq(1, 7, by = 2)]  <- coverage_mat   # hlen x 4 → cols 1,3,5,7
summary_mat[, seq(2, 8, by = 2)]  <- length_mat      # hlen x 4 → cols 2,4,6,8

col_names <- c("naive_nobc_cov", "naive_nobc_len",
               "naive_jk_cov",   "naive_jk_len",
               "main_nobc_cov",  "main_nobc_len",
               "main_jk_cov",    "main_jk_len")

# Count valid (non-NA) replications per (bandwidth, CI type): hlen x 4
n_valid <- apply(full_arr, c(2L, 4L), function(a) sum(!is.na(a[, 1L])))

# ── Print table ────────────────────────────────────────────────────────────────
# cov_col is an odd summary_mat column (1,3,5,7); CI-type index = (cov_col+1)/2
print_ci_block <- function(cov_col, len_col, label) {
  nv   <- n_valid[, (cov_col + 1L) / 2L]
  cov  <- summary_mat[, cov_col] / nv * 100   # empirical coverage (%)
  len  <- summary_mat[, len_col]

  cat(sprintf("\n  --- %s ---\n", label))
  cat(sprintf("  %-6s  %7s  %9s  %6s\n", "h", "Cov (%)", "Mean Len", "n_ok"))
  cat(strrep("-", 36), "\n")
  for (hh in seq_len(hlen)) {
    cat(sprintf("  %-6.3f  %7.2f  %9.4f  %6d\n",
                H_GRID[hh], cov[hh], len[hh], nv[hh]))
  }
  best <- which.min(abs(cov - 95))
  cat(sprintf("  => h closest to 95%% coverage: %.3f  (%.2f%%,  len = %.4f)\n",
              H_GRID[best], cov[best], len[best]))
}

cat("=================================================================\n")
cat("  Bootstrap CI Simulation Results\n")
cat(sprintf("  model = %d,  n = %d,  c1 = %g,  R = %d replications\n",
            model, N, c1, total_reps))
cat("=================================================================\n")

print_ci_block(1L, 2L, "Naive CI,  No Bias Correction  [hat_theta(h)]")
print_ci_block(3L, 4L, "Naive CI,  Jackknife L=2        [theta_tilde(h)]")
print_ci_block(5L, 6L, "Main CI,   No Bias Correction  [hat_theta(scale*h)]")
print_ci_block(7L, 8L, "Main CI,   Jackknife L=2        [theta_tilde(scale*h)]")

cat("\n")

# ── Save compiled results ──────────────────────────────────────────────────────
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
outfile <- file.path(outdir, sprintf("compiled_boot_model%d_n%d_c%g.rds", model, N, c1))
saveRDS(
  list(
    model       = model,
    N           = N,
    c1          = c1,
    R           = total_reps,
    H_GRID      = H_GRID,
    summary_mat = summary_mat,
    col_names   = col_names,
    n_valid     = n_valid,
    full_arr    = full_arr
  ),
  file = outfile
)
cat("Compiled results saved to:", outfile, "\n")
