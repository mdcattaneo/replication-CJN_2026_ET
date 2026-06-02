rm(list = ls())
source("function_plr_boot.R")


# -----------------------------------------------------------
# Command-line arguments
# Usage (single node):
#   Rscript boot_plr.R n=2000 B_boot=2000 R=2000 model=1 outdir=partials_plr jk_c=2
# Usage (distributed, node k of K):
#   Rscript boot_plr.R n=2000 B_boot=2000 R=2000 model=1 \
#           node_id=k n_nodes=K outdir=partials_plr jk_c=2
#
# Arguments:
#   n       : sample size                            [default 2000]
#   B_boot  : bootstrap draws per replication        [default 2000]
#   R       : total Monte Carlo replications         [default 2000]
#   model   : 1, 2, or 3                             [default 1]
#   node_id : this node's 1-based index              [default 1]
#   n_nodes : total number of nodes                  [default 1]
#   outdir  : directory for partial RDS files        [default "partials_plr"]
#   jk_c    : jackknife constant                     [default 2]
# -----------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(key, default) {
  hit <- grep(paste0("^", key, "="), args, value = TRUE)
  if (length(hit) == 0) return(default)
  sub(paste0("^", key, "="), "", hit[1])
}

n       <- as.integer(get_arg("n",       "2000"))
B_boot  <- as.integer(get_arg("B_boot", "2000"))
R       <- as.integer(get_arg("R",       "2000"))
model   <- as.integer(get_arg("model",   "1"))
node_id <- as.integer(get_arg("node_id", "1"))
n_nodes <- as.integer(get_arg("n_nodes", "1"))
outdir  <- get_arg("outdir", "partials_plr")


JK_C <- as.numeric(get_arg("jk_c","2"))
JK_LAM0 <-  JK_C^2 / (JK_C^2 - 1)
JK_LAM1 <- -1      / (JK_C^2 - 1)

if (is.na(node_id) || is.na(n_nodes) || node_id < 1 || node_id > n_nodes)
  stop("node_id must be an integer in 1..n_nodes")

stopifnot(model %in% c(1L, 2L, 3L))
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# -----------------------------------------------------------
# Bandwidth grids by model.
# -----------------------------------------------------------
if (model == 1){
  h0 <- 0.7
  h1 <- 1.0
}else if (model == 2){
  h0 <- 0.25
  h1 <- 0.75
}else {
  h0 <- 0.4
  h1 <- 1.0
}


theta0 <- 1.0    # true parameter
tmp    <- seq(0.5, 1.5, by = 0.1)
hvals  <- sort( unique( round(c(tmp * h0, tmp * h1),  10) ) )
hlen   <- length(hvals)


# -----------------------------------------------------------
# Split replications across nodes (disjoint, deterministic).
# -----------------------------------------------------------
rep_ids <- seq(from = node_id, to = R, by = n_nodes)
R_i     <- length(rep_ids)


# -----------------------------------------------------------
# Reproducible RNG: independent L'Ecuyer-CMRG streams.
# -----------------------------------------------------------
RNGkind("L'Ecuyer-CMRG")
set.seed(123)
if (node_id > 1) {
  for (k in seq(2L, node_id))
    .Random.seed <- parallel::nextRNGStream(.Random.seed)
}


# -----------------------------------------------------------
# Output arrays
#   out.est[r, i, j]        = t0[i, j]           point estimate at h_i
#   out.q[r, i, s, j]       = boot.quant          B=3^(1/d) quantiles
#   out.q.naive[r, i, s, j] = boot.quant.naive    B=1 quantiles
#
#   r : local replication  (1..R_i)
#   i : bandwidth index    (1..hlen)
#   j : estimator type     (1 = no bias correction, 2 = jackknife)
#   s : quantile side      (1 = lower alpha/2, 2 = upper 1-alpha/2)
# -----------------------------------------------------------
out.est     <- array(NA_real_, dim = c(R_i, hlen, 2))
out.q       <- array(NA_real_, dim = c(R_i, hlen, 2, 2))
out.q.naive <- array(NA_real_, dim = c(R_i, hlen, 2, 2))

model_str <- if (model == 1L) {
  "Model (6.1): y = x + (1+w^2)         + eps  [d=1]"
} else if (model == 2L) {
  "Model (new): y = x + (1+sum_j w_j^2) + eps  [d=2, q=2, equicorrelated N(1,3)]"
} else {
  "Model (new): y = x + (1+sum_j w_j^2) + eps  [d=3, q=3, equicorrelated N(1,3)]"
}

cat(sprintf(
  "PLR bootstrap simulation\n  %s\n  n=%d  B_boot=%d  R=%d  node=%d/%d\n",
  model_str, n, B_boot, R, node_id, n_nodes))
cat(sprintf(
  "  Jackknife: c=%g, lam0=%.4f, lam1=%.4f\n\n",
  JK_C, JK_LAM0, JK_LAM1))

current.time <- Sys.time()

for (r in seq_len(R_i)) {

  OBS <- generate_data_plr(n, theta0 = theta0, model = model)
  tmp <- bootstrap_once(OBS, B_boot, hvals)

  out.est[r, , ]       <- tmp$est
  out.q[r, , , ]       <- tmp$boot.quant
  out.q.naive[r, , , ] <- tmp$boot.quant.naive

  if (r %% 2 == 0) {
    cat(sprintf("  rep %3d / %3d   elapsed: %.2f min\n",
                r, R_i,
                as.numeric(difftime(Sys.time(), current.time, units = "mins"))))
    current.time <- Sys.time()
  }
}


# -----------------------------------------------------------
# Save partial results
# -----------------------------------------------------------
outfile <- file.path(outdir, sprintf(
  "plr_part_node%03d_of_%03d_n%04d_m%d_c%g.rds",
  node_id, n_nodes, n, model,JK_C
))

saveRDS(
  list(
    node_id  = as.integer(node_id),
    n_nodes  = as.integer(n_nodes),
    n        = as.integer(n),
    R_total  = as.integer(R),
    rep_ids  = rep_ids,
    hvals    = hvals,
    jk_c     = JK_C,
    jk_lam   = c(JK_LAM0, JK_LAM1),
    model    = as.integer(model),
    theta0   = theta0,
    est      = out.est,       # point estimates t0  [R_i x hlen x 2]
    ci       = out.q,         # B=3^(1/d) quantiles  [R_i x hlen x 2 x 2]
    ci.naive = out.q.naive    # B=1 quantiles  [R_i x hlen x 2 x 2]
  ),
  file = outfile
)

cat("Saved:", outfile, "\n")



