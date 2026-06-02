rm(list=ls())
source("function_logit_mse.R")

hL0   <- seq(0.2,0.6,by=0.01)
hL1  <- seq(0.4,0.8,by=0.01)
beta0 <- c(1.0,1.0)

## ---- parse args like R=2000 n=2000 c1=2 outdir=mse_partials node_id=1 n_nodes=10 ----
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(key, default) {
  hit <- grep(paste0("^", key, "="), args, value = TRUE)
  if (length(hit) == 0) return(default)
  sub(paste0("^", key, "="), "", hit[1])
}
n     <- as.integer(get_arg("n", "2000"))
R     <- as.integer(get_arg("R", "2000"))
node_id <- as.integer(get_arg("node_id", "1"))   # 1..n_nodes
n_nodes <- as.integer(get_arg("n_nodes", "10"))   # total number of nodes
outdir  <- get_arg("outdir", "mse_partials")
model   <- as.integer(get_arg("model", "1"))      # DGP model: 1, 2, or 3

## jackknife scaling constant
c1      <- as.numeric(get_arg("c1", "2"))
lambdas <- c(c1^2 / (c1^2 - 1), -1 / (c1^2 - 1))

if (is.na(node_id) || is.na(n_nodes) || node_id < 1 || node_id > n_nodes) {
  stop("node_id must be an integer in 1..n_nodes")
}

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

## ---- split replications across nodes ----
# Each node gets a disjoint set of replication indices
rep_ids <- seq(from = node_id, to = R, by = n_nodes)
R_i <- length(rep_ids)

## ---- reproducible RNG per node (independent streams) ----
# This ensures reproducibility regardless of how many nodes you use,
# and avoids overlapping random streams.
RNGkind("L'Ecuyer-CMRG")
set.seed(123)
# advance to node-specific stream deterministically
if (node_id > 1) {
  for (k in 2:node_id) .Random.seed <- parallel::nextRNGStream(.Random.seed)
}

## ---- run simulation for this node ----
out.L0 <- matrix(NA, nrow=R_i, ncol=length(hL0))
out.L1 <- matrix(NA, nrow=R_i, ncol=length(hL1))

current.time <- Sys.time()
for (r in seq_len(R_i)){
  OBS         <- generate_data(n, beta0, model)
  out.L0[r,] <- estL0(OBS, hL0)
  out.L1[r,] <- estL1(OBS, hL1, c1, lambdas)
  
  if (r == 1 || r%% 50 == 0 || r == R_i ){
  cat(sprintf("Iteration %3d / %3d finished in %.2f mins\n",
              r, R_i, as.numeric(difftime(Sys.time(), current.time, units = "mins"))))
  }
} 

## ---- save partial ----
outfile <- file.path(outdir, sprintf("mse_part_model%d_node%03d_of_%03d_n%d_c%s.rds",
                                     model, node_id, n_nodes, n,
                                     sprintf("%g", c1)))
saveRDS(
  list(
    model   = as.integer(model),
    node_id = as.integer(node_id),
    n_nodes = as.integer(n_nodes),
    n       = as.integer(n),
    R_total = as.integer(R),
    rep_ids = rep_ids,
    hvals   = list(hL0, hL1),
    c1      = c1,
    hatL0   = out.L0,
    hatL1   = out.L1
  ),
  file = outfile
)
