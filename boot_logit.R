################################################################################
# boot_logit.R
#
# Bootstrap confidence interval simulation for the pairwise logit estimator.
# Designed to run as a SLURM job-array task; each node handles a disjoint
# subset of Monte Carlo replications and writes a partial RDS file.
# Partial files are combined by combine_logit.R.
#
# Usage:
#   Rscript boot_logit.R n=2000 R=2000 B_boot=2000 c1=2 model=1 \
#           node_id=1 n_nodes=10 outdir=partials
#
# Arguments:
#   n       : sample size                        [default 2000]
#   R       : total Monte Carlo replications     [default 2000]
#   B_boot  : bootstrap replications per MC rep  [default 2000]
#   c1      : jackknife scaling constant         [default 2]
#   model   : DGP model index (1, 2, or 3)       [default 1]
#   node_id : this node's index (1..n_nodes)     [default 1]
#   n_nodes : total number of parallel nodes     [default 10]
#   outdir  : output directory for partial RDS   [default "partials"]
################################################################################
rm(list=ls())
source("function_logit_boot.R")

beta0  <- c(1.0,1.0)

## ---- parse args like R=2000 n=2000 B_boot=2000 c1=2 outdir=partials node_id=1 n_nodes=10----
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(key, default) {
  hit <- grep(paste0("^", key, "="), args, value = TRUE)
  if (length(hit) == 0) return(default)
  sub(paste0("^", key, "="), "", hit[1])
}
n     <- as.integer(get_arg("n", "2000"))
B_boot<- as.integer(get_arg("B_boot", "2000"))
R     <- as.integer(get_arg("R", "2000"))
node_id <- as.integer(get_arg("node_id", "1"))   # 1..n_nodes
n_nodes <- as.integer(get_arg("n_nodes", "10"))   # total number of nodes
outdir  <- get_arg("outdir", "partials")
model   <- as.integer(get_arg("model", "1"))      # DGP model: 1, 2, or 3

# MSE-optimal bandwidths (from separate MSE simulations via mse_combine_logit.R).
# hmseL0: min-MSE bandwidth for the L=0 estimator.
# hmseL1: min-MSE bandwidth for the jackknife (L=1) estimator with c1=2.
if (model==1){
  hmseL0 <- 0.2
  hmseL1 <- 0.36 # c1=2
} else if (model==2){
  hmseL0 <- 0.28
  hmseL1 <- 0.45 # c1=2
} else {
  hmseL0 <- 0.39
  hmseL1 <- 0.50 # c1=2
}
hvals  <- unique(sort(round(c(hmseL0*seq(0.5,1.5,by=0.1),hmseL1*seq(0.5,1.5,by=0.1)), 10) ) )
hlen   <- length(hvals)

## constant for jackknife estimator 

c1      <- as.numeric(get_arg("c1", "2"))
lambdas <- c(c1^2/(c1^2-1), -1/(c1^2-1) )



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
out.CI <- array(NA,dim=c(R_i,hlen,2,4))

cat(sprintf("Simulation start: n =%5d, boot iterations =%5d\n",n,B_boot))

current.time <- Sys.time()

for (r in seq_len(R_i)){
  OBS <- generate_data(n, beta0, model)
  orig.est <- thetahat(OBS,hvals,c1,lambdas)
  
  boot_store <- array(NA,dim=c(B_boot,hlen,4))
  for (b in seq_len(B_boot)){
    boot.OBS        <- OBS[sample.int(n, size = n, replace = TRUE),]
    boot_store[b,,] <- thetahat(boot.OBS,hvals,c1,lambdas)
  }
  
  out.CI[r, , , ] <- boot_result(boot_store,orig.est)
  
  cat(sprintf("Iteration %3d / %3d finished in %.2f hours\n",
              r, R_i, as.numeric(difftime(Sys.time(), current.time, units = "hours"))))
} 

## ---- save partial ----
outfile <- file.path(outdir, sprintf("boot_part_model%d_node%03d_of_%03d_n%d_c%g.rds",
                                     model, node_id, n_nodes, n, c1))
saveRDS(
  list(
    model   = as.integer(model),
    node_id = as.integer(node_id),
    n_nodes = as.integer(n_nodes),
    n       = as.integer(n),
    R_total = as.integer(R),
    rep_ids = rep_ids,
    hvals   = hvals,
    c1      = c1,
    cis     = out.CI
  ),
  file = outfile
)


cat("Saved:", outfile, "\n")
