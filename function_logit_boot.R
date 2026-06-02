################################################################################
# function_logit_boot.R
#
# Estimation and bootstrap helper functions for the pairwise logit simulation.
# Sourced by boot_logit.R.
#
# Functions:
#   thetahat_orig  — L=0 and jackknife (L=1) estimates for naive CIs (B=1 only)
#   thetahat       — L=0 and jackknife (L=1) estimates for all four CI variants
#   boot_result    — constructs 95% bootstrap CIs from a boot_store array
################################################################################
library(Rcpp)
sourceCpp("logit.cpp")


# --- Estimation functions -----------------------------------------------------

# Exact floating-point key for matching bandwidth values across sorted grids.
key     <- function(x) sprintf("%.17g", x)

# Compute the L=0 and jackknife-debiased (L=1) estimates for the naive
# CI variants only (bandwidth scale B=1).  Returns an (hlen x 2) matrix:
#   col 1: hat_theta(h)       — L=0 estimator
#   col 2: theta_tilde(h)     — jackknife estimator (lambdas[1]*h + lambdas[2]*c1*h)
thetahat_orig <- function(dat, hgrid, c1,lambdas) {
  n     <- length(dat$y)
  w.pos <- grepl("w",names(dat))
  B     <- 3^( 1/sum( w.pos ) )

  h_all    <- unique(sort( c(hgrid,c1*hgrid) ))
  hall_len <- length(h_all)

  idx_1   <- match(key(hgrid),      key(h_all))
  idx_c   <- match(key(c1*hgrid),      key(h_all))
  tmp     <- matrix(NA,nrow=hall_len,ncol=2)

  wmat      <- as.matrix(dat[, w.pos, drop=FALSE])
  prev_beta <- c(0, 0)
  for (k in seq_len(hall_len)) {
    prep      <- precompute_active_pairs(dat$y, cbind(dat$x1,dat$x2), wmat, h_all[k])
    prev_beta <- fast_wlogit(prep$dX, prep$yi, prep$kern_vals, start=prev_beta)
    tmp[k,]   <- prev_beta
  }

  # Construct estimators
  out     <- matrix(NA, nrow = length(hgrid),ncol=2)

  out[,1] <- tmp[idx_1,1]
  out[,2] <- lambdas[1]*tmp[idx_1,1] + lambdas[2]*tmp[idx_c,1]

  return(out)
}


# Compute all four estimators used in the bootstrap CI simulation.
# B = 3^(1/d) is the bandwidth scale factor for the main CI variants.
# Returns an (hlen x 4) matrix:
#   col 1: hat_theta(h)        — L=0,  bandwidth h
#   col 2: theta_tilde(h)      — L=1,  bandwidth h
#   col 3: hat_theta(B*h)      — L=0,  bandwidth B*h
#   col 4: theta_tilde(B*h)    — L=1,  bandwidth B*h
thetahat <- function(dat, hgrid, c1,lambdas) {
  n     <- length(dat$y)
  w.pos <- grepl("w",names(dat))
  B     <- 3^( 1/sum( w.pos ) )

  h_all    <- unique(sort( c(hgrid,c1*hgrid,B*hgrid,B*c1*hgrid) ))
  hall_len <- length(h_all)

  idx_1   <- match(key(hgrid),          key(h_all))
  idx_c   <- match(key(c1*hgrid),       key(h_all))
  idx_B   <- match(key(B*hgrid),        key(h_all))
  idx_Bc  <- match(key(B*c1*hgrid),     key(h_all))

  tmp     <- matrix(NA,nrow=hall_len,ncol=2)

  wmat      <- as.matrix(dat[, w.pos, drop=FALSE])
  prev_beta <- c(0, 0)
  for (k in seq_len(hall_len)) {
    prep      <- precompute_active_pairs(dat$y, cbind(dat$x1,dat$x2), wmat, h_all[k])
    prev_beta <- fast_wlogit(prep$dX, prep$yi, prep$kern_vals, start=prev_beta)
    tmp[k,]   <- prev_beta
  }

  # Construct estimators
  out <- matrix(NA, nrow = length(hgrid),ncol=4)

  out[,1] <- tmp[idx_1,1]
  out[,2] <- lambdas[1]*tmp[idx_1,1] + lambdas[2]*tmp[idx_c,1]
  out[,3] <- tmp[idx_B,1]
  out[,4] <- lambdas[1]*tmp[idx_B,1] + lambdas[2]*tmp[idx_Bc,1]

  return(out)
}

# Construct 95% percentile bootstrap CIs from a bootstrap-estimate array.
# boot_store: B_boot x hlen x 4 array of bootstrap estimates (from thetahat).
# thetahat:   hlen x 4 matrix of original-sample estimates.
# Returns a hlen x 2 x 4 array of (lower, upper) CI bounds.
# For the main CI variants (l=3,4) the CI is centred on the h-bandwidth
# estimate (cols 1–2) while the pivot is formed from the B*h estimates (cols 3–4).
boot_result <- function(boot_store,thetahat){
  hlen  <- ncol(boot_store[,,1,drop = FALSE])
  ps    <- c(2.5,97.5)/100
  CIs   <- array(NA, dim= c(hlen,2,4) )
  
  for (h in seq_len(hlen)){
    for (l in 1:4){
      qs        <- quantile(boot_store[,h,l] - thetahat[h,l],probs= ps)
      CIs[h,,l] <- c(thetahat[h,(l-1)%%2+1] -qs[2], thetahat[h,(l-1)%%2+1] -qs[1])
    }
  }
  return(CIs)
}