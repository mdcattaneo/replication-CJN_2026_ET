################################################################################
# function_logit_mse.R
#
# Estimation helper functions for the pairwise logit MSE simulation.
#
# Functions:
#   estL0 — L=0 estimate via weighted logistic regression
#   estL1 — jackknife-debiased (L=1) estimate via weighted logistic regression
################################################################################
library(Rcpp)
sourceCpp("logit.cpp")

# L=0 estimator.  Returns the first coefficient hat_theta(h) for each
# bandwidth in hgrid as a numeric vector of length hlen.
estL0 <- function(dat, hgrid) {
  n     <- length(dat$y)
  w.pos <- grepl("w",names(dat))
  wmat  <- as.matrix(dat[, w.pos, drop=FALSE])

  tmp     <- matrix(NA,nrow=length(hgrid),ncol=2)

  for (k in seq_along(hgrid)) {
    prep       <- precompute_active_pairs(dat$y, cbind(dat$x1,dat$x2), wmat, hgrid[k])
    
    dset       <- as.data.frame(prep$dX)
    colnames(dset) <- paste0("x",seq_len(ncol(dset)))
    dset$y     <- as.integer(prep$yi)
    
    tmp[k,]    <- glm(y~x1+x2-1,family=quasibinomial(link="logit"),data=dset,weights=prep$kern_vals)$coefficients
  }
  return(tmp[,1,drop=TRUE])
}

# Jackknife-debiased (L=1) estimator.  Returns theta_tilde(h) for each
# bandwidth in hgrid as a numeric vector of length hlen.
estL1 <- function(dat, hgrid, c1, lambdas) {
  w.pos   <- grepl("w",names(dat))
  wmat    <- as.matrix(dat[, w.pos, drop=FALSE])
  hc1     <- c1*hgrid
  out     <- rep(NA,length(hgrid))

  for (k in seq_along(hgrid)) {
    prep       <- precompute_active_pairs(dat$y, cbind(dat$x1,dat$x2), wmat, hgrid[k])
    dset       <- as.data.frame(prep$dX)
    colnames(dset) <- paste0("x",seq_len(ncol(dset)))
    dset$y     <- as.integer(prep$yi)
    est0       <- glm(y~x1+x2-1,family=quasibinomial(link="logit"),data=dset,weights=prep$kern_vals)$coefficients

    prep       <- precompute_active_pairs(dat$y, cbind(dat$x1,dat$x2), wmat, hc1[k])
    dset       <- as.data.frame(prep$dX)
    colnames(dset) <- paste0("x",seq_len(ncol(dset)))
    dset$y     <- as.integer(prep$yi)
    est1       <- glm(y~x1+x2-1,family=quasibinomial(link="logit"),data=dset,weights=prep$kern_vals)$coefficients
    
    out[k]     <- est0[1]*lambdas[1] + est1[1]*lambdas[2]
  }
  return(out)
}