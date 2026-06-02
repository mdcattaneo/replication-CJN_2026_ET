################################################################################
# plot_coverage_logit.R
#
# Reads a compiled RDS produced by combine_logit.R and plots coverage
# probabilities vs bandwidth h for all four CI variants, matching the style
# of the paper figure.
#
# Usage:
#   Rscript plot_coverage_logit.R model=1 n=2000 indir=results outdir=results c1=2
#
# Arguments:
#   model  : DGP model index                  [default 1]
#   n      : sample size                      [required]
#   indir  : directory with compiled RDS      [default "results"]
#   outdir : directory for PDF output         [default same as indir]
#   c1     : jackknife scaling constant       [default 2]
#   alpha  : nominal level                    [default 0.05]
#   yellow : h value for dashed vertical line (matches table highlight)  [optional]
#   orange : h value for dotted vertical line (matches table highlight)  [optional]
#   width  : plot width in inches             [default 6]
#   height : plot height in inches            [default 4]
################################################################################

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(key, default = NULL) {
  hit <- grep(paste0("^", key, "="), args, value = TRUE)
  if (length(hit) == 0) return(default)
  sub(paste0("^", key, "="), "", hit[1])
}

model  <- as.integer(get_arg("model", "1"))
N      <- as.integer(get_arg("n"))
indir  <- get_arg("indir", "results")
outdir <- get_arg("outdir", indir)
c1     <- as.numeric(get_arg("c1", "2"))
alpha  <- as.numeric(get_arg("alpha", "0.05"))
width  <- as.numeric(get_arg("width",  "6"))
height <- as.numeric(get_arg("height", "4"))

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

tmp    <- get_arg("yellow")
if (is.null(tmp)){
  yellow_h <- hmseL0
} else {
  yellow_h <- tmp
}

tmp    <- get_arg("orange")
if (is.null(tmp)){
  orange_h  <- hmseL1
} else {
  orange_h  <- tmp
}

if (is.na(N)) stop("Required argument missing: n=<integer>")

# -- Load compiled results -----------------------------------------------------
rds_file <- file.path(indir, sprintf("compiled_boot_model%d_n%d_c%g.rds", model, N, c1))
if (!file.exists(rds_file)) stop("Compiled RDS not found: ", rds_file)

obj <- readRDS(rds_file)

H_GRID      <- obj$H_GRID
summary_mat <- obj$summary_mat   # hlen x 8
n_valid     <- obj$n_valid       # hlen x 4
hlen        <- length(H_GRID)

cov <- matrix(NA, nrow = hlen, ncol = 4)
for (j in 1:4)
  cov[, j] <- summary_mat[, 2*j - 1] / n_valid[, j]

# col 1: naive no-bc  (B=1,   L=0)
# col 2: naive jk     (B=1,   L=1)
# col 3: main  no-bc  (B=B,   L=0)
# col 4: main  jk     (B=B,   L=1)

# -- Plot ----------------------------------------------------------------------
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
outfile <- file.path(outdir, sprintf("coverage_logit_model%d.pdf", model))

pdf(outfile, width = width, height = height)
par(mar = c(4, 4.5, 1.5, 1))

plot(H_GRID, cov[, 1],
     type = "n", ylim = c(0, 1),
     xlab = expression(italic(h)),
     ylab = "Coverage",
     las  = 1, bty = "l")

# 1-alpha reference line
abline(h = 1 - alpha, lty = "dashed", col = "grey50", lwd = 1.2)

# vertical marker lines (yellow = dashed, orange = dotted)
if (!is.na(yellow_h)) abline(v = yellow_h, lty = "dashed", col = "grey50", lwd = 1.2)
if (!is.na(orange_h)) abline(v = orange_h, lty = "dotted", col = "grey50", lwd = 1.2)

# four CI curves
lines(H_GRID, cov[, 1], col = "#E69F00",  lty = "dashed", lwd = 2)   # naive no-bc
lines(H_GRID, cov[, 2], col = "#56B4E9",  lty = "dashed", lwd = 2)   # naive jk
lines(H_GRID, cov[, 3], col = "#D55E00",  lty = "dotted", lwd = 2)   # main  no-bc
lines(H_GRID, cov[, 4], col = "#0072B2",  lty = "solid",  lwd = 2)   # main  jk

legend("right", bty = "n", cex = 1.0, lwd = 2,
       col    = c("#E69F00", "#56B4E9", "#D55E00", "#0072B2"),
       lty    = c("dashed", "dashed", "dotted", "solid"),
       legend = c(
         expression(CI[paste(n,","*1-alpha)](1*","*0)),
         expression(CI[paste(n,","*1-alpha)](1*","*1)),
         expression(CI[paste(n,","*1-alpha)](3^{1/d}*","*0)),
         expression(CI[paste(n,","*1-alpha)](3^{1/d}*","*1))
       ))

dev.off()
cat("Plot written to:", outfile, "\n")