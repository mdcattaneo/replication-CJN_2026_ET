################################################################################
# plot_coverage_plr.R  (PLR)
#
# Reads a combined RDS produced by combine_plr.R and plots coverage
# probabilities vs bandwidth h for all four CI variants, matching the style
# of the paper figure.
#
# Usage:
#   Rscript plot_coverage_plr.R model=1 n=2000 indir=. outdir=figures jk_c=2
#
# Arguments:
#   model  : DGP model index                  [default 1]
#   n      : sample size                      [required]
#   indir  : directory with combined RDS      [default "."]
#   outdir : directory for PDF output         [default same as indir]
#   jk_c   : jackknife scaling constant       [default 2]
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
indir  <- get_arg("indir", ".")
outdir <- get_arg("outdir", "figures")
jk_c   <- as.numeric(get_arg("jk_c", "2"))
alpha  <- as.numeric(get_arg("alpha", "0.05"))
width  <- as.numeric(get_arg("width",  "6"))
height <- as.numeric(get_arg("height", "4"))

if (model==1){
  hmseL0 <- 0.7
  hmseL1 <- 1.0 # c = 2
} else if (model==2){
  hmseL0 <- 0.25
  hmseL1 <- 0.75 # c = 2
} else {
  hmseL0 <- 0.4
  hmseL1 <- 1.0 # c = 2
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

# -- Load combined results -----------------------------------------------------
rds_file <- file.path(indir, sprintf("combined_plr_n%04d_m%d_c%g.rds", N, model, jk_c))
if (!file.exists(rds_file)) stop("Combined RDS not found: ", rds_file)

obj <- readRDS(rds_file)

hvals     <- obj$hvals
hlen      <- length(hvals)
res_naive <- obj$res_naive   # list(coverage, length): hlen x 2
res_main  <- obj$res_main    # list(coverage, length): hlen x 2

# col 1: naive no-bc  (B=1,   L=0)
# col 2: naive jk     (B=1,   L=1)
# col 3: main  no-bc  (B=B,   L=0)
# col 4: main  jk     (B=B,   L=1)
cov <- cbind(res_naive$coverage[, 1], res_naive$coverage[, 2],
             res_main$coverage[, 1],  res_main$coverage[, 2])

# -- Plot ----------------------------------------------------------------------
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
outfile <- file.path(outdir, sprintf("coverage_plr_model%d.pdf", model))

pdf(outfile, width = width, height = height)
par(mar = c(4, 4.5, 1.5, 1))

plot(hvals, cov[, 1],
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
lines(hvals, cov[, 1], col = "#E69F00",  lty = "dashed", lwd = 2)   # naive no-bc
lines(hvals, cov[, 2], col = "#56B4E9",  lty = "dashed", lwd = 2)   # naive jk
lines(hvals, cov[, 3], col = "#D55E00",  lty = "dotted", lwd = 2)   # main  no-bc
lines(hvals, cov[, 4], col = "#0072B2",  lty = "solid",  lwd = 2)   # main  jk

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