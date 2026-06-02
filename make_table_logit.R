################################################################################
# make_table_logit.R
#
# Reads a compiled RDS produced by combine_logit.R and writes a plain-text
# file containing a LaTeX tabular block (no table wrapper).
#
# Usage:
#   Rscript make_table_logit.R model=1 n=2000 indir=partials outdir=partials c1=2
#
# Arguments:
#   model  : DGP model index             [default 1]
#   n      : sample size                 [required]
#   indir  : directory with compiled RDS [default "partials"]
#   outdir : directory for .txt output   [default same as indir]
#   c1     : jackknife scaling constant  [default 2]
#   yellow : h value to highlight yellow [optional]
#   orange : h value to highlight orange [optional]
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

yellow_h <- suppressWarnings(as.numeric(get_arg("yellow")))
orange_h <- suppressWarnings(as.numeric(get_arg("orange")))

if (is.na(N)) stop("Required argument missing: n=<integer>")

# -- Load compiled results -----------------------------------------------------
rds_file <- file.path(indir, sprintf("compiled_boot_model%d_n%d_c%g.rds", model, N, c1))
if (!file.exists(rds_file)) stop("Compiled RDS not found: ", rds_file)

obj <- readRDS(rds_file)

H_GRID      <- obj$H_GRID
summary_mat <- obj$summary_mat   # hlen x 8
n_valid     <- obj$n_valid       # hlen x 4
hlen        <- length(H_GRID)
R           <- obj$R

# summary_mat column layout (from combine_logit.R):
#   1 B=1 L=0 coverage  2 B=1 L=0 length
#   3 B=1 L=1 coverage  4 B=1 L=1 length
#   5 B=3^(1/d) L=0 coverage  6 B=3^(1/d) L=0 length
#   7 B=3^(1/d) L=1 coverage  8 B=3^(1/d) L=1 length
#
# Table column order: B=1 L=0, B=1 L=1, B=B L=0, B=B L=1

cov <- matrix(NA, nrow = hlen, ncol = 4)
len <- matrix(NA, nrow = hlen, ncol = 4)
for (j in 1:4) {
  cov[, j] <- summary_mat[, 2*j - 1] / n_valid[, j]
  len[, j] <- summary_mat[, 2*j]
}

# -- Build tabular block -------------------------------------------------------
lines <- character(0)
add <- function(...) lines <<- c(lines, paste0(...))

add("\\begin{tabular}{r cccc cccc}")
add("\\toprule")
add("& \\multicolumn{4}{c}{$B = 1$} & \\multicolumn{4}{c}{$B = 3^{1/d}$} \\\\")
add("\\cmidrule(lr){2-5}\\cmidrule(lr){6-9}")
add("& \\multicolumn{2}{c}{$L = 0$} & \\multicolumn{2}{c}{$L = 1$} &",
    " \\multicolumn{2}{c}{$L = 0$} & \\multicolumn{2}{c}{$L = 1$} \\\\")
add("\\cmidrule(lr){2-3}\\cmidrule(lr){4-5}\\cmidrule(lr){6-7}\\cmidrule(lr){8-9}")
add("$h$ & Coverage & Length & Coverage & Length & Coverage & Length & Coverage & Length \\\\")
add("\\midrule")

for (i in seq_len(hlen)) {
  h <- H_GRID[i]
  row_prefix <- ""
  if (!is.na(yellow_h) && abs(h - yellow_h) < 1e-9) {
    row_prefix <- "\\rowcolor{yellow}  "
  } else if (!is.na(orange_h) && abs(h - orange_h) < 1e-9) {
    row_prefix <- "\\rowcolor{orange}  "
  }
  add(sprintf(
    "%s%.2f & %.3f & %.3f & %.3f & %.3f & %.3f & %.3f & %.3f & %.3f \\\\",
    row_prefix,
    h,
    cov[i,1], len[i,1],
    cov[i,2], len[i,2],
    cov[i,3], len[i,3],
    cov[i,4], len[i,4]
  ))
}

add("\\bottomrule")
add("\\end{tabular}")

# -- Write .txt file -----------------------------------------------------------
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
outfile <- file.path(outdir, sprintf("table_boot_model%d_n%d_c%g.txt", model, N, c1))
writeLines(lines, outfile)
cat("Table written to:", outfile, "\n")
cat(paste(lines, collapse = "\n"), "\n")