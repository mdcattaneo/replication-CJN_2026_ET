################################################################################
# make_table_plr.R
#
# Reads a combined RDS produced by combine_plr.R and writes a plain-text
# file containing a LaTeX tabular block (no table wrapper).
#
# Usage:
#   Rscript make_table_plr.R model=1 n=2000 indir=. outdir=. jk_c=2
#
# Arguments:
#   model  : DGP model index              [default 1]
#   n      : sample size                  [required]
#   indir  : directory with combined RDS  [default "."]
#   outdir : directory for .txt output    [default same as indir]
#   jk_c   : jackknife scaling constant   [default 2]
#   yellow : h value to highlight yellow  [optional]
#   orange : h value to highlight orange  [optional]
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
outdir <- get_arg("outdir", indir)
jk_c   <- as.numeric(get_arg("jk_c", "2"))

yellow_h <- suppressWarnings(as.numeric(get_arg("yellow", NA_character_)))
orange_h <- suppressWarnings(as.numeric(get_arg("orange", NA_character_)))

if (is.na(N)) stop("Required argument missing: n=<integer>")

# -- Load combined results -----------------------------------------------------
rds_file <- file.path(indir, sprintf("combined_plr_n%04d_m%d_c%g.rds", N, model, jk_c))
if (!file.exists(rds_file)) stop("Combined RDS not found: ", rds_file)

obj <- readRDS(rds_file)

hvals     <- obj$hvals
hlen      <- length(hvals)
R         <- obj$R
res_naive <- obj$res_naive   # list(coverage, length): hlen x 2
res_main  <- obj$res_main    # list(coverage, length): hlen x 2

# Column layout (j=1: no bias-correction, j=2: jackknife):
#   col 1: naive no-bc   col 2: naive jk
#   col 3: main  no-bc   col 4: main  jk

cov <- cbind(res_naive$coverage[, 1], res_naive$coverage[, 2],
             res_main$coverage[, 1],  res_main$coverage[, 2])
len <- cbind(res_naive$length[, 1],   res_naive$length[, 2],
             res_main$length[, 1],    res_main$length[, 2])

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
  h <- hvals[i]
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
outfile <- file.path(outdir, sprintf("table_plr_model%d.txt",  model))
writeLines(lines, outfile)
cat("Table written to:", outfile, "\n")
cat(paste(lines, collapse = "\n"), "\n")