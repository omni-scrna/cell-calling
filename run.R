#!/usr/bin/env Rscript
# 
# Reads a raw h5ad matrix (one pool/dataset) and runs emptyDrops, writing a TSV
# of called cell barcodes.

suppressPackageStartupMessages({
  library(DropletUtils)
  library(anndataR)
  library(SingleCellExperiment)
  library(data.table)
})

# arg parsing
source("src/common/cli.R")
p <- arg_parser("Cell calling module")
p <- add_base_args(p)
p <- add_stage_args(p, "cell_calling")
p <- add_argument(p, "--fdr_threshold", type = "numeric", default = 0.01, help = "FDR threshold for EmptyDrops")
args <- parse_args(p, argv = strip_output_args())

# logging
cat(sprintf("Full command: %s\n", paste(commandArgs(trailingOnly = FALSE), collapse = " ")))
cat(sprintf("LOG: command line args\n----------------------------------\n"))
for (i in 1:length(args)) {
  cat(sprintf("  %s: %s\n", names(args)[i], args[[i]]))
}
cat(sprintf("----------------------------------\n"))

main <- function() {
  dir.create(args$output_dir, showWarnings = FALSE, recursive = TRUE)

  message("  reading raw h5ad ..")
  # x_mapping: must match the x_mapping="counts" used on write, or the
  # reconstructed assay is named "X" instead of "counts" (DropletUtils/
  # decontX expect an assay literally named "counts").
  sce <- read_h5ad(args$rawdata_h5ad, as = "SingleCellExperiment", x_mapping = "counts")

  message("  running emptyDrops ..")
  edrops_out <- emptyDrops(counts(sce))

  is_cell <- edrops_out$FDR <= args$fdr_threshold
  is_cell[is.na(is_cell)] <- FALSE
  sce_filtered <- sce[, is_cell]
  message(sprintf("  %d cells after filtering", ncol(sce_filtered)))

  out_dt <- data.table(cell_id = colnames(sce_filtered))

  out_path <- file.path(args$output_dir, paste0(args$name, "_called_cells.tsv"))
  fwrite(out_dt, out_path, sep = "\t", quote = FALSE, row.names = FALSE)
  message(sprintf("  wrote: %s", out_path))
}

if (sys.nframe() == 0L) {
  main()
}
