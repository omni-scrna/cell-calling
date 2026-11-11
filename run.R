#!/usr/bin/env Rscript
# 
# Reads a raw h5ad matrix (one pool or sample /dataset) and runs emptyDrops, writing a TSV
# of called cell barcodes. For flex datasets (properties.yaml's tenx_assay == "flex"),
# emptyDrops is run separately per probe_id rather than once over the whole pool.

suppressPackageStartupMessages({
  library(DropletUtils)
  library(anndataR)
  library(SingleCellExperiment)
  library(data.table)
  library(yaml)
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
  sce <- read_h5ad(args$rawdata_h5ad, as = "SingleCellExperiment", x_mapping = "counts")

  message("  keeping only RNA (Gene Expression) features ..")
  sce <- sce[rowData(sce)$feature_type == "Gene Expression", ]

  properties <- read_yaml(args$properties_info)

  if (properties$tenx_assay == "flex") {
    # call cells separately for each probe_id
    message("  flex dataset: running emptyDrops per probe_id ..")
    if (!"probe_id" %in% colnames(colData(sce))) {
      stop("tenx_assay is 'flex' but colData(sce) has no 'probe_id' column")
    }

    called_ids <- unlist(lapply(unique(colData(sce)$probe_id), function(pid) {
      in_probe <- colData(sce)$probe_id == pid
      message(sprintf("    probe_id %s: %d barcodes", pid, sum(in_probe)))
      edrops_out <- emptyDrops(counts(sce)[, in_probe])
      is_cell <- edrops_out$FDR <= args$fdr_threshold
      is_cell[is.na(is_cell)] <- FALSE
      colnames(sce)[in_probe][is_cell]
    }))
  } else {
    message("  running emptyDrops ..")
    edrops_out <- emptyDrops(counts(sce))
    is_cell <- edrops_out$FDR <= args$fdr_threshold
    is_cell[is.na(is_cell)] <- FALSE
    called_ids <- colnames(sce)[is_cell]
  }

  message(sprintf("  %d cells after filtering", length(called_ids)))

  out_dt <- data.table(cell_id = called_ids)

  out_path <- file.path(args$output_dir, paste0(args$name, "_called_cells.tsv"))
  fwrite(out_dt, out_path, sep = "\t", quote = FALSE, row.names = FALSE)
  message(sprintf("  wrote: %s", out_path))
}

if (sys.nframe() == 0L) {
  main()
}
