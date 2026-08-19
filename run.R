#!/usr/bin/env Rscript
# cell calling module
#
# Reads a raw h5ad matrix runs emptyDrops per pool,
# assigns species (human/mouse/ambiguous),
# creates sample_id = "{pool_id}_{species}", and writes a TSV of called cell barcodes.

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
p <- add_argument(p, "--passthrough", flag = TRUE, help = "create empty placeholder instead of running cell calling")
args <- parse_args(p)

# logging
cat(sprintf("Full command: %s\n", paste(commandArgs(trailingOnly = FALSE), collapse = " ")))
cat(sprintf("LOG: command line args\n----------------------------------\n"))
for (i in 1:length(args)) {
  cat(sprintf("  %s: %s\n", names(args)[i], args[[i]]))
}
cat(sprintf("----------------------------------\n"))


assign_species <- function(sce, quantile = 5) {

  mat <- counts(sce)
  rn  <- rownames(mat)

  # get sum of counts for human and mouse genes per cell
  human_sums <- colSums(mat[startsWith(rn, "GRCh38_"), , drop = FALSE])
  mouse_sums <- colSums(mat[startsWith(rn, "GRCm39_"), , drop = FALSE])

  dt <- data.table(
    cell_id = colnames(sce),
    human_sums = human_sums,
    mouse_sums = mouse_sums
  )

  dt[, total_sums := human_sums + mouse_sums]
  dt[, human_ratios := human_sums / total_sums]
  dt[, mouse_ratios := mouse_sums / total_sums]
  dt[, species_majority_ratio := pmax(human_ratios, mouse_ratios)]

  # get lower bound for majority ratios based on the specified percentile.
  thresh <- quantile(dt$species_majority_ratio, probs = quantile / 100, na.rm = TRUE)

  # asign species
  dt[, species_raw := fifelse(human_ratios > mouse_ratios, "human", "mouse")]
  dt[, species := fifelse(species_majority_ratio > thresh, species_raw, "ambiguous")]

  return(dt[, .(cell_id, species, species_majority_ratio)])
}


main <- function() {
  dir.create(args$output_dir, showWarnings = FALSE, recursive = TRUE)

  if (args$passthrough) {
    out_path <- file.path(args$output_dir, paste0(args$name, "_called_cells.tsv"))
    file.create(out_path)
    cat(sprintf("  wrote empty file:: %s\n", out_path))
    quit(save = "no", status = 0)
  }

  message("  reading raw h5ad ..")
  sce <- read_h5ad(args$rawdata_raw_h5ad, as = "SingleCellExperiment")
  pools <- unique(colData(sce)$pool_id)

  cells_dt_ls <- lapply(pools, function(pool){

    message(sprintf("  --- processing pool: %s ---", pool))
    pool_idx <- which(colData(sce)$pool_id == pool)
    sce_pool <- sce[, pool_idx]

    message("  running emptyDrops ..")
    edrops_out <- emptyDrops(counts(sce_pool))

    is_cell <- edrops_out$FDR <= args$fdr_threshold
    is_cell[is.na(is_cell)] <- FALSE
    sce_filtered <- sce_pool[, is_cell]
    message(sprintf("  %d cells after filtering", ncol(sce_filtered)))

    message("  assigning species ..")
    cell_dt <- assign_species(sce_filtered)
    
    cell_dt[, pool_id := pool]
    cell_dt[, sample_id := paste0(pool, "_", species)]

    message(sprintf("  species counts: %s",
      paste(names(table(cell_dt$species)),
            table(cell_dt$species), sep = "=", collapse = ", ")))

    return(cell_dt)
  })

  out_dt <- rbindlist(cells_dt_ls)

  out_path <- file.path(args$output_dir, paste0(args$name, "_called_cells.tsv"))
  fwrite(out_dt, out_path,sep = "\t", quote = FALSE, row.names = FALSE)
  message(sprintf("  wrote: %s", out_path))
}

if (sys.nframe() == 0L) {
  main()
}
