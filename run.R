#!/usr/bin/env Rscript
# EmptyDrops cell calling module
#
# Reads a raw 10x H5 matrix, runs emptyDrops to identify real cells,
# assigns species (human/mouse/ambiguous), and writes a filtered h5ad matrix.

suppressPackageStartupMessages({
  library(DropletUtils)
  library(anndataR)
})

script_dir <- (function() {
  cargs <- commandArgs(trailingOnly = FALSE)
  m <- grep("^--file=", cargs)
  if (length(m) > 0) dirname(sub("^--file=", "", cargs[[m]])) else getwd()
})()
source(file.path(script_dir, "src", "cli.R"))

# Count human vs mouse reads per cell — genes are identified by genome prefix (e.g., hg19_, mm10_, GRCh38_)
# Compute ratios: human_ratio = human_counts / (human_counts + mouse_counts), same for mouse
# Assign species: whichever ratio is higher → that species. The max ratio is stored as majority_ratio
# Filter ambiguous cells: cells below the 5th percentile of majority_ratio (1st percentile for hgmm1k) are labeled "ambiguous" 

assign_species <- function(sce, quantile = 5) {
  genomes <- rowData(sce)$genome
  unique_genomes <- unique(genomes)

  is_human <- grepl("^(hg|GRCh)", unique_genomes)
  is_mouse <- grepl("^(mm|GRCm)", unique_genomes)
  human_genome <- unique_genomes[is_human]
  mouse_genome <- unique_genomes[is_mouse]

  mat <- counts(sce)
  human_counts <- colSums(mat[genomes == human_genome, , drop = FALSE])
  mouse_counts <- colSums(mat[genomes == mouse_genome, , drop = FALSE])
  total <- human_counts + mouse_counts

  human_ratio <- human_counts / total
  mouse_ratio <- mouse_counts / total

  species_raw <- ifelse(human_ratio > mouse_ratio, "human", "mouse")
  majority_ratio <- pmax(human_ratio, mouse_ratio)

  thresh <- quantile(majority_ratio, probs = quantile / 100)
  species <- ifelse(majority_ratio > thresh, species_raw, "ambiguous")

  colData(sce)$species <- species
  colData(sce)$majority_ratio <- majority_ratio
  sce
}


main <- function() {
  args <- parse_args_checked()
  message(sprintf("Full command: %s", paste(commandArgs(trailingOnly = FALSE), collapse = " ")))
  for (k in names(args)) message(sprintf("  %s: %s", k, args[[k]]))

  dir.create(args$output_dir, showWarnings = FALSE, recursive = TRUE)

  message("  reading raw matrix")
  sce <- read10xCounts(args$rawdata_h5, type = "HDF5")

  message("  running emptyDrops")
  e_out <- emptyDrops(counts(sce))

  is_cell <- e_out$FDR <= args$fdr_threshold
  is_cell[is.na(is_cell)] <- FALSE
  sce_filtered <- sce[, is_cell]

  message("  assigning species")
  sce_filtered <- assign_species(sce_filtered, quantile = args$species_quantile)
  message(sprintf("  species counts: %s",
    paste(names(table(colData(sce_filtered)$species)),
          table(colData(sce_filtered)$species), sep = "=", collapse = ", ")))

  out_path <- file.path(args$output_dir, paste0(args$name, "_filtered.h5ad"))
  write_h5ad(sce_filtered, out_path)
  message(sprintf("  wrote: %s", out_path))
}

if (sys.nframe() == 0L) {
  main()
}
