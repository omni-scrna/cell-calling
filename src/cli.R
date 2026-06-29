suppressPackageStartupMessages(library(optparse))

build_parser <- function() {
  option_list <- list(
    make_option("--output_dir", type = "character",
                help = "Output directory for results"),
    make_option("--name", type = "character",
                help = "Module name/identifier"),
    make_option("--rawdata.h5", type = "character",
                help = "Raw 10x H5 count matrix (unfiltered)"),
    make_option("--fdr_threshold", type = "double", default = 0.01,
                help = "FDR threshold for EmptyDrops cell calling"),
    make_option("--species_quantile", type = "double", default = 5,
                help = "Percentile of majority_ratio for ambiguous threshold [default: 5]")
  )
  OptionParser(
    option_list = option_list,
    description = "EmptyDrops cell calling module"
  )
}

parse_args_checked <- function() {
  parser <- build_parser()
  raw <- parse_args(parser)

  args <- list(
    output_dir    = raw$output_dir,
    name          = raw$name,
    rawdata_h5    = raw[["rawdata.h5"]],
    fdr_threshold = raw$fdr_threshold,
    species_quantile = raw$species_quantile
  )

  required <- c("output_dir", "name", "rawdata_h5")
  missing <- required[vapply(args[required], function(v) is.null(v) || is.na(v),
                             logical(1))]
  if (length(missing) > 0) {
    stop("Missing required argument(s): ", paste(missing, collapse = ", "))
  }

  args
}
