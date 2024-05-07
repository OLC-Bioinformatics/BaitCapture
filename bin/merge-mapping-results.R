#!/usr/bin/env Rscript

#===============================================================================
# merge-mapping-results.R
# Combine read mapping statistics into a single file for downstream analysis and
# visualization
# Author: Liam Brown
# License: MIT (2024)
# GitHub: https://github.com/OLC-Bioinformatics/BaitCapture
#===============================================================================

args = commandArgs(trailingOnly=TRUE)

suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(optparse))

#-------------------------------------------------------------------------------
# Load CLI args
#-------------------------------------------------------------------------------

option_list = list(
  # Input file args
  make_option(
    c("-a", "--aligncov"),
    type = "character",
    default = NULL,
    help = "Path to an AlignCov `*_stats.tsv` file",
    metavar = "path"
  ),
  make_option(
    c("-i", "--idxstats"),
    type = "character",
    default = NULL,
    help = "Path to a TSV file produced by `samtools idxstats`",
    metavar = "path"
  ),
  make_option(
    c("-k", "--kma"),
    type = "character",
    default = NULL,
    help = "Path to a KMA `*.res` file [optional]",
    metavar = "path"
  ),
  # Output options
  make_option(
    c("-o", "--outdir"),
    type = "character",
    default = "merged-results/",
    help = "Path to an output directory to save files [default: %default]",
    metavar = "path"
  ),
  make_option(
    c("-p", "--output_prefix"),
    type = "character",
    default = NULL,
    help = "A prefix to add to the output file names and as a column 'sampleid'",
    metavar = "string"
  ),  
  # Adjustable reporting thresholds
  make_option(
    c("--len_cov_threshold"),
    type = "numeric",
    default = 0,
    help = "The minimum length (in bp) that a target must be covered by to call a positive detection [default: %default] [optional]",
    metavar = "number"
  ),
  make_option(
    c("--prop_cov_threshold"),
    type = "numeric",
    default = 0.9,
    help = "The minimum percentage of length (in bp) that a target must be covered by to call a positive detection [default: %default] [optional]",
    metavar = "number"
  ),
  make_option(
    c("--fold_cov_threshold"),
    type = "numeric",
    default = 0.9,
    help = "The minimum fold-coverage of a target that must be achieved to call a positive detection [default: %default] [optional]",
    metavar = "number"
  ),
  make_option(
    c("--mapped_reads_threshold"),
    type = "numeric",
    default = 2,
    help = "The minimum number of reads that must be mapped to a target to call a positive detection [default: %default] [optional]",
    metavar = "number"
  ),
  make_option(
    c("--percent_identity_threshold"),
    type = "numeric",
    default = NULL,
    help = "The minimum percentage identity match to a target that must be achieved to call a positive detection (only available for KMA results) [optional]",
    metavar = "number"
  ),
  # Workflow args
  make_option(
    c("-r", "--report_all"),
    action = "store_true",
    type = "logical",
    default = FALSE,
    help = "Show non-detected targets in merged result files [optional]",
    metavar = "logical"
  ),
  make_option(
    c("-f", "--force"),
    action = "store_true",
    type = "logical",
    default = FALSE,
    help = "Overwrite existing output file if one already exists within `--outdir` [optional]",
    metavar = "logical"
  )
)

#-------------------------------------------------------------------------------
# Parse CLI args
#-------------------------------------------------------------------------------

# # Debugging
# opt = parse_args(
#   OptionParser(option_list = option_list),
#   args = c(
#     "--aligncov=results/aligncov/kma/Chicken-10-S1-sub_stats.tsv",
#     "--idxstats=results/samtools_stats/kma/Chicken-10-S1-sub.idxstats",
#     # "--kma=results/kma/Chicken-10-S1-sub.res",
#     "--outdir=merged-results/",
#     "--output_prefix=Chicken-10-S1-sub",
#     "--len_cov_threshold=0",
#     "--prop_cov_threshold=0.9",
#     "--fold_cov_threshold=0.9",
#     "--mapped_reads_threshold=2",
#     "--percent_identity_threshold=0.9",
#     "--report_all",
#     "--force"
#     )
#   )

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

# File/dir paths
if (is.null(opt$aligncov)) {
  stop("An AlignCov `*_stats.tsv` file must be provided with `--aligncov`.")
} else {
  aligncov_stats_file = opt$aligncov
}
if (is.null(opt$idxstats)) {
  stop("A TSV file produced by `samtools idxstats` must be provided with `--idxstats`.")
} else {
  idxstats_file = opt$idxstats
}
kma_res_file = opt$kma
outdir = opt$outdir

# Thresholds
len_cov_threshold = opt$len_cov_threshold
prop_cov_threshold = opt$prop_cov_threshold
fold_cov_threshold = opt$fold_cov_threshold
mapped_reads_threshold = opt$mapped_reads_threshold
if (!is.null(opt$percent_identity_threshold) & is.null(kma_res_file)) {
  stop("`--percent_identity_threshold` cannot be used without providing a KMA `*.res` file with `--kma`.")
} else {
  percent_identity_threshold = opt$percent_identity_threshold
}

# Other variables
output_prefix = opt$output_prefix
report_all_targets = opt$report_all
force_overwrite = opt$force

# Get output file names
if (is.null(output_prefix)) {
  merged_out = file.path(outdir, "mapstats.tsv")
  presence_absence_out = file.path(outdir, "presence_absence.tsv")
} else {
  merged_out = file.path(outdir, paste(output_prefix, "mapstats.tsv", sep = "."))
  presence_absence_out = file.path(outdir, paste(output_prefix, "presence_absence.tsv", sep = "."))
}

#-------------------------------------------------------------------------------
# Read and parse files
#-------------------------------------------------------------------------------

aligncov_stats = read_tsv(aligncov_stats_file)

idxstats = read_tsv(idxstats_file,
                    col_select = 1:3,
                    col_names = c(
                      'target',
                      'seqlen',
                      'mapped_reads'
                    )
) |> 
  filter(target != "*")

if (!is.null(kma_res_file)) {
  kma_res = read_tsv(kma_res_file) |> 
    rename_all(~ tolower(.)) |> 
    rename(target = `#template`,
           seqlen = template_length) |>
    relocate(seqlen, .after = target) |> 
    rename_with(.cols = score:p_value,
                .fn = ~ paste("kma", .x, sep = "_"))
}

#-------------------------------------------------------------------------------
# Merge read mapping results
#-------------------------------------------------------------------------------

merged = aligncov_stats |> 
  left_join(idxstats,
            by = c('target', 'seqlen'))

if (!is.null(kma_res_file)) {
  merged = merged |> 
    left_join(kma_res,
              by = c('target', 'seqlen'))
}

if (report_all_targets == FALSE) {
  merged = merged |> 
    filter(depth != 0)
}

#-------------------------------------------------------------------------------
# Create presence-absence matrices based upon thresholds
#-------------------------------------------------------------------------------

if (!is.null(kma_res_file) & !is.null(percent_identity_threshold)) {
  presence_absence = merged |> 
    mutate(presence_absence = if_else(
      len_cov >= len_cov_threshold &
        prop_cov >= prop_cov_threshold &
        fold_cov >= fold_cov_threshold &
        mapped_reads >= mapped_reads_threshold &
        kma_template_identity >= percent_identity_threshold,
      1,
      0)) |> 
    select(target, presence_absence)
} else {
  presence_absence = merged |> 
    mutate(presence_absence = if_else(
      len_cov >= len_cov_threshold &
        prop_cov >= prop_cov_threshold &
        fold_cov >= fold_cov_threshold &
        mapped_reads >= mapped_reads_threshold,
      1,
      0))
}

#-------------------------------------------------------------------------------
# Add 'sampleid' column if --output_prefix is provided
#-------------------------------------------------------------------------------

if (!is.null(output_prefix)) {
  merged = merged |> 
    mutate(sampleid = output_prefix) |> 
    relocate(sampleid, .before = everything())
  presence_absence = presence_absence |> 
    mutate(sampleid = output_prefix) |> 
    relocate(sampleid, .before = everything())
}

#-------------------------------------------------------------------------------
# Save merged results, presence-absence, and package versions
#-------------------------------------------------------------------------------

# Collect package versions
package_versions <- c(
  "readr" = packageVersion("readr"),
  "dplyr" = packageVersion("dplyr"),
  "tidyr" = packageVersion("tidyr"),
  "optparse" = packageVersion("optparse")
)

if (force_overwrite == FALSE && file.exists(merged_out)) {
  # Raise an error if the merged results file already exists and `--force` isn't used
  stop(
    paste(
      "Output file",
      merged_out,
      "already exists. Either remove this file or re-run script with the",
      "`--force` argument to overwrite."
    )
  )
} else if (force_overwrite == FALSE && file.exists(presence_absence_out)) {
  # Raise an error if the presence-absence matrix file already exists and `--force` isn't used
    stop(
      paste(
        "Output file",
        presence_absence_out,
        "already exists. Either remove this file or re-run script with the",
        "`--force` argument to overwrite."
        )
      )
  } else {
  # Create output directory if it doesn't already exist
  if (!dir.exists(outdir)) {
    dir.create(outdir)
  }
  # Write merged results table
  merged |>
    write_tsv(merged_out)
  # Write presence-absence table
  presence_absence |> 
    write_tsv(presence_absence_out)  
  # Write package versions
  writeLines(paste(names(package_versions), package_versions, sep = " = "), 
             file.path(outdir, "package-versions.txt"))
}