#!/usr/bin/env Rscript

#===============================================================================
# summarize-statistics.R
# Summarize reads filtered/lost/mapped at each step
# Author: Liam Brown
# License: MIT (2024)
# GitHub: https://github.com/OLC-Bioinformatics/BaitCapture
#===============================================================================

args = commandArgs(trailingOnly=TRUE)

suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(jsonlite))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(optparse))

#-------------------------------------------------------------------------------
# Load CLI args
#-------------------------------------------------------------------------------

option_list = list(
  # Input file args
  make_option(
    c("-r", "--raw"),
    type = "character",
    default = NULL,
    help = "Path to a `*.json` file produced by fastqscan for raw sequence data",
    metavar = "path"
  ),
  make_option(
    c("-p", "--preprocessed"),
    type = "character",
    default = NULL,
    help = "Path to a `*.json` file produced by fastqscan for preprocessed sequence data [optional]",
    metavar = "path"
  ),  
  make_option(
    c("-t", "--trimmed"),
    type = "character",
    default = NULL,
    help = "Path to a `*.log` file produced by fastp for trimmed sequence data [optional]",
    metavar = "path"
  ),
  make_option(
    c("-s", "--stats"),
    type = "character",
    default = NULL,
    help = "Path to a `*.stats` file produced by `samtools stats` for target alignment",
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
    c("--output_prefix"),
    type = "character",
    default = NULL,
    help = "A prefix to add to the output file names and as a column 'sampleid'",
    metavar = "string"
  ),
  # Workflow args
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

# # Uncomment for debugging
# opt = parse_args(
#   OptionParser(option_list = option_list),
#   args = c(
#     "--raw=results/fastqscan_raw/Chicken-10-S1-sub_raw.json",
#     # "--preprocessed=results/fastqscan_preprocessed/Chicken-10-S1-sub_preprocessed.json",
#     "--trimmed=results/fastp/Chicken-10-S1-sub.fastp.log",
#     "--stats=results/samtools_stats_targets/kma/Chicken-10-S1-sub.stats",
#     "--outdir=summarized-stats/",
#     "--output_prefix=Chicken-10-S1-sub",
#     "--force"
#     )
#   )

# Comment for debugging
opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

# File/dir paths
preprocessed_fastqscan_file = opt$preprocessed
trimmed_fastp_file = opt$trimmed
outdir = opt$outdir

if (is.null(opt$raw)) {
  stop("A fastqscan `*.json` file for raw sequence data must be provided with `--raw`.")
} else {
  raw_fastqscan_file = opt$raw
}
if (is.null(opt$stats)) {
  stop("A `*.stats` file produced by `samtools stats` for target alignment must be provided with `--stats`.")
} else {
  stats_file = opt$stats
}

# Other variables
output_prefix = opt$output_prefix
force_overwrite = opt$force

# Get output file names
sumstats_out = ifelse(is.null(output_prefix),
                    file.path(outdir, "sumstats.tsv"),
                    file.path(outdir, paste(output_prefix, "sumstats.tsv", sep = ".")))

#-------------------------------------------------------------------------------
# Read and parse files
#-------------------------------------------------------------------------------

stats = read_lines(stats_file,
           skip = 7,
           n_max = 39) |> 
  as_tibble() |> 
  mutate(value = gsub("SN\t", "", value)) |> 
  separate(value, 
           into = c("name", "value"), 
           extra = "drop",
           sep = "\t") |> 
  mutate(name = gsub(":", "", name)) |> 
  filter(name %in% c("reads mapped", "bases mapped")) |> 
  pivot_wider(names_from = name,
              values_from = value) |> 
  rename(mapped_total_reads = `reads mapped`,
         mapped_total_bp = `bases mapped`) |> 
  mutate(across(everything(),
                ~ as.integer(.x)))

raw_fastqscan = fromJSON(raw_fastqscan_file, simplifyDataFrame = TRUE)$qc_stats |> 
  as_tibble() |> 
  select(raw_total_reads = read_total, raw_total_bp = total_bp)

if (!is.null(preprocessed_fastqscan_file)) {
  preprocessed_fastqscan = fromJSON(preprocessed_fastqscan_file, simplifyDataFrame = TRUE)$qc_stats |> 
    as_tibble() |> 
    select(preprocessed_total_reads = read_total, preprocessed_total_bp = total_bp)
}

if (!is.null(trimmed_fastp_file)) {
  trimmed_fastp = read_file(trimmed_fastp_file) |> 
    str_split("\n\n") |> 
    unlist() |> 
    as_tibble() |> 
    filter(grepl("filtering", value)) |> 
    extract(
      col = value,
      into = c("r1_r2", "before_after_filtering", "total_reads", "total_bp"),
      regex = "Read([1,2]) ([a-z]+) filtering:\ntotal reads: ([0-9]+)\ntotal bases: ([0-9]+).*"
    ) |> 
    mutate(across(c(total_reads, total_bp),
           ~ as.numeric(.x))) |> 
    group_by(before_after_filtering) |> 
    summarize(total_reads = sum(total_reads),
              total_bp = sum(total_bp)) |> 
    ungroup() |> 
    filter(before_after_filtering == "after") |> 
    select(-before_after_filtering) |> 
    rename(fastp_total_reads = total_reads, fastp_total_bp = total_bp)
}

#-------------------------------------------------------------------------------
# Create a table of read mapping/filtering/loss statistics
#-------------------------------------------------------------------------------

# Using samtools stats output for obtaining reads and bases mapped, because
# the sum of the `mapped_reads` column in samtools idxstats output will count
# bases from reads that mapped to multiple targets if not using KMA

# Case 1: Trimming and dehosting
if (!is.null(trimmed_fastp_file) && !is.null(preprocessed_fastqscan_file)) {
  sumstats = bind_cols(raw_fastqscan, trimmed_fastp, preprocessed_fastqscan, stats) |> 
    mutate(percent_reads_lost_fastp = (1 - fastp_total_reads/raw_total_reads) * 100) |> 
    mutate(percent_reads_lost_dehosting = (1 - preprocessed_total_reads/fastp_total_reads) * 100) |> 
    mutate(percent_reads_on_target = (1 - mapped_total_reads/preprocessed_total_reads) * 100)
# Case 2: Trimming only
} else if (!is.null(trimmed_fastp_file) && is.null(preprocessed_fastqscan_file)) {
  sumstats = bind_cols(raw_fastqscan, trimmed_fastp, stats) |> 
    mutate(percent_reads_lost_fastp = (1 - fastp_total_reads/raw_total_reads) * 100) |> 
    mutate(percent_reads_on_target = (1 - mapped_total_reads/fastp_total_reads) * 100)
# Case 3: Dehosting only
} else if (is.null(trimmed_fastp_file) && !is.null(preprocessed_fastqscan_file)) {
  sumstats = bind_cols(raw_fastqscan, preprocessed_fastqscan, stats) |> 
    mutate(percent_reads_lost_dehosting = (1 - preprocessed_total_reads/raw_total_reads) * 100) |> 
    mutate(percent_reads_on_target = (1 - mapped_total_reads/preprocessed_total_reads) * 100)
# Case 4: No trimming or dehosting
} else if (is.null(trimmed_fastp_file) && is.null(preprocessed_fastqscan_file)) {
  sumstats = bind_cols(raw_fastqscan, stats) |> 
    mutate(percent_reads_on_target = (1 - mapped_total_reads/raw_total_reads) * 100)
}

#-------------------------------------------------------------------------------
# Add 'sampleid' column if --output_prefix is provided
#-------------------------------------------------------------------------------

if (!is.null(output_prefix)) {
  sumstats = sumstats |> 
    mutate(sampleid = output_prefix) |> 
    relocate(sampleid, .before = everything())
}

#-------------------------------------------------------------------------------
# Write summarized statistics to file
#-------------------------------------------------------------------------------

# Collect package versions
package_versions <- c(
  "readr" = packageVersion("readr"),
  "dplyr" = packageVersion("dplyr"),
  "jsonlite" = packageVersion("jsonlite"),
  "stringr" = packageVersion("stringr"),
  "tidyr" = packageVersion("tidyr"),
  "optparse" = packageVersion("optparse")
)

# Check for existing files and and raise error if force_overwrite is disabled
if (force_overwrite == FALSE && file.exists(sumstats_out)) {
  stop(paste("Output file", sumstats_out, "already exists. Either remove this file or re-run script with the `--force` argument to overwrite."))
}

# Create outdir if it doesn't already exist
if (!dir.exists(outdir)) {
  dir.create(outdir)
}

# Write summarized statistics table
sumstats |>
  write_tsv(sumstats_out)

# Write package versions
writeLines(paste(names(package_versions), package_versions, sep = " = "), 
           file.path(outdir, "package-versions.txt"))
