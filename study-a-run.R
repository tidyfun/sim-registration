# Study A production run
# Usage: Rscript study-a-run.R [n_cores]

library(tf)
library(parallel)

base_dir <- here::here()
results_dir <- file.path(base_dir, "results")
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

source(file.path(base_dir, "sim-dgp.R"))
source(file.path(base_dir, "sim-methods.R"))
source(file.path(base_dir, "sim-metrics.R"))
source(file.path(base_dir, "sim-config.R"))
source(file.path(base_dir, "sim-run.R"))

args <- commandArgs(trailingOnly = TRUE)
n_cores <- if (length(args) >= 1) as.integer(args[1]) else 4L

run_benchmark(study = "A", n_cores = n_cores)
