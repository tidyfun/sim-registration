# Study A cleanup: finish the 4 remaining FDA batches for D12/D13
# Usage: Rscript study-a-cleanup.R [n_cores]

library(tf)
library(parallel)

base_dir <- here::here()
results_dir <- file.path(base_dir, "results")

source(file.path(base_dir, "sim-dgp.R"))
source(file.path(base_dir, "sim-methods.R"))
source(file.path(base_dir, "sim-metrics.R"))
source(file.path(base_dir, "sim-config.R"))
source(file.path(base_dir, "sim-run.R"))

args <- commandArgs(trailingOnly = TRUE)
n_cores <- if (length(args) >= 1) as.integer(args[1]) else 4L

design <- full_design("A")
all_tasks <- create_tasks(design)

task_groups <- split(
  all_tasks,
  sapply(all_tasks, function(t) {
    paste(t$dgp, t$method, t$study, sep = "_")
  })
)

cat(sprintf("=== Study A Cleanup ===\n"))
cat(sprintf("Cores: %d | Start: %s\n\n", n_cores, format(Sys.time())))

t0 <- proc.time()
n_done <- 0L

for (nm in names(task_groups)) {
  out_file <- file.path(results_dir, sprintf("results_%s.rds", nm))
  if (file.exists(out_file)) next

  cat(sprintf(
    "[%s] Running %d tasks for %s...\n",
    format(Sys.time(), "%H:%M:%S"),
    length(task_groups[[nm]]),
    nm
  ))
  run_batch(task_groups[[nm]], n_cores = n_cores)
  n_done <- n_done + 1L
}

elapsed <- (proc.time() - t0)["elapsed"]
cat(sprintf(
  "\n=== Cleanup Complete: %d batches in %.1f min ===\n",
  n_done,
  elapsed / 60
))
