# Study A reverse-order run (complements the forward-order job)
# Skips batches that already have result files to avoid collisions.
# Usage: Rscript study-a-reverse.R [n_cores]

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

# --- Build task groups in REVERSE order ---------------------------------------

design <- full_design("A")
all_tasks <- create_tasks(design)

task_groups <- split(
  all_tasks,
  sapply(all_tasks, function(t) {
    paste(t$dgp, t$method, t$study, sep = "_")
  })
)

# Reverse the method order (fast methods first = reverse of slow-first)
method_order <- rev(c(
  "fda_default",
  "fda_crit1",
  "srvf",
  "affine_ss",
  "landmark_auto"
))
batch_method <- sapply(task_groups, function(tl) tl[[1]]$method)
batch_order <- order(match(batch_method, method_order))

cat(sprintf("=== Registration Benchmark v2 (REVERSE) ===\n"))
cat(sprintf("Study: A | Cores: %d\n", n_cores))
cat(sprintf("Start: %s\n\n", format(Sys.time())))
cat(sprintf("Total batches: %d\n\n", length(task_groups)))

t0 <- proc.time()
all_results <- list()

for (idx in batch_order) {
  group_name <- names(task_groups)[idx]
  out_file <- file.path(results_dir, sprintf("results_%s.rds", group_name))

  if (file.exists(out_file)) {
    cat(sprintf(
      "[%s] SKIP %s (already exists)\n",
      format(Sys.time(), "%H:%M:%S"),
      group_name
    ))
    next
  }

  all_results[[group_name]] <- run_batch(
    task_groups[[group_name]],
    n_cores = n_cores
  )
}

elapsed <- (proc.time() - t0)["elapsed"]
cat(sprintf("\n=== Reverse Run Complete ===\n"))
cat(sprintf("Total time: %.1f minutes\n", elapsed / 60))
cat(sprintf("Batches completed: %d\n", length(all_results)))
