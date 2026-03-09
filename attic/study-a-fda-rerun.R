# Study A FDA rerun -- uses new tf-native FDA backend
# Reruns fda_default (crit=2) and fda_crit1 (crit=1) for all 13 DGPs.
# Usage: Rscript study-a-fda-rerun.R [n_cores]

library(tf)
library(parallel)

base_dir <- here::here("attic", "sim-registration")
results_dir <- file.path(base_dir, "results")
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

source(file.path(base_dir, "sim-dgp.R"))
source(file.path(base_dir, "sim-methods.R"))
source(file.path(base_dir, "sim-metrics.R"))
source(file.path(base_dir, "sim-config.R"))
source(file.path(base_dir, "sim-run.R"))

args <- commandArgs(trailingOnly = TRUE)
n_cores <- if (length(args) >= 1) as.integer(args[1]) else 4L

cat("=== FDA Rerun (tf-native backend) ===\n")
cat(sprintf("tf version: %s\n", packageVersion("tf")))
cat(sprintf("Cores: %d\n\n", n_cores))

# Build full Study A design, then filter to FDA methods only
design <- full_design("A")
all_tasks <- create_tasks(design)
fda_tasks <- Filter(
  function(t) t$method %in% c("fda_default", "fda_crit1"),
  all_tasks
)

cat(sprintf(
  "FDA tasks: %d (of %d total)\n\n",
  length(fda_tasks),
  length(all_tasks)
))

# Group by DGP + method for incremental saves
task_groups <- split(
  fda_tasks,
  sapply(fda_tasks, function(t) {
    paste(t$dgp, t$method, t$study, sep = "_")
  })
)

t0 <- proc.time()
for (group_name in sort(names(task_groups))) {
  run_batch(task_groups[[group_name]], n_cores = n_cores)
}

elapsed <- (proc.time() - t0)["elapsed"]
cat(sprintf("\n=== FDA Rerun Complete: %.1f minutes ===\n", elapsed / 60))
