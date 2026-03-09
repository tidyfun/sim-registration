# Smoke test: 2 reps of D01, single core
# Run from tf repo root

library(tf)
library(parallel)

base_dir <- here::here()
results_dir <- file.path(base_dir, "results")
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

source(file.path(base_dir, "sim-dgp.R"))
source(file.path(base_dir, "sim-methods.R"))
source(file.path(base_dir, "sim-config.R"))
source(file.path(base_dir, "sim-metrics.R"))
source(file.path(base_dir, "sim-run.R"))

# Override design: just D01, 2 reps, all 5 methods
design <- study_a_design()
design <- design[
  design$dgp == "D01" & design$noise_sd == 0 & design$severity == 0.5,
]
design$reps <- 2

cat(sprintf(
  "Smoke test: %d cells, %d total runs\n",
  nrow(design),
  sum(design$reps)
))

tasks <- create_tasks(design)
cat(sprintf("Tasks created: %d\n", length(tasks)))

t0 <- proc.time()
results <- lapply(tasks, run_one_task)
results_df <- do.call(rbind, results)
elapsed <- (proc.time() - t0)["elapsed"]

cat(sprintf("\nCompleted in %.1f seconds\n", elapsed))
cat(sprintf("Failures: %d/%d\n", sum(results_df$failure), nrow(results_df)))
print(results_df[, c(
  "method",
  "rep",
  "warp_mise",
  "alignment_error",
  "failure"
)])

saveRDS(results_df, file.path(results_dir, "smoke_test.rds"))
cat("Smoke test done.\n")
