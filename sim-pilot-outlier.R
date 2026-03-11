# sim-pilot-outlier.R -- Study E pilot: outlier contamination
#
# Smoke test for Study E with reduced fractions and reps.
# 3 DGPs x 2 outlier types x 2 fractions (10%, 30%) x 5 methods
# x 2 noise x severity=0.5 x 5 reps = 600 runs
#
# Verifies:
#   - contaminate_data() runs without error for all DGP x outlier combos
#   - Seed nesting: 10% outlier set is subset of 30% set (paired fractions)
#   - Failure rates are reasonable (not inflated by contamination)
#   - Warp MISE increases monotonically with contamination fraction
#   - Shape vs phase outliers have differential effects
#
# Usage:
#   Rscript sim-pilot-outlier.R [n_cores]

library(parallel)

base_dir <- here::here()
results_dir <- file.path(base_dir, "results")
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

source(file.path(base_dir, "sim-dgp.R"))
source(file.path(base_dir, "sim-methods.R"))
source(file.path(base_dir, "sim-metrics.R"))
source(file.path(base_dir, "sim-config.R"))
source(file.path(base_dir, "sim-run.R"))

# --- Design -------------------------------------------------------------------

pilot_outlier_design <- function() {
  grid <- expand.grid(
    dgp = c("D02", "D09", "D12"),
    n_curves = 50,
    n_grid = 101,
    noise_sd = c(0.1, 0.3),
    severity = 0.5,
    contam_frac = c(0.10, 0.30),
    outlier_type = c("shape", "phase"),
    stringsAsFactors = FALSE
  )
  grid <- expand_methods(grid)
  grid$reps <- 5
  grid$study <- "E_pilot"
  grid$use_true_template <- FALSE
  grid$lambda <- NA_real_
  grid
}

# --- Runner -------------------------------------------------------------------

run_outlier_pilot <- function(n_cores = 4) {
  cat("=== Study E Pilot: Outlier Contamination ===\n")
  cat(sprintf("Cores: %d\n", n_cores))
  cat(sprintf("Start: %s\n\n", format(Sys.time())))

  design <- pilot_outlier_design()
  tasks <- create_tasks(design)
  cat(sprintf("Total tasks: %d\n\n", length(tasks)))

  # Group by DGP + method + outlier_type for incremental saves
  task_groups <- split(
    tasks,
    sapply(tasks, function(t) {
      paste(t$dgp, t$method, t$outlier_type, t$study, sep = "_")
    })
  )

  all_results <- list()
  t0 <- proc.time()

  for (group_name in names(task_groups)) {
    all_results[[group_name]] <- run_batch(
      task_groups[[group_name]],
      n_cores = n_cores
    )
  }

  results_df <- do.call(rbind, all_results)
  rownames(results_df) <- NULL
  out_file <- file.path(results_dir, "results_pilot_outlier.rds")
  saveRDS(results_df, out_file)

  elapsed <- (proc.time() - t0)["elapsed"]
  cat(sprintf("\n=== Pilot Complete ===\n"))
  cat(sprintf("Time: %.1f minutes\n", elapsed / 60))
  cat(sprintf(
    "Rows: %d | Failures: %d (%.1f%%)\n",
    nrow(results_df),
    sum(results_df$failure),
    100 * mean(results_df$failure)
  ))
  cat(sprintf("Results: %s\n", out_file))

  # --- Quick diagnostics -----------------------------------------------------
  library(data.table)
  dt <- as.data.table(results_df)

  # Failure rates by contamination and method
  cat("\n--- Failure rates by method x outlier_type x contam_frac ---\n")
  fail_summary <- dt[,
    .(
      n_total = .N,
      n_fail = sum(failure),
      rate_pct = sprintf("%.1f%%", 100 * mean(failure))
    ),
    by = .(method, outlier_type, contam_frac)
  ]
  print(fail_summary[order(method, outlier_type, contam_frac)])

  # Warp MISE by contamination fraction (should increase with fraction)
  cat("\n--- Median warp MISE by method x outlier_type x contam_frac ---\n")
  dt_ok <- dt[failure == FALSE]
  wmise <- dt_ok[,
    .(warp_mise = median(warp_mise, na.rm = TRUE)),
    by = .(method, outlier_type, contam_frac, dgp)
  ]
  wmise_wide <- dcast(
    wmise,
    method + outlier_type + dgp ~ contam_frac,
    value.var = "warp_mise"
  )
  setnames(wmise_wide, c("0.1", "0.3"), c("frac_10", "frac_30"))
  wmise_wide[, ratio_30_10 := frac_30 / frac_10]
  print(wmise_wide[order(method, outlier_type, dgp)])

  # Shape vs phase comparison
  cat("\n--- Shape vs phase: median warp MISE ratio ---\n")
  type_comparison <- dt_ok[,
    .(warp_mise = median(warp_mise, na.rm = TRUE)),
    by = .(method, outlier_type, dgp, noise_sd)
  ]
  type_wide <- dcast(
    type_comparison,
    method + dgp + noise_sd ~ outlier_type,
    value.var = "warp_mise"
  )
  type_wide[, shape_phase_ratio := shape / phase]
  cat("Median shape/phase ratio by method:\n")
  print(type_wide[,
    .(median_ratio = median(shape_phase_ratio, na.rm = TRUE)),
    by = method
  ])

  # Template bias: does contamination bias template estimation?
  cat("\n--- Template MISE by outlier_type x contam_frac ---\n")
  tmise <- dt_ok[,
    .(template_mise = median(template_mise, na.rm = TRUE)),
    by = .(method, outlier_type, contam_frac)
  ]
  tmise_wide <- dcast(
    tmise,
    method + outlier_type ~ contam_frac,
    value.var = "template_mise"
  )
  print(tmise_wide[order(method, outlier_type)])

  invisible(results_df)
}

# --- CLI ----------------------------------------------------------------------

if (sys.nframe() == 0) {
  args <- commandArgs(trailingOnly = TRUE)
  n_cores <- if (length(args) >= 1) as.integer(args[1]) else 4
  run_outlier_pilot(n_cores = n_cores)
}
