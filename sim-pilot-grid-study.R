# sim-pilot-grid-study.R -- Study C pilot: grid resolution sensitivity
#
# Smoke test for Study C with reduced DGPs and reps.
# 3 DGPs (D01 easy, D04 affine, D12 hard) x 3 grids x 2 noise x 5 methods x 5 reps
# = 450 runs
#
# Usage:
#   Rscript sim-pilot-grid-study.R [n_cores]

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

pilot_grid_study_design <- function() {
  # Subset of Study C: 3 DGPs, 2 noise levels, 5 reps
  grid <- expand.grid(
    dgp = c("D01", "D04", "D12"),
    n_curves = 50,
    n_grid = c(51L, 101L, 201L),
    noise_sd = c(0, 0.3),
    severity = 1.0,
    stringsAsFactors = FALSE
  )
  grid <- expand_methods(grid)
  grid$reps <- 5
  grid$study <- "C_pilot"
  grid$use_true_template <- FALSE
  grid$contam_frac <- NA_real_
  grid$outlier_type <- NA_character_
  grid$lambda <- NA_real_
  grid
}

# --- Runner -------------------------------------------------------------------

run_grid_study_pilot <- function(n_cores = 4) {
  cat("=== Study C Pilot: Grid Resolution Sensitivity ===\n")
  cat(sprintf("Cores: %d\n", n_cores))
  cat(sprintf("Start: %s\n\n", format(Sys.time())))

  design <- pilot_grid_study_design()
  tasks <- create_tasks(design)
  cat(sprintf("Total tasks: %d\n\n", length(tasks)))

  # Group by DGP + method + grid for incremental saves
  task_groups <- split(
    tasks,
    sapply(tasks, function(t) {
      paste(t$dgp, t$method, t$n_grid, t$study, sep = "_")
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
  out_file <- file.path(results_dir, "results_pilot_grid_study.rds")
  saveRDS(results_df, out_file)

  elapsed <- (proc.time() - t0)["elapsed"]
  cat(sprintf("\n=== Pilot Complete ===\n"))
  cat(sprintf("Time: %.1f minutes\n", elapsed / 60))
  cat(sprintf(
    "Rows: %d | Failures: %d\n",
    nrow(results_df),
    sum(results_df$failure)
  ))
  cat(sprintf("Results: %s\n", out_file))

  # --- Quick diagnostics -----------------------------------------------------
  library(data.table)
  dt <- as.data.table(results_df)
  dt <- dt[failure == FALSE]

  # Grid comparison by method (pooled over DGPs)
  comp <- dt[,
    .(
      warp_mise = median(warp_mise, na.rm = TRUE),
      template_mise = median(template_mise, na.rm = TRUE),
      alignment_error = median(alignment_error, na.rm = TRUE),
      elastic_dist = median(template_elastic_dist, na.rm = TRUE)
    ),
    by = .(n_grid, method, noise_sd)
  ]

  cat("\n--- Median metrics by grid size x method (noise=0) ---\n")
  print(comp[
    noise_sd == 0,
    .(method, n_grid, warp_mise, template_mise, elastic_dist)
  ])

  cat("\n--- Median metrics by grid size x method (noise=0.3) ---\n")
  print(comp[
    noise_sd == 0.3,
    .(method, n_grid, warp_mise, template_mise, elastic_dist)
  ])

  # Grid ratios relative to grid=101
  baseline <- comp[n_grid == 101]
  setnames(
    baseline,
    c("warp_mise", "template_mise"),
    c("base_warp", "base_template"),
    skip_absent = TRUE
  )
  ratio_dt <- merge(
    comp[, .(method, noise_sd, n_grid, warp_mise, template_mise)],
    baseline[, .(method, noise_sd, base_warp, base_template)],
    by = c("method", "noise_sd")
  )
  ratio_dt[, `:=`(
    warp_ratio = warp_mise / base_warp,
    template_ratio = template_mise / base_template
  )]

  cat("\n--- Warp MISE ratio vs grid=101 ---\n")
  print(dcast(ratio_dt, method + noise_sd ~ n_grid, value.var = "warp_ratio"))

  cat("\n--- Template MISE ratio vs grid=101 ---\n")
  print(dcast(
    ratio_dt,
    method + noise_sd ~ n_grid,
    value.var = "template_ratio"
  ))

  # Timing
  cat("\n--- Median time per curve (seconds) by method x grid ---\n")
  timing <- dt[,
    .(time_per_curve = median(time_per_curve, na.rm = TRUE)),
    by = .(method, n_grid)
  ]
  print(dcast(timing, method ~ n_grid, value.var = "time_per_curve"))

  invisible(results_df)
}

# --- CLI ----------------------------------------------------------------------

if (sys.nframe() == 0) {
  args <- commandArgs(trailingOnly = TRUE)
  n_cores <- if (length(args) >= 1) as.integer(args[1]) else 4
  run_grid_study_pilot(n_cores = n_cores)
}
