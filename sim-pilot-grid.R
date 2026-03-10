# sim-pilot-grid.R -- Grid-size sensitivity pilot
#
# Compares n_grid=101 (current) vs n_grid=201 for all DGPs/methods.
# 1 paired rep at severity=1.0, noise=0.1.
# Total: 15 DGPs x 5 methods x 2 grids x 1 rep = 150 runs.
#
# Usage:
#   Rscript sim-pilot-grid.R [n_cores]

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

pilot_grid_design <- function() {
  grid <- expand.grid(
    dgp = paste0("D", sprintf("%02d", 1:15)),
    n_curves = 50,
    n_grid = c(101, 201),
    noise_sd = 0.1,
    severity = 1.0,
    stringsAsFactors = FALSE
  )
  grid <- expand_methods(grid)
  grid$reps <- 1
  grid$study <- "grid_pilot"
  grid$use_true_template <- FALSE
  grid$contam_frac <- NA_real_
  grid$outlier_type <- NA_character_
  grid$lambda <- NA_real_
  grid
}

# --- Runner -------------------------------------------------------------------

run_grid_pilot <- function(n_cores = 4) {
  cat("=== Grid-Size Sensitivity Pilot ===\n")
  cat(sprintf("Cores: %d\n", n_cores))
  cat(sprintf("Start: %s\n\n", format(Sys.time())))

  design <- pilot_grid_design()
  tasks <- create_tasks(design)
  cat(sprintf("Total tasks: %d\n\n", length(tasks)))

  t0 <- proc.time()

  if (n_cores > 1) {
    results <- mclapply(tasks, run_one_task, mc.cores = n_cores)
  } else {
    results <- lapply(tasks, run_one_task)
  }

  # Handle crashes
  failed <- vapply(
    results,
    function(r) {
      inherits(r, "try-error") || is.null(r) || !is.data.frame(r)
    },
    logical(1)
  )
  if (any(failed)) {
    warning(sprintf("%d/%d tasks crashed", sum(failed), length(tasks)))
    for (i in which(failed)) {
      results[[i]] <- make_result_row(
        tasks[[i]],
        failure_metrics(list()),
        "Worker crash"
      )
    }
  }

  results_df <- do.call(rbind, results)
  out_file <- file.path(results_dir, "results_grid_pilot.rds")
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

  # Quick comparison
  library(data.table)
  dt <- as.data.table(results_df)
  dt <- dt[failure == FALSE]
  comp <- dt[,
    .(
      warp_mise = median(warp_mise, na.rm = TRUE),
      template_mise = median(template_mise, na.rm = TRUE),
      alignment_error = median(alignment_error, na.rm = TRUE)
    ),
    by = .(n_grid, method)
  ]
  comp_wide <- dcast(
    melt(comp, id.vars = c("n_grid", "method")),
    method + variable ~ n_grid,
    value.var = "value"
  )
  setnames(comp_wide, c("101", "201"), c("grid_101", "grid_201"))
  comp_wide[, ratio := grid_201 / grid_101]
  cat("\n--- Grid 201 / Grid 101 ratio (per method, pooled over DGPs) ---\n")
  print(comp_wide[order(variable, method)])

  # Per-DGP paired comparison
  paired <- dt[,
    .(warp_mise = median(warp_mise, na.rm = TRUE)),
    by = .(dgp, method, n_grid)
  ]
  paired_wide <- dcast(paired, dgp + method ~ n_grid, value.var = "warp_mise")
  setnames(paired_wide, c("101", "201"), c("grid_101", "grid_201"))
  paired_wide[, ratio := grid_201 / grid_101]
  cat("\n--- Per-DGP warp MISE ratio (grid_201 / grid_101) ---\n")
  cat(sprintf(
    "Median ratio: %.4f | Range: [%.4f, %.4f]\n",
    median(paired_wide$ratio, na.rm = TRUE),
    min(paired_wide$ratio, na.rm = TRUE),
    max(paired_wide$ratio, na.rm = TRUE)
  ))
  cat(sprintf(
    "Wilcoxon signed-rank test p-value: %.4f\n",
    wilcox.test(
      paired_wide$grid_101,
      paired_wide$grid_201,
      paired = TRUE
    )$p.value
  ))

  invisible(results_df)
}

# --- CLI ----------------------------------------------------------------------

if (sys.nframe() == 0) {
  args <- commandArgs(trailingOnly = TRUE)
  n_cores <- if (length(args) >= 1) as.integer(args[1]) else 4
  run_grid_pilot(n_cores = n_cores)
}
