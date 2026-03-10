# sim-pilot-lambda.R -- Lambda grid pilot for Study B
#
# Runs a fine lambda grid for srvf, fda_default, fda_crit1 on 2 DGPs
# to identify where the MISE-vs-lambda U-shape minimum lies for each method.
#
# Design: 2 DGPs x 10 lambdas x 3 methods x 3 noise x 10 reps = 1,800 runs
#
# Usage:
#   Rscript sim-pilot-lambda.R [n_cores]

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

pilot_lambda_design <- function() {
  lambdas <- c(0, 1e-4, 1e-3, 5e-3, 0.01, 0.05, 0.1, 0.5, 1, 10)
  methods <- c("srvf", "fda_default", "fda_crit1")

  grid <- expand.grid(
    dgp = c("D02", "D12"),
    n_curves = 50,
    n_grid = 101,
    noise_sd = c(0, 0.1, 0.3),
    severity = 1.0,
    lambda = lambdas,
    method = methods,
    stringsAsFactors = FALSE
  )
  grid$reps <- 10
  grid$study <- "lambda_pilot"
  grid$use_true_template <- FALSE
  grid$contam_frac <- NA_real_
  grid$outlier_type <- NA_character_
  grid
}

# --- Runner -------------------------------------------------------------------

run_lambda_pilot <- function(n_cores = 4) {
  cat("=== Lambda Grid Pilot (Study B) ===\n")
  cat(sprintf("Cores: %d\n", n_cores))
  cat(sprintf("Start: %s\n\n", format(Sys.time())))

  design <- pilot_lambda_design()
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
  out_file <- file.path(results_dir, "results_lambda_pilot.rds")
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

  # --- Analysis ---------------------------------------------------------------
  library(data.table)
  dt <- as.data.table(results_df)
  dt <- dt[failure == FALSE]

  # Per-method, per-DGP, per-noise optimal lambda
  cat(
    "\n--- Optimal lambda (lowest median warp MISE) per method × DGP × noise ---\n"
  )
  best <- dt[,
    .(warp_mise = median(warp_mise, na.rm = TRUE)),
    by = .(method, dgp, noise_sd, lambda)
  ]
  optimal <- best[, .SD[which.min(warp_mise)], by = .(method, dgp, noise_sd)]
  print(optimal[
    order(method, dgp, noise_sd),
    .(method, dgp, noise_sd, lambda, warp_mise)
  ])

  # MISE vs lambda profiles (wide format for quick comparison)
  cat("\n--- Warp MISE vs lambda (median, D02, noise=0.3) ---\n")
  profile <- best[dgp == "D02" & noise_sd == 0.3]
  profile_wide <- dcast(profile, lambda ~ method, value.var = "warp_mise")
  print(profile_wide)

  cat("\n--- Warp MISE vs lambda (median, D12, noise=0.3) ---\n")
  profile2 <- best[dgp == "D12" & noise_sd == 0.3]
  profile_wide2 <- dcast(profile2, lambda ~ method, value.var = "warp_mise")
  print(profile_wide2)

  # Ratio vs lambda=0 baseline
  cat("\n--- MISE ratio vs lambda=0 baseline (D02, noise=0.3) ---\n")
  baseline <- best[lambda == 0, .(method, dgp, noise_sd, base_mise = warp_mise)]
  ratio_dt <- merge(best, baseline, by = c("method", "dgp", "noise_sd"))
  ratio_dt[, ratio := warp_mise / base_mise]
  ratio_sub <- ratio_dt[dgp == "D02" & noise_sd == 0.3]
  ratio_wide <- dcast(ratio_sub, lambda ~ method, value.var = "ratio")
  print(ratio_wide)

  invisible(results_df)
}

# --- CLI ----------------------------------------------------------------------

if (sys.nframe() == 0) {
  args <- commandArgs(trailingOnly = TRUE)
  n_cores <- if (length(args) >= 1) as.integer(args[1]) else 4
  run_lambda_pilot(n_cores = n_cores)
}
