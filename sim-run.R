# sim-run.R -- Main Runner for Registration Benchmark v2
#
# Usage:
#   Rscript sim-run.R [study] [n_cores]
#
# Arguments:
#   study:   "A", "B", "C", "D", "pilot", "all"
#   n_cores: number of parallel cores (default: 4)
#
# Results saved incrementally per DGP to results/

library(parallel)

`%||%` <- function(x, y) if (is.null(x)) y else x

base_dir <- here::here()
results_dir <- file.path(base_dir, "results")
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

source(file.path(base_dir, "sim-dgp.R"))
source(file.path(base_dir, "sim-methods.R"))
source(file.path(base_dir, "sim-metrics.R"))
source(file.path(base_dir, "sim-config.R"))

# --- Task Definition ----------------------------------------------------------

#' Create task list from design table
create_tasks <- function(design) {
  tasks <- list()
  for (i in seq_len(nrow(design))) {
    row <- design[i, ]
    for (rep in seq_len(row$reps)) {
      task <- list(
        dgp = row$dgp,
        n_curves = row$n_curves,
        n_grid = row$n_grid,
        noise_sd = row$noise_sd,
        severity = row$severity,
        method = row$method,
        rep = rep,
        study = row$study,
        seed = make_seed(row$dgp, rep),
        use_true_template = row$use_true_template %||% FALSE
      )
      # Study B: lambda
      if (!is.na(row$lambda %||% NA)) {
        task$lambda <- row$lambda
      }
      # Study D: contamination
      if (!is.na(row$contam_frac %||% NA)) {
        task$contam_frac <- row$contam_frac
        task$outlier_type <- row$outlier_type
      }
      tasks[[length(tasks) + 1]] <- task
    }
  }
  tasks
}

# --- Single Task Runner -------------------------------------------------------

#' Metric column names (flat)
metric_cols <- c(
  "warp_mise",
  "alignment_error",
  "template_mise",
  "alignment_cc",
  "time_per_curve",
  "registration_ratio_median",
  "registration_ratio_iqr",
  "amp_variance_ratio",
  "warp_slope_ratio_median",
  "warp_slope_ratio_iqr",
  "time",
  "failure"
)

#' Run one task
run_one_task <- function(task) {
  # Generate data
  data <- tryCatch(
    generate_data(
      dgp = task$dgp,
      n = task$n_curves,
      n_grid = task$n_grid,
      severity = task$severity,
      noise_sd = task$noise_sd,
      seed = task$seed
    ),
    error = function(e) list(error = conditionMessage(e))
  )

  # Apply contamination (Study D)
  if (!is.null(data$x) && !is.null(task$contam_frac)) {
    data <- tryCatch(
      contaminate_data(
        data,
        contam_frac = task$contam_frac,
        outlier_type = task$outlier_type,
        seed = task$seed
      ),
      error = function(e) {
        data$error <- conditionMessage(e)
        data
      }
    )
  }

  if (!is.null(data$error)) {
    return(make_result_row(task, failure_metrics(list()), data$error))
  }

  # Fit method + extract metrics
  metrics <- tryCatch(
    {
      result <- fit_method(
        data,
        task$method,
        use_true_template = task$use_true_template %||% FALSE,
        lambda = task$lambda
      )
      m <- extract_metrics(data, result)
      m$error_msg <- result$error %||% ""
      m
    },
    error = function(e) {
      m <- failure_metrics(list())
      m$error_msg <- conditionMessage(e)
      m
    }
  )

  make_result_row(task, metrics, metrics$error_msg %||% "")
}

#' Assemble result row
#' @keywords internal
make_result_row <- function(task, metrics, error_msg = "") {
  flat <- flatten_metrics(metrics)
  data.frame(
    dgp = task$dgp,
    n_curves = task$n_curves,
    n_grid = task$n_grid,
    noise_sd = task$noise_sd,
    severity = task$severity,
    method = task$method,
    rep = task$rep,
    seed = task$seed,
    study = task$study,
    use_true_template = task$use_true_template %||% FALSE,
    lambda = task$lambda %||% NA_real_,
    contam_frac = task$contam_frac %||% NA_real_,
    outlier_type = task$outlier_type %||% NA_character_,
    warp_mise = flat$warp_mise %||% NA_real_,
    alignment_error = flat$alignment_error %||% NA_real_,
    template_mise = flat$template_mise %||% NA_real_,
    alignment_cc = flat$alignment_cc %||% NA_real_,
    time_per_curve = flat$time_per_curve %||% NA_real_,
    registration_ratio_median = flat$registration_ratio_median %||% NA_real_,
    registration_ratio_iqr = flat$registration_ratio_iqr %||% NA_real_,
    amp_variance_ratio = flat$amp_variance_ratio %||% NA_real_,
    warp_slope_ratio_median = flat$warp_slope_ratio_median %||% NA_real_,
    warp_slope_ratio_iqr = flat$warp_slope_ratio_iqr %||% NA_real_,
    time = flat$time %||% NA_real_,
    failure = flat$failure %||% TRUE,
    error_msg = error_msg,
    stringsAsFactors = FALSE
  )
}

# --- Batch Runner with Incremental Saves --------------------------------------

#' Run all tasks for one group, save results
run_batch <- function(tasks, n_cores = 4) {
  group_key <- paste(
    tasks[[1]]$dgp,
    tasks[[1]]$method,
    tasks[[1]]$study,
    sep = "_"
  )

  out_file <- file.path(
    results_dir,
    sprintf("results_%s.rds", group_key)
  )

  cat(sprintf(
    "[%s] Running %d tasks for %s...\n",
    format(Sys.time(), "%H:%M:%S"),
    length(tasks),
    group_key
  ))

  if (n_cores > 1) {
    results <- mclapply(tasks, run_one_task, mc.cores = n_cores)
  } else {
    results <- lapply(tasks, run_one_task)
  }

  # Handle worker crashes
  failed <- vapply(
    results,
    function(r) {
      inherits(r, "try-error") || is.null(r) || !is.data.frame(r)
    },
    logical(1)
  )

  if (any(failed)) {
    warning(sprintf(
      "%d/%d tasks crashed in %s",
      sum(failed),
      length(tasks),
      group_key
    ))
    for (i in which(failed)) {
      results[[i]] <- make_result_row(
        tasks[[i]],
        failure_metrics(list()),
        "Worker crash"
      )
    }
  }

  results_df <- do.call(rbind, results)
  saveRDS(results_df, out_file)

  cat(sprintf(
    "[%s] Saved %d rows to %s\n",
    format(Sys.time(), "%H:%M:%S"),
    nrow(results_df),
    basename(out_file)
  ))
  results_df
}

# --- Main Runner --------------------------------------------------------------

#' Run the benchmark
#'
#' @param study "A", "B", "C", "D", "pilot", "all"
#' @param n_cores number of parallel cores
run_benchmark <- function(study = "pilot", n_cores = 4) {
  cat(sprintf("=== Registration Benchmark v2 ===\n"))
  cat(sprintf("Study: %s | Cores: %d\n", study, n_cores))
  cat(sprintf("Start: %s\n\n", format(Sys.time())))

  design <- switch(
    study,
    A = full_design("A"),
    B = full_design("B"),
    C = full_design("C"),
    D = full_design("D"),
    pilot = {
      d <- full_design("A")
      d$reps <- 10
      d
    },
    all = full_design(c("A", "B", "C", "D")),
    cli::cli_abort("Unknown study: {study}")
  )

  all_tasks <- create_tasks(design)
  cat(sprintf("Total tasks: %d\n\n", length(all_tasks)))

  # Group by DGP + method + study for incremental saves and
  # homogeneous task runtimes within each batch
  task_groups <- split(
    all_tasks,
    sapply(all_tasks, function(t) {
      paste(t$dgp, t$method, t$study, sep = "_")
    })
  )

  all_results <- list()
  t0 <- proc.time()

  # Sort batches so slow methods (fda_default) start first for better
  # overall load balancing across the sequential batch loop
  method_order <- c(
    "fda_default",
    "fda_crit1",
    "srvf",
    "affine_ss",
    "landmark_auto"
  )
  batch_method <- sapply(task_groups, function(tl) tl[[1]]$method)
  batch_order <- order(match(batch_method, method_order))

  for (idx in batch_order) {
    group_name <- names(task_groups)[idx]
    all_results[[group_name]] <- run_batch(
      task_groups[[group_name]],
      n_cores = n_cores
    )
  }

  final_results <- do.call(rbind, all_results)
  rownames(final_results) <- NULL

  combined_file <- file.path(
    results_dir,
    sprintf("results_combined_%s.rds", study)
  )
  saveRDS(final_results, combined_file)

  elapsed <- (proc.time() - t0)["elapsed"]
  cat(sprintf("\n=== Benchmark Complete ===\n"))
  cat(sprintf("Total time: %.1f minutes\n", elapsed / 60))
  cat(sprintf("Total rows: %d\n", nrow(final_results)))
  cat(sprintf("Failure rate: %.1f%%\n", 100 * mean(final_results$failure)))
  cat(sprintf("Results: %s\n", combined_file))

  writeLines(
    capture.output(sessionInfo()),
    file.path(results_dir, "sessionInfo.txt")
  )

  invisible(final_results)
}

# --- CLI Entry Point ----------------------------------------------------------

if (sys.nframe() == 0) {
  args <- commandArgs(trailingOnly = TRUE)
  study <- if (length(args) >= 1) args[1] else "pilot"
  n_cores <- if (length(args) >= 2) as.integer(args[2]) else 4
  run_benchmark(study = study, n_cores = n_cores)
}
