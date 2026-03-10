# sim-pilot-oracle.R -- Study D pilot: oracle template comparison
#
# Smoke test for Study D with reduced DGPs and reps.
# 3 DGPs (D01 easy, D03 wiggly, D10 complex+highrank) x 2 template modes
# x 2 conditions x 5 methods x 5 reps = 300 runs
#
# Usage:
#   Rscript sim-pilot-oracle.R [n_cores]

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

pilot_oracle_design <- function() {
  # Subset of Study D: 3 DGPs, both conditions, both template modes
  dgps <- c("D01", "D03", "D10")

  # Two conditions: easy + hard
  conditions <- data.frame(
    noise_sd = c(0.1, 0.3),
    severity = c(0.5, 1.0),
    stringsAsFactors = FALSE
  )

  grids <- list()
  for (i in seq_len(nrow(conditions))) {
    for (oracle in c(TRUE, FALSE)) {
      g <- expand.grid(
        dgp = dgps,
        n_curves = 50,
        n_grid = 101,
        noise_sd = conditions$noise_sd[i],
        severity = conditions$severity[i],
        use_true_template = oracle,
        stringsAsFactors = FALSE
      )
      grids[[length(grids) + 1]] <- g
    }
  }

  grid <- do.call(rbind, grids)
  grid <- expand_methods(grid)
  grid$reps <- 5
  grid$study <- "D_pilot"
  grid$contam_frac <- NA_real_
  grid$outlier_type <- NA_character_
  grid$lambda <- NA_real_
  grid
}

# --- Runner -------------------------------------------------------------------

run_oracle_pilot <- function(n_cores = 4) {
  cat("=== Study D Pilot: Oracle Template Comparison ===\n")
  cat(sprintf("Cores: %d\n", n_cores))
  cat(sprintf("Start: %s\n\n", format(Sys.time())))

  design <- pilot_oracle_design()
  tasks <- create_tasks(design)
  cat(sprintf("Total tasks: %d\n\n", length(tasks)))

  # Group by DGP + method + oracle mode for incremental saves
  task_groups <- split(
    tasks,
    sapply(tasks, function(t) {
      paste(t$dgp, t$method, t$use_true_template, t$study, sep = "_")
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
  out_file <- file.path(results_dir, "results_pilot_oracle.rds")
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

  # Sanity check: template MISE should be ~0 in oracle mode
  cat("\n--- Template MISE: oracle vs estimated ---\n")
  tmise <- dt[,
    .(template_mise = median(template_mise, na.rm = TRUE)),
    by = .(use_true_template, method, dgp)
  ]
  tmise_wide <- dcast(
    tmise,
    method + dgp ~ use_true_template,
    value.var = "template_mise"
  )
  setnames(tmise_wide, c("FALSE", "TRUE"), c("estimated", "oracle"))
  print(tmise_wide[order(dgp, method)])

  # Oracle benefit: warp MISE reduction
  cat("\n--- Oracle benefit: warp MISE (estimated - oracle) ---\n")
  wmise <- dt[,
    .(warp_mise = median(warp_mise, na.rm = TRUE)),
    by = .(use_true_template, method, dgp, noise_sd, severity)
  ]
  wmise_wide <- dcast(
    wmise,
    method + dgp + noise_sd + severity ~ use_true_template,
    value.var = "warp_mise"
  )
  setnames(wmise_wide, c("FALSE", "TRUE"), c("estimated", "oracle"))
  wmise_wide[, oracle_benefit := estimated - oracle]
  wmise_wide[, pct_improvement := 100 * oracle_benefit / estimated]
  print(wmise_wide[
    order(dgp, noise_sd, method),
    .(method, dgp, noise_sd, severity, estimated, oracle, pct_improvement)
  ])

  # Per-method oracle benefit summary
  cat("\n--- Median oracle benefit (% improvement) by method ---\n")
  method_summary <- wmise_wide[,
    .(
      median_pct = median(pct_improvement, na.rm = TRUE),
      min_pct = min(pct_improvement, na.rm = TRUE),
      max_pct = max(pct_improvement, na.rm = TRUE)
    ),
    by = .(method)
  ]
  print(method_summary[order(-median_pct)])

  # Landmark sanity: oracle benefit should be ~0
  cat("\n--- Landmark oracle benefit (expect ~0): ---\n")
  lm_benefit <- wmise_wide[method == "landmark_auto"]
  if (nrow(lm_benefit) > 0) {
    cat(sprintf(
      "  Median pct improvement: %.2f%%\n",
      median(lm_benefit$pct_improvement, na.rm = TRUE)
    ))
  }

  invisible(results_df)
}

# --- CLI ----------------------------------------------------------------------

if (sys.nframe() == 0) {
  args <- commandArgs(trailingOnly = TRUE)
  n_cores <- if (length(args) >= 1) as.integer(args[1]) else 4
  run_oracle_pilot(n_cores = n_cores)
}
