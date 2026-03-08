#!/usr/bin/env Rscript

# Informal sandbox for registration identifiability intuition:
# - compares SRVF vs FDA (crit = 1/2)
# - varies amplitude rank and phase severity
# - reports warp/alignment RMSE against known truth
#
# Important:
# template_mode = "method_default" uses each method's own default template rule:
# - SRVF: internal Karcher mean (fdasrvf::time_warping, template = NULL)
# - FDA: arithmetic mean
# template_mode = "true_euclidean_mean" fixes template to the true Euclidean
# mean of aligned curves for all methods.
#
# Run from repository root:
#   Rscript informal-identifiability-check.R
# Optional:
#   Rscript informal-identifiability-check.R 5
# where "5" is the number of replicates per scenario.

if (!file.exists("DESCRIPTION")) {
  stop(
    "Please run this script from the repository root (where DESCRIPTION lives)."
  )
}

required_pkgs <- c("fda", "fdasrvf", "devtools")
missing_pkgs <- required_pkgs[
  !vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_pkgs) > 0) {
  stop("Missing required packages: ", paste(missing_pkgs, collapse = ", "))
}

devtools::load_all(".", quiet = TRUE)

args <- commandArgs(trailingOnly = TRUE)
n_rep <- if (length(args) >= 1L) as.integer(args[1]) else 2L
if (!is.finite(n_rep) || n_rep < 1L) {
  stop("n_rep must be a positive integer.")
}

set.seed(20260303)

n_curve <- 12L
n_grid <- 81L
t_grid <- seq(0, 1, length.out = n_grid)
amp_ranks <- c(1L, 3L)
phase_sds <- c(0.20, 0.65)
template_modes <- c("method_default", "true_euclidean_mean")
dgp_families <- c("sinusoidal", "asymmetric")
methods <- c("srvf", "fda_crit1", "fda_crit2")
method_labels <- c(
  srvf = "SRVF",
  fda_crit1 = "FDA crit = 1",
  fda_crit2 = "FDA crit = 2"
)
phase_ratio_min_fraction <- 0.05
k_eig <- 5L

basis_library_sinusoidal <- function(t) {
  rbind(
    sin(2 * pi * t),
    cos(2 * pi * t),
    sin(4 * pi * t),
    cos(4 * pi * t)
  )
}

basis_library_asymmetric <- function(t) {
  rbind(
    (t - 0.5) * dnorm(t, mean = 0.6, sd = 0.2),
    sin(pi * t) * (t - 0.2),
    cos(2 * pi * t) * t * (1 - t),
    (t - 0.7)^2 - mean((t - 0.7)^2)
  )
}

mean_curve <- function(t, family) {
  if (family == "sinusoidal") {
    return(0.8 * sin(2 * pi * t))
  }
  if (family == "asymmetric") {
    return(dnorm(t, 0.25, 0.08) + 0.7 * dnorm(t, 0.65, 0.12) + 0.4 * t)
  }
  stop("Unknown family: ", family)
}

simulate_aligned_curves <- function(n, t, amp_rank, family) {
  basis_mat <- switch(
    family,
    sinusoidal = basis_library_sinusoidal(t),
    asymmetric = basis_library_asymmetric(t)
  )

  coefs <- matrix(rnorm(n * amp_rank), nrow = n, ncol = amp_rank)
  scales <- 1 / sqrt(seq_len(amp_rank))
  x_var <- coefs %*%
    (diag(scales, amp_rank) %*% basis_mat[seq_len(amp_rank), , drop = FALSE])
  x_mean <- matrix(
    mean_curve(t, family),
    nrow = n,
    ncol = length(t),
    byrow = TRUE
  )
  x_mean + x_var
}

simulate_warps <- function(n, t, phase_sd) {
  # Smooth monotone perturbation with fixed endpoints:
  # h(t) = t + alpha * t(1-t)(2t-1), alpha in [-0.9, 0.9]
  shape <- t * (1 - t) * (2 * t - 1)
  alpha <- pmax(pmin(rnorm(n, sd = phase_sd), 0.9), -0.9)
  warp_cols <- vapply(alpha, \(a) t + a * shape, numeric(length(t)))
  t(warp_cols)
}

row_trapz <- function(mat, t) {
  dt <- matrix(diff(t), nrow = nrow(mat), ncol = ncol(mat) - 1, byrow = TRUE)
  rowSums((mat[, -1, drop = FALSE] + mat[, -ncol(mat), drop = FALSE]) * dt / 2)
}

functional_variance_mass <- function(mat, t) {
  var_t <- apply(mat, 2, var)
  sum((var_t[-1] + var_t[-length(var_t)]) * diff(t) / 2)
}

top_pve <- function(mat, k = 5L) {
  eig <- eigen(stats::cov(mat), symmetric = TRUE, only.values = TRUE)$values
  eig <- pmax(eig, 0)
  tot <- sum(eig)
  if (tot <= 0) {
    return(rep(NA_real_, k))
  }
  pve <- eig / tot
  if (length(pve) >= k) {
    pve[seq_len(k)]
  } else {
    c(pve, rep(0, k - length(pve)))
  }
}

estimate_warp <- function(y_obs, method, template_mode, true_template) {
  template_arg <- if (template_mode == "true_euclidean_mean") true_template else
    NULL

  if (method == "srvf") {
    return(suppressWarnings(tf_register(
      y_obs,
      method = "srvf",
      template = template_arg
    )))
  }
  if (method == "fda_crit1") {
    return(suppressWarnings(tf_register(
      y_obs,
      method = "fda",
      crit = 1,
      template = template_arg
    )))
  }
  if (method == "fda_crit2") {
    return(suppressWarnings(tf_register(
      y_obs,
      method = "fda",
      crit = 2,
      template = template_arg
    )))
  }

  stop("Unknown method: ", method)
}

scenario_run <- function(
  amp_rank,
  phase_sd,
  template_mode,
  family,
  rep_id,
  capture_viz = FALSE
) {
  x_true_mat <- simulate_aligned_curves(n_curve, t_grid, amp_rank, family)
  w_true_mat <- simulate_warps(n_curve, t_grid, phase_sd)

  x_true <- tfd(x_true_mat, arg = t_grid)
  w_true <- tfd(w_true_mat, arg = t_grid)
  y_obs <- tf_warp(x_true, w_true)
  y_obs_mat <- as.matrix(tfd(y_obs, arg = t_grid))

  true_template <- mean(x_true)
  observed_template_rmse <- sqrt(mean(
    (as.matrix(mean(y_obs)) - as.matrix(true_template))^2
  ))

  arg_mat <- matrix(t_grid, nrow = n_curve, ncol = n_grid, byrow = TRUE)
  true_warp_strength <- mean(row_trapz(abs(w_true_mat - arg_mat), t_grid))
  amp_var_true <- functional_variance_mass(x_true_mat, t_grid)
  amp_var_obs <- functional_variance_mass(y_obs_mat, t_grid)
  phase_den <- amp_var_obs - amp_var_true
  phase_ratio_defined <- phase_den > phase_ratio_min_fraction * amp_var_true
  pve_true <- top_pve(x_true_mat, k = k_eig)
  pve_obs <- top_pve(y_obs_mat, k = k_eig)
  pve_obs_to_true_l1 <- sum(abs(pve_obs - pve_true))

  rows <- vector("list", length(methods))
  names(rows) <- methods
  method_viz <- if (capture_viz) vector("list", length(methods)) else NULL
  if (capture_viz) {
    names(method_viz) <- methods
  }

  for (method_name in methods) {
    w_est <- estimate_warp(y_obs, method_name, template_mode, true_template)
    x_est <- suppressWarnings(tf_unwarp(y_obs, w_est))

    w_est_mat <- as.matrix(tfd(w_est, arg = t_grid))
    x_est_mat <- suppressWarnings(as.matrix(tfd(x_est, arg = t_grid)))

    warp_rmse <- sqrt(mean((w_est_mat - w_true_mat)^2))
    align_rmse <- sqrt(mean((x_est_mat - x_true_mat)^2, na.rm = TRUE))
    warp_strength_est <- mean(row_trapz(abs(w_est_mat - arg_mat), t_grid))
    amp_var_est <- functional_variance_mass(x_est_mat, t_grid)

    warp_strength_ratio <- warp_strength_est / true_warp_strength
    amp_preservation_ratio <- amp_var_est / amp_var_true
    phase_removal_ratio <- if (!phase_ratio_defined || abs(phase_den) < 1e-12) {
      NA_real_
    } else {
      (amp_var_obs - amp_var_est) / phase_den
    }
    pve_est <- top_pve(x_est_mat, k = k_eig)
    pve_to_true_l1 <- sum(abs(pve_est - pve_true))
    pve1_abs_error <- abs(pve_est[1] - pve_true[1])

    # Over/under alignment diagnostic from amplitude preservation:
    # - over: aligned curves too homogeneous (amplitude variance too low)
    # - under: residual phase variation left (amplitude variance too high)
    alignment_state <- if (amp_preservation_ratio < 0.90) {
      "over"
    } else if (amp_preservation_ratio > 1.10) {
      "under"
    } else {
      "balanced"
    }

    rows[[method_name]] <- data.frame(
      rep = rep_id,
      dgp_family = family,
      amp_rank = amp_rank,
      phase_sd = phase_sd,
      template_mode = template_mode,
      method = method_name,
      warp_rmse = warp_rmse,
      align_rmse = align_rmse,
      observed_template_rmse = observed_template_rmse,
      warp_strength_true = true_warp_strength,
      warp_strength_est = warp_strength_est,
      warp_strength_ratio = warp_strength_ratio,
      amp_var_true = amp_var_true,
      amp_var_obs = amp_var_obs,
      amp_var_est = amp_var_est,
      amp_preservation_ratio = amp_preservation_ratio,
      phase_removal_ratio = phase_removal_ratio,
      phase_ratio_defined = phase_ratio_defined,
      alignment_state = alignment_state,
      pve_to_true_l1 = pve_to_true_l1,
      pve1_abs_error = pve1_abs_error,
      pve_obs_to_true_l1 = pve_obs_to_true_l1
    )
    if (capture_viz) {
      method_viz[[method_name]] <- list(
        warp_mat = w_est_mat,
        aligned_mat = x_est_mat,
        warp_rmse = warp_rmse,
        align_rmse = align_rmse,
        warp_strength_ratio = warp_strength_ratio,
        amp_preservation_ratio = amp_preservation_ratio,
        alignment_state = alignment_state
      )
    }
  }

  metrics <- do.call(rbind, rows)
  if (!capture_viz) {
    return(metrics)
  }

  list(
    metrics = metrics,
    viz = list(
      dgp_family = family,
      amp_rank = amp_rank,
      phase_sd = phase_sd,
      template_mode = template_mode,
      rep = rep_id,
      y_obs_mat = y_obs_mat,
      x_true_mat = x_true_mat,
      w_true_mat = w_true_mat,
      method_results = method_viz
    )
  )
}

plot_scenario_viz <- function(viz_obj, curve_cols, lasagna_cols) {
  layout(t(matrix(seq_len(16L), nrow = 4L, ncol = 4L)))
  par(
    mar = c(3.2, 3.2, 2.2, 0.8),
    oma = c(0, 0, 2.6, 0),
    mgp = c(2.0, 0.6, 0)
  )

  plot(
    tfd(viz_obj$y_obs_mat, arg = t_grid),
    main = "Observed (Generated)",
    col = curve_cols,
    lwd = 1.2
  )
  plot(
    tfd(viz_obj$w_true_mat, arg = t_grid),
    main = "True Warps",
    ylab = "",
    points = FALSE,
    col = curve_cols,
    lwd = 1.2
  )
  abline(0, 1, lty = 3)
  plot(
    tfd(viz_obj$x_true_mat, arg = t_grid),
    main = "True Aligned",
    col = curve_cols,
    lwd = 1.2,
    points = FALSE
  )
  plot(
    tfd(viz_obj$x_true_mat, arg = t_grid),
    main = "True Aligned (Lasagna)",
    type = "lasagna",
    col = lasagna_cols
  )

  for (method_name in methods) {
    method_obj <- viz_obj$method_results[[method_name]]
    label <- method_labels[[method_name]]
    plot.new()
    text(
      x = 0.02,
      y = 0.98,
      adj = c(0, 1),
      cex = 0.88,
      labels = paste(
        label,
        paste0("warp RMSE: ", sprintf("%.3f", method_obj$warp_rmse)),
        paste0("align RMSE: ", sprintf("%.3f", method_obj$align_rmse)),
        paste0(
          "warp strength ratio: ",
          sprintf("%.2f", method_obj$warp_strength_ratio)
        ),
        paste0(
          "amp preservation ratio: ",
          sprintf("%.2f", method_obj$amp_preservation_ratio)
        ),
        paste0("state: ", method_obj$alignment_state),
        sep = "\n"
      )
    )
    plot(
      tfd(method_obj$warp_mat, arg = t_grid),
      main = paste(label, "Warps"),
      ylab = "",
      points = FALSE,
      col = curve_cols,
      lwd = 1.2
    )
    abline(0, 1, lty = 3)
    plot(
      tfd(method_obj$aligned_mat, arg = t_grid),
      main = paste(label, "Aligned"),
      col = curve_cols,
      lwd = 1.2,
      points = FALSE
    )
    plot(
      tfd(method_obj$aligned_mat, arg = t_grid),
      main = paste(label, "Lasagna"),
      type = "lasagna",
      col = lasagna_cols
    )
  }

  mtext(
    paste0(
      "family = ",
      viz_obj$dgp_family,
      ", amp_rank = ",
      viz_obj$amp_rank,
      ", phase_sd = ",
      viz_obj$phase_sd,
      ", template_mode = ",
      viz_obj$template_mode,
      ", rep = ",
      viz_obj$rep
    ),
    outer = TRUE,
    cex = 1
  )
}

cat("\n=== Informal registration sandbox ===\n")
cat("replicates per scenario:", n_rep, "\n")
cat("n_curve:", n_curve, " n_grid:", n_grid, "\n\n")

res_list <- list()
viz_list <- list()
idx <- 1L
viz_idx <- 1L
for (amp_rank in amp_ranks) {
  for (phase_sd in phase_sds) {
    for (family in dgp_families) {
      for (template_mode in template_modes) {
        for (rep_id in seq_len(n_rep)) {
          run_obj <- scenario_run(
            amp_rank,
            phase_sd,
            template_mode,
            family,
            rep_id,
            capture_viz = rep_id == 1L
          )
          if (rep_id == 1L) {
            res_list[[idx]] <- run_obj$metrics
            viz_list[[viz_idx]] <- run_obj$viz
            viz_idx <- viz_idx + 1L
          } else {
            res_list[[idx]] <- run_obj
          }
          idx <- idx + 1L
        }
        cat(
          "completed: family = ",
          family,
          ", amp_rank = ",
          amp_rank,
          ", phase_sd = ",
          phase_sd,
          ", template_mode = ",
          template_mode,
          "\n",
          sep = ""
        )
      }
    }
  }
}

results <- do.call(rbind, res_list)

median_na <- function(x) {
  if (all(is.na(x))) NA_real_ else median(x, na.rm = TRUE)
}

summary_tbl <- aggregate(
  cbind(
    warp_rmse,
    align_rmse,
    observed_template_rmse,
    warp_strength_ratio,
    amp_preservation_ratio,
    phase_removal_ratio,
    pve_to_true_l1,
    pve1_abs_error,
    pve_obs_to_true_l1
  ) ~
    dgp_family + amp_rank + phase_sd + template_mode + method,
  data = results,
  FUN = median_na
)

summary_tbl <- summary_tbl[
  order(
    summary_tbl$dgp_family,
    summary_tbl$amp_rank,
    summary_tbl$phase_sd,
    summary_tbl$template_mode,
    summary_tbl$method
  ),
]

cat("\n=== Median metrics by scenario/method ===\n")
print(summary_tbl, row.names = FALSE)

phase_ratio_coverage <- aggregate(
  phase_ratio_defined ~
    dgp_family + amp_rank + phase_sd + template_mode + method,
  data = results,
  FUN = mean
)

cat("\n=== phase_removal_ratio availability (fraction defined) ===\n")
print(phase_ratio_coverage, row.names = FALSE)

state_counts <- as.data.frame(table(
  results$dgp_family,
  results$amp_rank,
  results$phase_sd,
  results$template_mode,
  results$method,
  results$alignment_state
))
names(state_counts) <- c(
  "dgp_family",
  "amp_rank",
  "phase_sd",
  "template_mode",
  "method",
  "alignment_state",
  "count"
)
state_counts <- state_counts[state_counts$count > 0, ]

cat("\n=== Over/Under/Balanced counts (from amp_preservation_ratio) ===\n")
print(state_counts, row.names = FALSE)

cat("\nInterpretation guide:\n")
cat(
  "- template_mode = method_default is not a shared template: SRVF uses Karcher mean; FDA uses arithmetic mean.\n"
)
cat(
  "- warp_strength_ratio near 1 means estimated total warping magnitude matches generated magnitude.\n"
)
cat(
  "- amp_preservation_ratio near 1 means aligned functions preserve true amplitude variance.\n"
)
cat(
  "- alignment_state is based on amp_preservation_ratio: < 0.90 over-align, > 1.10 under-align.\n"
)
cat(
  "- phase_removal_ratio is only reported when phase variance is large enough (see availability table).\n"
)
cat(
  "- pve_to_true_l1 (top-",
  k_eig,
  " PVE L1 distance) quantifies amplitude eigenspectrum distortion after alignment.\n",
  sep = ""
)
cat(
  "- pve1_abs_error is absolute error in leading PVE proportion (rank-1 dominance mismatch).\n"
)
cat(
  "- Compare sinusoidal vs asymmetric families to check whether conclusions are DGP-specific.\n"
)

output_dir <- file.path("results")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
}

stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
latest_path <- file.path(output_dir, "informal-identifiability-latest.rds")
stamped_path <- file.path(
  output_dir,
  paste0("informal-identifiability-", stamp, ".rds")
)
viz_latest_path <- file.path(
  output_dir,
  "informal-identifiability-viz-latest.pdf"
)
viz_stamped_path <- file.path(
  output_dir,
  paste0("informal-identifiability-viz-", stamp, ".pdf")
)

curve_cols <- grDevices::adjustcolor(
  grDevices::hcl.colors(n_curve, palette = "Dark 3"),
  alpha.f = 0.65
)
lasagna_cols <- grDevices::hcl.colors(12, palette = "YlOrRd", rev = TRUE)

grDevices::pdf(viz_stamped_path, width = 13.5, height = 10.5, onefile = TRUE)
for (viz_obj in viz_list) {
  plot_scenario_viz(
    viz_obj,
    curve_cols = curve_cols,
    lasagna_cols = lasagna_cols
  )
}
grDevices::dev.off()
file.copy(viz_stamped_path, viz_latest_path, overwrite = TRUE)

out_obj <- list(
  meta = list(
    timestamp = stamp,
    n_rep = n_rep,
    n_curve = n_curve,
    n_grid = n_grid,
    amp_ranks = amp_ranks,
    phase_sds = phase_sds,
    template_modes = template_modes,
    dgp_families = dgp_families,
    methods = methods,
    method_labels = method_labels,
    k_eig = k_eig,
    phase_ratio_min_fraction = phase_ratio_min_fraction,
    viz_latest_path = viz_latest_path,
    viz_stamped_path = viz_stamped_path
  ),
  results = results,
  summary_tbl = summary_tbl,
  phase_ratio_coverage = phase_ratio_coverage,
  state_counts = state_counts,
  first_rep_scenarios = viz_list
)

saveRDS(out_obj, latest_path)
saveRDS(out_obj, stamped_path)

cat("\nSaved RDS objects:\n")
cat("- ", latest_path, "\n", sep = "")
cat("- ", stamped_path, "\n", sep = "")
cat("\nSaved visualization PDFs:\n")
cat("- ", viz_latest_path, "\n", sep = "")
cat("- ", viz_stamped_path, "\n", sep = "")
