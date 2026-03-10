# sim-metrics.R -- Performance Measures for Registration Benchmark v2
#
# All metrics use the tf API (tf_integrate, tf_derive, tf_crosscor).
# No manual trapez_int.
#
# Provides:
#   extract_metrics()       - main entry: all metrics for one result
#   Primary: warp_mise, alignment_error, template_mise
#   Secondary: alignment_cc, template_elastic_dist, time_per_curve, failure_rate
#   Diagnostics: registration_ratio, amplitude_variance_ratio, warp_slope_ratio

if (requireNamespace("tf", quietly = TRUE)) library(tf) else
  devtools::load_all()

# --- Main Entry Point ---------------------------------------------------------

#' Extract all performance metrics for one method result
#'
#' @param data list from generate_data()
#' @param result list from fit_method() with: registration, time, error
#' @return named list of metric values
extract_metrics <- function(data, result) {
  if (!is.null(result$error)) {
    return(failure_metrics(result))
  }

  reg <- result$registration
  aligned <- tf_aligned(reg)
  inv_warps_est <- tf_inv_warps(reg)
  template_est <- tf_template(reg)
  arg <- data$arg
  template_true <- data$template
  warps_true <- data$warps
  true_prewarp <- data$true_prewarp
  n <- length(aligned)

  # Compare forward warps for all methods and DGPs.
  # tf_registration stores inverse warps, so invert once here to recover the
  # estimated forward warp h(s) = t on the benchmark grid.
  warp_cmp_est <- safe_metric(function() tf_invert(inv_warps_est))
  warp_cmp_true <- warps_true

  metrics <- list(
    # Primary
    warp_mise = safe_metric(
      compute_warp_mise,
      warp_cmp_est,
      warp_cmp_true,
      arg
    ),
    alignment_error = safe_metric(
      compute_alignment_error,
      aligned,
      true_prewarp
    ),
    template_mise = safe_metric(
      compute_template_mise,
      template_est,
      template_true
    ),
    template_elastic_dist = safe_metric(
      compute_template_elastic_dist,
      template_est,
      template_true,
      arg
    ),
    # Secondary
    alignment_cc = safe_metric(
      compute_alignment_cc,
      aligned,
      template_true
    ),
    time_per_curve = result$time / n,
    # Diagnostics
    registration_ratio = safe_metric(
      compute_registration_ratio,
      warp_cmp_est,
      warp_cmp_true,
      arg
    ),
    amp_variance_ratio = safe_metric(
      compute_amp_variance_ratio,
      aligned,
      template_est,
      true_prewarp,
      template_true
    ),
    warp_slope_ratio = safe_metric(
      compute_warp_slope_ratio,
      warp_cmp_est,
      warp_cmp_true
    ),
    # Meta
    time = result$time,
    failure = FALSE
  )
  metrics
}

#' Return all-NA metrics for failed runs
#' @keywords internal
failure_metrics <- function(result) {
  list(
    warp_mise = NA_real_,
    alignment_error = NA_real_,
    template_mise = NA_real_,
    template_elastic_dist = NA_real_,
    alignment_cc = NA_real_,
    time_per_curve = NA_real_,
    registration_ratio_median = NA_real_,
    registration_ratio_iqr = NA_real_,
    amp_variance_ratio = NA_real_,
    warp_slope_ratio_median = NA_real_,
    warp_slope_ratio_iqr = NA_real_,
    time = result$time %||% NA_real_,
    failure = TRUE
  )
}

#' Safely evaluate a metric function (returns NA on error)
#' @keywords internal
safe_metric <- function(f, ...) {
  tryCatch(f(...), error = function(e) NA_real_)
}

`%||%` <- function(x, y) if (is.null(x)) y else x

# --- Primary Metrics ----------------------------------------------------------

#' Warp Mean Integrated Squared Error
#'
#' Compares forward warps: mean(tf_integrate((h_est - h_true)^2))
compute_warp_mise <- function(warp_est, warp_true, arg) {
  if (is.null(warp_est) || inherits(warp_est, "numeric")) return(NA_real_)
  # Ensure both on same grid
  warp_est <- tfd(warp_est, arg = arg)
  warp_true <- tfd(warp_true, arg = arg)
  sq_diff <- (warp_est - warp_true)^2
  ise <- tf_integrate(sq_diff)
  mean(ise, na.rm = TRUE)
}

#' Alignment Error: how well aligned curves recover true pre-warp curves
#'
#' mean(tf_integrate((aligned_i - true_prewarp_i)^2))
compute_alignment_error <- function(aligned, true_prewarp) {
  sq_diff <- (aligned - true_prewarp)^2
  ise <- tf_integrate(sq_diff)
  mean(ise, na.rm = TRUE)
}

#' Template MISE: quality of template estimation
#'
#' tf_integrate((template_est - template_true)^2)
compute_template_mise <- function(template_est, template_true) {
  sq_diff <- (template_est - template_true)^2
  as.numeric(tf_integrate(sq_diff))
}

#' Template elastic distance (Fisher-Rao amplitude distance)
#'
#' Computes the reparameterization-invariant amplitude distance between
#' estimated and true templates using the Fisher-Rao framework:
#'   d_a(f,g) = min_gamma ||q_f - (q_g o gamma) sqrt(gamma')||_L2
#' where q_f is the SRSF of f. Uses fdasrvf's DP alignment to find
#' the optimal reparameterization, then computes aligned SRSF L2 distance.
#'
#' Unlike L2 MISE, this is invariant to phase shifts in the estimated template.
compute_template_elastic_dist <- function(template_est, template_true, arg) {
  f_est <- as.numeric(tf_evaluations(template_est)[[1]])
  f_true <- as.numeric(tf_evaluations(template_true)[[1]])

  # elastic.distance returns Dy (amplitude) and Dx (phase) distances
  ed <- fdasrvf::elastic.distance(f_est, f_true, time = arg)
  ed$Dy
}

# --- Secondary Metrics --------------------------------------------------------

#' Alignment cross-correlation
compute_alignment_cc <- function(aligned, template_true) {
  template_rep <- rep(template_true, length(aligned))
  cc_vals <- tf_crosscor(aligned, template_rep)
  mean(cc_vals, na.rm = TRUE)
}

# --- Diagnostic Metrics -------------------------------------------------------

#' Registration ratio per curve
#'
#' ratio_i = tf_integrate((h_est_i - id)^2) / tf_integrate((h_true_i - id)^2)
#' < 1 means under-registration, > 1 means over-registration
#'
#' @return list with median and IQR
compute_registration_ratio <- function(warp_est, warp_true, arg) {
  if (is.null(warp_est) || inherits(warp_est, "numeric")) return(NA_real_)
  warp_est <- tfd(warp_est, arg = arg)
  warp_true <- tfd(warp_true, arg = arg)

  identity <- tfd(matrix(arg, nrow = 1), arg = arg)
  identity_n <- rep(identity, length(warp_est))

  est_dev <- tf_integrate((warp_est - identity_n)^2)
  true_dev <- tf_integrate((warp_true - identity_n)^2)

  # Avoid division by near-zero
  valid <- true_dev > 1e-10
  if (sum(valid) < 2) return(list(median = NA_real_, iqr = NA_real_))

  ratios <- est_dev[valid] / true_dev[valid]
  list(
    median = median(ratios, na.rm = TRUE),
    iqr = IQR(ratios, na.rm = TRUE)
  )
}

#' Amplitude variance ratio
#'
#' mean_t(Var_i(aligned - template_est)) / mean_t(Var_i(true_amp_component))
#' < 1 means over-registration absorbed amplitude into phase
#' > 1 means under-registration left phase in amplitude
compute_amp_variance_ratio <- function(
  aligned,
  template_est,
  true_prewarp,
  template_true
) {
  # Residuals from estimated template
  resid_est <- aligned - rep(template_est, length(aligned))
  var_est <- var(resid_est)
  mean_var_est <- mean(tf_evaluations(var_est)[[1]], na.rm = TRUE)

  # True amplitude component
  resid_true <- true_prewarp - rep(template_true, length(true_prewarp))
  var_true <- var(resid_true)
  mean_var_true <- mean(tf_evaluations(var_true)[[1]], na.rm = TRUE)

  if (mean_var_true < 1e-10) return(NA_real_)
  mean_var_est / mean_var_true
}

#' Warp slope ratio
#'
#' Compare range(h'_est(t)) vs range(h'_true(t)) per curve.
#' Detects whether estimated warps are too smooth or too wiggly.
#'
#' @return list with median and IQR of per-curve slope range ratios
compute_warp_slope_ratio <- function(warp_est, warp_true) {
  if (is.null(warp_est) || inherits(warp_est, "numeric")) return(NA_real_)

  # Derivatives
  d_est <- tf_derive(warp_est)
  d_true <- tf_derive(warp_true)

  slope_range <- function(d) {
    sapply(tf_evaluations(d), function(v) {
      v <- v[is.finite(v)]
      if (length(v) < 2) return(NA_real_)
      diff(range(v))
    })
  }

  range_est <- slope_range(d_est)
  range_true <- slope_range(d_true)

  valid <- is.finite(range_est) & is.finite(range_true) & range_true > 1e-10
  if (sum(valid) < 2) return(list(median = NA_real_, iqr = NA_real_))

  ratios <- range_est[valid] / range_true[valid]
  list(
    median = median(ratios, na.rm = TRUE),
    iqr = IQR(ratios, na.rm = TRUE)
  )
}

# --- Flatten Metrics for Data Frame -------------------------------------------

#' Flatten nested metric list to single-level named list
#' @keywords internal
flatten_metrics <- function(metrics) {
  flat <- list()
  for (nm in names(metrics)) {
    val <- metrics[[nm]]
    if (is.list(val) && !is.null(names(val))) {
      for (sub_nm in names(val)) {
        flat[[paste0(nm, "_", sub_nm)]] <- val[[sub_nm]]
      }
    } else {
      flat[[nm]] <- val
    }
  }
  flat
}
