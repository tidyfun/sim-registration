# sim-methods.R -- Method Wrappers for Registration Benchmark v2
#
# 5 methods, all return tf_registration objects.
#
# Provides:
#   fit_method()      - main entry: run one method on one dataset
#   method_configs()  - table of 5 method configurations

if (requireNamespace("tf", quietly = TRUE)) library(tf) else
  devtools::load_all()

# --- Method Configuration Table -----------------------------------------------

method_configs <- function() {
  list(
    srvf = list(method = "srvf", args = list()),
    cc_default = list(method = "fda", args = list()),
    cc_crit1 = list(method = "fda", args = list(crit = 1)),
    affine_ss = list(method = "affine", args = list(type = "shift_scale")),
    landmark_auto = list(method = "landmark", args = list())
  )
}

study_f_preproc_configs <- function() {
  list(
    none = list(
      family = "none",
      primary = FALSE,
      transfer_raw = TRUE,
      transform = function(x) x
    ),
    lowess_f010 = list(
      family = "tf_smooth",
      primary = FALSE,
      transfer_raw = TRUE,
      transform = function(x) {
        tf_smooth(x, method = "lowess", f = 0.10, verbose = FALSE)
      }
    ),
    lowess_f015 = list(
      family = "tf_smooth",
      primary = TRUE,
      transfer_raw = TRUE,
      transform = function(x) {
        tf_smooth(x, method = "lowess", f = 0.15, verbose = FALSE)
      }
    ),
    spline_local_k15 = list(
      family = "tfb_spline",
      primary = FALSE,
      transfer_raw = TRUE,
      transform = function(x) {
        tfb_spline(
          x,
          penalized = TRUE,
          global = FALSE,
          k = 15,
          verbose = FALSE
        )
      }
    ),
    spline_local_k25 = list(
      family = "tfb_spline",
      primary = TRUE,
      transfer_raw = TRUE,
      transform = function(x) {
        tfb_spline(
          x,
          penalized = TRUE,
          global = FALSE,
          k = 25,
          verbose = FALSE
        )
      }
    ),
    spline_global_k25 = list(
      family = "tfb_spline",
      primary = FALSE,
      transfer_raw = TRUE,
      transform = function(x) {
        tfb_spline(
          x,
          penalized = TRUE,
          global = TRUE,
          k = 25,
          verbose = FALSE
        )
      }
    )
  )
}

prepare_estimation_input <- function(data, preproc_id = NULL) {
  if (is.null(preproc_id) || is.na(preproc_id) || identical(preproc_id, "")) {
    preproc_id <- "none"
  }
  configs <- study_f_preproc_configs()
  if (!preproc_id %in% names(configs)) {
    cli::cli_abort(
      "Unknown preproc_id: {preproc_id}. Must be one of: {paste(names(configs), collapse = ', ')}"
    )
  }
  config <- configs[[preproc_id]]
  x_fit <- config$transform(data$x)
  x_eval <- as.tfd(x_fit)
  if (inherits(x_eval, "tfd_irreg")) {
    x_eval <- tfd(x_eval, arg = data$arg)
  }
  list(
    id = preproc_id,
    family = config$family,
    primary = isTRUE(config$primary),
    transfer_raw = isTRUE(config$transfer_raw),
    x_fit = x_fit,
    x_eval = x_eval
  )
}

# --- Main Method Wrapper ------------------------------------------------------

#' Fit one method on one dataset
#'
#' @param data list from generate_data()
#' @param config_name character: name from method_configs()
#' @param use_true_template logical: if TRUE, pass true template to method
#' @param lambda numeric or NULL: regularization parameter (Study B)
#' @return list with:
#'   registration: tf_registration object (or NULL on failure)
#'   time: elapsed time in seconds
#'   error: error message (or NULL on success)
fit_method <- function(
  data,
  config_name,
  use_true_template = FALSE,
  lambda = NULL,
  preproc_id = NULL
) {
  configs <- method_configs()
  if (!config_name %in% names(configs)) {
    cli::cli_abort(
      "Unknown config: {config_name}. Must be one of: {paste(names(configs), collapse = ', ')}"
    )
  }
  config <- configs[[config_name]]

  prep <- prepare_estimation_input(data, preproc_id = preproc_id)
  x <- prep$x_fit
  template <- if (use_true_template) data$template else NULL

  # Build call arguments
  call_args <- list(
    x = x,
    method = config$method,
    template = template
  )
  extra <- config$args

  # Method-specific setup
  if (config_name == "landmark_auto") {
    # Auto-detect landmarks
    lm_detected <- tryCatch(
      tf_landmarks_extrema(tf_smooth(x, verbose = FALSE), "both"),
      error = function(e) NULL
    )
    if (is.null(lm_detected)) {
      return(list(
        registration = NULL,
        time = NA_real_,
        error = "Landmark detection failed"
      ))
    }
    extra$landmarks <- lm_detected
    call_args$template <- NULL
  }

  # Template-based methods: increase Procrustes iterations when estimating
  if (
    !use_true_template &&
      config$method %in% c("fda", "affine") &&
      config_name != "landmark_auto"
  ) {
    call_args$max_iter <- 10L
  }

  # Study B: pass lambda if provided
  if (!is.null(lambda)) {
    extra$lambda <- lambda
  }

  call_args <- c(call_args, extra)

  # Run with timing and error handling
  t0 <- proc.time()
  result <- tryCatch(
    {
      reg <- do.call(tf_register, call_args)
      list(registration = reg, error = NULL)
    },
    error = function(e) {
      list(registration = NULL, error = conditionMessage(e))
    }
  )
  elapsed <- (proc.time() - t0)["elapsed"]

  list(
    registration = result$registration,
    time = as.numeric(elapsed),
    error = result$error,
    x_eval = prep$x_eval,
    preproc_id = prep$id,
    preproc_family = prep$family,
    preproc_primary = prep$primary,
    transfer_raw = prep$transfer_raw
  )
}
