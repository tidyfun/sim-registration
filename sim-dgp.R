# sim-dgp.R -- Data-Generating Processes for Registration Benchmark v2
#
# Provides:
#   make_template()       - 2 template types (harmonic, wiggly)
#   generate_warps()      - 3 warp types (simple, complex, affine)
#   generate_amplitude()  - 4 amplitude types (none, rank1_mult, rank2, highrank)
#   generate_data()       - main entry point
#   contaminate_data()    - outlier injection for Study D

if (requireNamespace("tf", quietly = TRUE)) library(tf) else
  devtools::load_all()

# --- Templates ----------------------------------------------------------------

#' Create z-scored template function on [0,1]
#'
#' @param type "harmonic" or "wiggly"
#' @param arg evaluation grid (numeric vector in [0,1])
#' @return tfd of length 1
make_template <- function(type = c("harmonic", "wiggly"), arg) {
  type <- match.arg(type)
  vals <- switch(
    type,
    harmonic = sin(2 * pi * arg) + 0.5 * sin(4 * pi * arg),
    wiggly = 0.6 *
      sin(2 * pi * arg) +
      0.5 * sin(5 * pi * arg + 0.3) +
      0.45 * cos(8 * pi * arg - 0.7) +
      0.4 * sin(11 * pi * arg + 1.2) +
      0.3 * cos(14 * pi * arg)
  )
  # Z-score for comparability
  vals <- (vals - mean(vals)) / sd(vals)
  tfd(matrix(vals, nrow = 1), arg = arg)
}

# --- Warp Generators ----------------------------------------------------------

#' Generate warping functions
#'
#' @param n number of curves
#' @param arg evaluation grid
#' @param type "simple", "complex", or "affine"
#' @param severity controls warp magnitude (default 0.5)
#' @return tfd of warping functions (forward warps h(s) = t)
generate_warps <- function(
  n,
  arg,
  type = c("simple", "complex", "affine"),
  severity = 0.5
) {
  type <- match.arg(type)
  switch(
    type,
    simple = generate_smooth_warps(n, arg, severity, gp_type = "simple"),
    complex = generate_smooth_warps(n, arg, severity, gp_type = "complex"),
    affine = generate_affine_warps(n, arg, severity)
  )
}

#' Smooth elastic warps via GP -> positivity -> integrate -> normalize
#'
#' @param n number of curves
#' @param arg evaluation grid
#' @param severity GP amplitude multiplier
#' @param gp_type "simple" (squareexp, long length-scale) or
#'   "complex" (matern, short length-scale, wiggly)
#' @return tfd of domain-preserving warps (h(0)=0, h(1)=1)
generate_smooth_warps <- function(
  n,
  arg,
  severity = 0.5,
  gp_type = c("simple", "complex")
) {
  gp_type <- match.arg(gp_type)
  gp <- switch(
    gp_type,
    simple = tf_rgp(
      n,
      arg = arg,
      cov = "squareexp",
      scale = 0.3,
      nugget = 0.001
    ),
    complex = tf_rgp(
      n,
      arg = arg,
      cov = "matern",
      order = 1.5,
      scale = 0.05,
      nugget = 0.001
    )
  )
  # Positivity: exp(severity * centered_gp)
  positive_funcs <- exp(severity * (gp - mean(gp)))
  # Integrate and normalize to [0,1] -> [0,1]
  warps <- tf_integrate(positive_funcs, definite = FALSE)
  # Use endpoint evaluation h(1) for normalization (monotone => max at endpoint)
  endpoint_vals <- as.numeric(warps[, max(arg)])
  warps / endpoint_vals
}

#' Affine warps: h_i(s) = a_i * s + b_i
#'
#' Non-domain-preserving. Post-hoc re-centered.
#'
#' @param n number of curves
#' @param arg evaluation grid
#' @param severity controls warp magnitude
#' @return tfd of forward affine warps
generate_affine_warps <- function(n, arg, severity = 0.5) {
  sd_a <- 0.10 * severity
  sd_b <- 0.06 * severity
  lower_a <- 1 - 0.30 * severity
  upper_a <- 1 + 0.30 * severity
  lower_b <- -0.15 * severity
  upper_b <- 0.15 * severity

  a <- rtruncnorm(n, mean = 1, sd = sd_a, lower = lower_a, upper = upper_a)
  b <- rtruncnorm(n, mean = 0, sd = sd_b, lower = lower_b, upper = upper_b)

  # Re-center for identifiability
  a <- a - mean(a) + 1
  b <- b - mean(b)

  warp_mat <- t(sapply(seq_len(n), function(i) a[i] * arg + b[i]))
  tfd(warp_mat, arg = arg)
}

#' Evaluate inverse affine warps on an observed-time grid
#'
#' Converts forward affine warps h_i(s) = a_i * s + b_i to inverse warp values
#' h_i^{-1}(t) on the supplied grid t.
#'
#' @param warps tfd of forward affine warps on aligned-time grid
#' @param arg observed-time grid for evaluating h^{-1}(t)
#' @return numeric matrix with one row per curve and one column per arg value
evaluate_inverse_affine_warps <- function(warps, arg) {
  warp_mat <- as.matrix(tfd(warps, arg = tf_arg(warps)))
  slope <- warp_mat[, ncol(warp_mat)] - warp_mat[, 1]
  intercept <- warp_mat[, 1]

  t(vapply(
    seq_len(nrow(warp_mat)),
    function(i) (arg - intercept[i]) / slope[i],
    numeric(length(arg))
  ))
}

#' Truncated normal distribution (rejection sampling)
#' @keywords internal
rtruncnorm <- function(
  n,
  mean = 0,
  sd = 1,
  lower = -Inf,
  upper = Inf,
  max_iter = 10000L
) {
  out <- numeric(n)
  for (i in seq_len(n)) {
    iter <- 0L
    repeat {
      iter <- iter + 1L
      if (iter > max_iter) {
        cli::cli_abort("rtruncnorm: exceeded {max_iter} iterations")
      }
      x <- rnorm(1, mean, sd)
      if (x >= lower && x <= upper) {
        out[i] <- x
        break
      }
    }
  }
  out
}

# --- Amplitude Variation ------------------------------------------------------

#' Generate amplitude variation components
#'
#' @param type "none", "rank1_mult", "rank2", or "highrank"
#' @param n number of curves
#' @param arg evaluation grid (needed for highrank interpolation)
#' @param amp_sd amplitude effect size (default 0.3)
#' @return list with $type and type-specific components
generate_amplitude <- function(
  type = c("none", "rank1_mult", "rank2", "highrank"),
  n,
  arg = NULL,
  amp_sd = 0.3
) {
  type <- match.arg(type)

  switch(
    type,
    none = list(type = "none"),

    rank1_mult = {
      # LogNormal(-s^2/2, s^2) has E[a] = 1
      a <- rlnorm(n, meanlog = -amp_sd^2 / 2, sdlog = amp_sd)
      list(type = "rank1_mult", a = a)
    },

    rank2 = {
      # Split total target variance amp_sd^2 equally between scale and shift
      # so that total Var(a_i * m + c_i) ~ amp_sd^2 (comparable to other types)
      scale_sd <- amp_sd / sqrt(2)
      shift_sd <- amp_sd / sqrt(2)
      a <- rlnorm(n, meanlog = -scale_sd^2 / 2, sdlog = scale_sd)
      c_shift <- rnorm(n, 0, shift_sd)
      list(type = "rank2", a = a, c = c_shift)
    },

    highrank = {
      checkmate::assert_numeric(arg, min.len = 10)
      fpc_data <- load_gait_fpcs()

      # Interpolate 3 eigenfunctions from gait's 20-point grid to benchmark grid
      phi_interp <- apply(fpc_data$phi, 2, function(col) {
        approx(fpc_data$arg, col, xout = arg, rule = 2)$y
      })

      # Re-orthonormalize after interpolation (QR decomposition)
      dt <- arg[2] - arg[1]
      qr_phi <- qr(phi_interp * sqrt(dt))
      phi_interp <- qr.Q(qr_phi) / sqrt(dt)

      # Flatten eigenvalue decay: target ratio 3:2:1
      lambda_raw <- fpc_data$lambda
      lambda_flat <- c(3, 2, 1) * mean(lambda_raw) / 2

      # Scale total amplitude variance to be comparable with rank1/rank2
      # Target: total amp SD ~ amp_sd (same as multiplicative LogNormal sd)
      total_var <- sum(lambda_flat)
      target_var <- amp_sd^2 # comparable to multiplicative a_i variance
      scale_factor <- sqrt(target_var / total_var)
      lambda_scaled <- lambda_flat * scale_factor^2

      # Generate scores c_ij ~ N(0, lambda_j)
      scores <- matrix(nrow = n, ncol = 3)
      for (j in 1:3) {
        scores[, j] <- rnorm(n, 0, sqrt(lambda_scaled[j]))
      }

      # amplitude_i(t) = sum_j c_ij * phi_j(t)
      amp_mat <- scores %*% t(phi_interp)

      list(
        type = "highrank",
        amp_mat = amp_mat,
        scores = scores,
        phi = phi_interp,
        lambda = lambda_scaled
      )
    }
  )
}

#' Load cached gait knee FPCs
#' @keywords internal
load_gait_fpcs <- function() {
  fpc_file <- file.path(
    here::here(),
    "gait_knee_fpcs.rds"
  )
  if (!file.exists(fpc_file)) {
    cli::cli_abort(c(
      "Gait FPC cache not found at {.file {fpc_file}}.",
      "i" = "Run the FPC extraction script first."
    ))
  }
  readRDS(fpc_file)
}

# --- Apply Amplitude Variation ------------------------------------------------

#' Apply amplitude variation to warped curves
#'
#' @param warped tfd of warped curves (m(h_i(t)))
#' @param amp list from generate_amplitude()
#' @param template tfd of length 1
#' @param arg evaluation grid
#' @return tfd of observed curves (before noise)
apply_amplitude <- function(warped, amp, template, arg) {
  switch(
    amp$type,
    none = warped,
    rank1_mult = warped * amp$a,
    rank2 = warped * amp$a + amp$c,
    highrank = warped + tfd(amp$amp_mat, arg = arg)
  )
}

#' Compute the true amplitude component for each curve
#'
#' Returns what the aligned curve SHOULD look like: template + amplitude.
#' Used for alignment error metric.
#'
#' @param amp list from generate_amplitude()
#' @param template tfd of length 1
#' @param n number of curves
#' @param arg evaluation grid
#' @return tfd of true pre-warp curves (template + amplitude component)
true_prewarp_curves <- function(amp, template, n, arg) {
  m <- rep(template, n)
  switch(
    amp$type,
    none = m,
    rank1_mult = m * amp$a,
    rank2 = m * amp$a + amp$c,
    highrank = m + tfd(amp$amp_mat, arg = arg)
  )
}

# --- DGP Specification Table --------------------------------------------------

#' Get DGP specification
#' @param dgp character: D01-D15
#' @return list with template, phase, amplitude
dgp_spec <- function(dgp) {
  specs <- list(
    # Baseline (phase only)
    D01 = list(template = "harmonic", phase = "simple", amplitude = "none"),
    D02 = list(template = "harmonic", phase = "complex", amplitude = "none"),
    D03 = list(template = "wiggly", phase = "complex", amplitude = "none"),
    D04 = list(template = "harmonic", phase = "affine", amplitude = "none"),
    # Low-rank amplitude
    D05 = list(
      template = "harmonic",
      phase = "simple",
      amplitude = "rank1_mult"
    ),
    D06 = list(
      template = "wiggly",
      phase = "simple",
      amplitude = "rank1_mult"
    ),
    D07 = list(template = "harmonic", phase = "complex", amplitude = "rank2"),
    D08 = list(template = "wiggly", phase = "complex", amplitude = "rank2"),
    # High-rank amplitude
    D09 = list(
      template = "harmonic",
      phase = "simple",
      amplitude = "highrank"
    ),
    D10 = list(
      template = "harmonic",
      phase = "complex",
      amplitude = "highrank"
    ),
    D11 = list(
      template = "wiggly",
      phase = "simple",
      amplitude = "highrank"
    ),
    D12 = list(template = "wiggly", phase = "complex", amplitude = "highrank"),
    D13 = list(template = "wiggly", phase = "affine", amplitude = "highrank"),
    # Added to complete key ladders (council review)
    D14 = list(template = "wiggly", phase = "simple", amplitude = "none"),
    D15 = list(template = "harmonic", phase = "affine", amplitude = "highrank")
  )
  if (!dgp %in% names(specs)) {
    cli::cli_abort(
      "Unknown DGP: {dgp}. Must be one of: {paste(names(specs), collapse = ', ')}"
    )
  }
  specs[[dgp]]
}

# --- Main Data Generator ------------------------------------------------------

#' Generate a complete dataset for one DGP
#'
#' @param dgp character DGP identifier (D01-D13)
#' @param n number of curves
#' @param n_grid number of grid points
#' @param severity warp severity
#' @param noise_sd noise standard deviation
#' @param seed random seed
#' @param amp_sd amplitude effect size (default 0.3)
#' @return list with components:
#'   x (tfd): observed curves
#'   template (tfd): true template
#'   warps (tfd): true warping functions (forward: h(s) = t)
#'   true_prewarp (tfd): template + amplitude component (before warping)
#'   arg (numeric): evaluation grid
#'   amp_coefs (list): amplitude coefficients used
#'   seed (integer): seed used
generate_data <- function(
  dgp,
  n = 50,
  n_grid = 101,
  severity = 0.5,
  noise_sd = 0,
  seed = NULL,
  amp_sd = 0.3
) {
  if (!is.null(seed)) set.seed(seed)
  arg <- seq(0, 1, length.out = n_grid)

  spec <- dgp_spec(dgp)

  # 1. Template
  template <- make_template(spec$template, arg)

  # 2. Warps
  warps <- generate_warps(n, arg, type = spec$phase, severity = severity)

  # 3. Apply warps to template
  is_non_dp <- spec$phase == "affine"
  if (is_non_dp) {
    # Non-domain-preserving: evaluate m at h_i^{-1}(t) so data$warps
    # remains a forward warp, consistent with the other DGPs.
    template_vals <- tf_evaluations(template)[[1]]
    inv_warp_mat <- evaluate_inverse_affine_warps(warps, arg)
    warped_mat <- t(vapply(
      seq_len(nrow(inv_warp_mat)),
      function(i) {
        approx(arg, template_vals, xout = inv_warp_mat[i, ], rule = 2)$y
      },
      numeric(length(arg))
    ))
    warped <- tfd(warped_mat, arg = arg)
  } else {
    warped <- tf_warp(rep(template, n), warps)
  }

  # 4. Amplitude variation
  amp <- generate_amplitude(
    type = spec$amplitude,
    n = n,
    arg = arg,
    amp_sd = amp_sd
  )
  prewarp <- true_prewarp_curves(amp, template, n, arg)
  x <- apply_amplitude(warped, amp, template, arg)

  # 5. Noise
  if (noise_sd > 0) {
    noise_mat <- matrix(rnorm(n * n_grid, 0, noise_sd), nrow = n)
    x <- x + tfd(noise_mat, arg = arg)
  }

  # Ensure regular tfd
  if (inherits(x, "tfd_irreg")) {
    x <- tfd(x, arg = arg)
  }

  list(
    x = x,
    template = template,
    warps = warps,
    true_prewarp = prewarp,
    arg = arg,
    amp_coefs = amp,
    seed = seed,
    dgp = dgp,
    spec = spec
  )
}

# --- Outlier Contamination ----------------------------------------------------

#' Contaminate a generated dataset with outliers
#'
#' @param data list from generate_data()
#' @param contam_frac fraction of curves to contaminate
#' @param outlier_type "shape" or "phase"
#' @param seed seed offset for deterministic outlier selection
#' @return modified data list with $outlier_mask
contaminate_data <- function(
  data,
  contam_frac = 0.10,
  outlier_type = c("shape", "phase"),
  seed = NULL
) {
  outlier_type <- match.arg(outlier_type)
  checkmate::assert_number(contam_frac, lower = 0, upper = 1)
  n <- length(data$x)
  n_outliers <- floor(n * contam_frac)
  if (n_outliers == 0) {
    data$outlier_mask <- rep(FALSE, n)
    return(data)
  }

  if (!is.null(seed)) {
    old_seed <- if (exists(".Random.seed", envir = .GlobalEnv)) {
      .Random.seed
    } else {
      NULL
    }
    on.exit(
      {
        if (!is.null(old_seed)) {
          .Random.seed <<- old_seed
        }
      },
      add = TRUE
    )
    set.seed(seed + 10000L)
  }
  outlier_idx <- sample.int(n, n_outliers)
  outlier_mask <- seq_len(n) %in% outlier_idx

  arg <- data$arg
  x_mat <- as.matrix(data$x)

  switch(
    outlier_type,

    shape = {
      # Add random sinusoids (scaled to ~2x template range)
      template_range <- diff(range(tf_evaluations(data$template)[[1]]))
      for (j in seq_along(outlier_idx)) {
        freq <- runif(1, 3, 8)
        phase <- runif(1, 0, 2 * pi)
        amp <- template_range * runif(1, 0.5, 1.5)
        x_mat[outlier_idx[j], ] <- x_mat[outlier_idx[j], ] +
          amp * sin(freq * pi * arg + phase)
      }
    },

    phase = {
      # Replace with extreme warps (severity = 3.0)
      extreme_warps <- generate_smooth_warps(
        n_outliers,
        arg,
        severity = 3.0,
        gp_type = "simple"
      )
      extreme_warp_evals <- tf_evaluations(extreme_warps)
      template_vals <- tf_evaluations(data$template)[[1]]
      for (j in seq_along(outlier_idx)) {
        w <- extreme_warp_evals[[j]]
        x_mat[outlier_idx[j], ] <- approx(
          arg,
          template_vals,
          xout = w,
          rule = 2
        )$y
      }
    }
  )

  data$x <- tfd(x_mat, arg = arg)
  data$outlier_mask <- outlier_mask
  data
}
