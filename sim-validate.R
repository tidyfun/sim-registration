# sim-validate.R -- Known-Answer Tests for Registration Benchmark v2
#
# Provides:
#   run_known_answer_tests()  - Sanity checks with known correct answers
#
# Run this before the full benchmark to catch bugs early.

base_dir <- here::here()
source(file.path(base_dir, "sim-dgp.R"))
source(file.path(base_dir, "sim-methods.R"))
source(file.path(base_dir, "sim-metrics.R"))
source(file.path(base_dir, "sim-config.R"))

# ==============================================================================
# Known-Answer Tests
# ==============================================================================

run_known_answer_tests <- function(verbose = TRUE) {
  n_pass <- 0
  n_fail <- 0

  report <- function(name, passed, detail = "") {
    status <- if (passed) "PASS" else "FAIL"
    if (verbose) cat(sprintf("  [%s] %s %s\n", status, name, detail))
    if (passed) n_pass <<- n_pass + 1 else n_fail <<- n_fail + 1
  }

  cat("=== Known-Answer Tests (v2) ===\n\n")

  # --- Test 1: Identity warp -> near-zero MISE --------------------------------
  cat("Test 1: Identity warp recovery\n")
  tryCatch(
    {
      arg <- seq(0, 1, length.out = 101)
      template <- make_template("harmonic", arg)
      x <- rep(template, 10) +
        tfd(matrix(rnorm(10 * 101, 0, 0.001), nrow = 10), arg = arg)
      identity_warps <- tfd(matrix(rep(arg, each = 10), nrow = 10), arg = arg)
      prewarp <- rep(template, 10)

      data <- list(
        x = x,
        template = template,
        warps = identity_warps,
        true_prewarp = prewarp,
        arg = arg,
        amp_coefs = list(type = "none"),
        spec = list(phase = "simple")
      )
      result <- fit_method(data, "srvf")
      if (!is.null(result$error)) {
        report("identity_srvf", FALSE, paste("Error:", result$error))
      } else {
        metrics <- extract_metrics(data, result)
        report(
          "identity_srvf",
          metrics$warp_mise < 0.001,
          sprintf("MISE = %.2e", metrics$warp_mise)
        )
      }
    },
    error = function(e) report("identity_warp", FALSE, e$message)
  )

  # --- Test 2: Affine method on affine DGP -> good recovery -------------------
  cat("\nTest 2: Affine method on affine DGP (D04)\n")
  tryCatch(
    {
      data <- generate_data(
        "D04",
        n = 30,
        severity = 0.5,
        noise_sd = 0,
        seed = 42
      )
      result <- fit_method(data, "affine_ss")
      if (!is.null(result$error)) {
        report("affine_on_D04", FALSE, paste("Error:", result$error))
      } else {
        metrics <- extract_metrics(data, result)
        report(
          "affine_on_D04",
          metrics$warp_mise < 0.05,
          sprintf("MISE = %.2e", metrics$warp_mise)
        )
      }
    },
    error = function(e) report("affine_on_D04", FALSE, e$message)
  )

  # --- Test 3: SRVF on easy DGP (D01) ----------------------------------------
  cat("\nTest 3: SRVF on D01 (easy baseline)\n")
  tryCatch(
    {
      data <- generate_data(
        "D01",
        n = 30,
        severity = 0.5,
        noise_sd = 0,
        seed = 123
      )
      result <- fit_method(data, "srvf")
      if (!is.null(result$error)) {
        report("srvf_D01", FALSE, paste("Error:", result$error))
      } else {
        metrics <- extract_metrics(data, result)
        report(
          "srvf_D01_mise",
          metrics$warp_mise < 0.001,
          sprintf("MISE = %.2e", metrics$warp_mise)
        )
        report(
          "srvf_D01_cc",
          metrics$alignment_cc > 0.99,
          sprintf("CC = %.4f", metrics$alignment_cc)
        )
      }
    },
    error = function(e) report("srvf_D01", FALSE, e$message)
  )

  # --- Test 4: Metric sanity checks -------------------------------------------
  cat("\nTest 4: Metric sanity\n")
  tryCatch(
    {
      arg <- seq(0, 1, length.out = 101)
      perfect_warps <- tfd(matrix(rep(arg, each = 5), nrow = 5), arg = arg)

      # MISE = 0 when identical
      mise <- compute_warp_mise(perfect_warps, perfect_warps, arg)
      report("mise_zero", abs(mise) < 1e-10, sprintf("%.2e", mise))

      # Alignment CC = 1 when perfectly aligned
      template <- make_template("harmonic", arg)
      aligned_perfect <- rep(template, 5)
      cc <- compute_alignment_cc(aligned_perfect, template)
      report("cc_one", abs(cc - 1) < 0.01, sprintf("CC = %.4f", cc))
    },
    error = function(e) report("metric_sanity", FALSE, e$message)
  )

  # --- Test 5: DGP generation correctness ------------------------------------
  cat("\nTest 5: DGP generation\n")
  tryCatch(
    {
      # All 13 DGPs should generate without error
      all_ok <- TRUE
      for (d in paste0("D", sprintf("%02d", 1:13))) {
        r <- tryCatch(
          generate_data(d, n = 10, n_grid = 51, seed = 1),
          error = function(e) e
        )
        if (inherits(r, "error")) {
          report(paste("generate", d), FALSE, r$message)
          all_ok <- FALSE
        }
      }
      if (all_ok) report("all_13_dgps_generate", TRUE)

      # Check output structure
      data <- generate_data("D09", n = 10, n_grid = 51, seed = 1)
      report("has_x", inherits(data$x, "tfd"))
      report("has_template", inherits(data$template, "tfd"))
      report("has_warps", inherits(data$warps, "tfd"))
      report("has_prewarp", inherits(data$true_prewarp, "tfd"))
      report("correct_n", length(data$x) == 10)
    },
    error = function(e) report("dgp_generation", FALSE, e$message)
  )

  # --- Test 6: Seed reproducibility -------------------------------------------
  cat("\nTest 6: Seed reproducibility\n")
  tryCatch(
    {
      d1 <- generate_data("D05", n = 10, seed = 42)
      d2 <- generate_data("D05", n = 10, seed = 42)
      x1 <- as.matrix(d1$x)
      x2 <- as.matrix(d2$x)
      report("seed_repro", max(abs(x1 - x2)) < 1e-12)
    },
    error = function(e) report("seed_repro", FALSE, e$message)
  )

  # --- Summary ----------------------------------------------------------------
  cat(sprintf("\n=== Results: %d passed, %d failed ===\n", n_pass, n_fail))
  invisible(list(passed = n_pass, failed = n_fail))
}

# ==============================================================================
# Run if sourced directly
# ==============================================================================

if (sys.nframe() == 0) {
  run_known_answer_tests()
}
