# sim-config.R -- All Study Designs for Registration Benchmark v2
#
# Studies A-E in one file.
#
# Provides:
#   study_a_design()  - Main study (15 DGPs x severity x noise x 5 methods)
#   study_b_design()  - Penalization (lambda grid, 3 methods)
#   study_c_design()  - Grid resolution sensitivity
#   study_d_design()  - Oracle template comparison
#   study_e_design()  - Outlier contamination
#   full_design()     - Combined design

# --- Study A: Main Study -----------------------------------------------------

#' Study A design: 15 DGPs, always estimate template
#'
#' 15 DGPs x 2 severity x 3 noise x 5 methods = 450 cells x 100 reps
study_a_design <- function() {
  grid <- expand.grid(
    dgp = paste0("D", sprintf("%02d", 1:15)),
    n_curves = 50,
    n_grid = 101,
    noise_sd = c(0, 0.1, 0.3),
    severity = c(0.5, 1.0),
    stringsAsFactors = FALSE
  )
  grid <- expand_methods(grid)
  grid$reps <- 100
  grid$study <- "A"
  grid$use_true_template <- FALSE
  grid$contam_frac <- NA_real_
  grid$outlier_type <- NA_character_
  grid$lambda <- NA_real_
  grid
}

# --- Study B: Penalization Arm ------------------------------------------------

#' Study B design: penalization study
#'
#' Lambda grid calibrated by pilot (job 5131464): methods need very different
#' lambda ranges but a shared 8-value grid covers all optima.
#'
#' 4 DGPs x 8 lambdas x 3 methods x 3 noise x 2 severity = 576 cells x 33 reps
#'
#' @param dgps character vector of DGP names
#' @param lambdas numeric vector of lambda values (shared across methods)
study_b_design <- function(
  dgps = c("D01", "D02", "D09", "D12"),
  lambdas = c(0, 1e-4, 1e-3, 0.01, 0.05, 0.1, 1, 10)
) {
  methods <- c("srvf", "fda_default", "fda_crit1")

  grid <- expand.grid(
    dgp = dgps,
    n_curves = 50,
    n_grid = 101,
    noise_sd = c(0, 0.1, 0.3),
    severity = c(0.5, 1.0),
    lambda = lambdas,
    method = methods,
    stringsAsFactors = FALSE
  )

  grid$reps <- 33
  grid$study <- "B"
  grid$use_true_template <- FALSE
  grid$contam_frac <- NA_real_
  grid$outlier_type <- NA_character_
  grid
}

# --- Study C: Grid Resolution Sensitivity ------------------------------------

#' Study C design: grid resolution sensitivity
#'
#' Tests whether grid resolution is a confound for method comparison.
#' SRVF operates on SRSFs (numerical derivatives) which amplify discretization
#' artifacts at finer grids. Pilot (job 5131384) showed 1.82x worse warp MISE
#' at n_grid=201 vs 101 for SRVF.
#'
#' 5 DGPs x 3 grid sizes x 3 noise x 5 methods x severity=1.0
#' = 225 cells x 50 reps = 11,250 runs
#'
#' @param dgps character vector of DGP names
#' @param grids integer vector of grid sizes
study_c_design <- function(
  dgps = c("D01", "D03", "D04", "D09", "D12"),
  grids = c(51L, 101L, 201L)
) {
  grid <- expand.grid(
    dgp = dgps,
    n_curves = 50,
    n_grid = grids,
    noise_sd = c(0, 0.1, 0.3),
    severity = 1.0,
    stringsAsFactors = FALSE
  )
  grid <- expand_methods(grid)
  grid$reps <- 50
  grid$study <- "C"
  grid$use_true_template <- FALSE
  grid$contam_frac <- NA_real_
  grid$outlier_type <- NA_character_
  grid$lambda <- NA_real_
  grid
}

# --- Study D: Oracle Template Comparison --------------------------------------

#' Study D design: oracle template comparison
#'
#' Decomposes registration error into template estimation error vs warp
#' estimation error. When given the true template, template MISE = 0 and any
#' remaining warp MISE is purely from the warp estimation algorithm.
#'
#' 6 DGPs x 2 template modes x 2 conditions x 5 methods
#' = 120 cells x 50 reps = 6,000 runs
#'
#' Two conditions:
#'   Easy: noise=0.1, severity=0.5 (anchor)
#'   Hard: noise=0.3, severity=1.0 (stress test)
#'
#' @param dgps character vector of DGP names
study_d_design <- function(
  dgps = c("D01", "D03", "D09", "D10", "D13", "D14")
) {
  # Two conditions: easy anchor + hard stress test
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
  grid$reps <- 50
  grid$study <- "D"
  grid$contam_frac <- NA_real_
  grid$outlier_type <- NA_character_
  grid$lambda <- NA_real_
  grid
}

# --- Study E: Outlier Contamination -------------------------------------------

#' Study E design: outlier contamination
#'
#' Tests method robustness when a fraction of curves are outliers.
#' Contamination fractions start at 10% (at n=50, 5% = 2-3 curves has
#' discrete rounding artifacts).
#'
#' 3 DGPs x 2 outlier types x 3 fractions x 5 methods x 2 noise x severity=0.5
#' = 180 cells x 50 reps = 9,000 runs
#'
#' @param dgps character vector of 3 representative DGP names
study_e_design <- function(dgps = c("D02", "D09", "D12")) {
  grid <- expand.grid(
    dgp = dgps,
    n_curves = 50,
    n_grid = 101,
    noise_sd = c(0.1, 0.3),
    severity = 0.5,
    contam_frac = c(0.10, 0.20, 0.30),
    outlier_type = c("shape", "phase"),
    stringsAsFactors = FALSE
  )
  grid <- expand_methods(grid)
  grid$reps <- 50
  grid$study <- "E"
  grid$use_true_template <- FALSE
  grid$lambda <- NA_real_
  grid
}

# --- Helpers ------------------------------------------------------------------

#' Expand grid to include all 5 method configs
#' @keywords internal
expand_methods <- function(grid) {
  methods <- all_method_names()
  rows <- lapply(seq_len(nrow(grid)), function(i) {
    cbind(
      grid[rep(i, length(methods)), , drop = FALSE],
      method = methods,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

#' All method config names (must match sim-methods.R::method_configs())
#' @keywords internal
all_method_names <- function() {
  c("srvf", "fda_default", "fda_crit1", "affine_ss", "landmark_auto")
}

# --- Seed Generation ----------------------------------------------------------

#' Generate deterministic seed for a (dgp, rep) combination
make_seed <- function(dgp, rep) {
  digest::digest2int(paste(dgp, rep, sep = "_"))
}

# --- Full Design --------------------------------------------------------------

#' Combine studies into one design table
#'
#' @param studies character vector: subset of c("A", "B", "C", "D", "E")
#' @return data.frame with all design rows
full_design <- function(studies = "A") {
  designs <- list()
  if ("A" %in% studies) designs$A <- study_a_design()
  if ("B" %in% studies) designs$B <- study_b_design()
  if ("C" %in% studies) designs$C <- study_c_design()
  if ("D" %in% studies) designs$D <- study_d_design()
  if ("E" %in% studies) designs$E <- study_e_design()

  design <- do.call(rbind, designs)
  rownames(design) <- NULL
  design$cell_id <- seq_len(nrow(design))

  message(sprintf(
    "Design: %d cells, %d total runs (studies: %s)",
    nrow(design),
    sum(design$reps),
    paste(studies, collapse = ", ")
  ))
  design
}

# --- Design Summary -----------------------------------------------------------

summarize_design <- function(design = NULL) {
  if (is.null(design)) design <- full_design("A")

  cat("=== Registration Benchmark v2 Design ===\n\n")
  for (s in unique(design$study)) {
    d <- design[design$study == s, ]
    cat(sprintf(
      "Study %s: %d cells, %d runs\n",
      s,
      nrow(d),
      sum(d$reps)
    ))
    cat(sprintf("  DGPs: %s\n", paste(sort(unique(d$dgp)), collapse = ", ")))
    cat(sprintf(
      "  Methods: %s\n",
      paste(sort(unique(d$method)), collapse = ", ")
    ))
    if ("n_grid" %in% names(d) && length(unique(d$n_grid)) > 1) {
      cat(sprintf(
        "  Grids: %s\n",
        paste(sort(unique(d$n_grid)), collapse = ", ")
      ))
    }
    if (any(!is.na(d$use_true_template)) && any(d$use_true_template)) {
      cat("  Template modes: oracle + estimated\n")
    }
  }
  cat(sprintf("\nTotal cells: %d\n", nrow(design)))
  cat(sprintf("Total runs: %d\n", sum(design$reps)))
}
