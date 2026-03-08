# sim-config.R -- All Study Designs for Registration Benchmark v2
#
# Studies A-D in one file.
#
# Provides:
#   study_a_design()  - Main study (13 DGPs x severity x noise x 5 methods)
#   study_b_design()  - Penalization arm (data-driven DGP selection)
#   study_c_design()  - Oracle paired sub-study
#   study_d_design()  - Outlier contamination
#   full_design()     - Combined design

# --- Study A: Main Study -----------------------------------------------------

#' Study A design: 13 DGPs, always estimate template
#'
#' 13 DGPs x 2 severity x 2 noise x 5 methods = 260 cells x 100 reps
study_a_design <- function() {
  grid <- expand.grid(
    dgp = paste0("D", sprintf("%02d", 1:13)),
    n_curves = 50,
    n_grid = 101,
    noise_sd = c(0, 0.1),
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
#' DGPs selected post-pilot based on over/under-registration signals.
#' 3 DGPs x 6 lambdas x 2 methods x 2 noise x 2 severity = 144 cells x 50 reps
#'
#' @param dgps character vector of 3 DGP names (data-driven from Study A)
#' @param lambda_srvf numeric vector of lambda values for SRVF
#' @param lambda_fda numeric vector of lambda values for FDA
study_b_design <- function(
  dgps = c("D02", "D10", "D12"),
  lambda_srvf = c(0, 0.001, 0.01, 0.1, 1, 10),
  lambda_fda = c(0, 0.001, 0.01, 0.1, 1, 10)
) {
  # SRVF arm
  srvf_grid <- expand.grid(
    dgp = dgps,
    n_curves = 50,
    n_grid = 101,
    noise_sd = c(0, 0.1),
    severity = c(0.5, 1.0),
    lambda = lambda_srvf,
    method = "srvf",
    stringsAsFactors = FALSE
  )

  # FDA arm
  fda_grid <- expand.grid(
    dgp = dgps,
    n_curves = 50,
    n_grid = 101,
    noise_sd = c(0, 0.1),
    severity = c(0.5, 1.0),
    lambda = lambda_fda,
    method = "fda_default",
    stringsAsFactors = FALSE
  )

  grid <- rbind(srvf_grid, fda_grid)
  grid$reps <- 50
  grid$study <- "B"
  grid$use_true_template <- FALSE
  grid$contam_frac <- NA_real_
  grid$outlier_type <- NA_character_
  grid
}

# --- Study C: Oracle Paired Sub-Study -----------------------------------------

#' Study C design: oracle template comparison
#'
#' 4 DGPs x 2 template modes x 5 methods x noise=0.1 x severity=0.5
#' = 40 cells x 50 reps
#'
#' @param dgps character vector of 4 DGP names (data-driven from Study A)
study_c_design <- function(dgps = c("D03", "D08", "D10", "D12")) {
  grid_oracle <- expand.grid(
    dgp = dgps,
    n_curves = 50,
    n_grid = 101,
    noise_sd = 0.1,
    severity = 0.5,
    use_true_template = TRUE,
    stringsAsFactors = FALSE
  )
  grid_estimated <- grid_oracle
  grid_estimated$use_true_template <- FALSE

  grid <- rbind(grid_oracle, grid_estimated)
  grid <- expand_methods(grid)
  grid$reps <- 50
  grid$study <- "C"
  grid$contam_frac <- NA_real_
  grid$outlier_type <- NA_character_
  grid$lambda <- NA_real_
  grid
}

# --- Study D: Outlier Contamination -------------------------------------------

#' Study D design: outlier contamination
#'
#' 3 DGPs x 2 outlier types x 3 fractions x 5 methods x noise=0.1 x severity=0.5
#' = 90 cells x 50 reps
#'
#' @param dgps character vector of 3 representative DGP names
study_d_design <- function(dgps = c("D01", "D03", "D09")) {
  grid <- expand.grid(
    dgp = dgps,
    n_curves = 50,
    n_grid = 101,
    noise_sd = 0.1,
    severity = 0.5,
    contam_frac = c(0.05, 0.10, 0.20),
    outlier_type = c("shape", "phase"),
    stringsAsFactors = FALSE
  )
  grid <- expand_methods(grid)
  grid$reps <- 50
  grid$study <- "D"
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
#' @param studies character vector: subset of c("A", "B", "C", "D")
#' @return data.frame with all design rows
full_design <- function(studies = "A") {
  designs <- list()
  if ("A" %in% studies) designs$A <- study_a_design()
  if ("B" %in% studies) designs$B <- study_b_design()
  if ("C" %in% studies) designs$C <- study_c_design()
  if ("D" %in% studies) designs$D <- study_d_design()

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
  }
  cat(sprintf("\nTotal cells: %d\n", nrow(design)))
  cat(sprintf("Total runs: %d\n", sum(design$reps)))
}
