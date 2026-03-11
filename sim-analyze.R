# sim-analyze.R -- Analysis Helpers for Registration Benchmark v2
#
# Provides:
#   load_study_a()        - Load and annotate Study A results
#   dgp_factors()         - DGP factor lookup table
#   dgp_desc_labels()     - Descriptive DGP labels for plots
#   aggregate_metric()    - Aggregate one metric by cells
#   aggregate_median_iqr() - Median + IQR summaries by cells
#   ratio_to_baseline()   - Ratios relative to a within-table baseline
#   find_optimal_by_metric() - Optimal setting per cell
#   compare_optima()      - Compare optimal settings across metrics
#   rank_summary_table()  - Rank summaries with uncertainty
#   mean_se_table()       - Mean/SE summaries for derived quantities
#   mc_summary()          - Aggregation with MC standard errors
#   amplitude_ladders()   - Paired DGP comparison ladders for Q2
#   warp_ladders()        - Paired DGP comparison ladders for Q3
#   degradation_ratios()  - Noise/severity sensitivity ratios for Q6
#   make_heatmap()        - DGP x Method heatmap (rank-based coloring)
#   make_ladder_plot()    - Amplitude/warp ladder line plot
#   make_reg_ratio_plot() - Registration ratio dot plot
#   make_dgp_example()    - Example data visualizations per DGP
#
# Plotting conventions:
#   - Heatmaps: rank-based coloring per row (yellow = rank 1, purple = rank 5)
#   - No scientific notation; round to 3 significant digits
#   - All ratios on log2 scale
#   - Horizontal boxplots where method is on the axis (text labels on y)

library(data.table)

base_dir <- here::here()
results_dir <- file.path(base_dir, "results")

# --- Formatting ---------------------------------------------------------------

#' Format numbers: avoid scientific notation, 3 significant digits
fmt <- function(x, digits = 3) {
  ifelse(
    is.na(x),
    "NA",
    ifelse(
      abs(x) < 0.001 & x != 0,
      formatC(signif(x, digits), format = "e", digits = 1),
      formatC(signif(x, digits), format = "fg", flag = "")
    )
  )
}

# --- DGP Factor Table --------------------------------------------------------

dgp_factors <- function() {
  data.table(
    dgp = paste0("D", sprintf("%02d", 1:15)),
    template = c(
      "harmonic",
      "harmonic",
      "wiggly",
      "harmonic",
      "harmonic",
      "wiggly",
      "harmonic",
      "wiggly",
      "harmonic",
      "harmonic",
      "wiggly",
      "wiggly",
      "wiggly",
      "wiggly",
      "harmonic"
    ),
    warp_type = c(
      "simple",
      "complex",
      "complex",
      "affine",
      "simple",
      "simple",
      "complex",
      "complex",
      "simple",
      "complex",
      "simple",
      "complex",
      "affine",
      "simple",
      "affine"
    ),
    amplitude = c(
      "none",
      "none",
      "none",
      "none",
      "rank1",
      "rank1",
      "rank2",
      "rank2",
      "highrank",
      "highrank",
      "highrank",
      "highrank",
      "highrank",
      "none",
      "highrank"
    )
  )
}

#' Descriptive DGP labels: "D01 (harmonic-simple-none)"
dgp_desc_labels <- function() {
  f <- dgp_factors()
  setNames(
    paste0(f$dgp, " (", f$template, "-", f$warp_type, "-", f$amplitude, ")"),
    f$dgp
  )
}

#' Short descriptive DGP labels for axis/facet use: "D01\nharm-simp-none"
dgp_short_labels <- function() {
  f <- dgp_factors()
  short_tpl <- ifelse(f$template == "harmonic", "harm", "wigg")
  short_warp <- c(
    simple = "simp",
    complex = "comp",
    affine = "affi"
  )[f$warp_type]
  setNames(
    paste0(f$dgp, "\n", short_tpl, "-", short_warp, "-", f$amplitude),
    f$dgp
  )
}

# --- Load Results -------------------------------------------------------------

load_study_a <- function(results_dir = file.path(base_dir, "results")) {
  files <- list.files(
    results_dir,
    pattern = "^results_D\\d+_\\w+_A\\.rds$",
    full.names = TRUE
  )
  if (length(files) == 0) stop("No Study A result files found in ", results_dir)

  dt <- rbindlist(lapply(files, readRDS))

  factors <- dgp_factors()
  dt <- merge(dt, factors, by = "dgp", all.x = TRUE)

  dt[, dgp := factor(dgp, levels = paste0("D", sprintf("%02d", 1:15)))]
  dt[,
    method := factor(
      method,
      levels = c(
        "srvf",
        "fda_default",
        "fda_crit1",
        "affine_ss",
        "landmark_auto"
      )
    )
  ]
  dt[,
    warp_type := factor(warp_type, levels = c("affine", "simple", "complex"))
  ]
  dt[,
    amplitude := factor(
      amplitude,
      levels = c("none", "rank1", "rank2", "highrank")
    )
  ]
  dt[, severity := factor(severity)]
  dt[, noise_sd := factor(noise_sd)]

  message(sprintf(
    "Loaded %d rows | %d DGPs | %d methods | %.1f%% failures",
    nrow(dt),
    uniqueN(dt$dgp),
    uniqueN(dt$method),
    100 * mean(dt$failure, na.rm = TRUE)
  ))
  dt
}

load_study_b <- function(results_dir = file.path(base_dir, "results")) {
  files <- list.files(
    results_dir,
    pattern = "^results_D\\d+_\\w+_B\\.rds$",
    full.names = TRUE
  )
  if (length(files) == 0) stop("No Study B result files found in ", results_dir)

  dt <- rbindlist(lapply(files, readRDS))

  factors <- dgp_factors()
  dt <- merge(dt, factors, by = "dgp", all.x = TRUE)

  dt[, dgp := factor(dgp, levels = sort(unique(dgp)))]
  dt[,
    method := factor(
      method,
      levels = c("srvf", "fda_default", "fda_crit1")
    )
  ]
  dt[, severity := factor(severity)]
  dt[, noise_sd := factor(noise_sd)]
  dt[, lambda := as.numeric(lambda)]

  message(sprintf(
    "Loaded %d rows | %d DGPs | %d methods | %d lambdas | %.1f%% failures",
    nrow(dt),
    uniqueN(dt$dgp),
    uniqueN(dt$method),
    uniqueN(dt$lambda),
    100 * mean(dt$failure, na.rm = TRUE)
  ))
  dt
}

load_study_c <- function(results_dir = file.path(base_dir, "results")) {
  files <- list.files(
    results_dir,
    pattern = "^results_D\\d+_\\w+_C\\.rds$",
    full.names = TRUE
  )
  if (length(files) == 0) stop("No Study C result files found in ", results_dir)

  dt <- rbindlist(lapply(files, readRDS))

  factors <- dgp_factors()
  dt <- merge(dt, factors, by = "dgp", all.x = TRUE)

  dt[, dgp := factor(dgp, levels = sort(unique(dgp)))]
  dt[,
    method := factor(
      method,
      levels = c(
        "srvf",
        "fda_default",
        "fda_crit1",
        "affine_ss",
        "landmark_auto"
      )
    )
  ]
  dt[,
    warp_type := factor(warp_type, levels = c("affine", "simple", "complex"))
  ]
  dt[, n_grid := factor(n_grid)]
  dt[, noise_sd := factor(noise_sd)]

  message(sprintf(
    "Loaded %d rows | %d DGPs | %d methods | %d grids | %.1f%% failures",
    nrow(dt),
    uniqueN(dt$dgp),
    uniqueN(dt$method),
    uniqueN(dt$n_grid),
    100 * mean(dt$failure, na.rm = TRUE)
  ))
  dt
}

load_study_d <- function(results_dir = file.path(base_dir, "results")) {
  files <- list.files(
    results_dir,
    pattern = "^results_D\\d+_\\w+_D\\.rds$",
    full.names = TRUE
  )
  if (length(files) == 0) stop("No Study D result files found in ", results_dir)

  dt <- rbindlist(lapply(files, readRDS))

  factors <- dgp_factors()
  dt <- merge(dt, factors, by = "dgp", all.x = TRUE)

  dt[, dgp := factor(dgp, levels = sort(unique(dgp)))]
  dt[,
    method := factor(
      method,
      levels = c(
        "srvf",
        "fda_default",
        "fda_crit1",
        "affine_ss",
        "landmark_auto"
      )
    )
  ]
  dt[,
    warp_type := factor(warp_type, levels = c("affine", "simple", "complex"))
  ]
  dt[, severity := factor(severity)]
  dt[, noise_sd := factor(noise_sd)]
  dt[,
    template_mode := factor(
      ifelse(use_true_template, "oracle", "estimated"),
      levels = c("estimated", "oracle")
    )
  ]
  dt[,
    condition := factor(
      ifelse(
        noise_sd == "0.1" & severity == "0.5",
        "easy",
        "hard"
      ),
      levels = c("easy", "hard")
    )
  ]

  message(sprintf(
    "Loaded %d rows | %d DGPs | %d methods | %.1f%% failures",
    nrow(dt),
    uniqueN(dt$dgp),
    uniqueN(dt$method),
    100 * mean(dt$failure, na.rm = TRUE)
  ))
  dt
}

load_study_e <- function(results_dir = file.path(base_dir, "results")) {
  files <- list.files(
    results_dir,
    pattern = "^results_D\\d+_\\w+_E\\.rds$",
    full.names = TRUE
  )
  if (length(files) == 0) stop("No Study E result files found in ", results_dir)

  dt <- rbindlist(lapply(files, readRDS))

  factors <- dgp_factors()
  dt <- merge(dt, factors, by = "dgp", all.x = TRUE)

  dt[, dgp := factor(dgp, levels = sort(unique(dgp)))]
  dt[,
    method := factor(
      method,
      levels = c(
        "srvf",
        "fda_default",
        "fda_crit1",
        "affine_ss",
        "landmark_auto"
      )
    )
  ]
  dt[,
    warp_type := factor(warp_type, levels = c("affine", "simple", "complex"))
  ]
  dt[, noise_sd := factor(noise_sd)]
  dt[, contam_frac := factor(contam_frac)]
  dt[, outlier_type := factor(outlier_type, levels = c("shape", "phase"))]

  message(sprintf(
    "Loaded %d rows | %d DGPs | %d methods | %d fractions | %d outlier types | %.1f%% failures",
    nrow(dt),
    uniqueN(dt$dgp),
    uniqueN(dt$method),
    uniqueN(dt$contam_frac),
    uniqueN(dt$outlier_type),
    100 * mean(dt$failure, na.rm = TRUE)
  ))
  dt
}

# --- Aggregation with MC SEs -------------------------------------------------

aggregate_metric <- function(
  dt,
  metric,
  by,
  stat = c("median", "mean"),
  value_name = NULL
) {
  stopifnot(
    "dt must be data.table" = data.table::is.data.table(dt),
    "metric must be length-1 character" = is.character(metric) &&
      length(metric) == 1L,
    "by must be non-empty character" = is.character(by) && length(by) >= 1L
  )
  stat <- match.arg(stat)
  if (is.null(value_name)) {
    value_name <- metric
  }

  dt_ok <- dt[failure == FALSE & !is.na(get(metric))]
  agg <- switch(
    stat,
    median = dt_ok[, .(value = median(get(metric), na.rm = TRUE)), by = by],
    mean = dt_ok[, .(value = mean(get(metric), na.rm = TRUE)), by = by]
  )
  data.table::setnames(agg, "value", value_name)
  agg[]
}

aggregate_median_iqr <- function(dt, metric, by, value_name = NULL) {
  stopifnot(
    "dt must be data.table" = data.table::is.data.table(dt),
    "metric must be length-1 character" = is.character(metric) &&
      length(metric) == 1L,
    "by must be non-empty character" = is.character(by) && length(by) >= 1L
  )
  if (is.null(value_name)) {
    value_name <- metric
  }

  out <- dt[
    failure == FALSE & !is.na(get(metric)),
    .(
      value = median(get(metric), na.rm = TRUE),
      q25 = quantile(get(metric), 0.25, na.rm = TRUE),
      q75 = quantile(get(metric), 0.75, na.rm = TRUE)
    ),
    by = by
  ]
  data.table::setnames(out, "value", value_name)
  out[]
}

ratio_to_baseline <- function(
  dt,
  metric,
  by,
  baseline_col,
  baseline_level,
  stat = c("median", "mean"),
  value_name = NULL,
  baseline_name = "baseline_value",
  ratio_name = "ratio",
  keep_baseline = FALSE
) {
  stopifnot(
    "baseline_col must be length-1 character" = is.character(baseline_col) &&
      length(baseline_col) == 1L,
    "baseline_col must appear in by" = baseline_col %in% by
  )
  if (is.null(value_name)) {
    value_name <- metric
  }

  agg <- aggregate_metric(
    dt = dt,
    metric = metric,
    by = by,
    stat = stat,
    value_name = value_name
  )
  key_cols <- setdiff(by, baseline_col)
  baseline_dt <- agg[
    get(baseline_col) == baseline_level,
    c(key_cols, value_name),
    with = FALSE
  ]
  data.table::setnames(baseline_dt, value_name, baseline_name)

  out <- merge(agg, baseline_dt, by = key_cols)
  out[, (ratio_name) := get(value_name) / get(baseline_name)]
  if (!keep_baseline) {
    out <- out[get(baseline_col) != baseline_level]
  }
  out[]
}

find_optimal_by_metric <- function(
  dt,
  metric,
  by,
  option_col,
  stat = c("median", "mean"),
  lower_is_better = TRUE,
  value_name = NULL
) {
  stopifnot(
    "option_col must be length-1 character" = is.character(option_col) &&
      length(option_col) == 1L,
    "option_col must not appear in by" = !option_col %in% by
  )
  if (is.null(value_name)) {
    value_name <- metric
  }

  agg <- aggregate_metric(
    dt = dt,
    metric = metric,
    by = c(by, option_col),
    stat = match.arg(stat),
    value_name = value_name
  )
  optimal <- if (lower_is_better) {
    agg[, .SD[which.min(get(value_name))], by = by]
  } else {
    agg[, .SD[which.max(get(value_name))], by = by]
  }
  optimal[]
}

compare_optima <- function(
  opt_a,
  opt_b,
  key_cols,
  value_cols,
  value_names = c("value_a", "value_b")
) {
  stopifnot(
    "opt_a must be data.table" = data.table::is.data.table(opt_a),
    "opt_b must be data.table" = data.table::is.data.table(opt_b),
    "value_cols must have length 2" = length(value_cols) == 2L,
    "value_names must have length 2" = length(value_names) == 2L
  )

  a_dt <- opt_a[, c(key_cols, value_cols[1]), with = FALSE]
  b_dt <- opt_b[, c(key_cols, value_cols[2]), with = FALSE]
  data.table::setnames(a_dt, value_cols[1], value_names[1])
  data.table::setnames(b_dt, value_cols[2], value_names[2])

  out <- merge(a_dt, b_dt, by = key_cols)
  out[, agree := get(value_names[1]) == get(value_names[2])]
  out[]
}

rank_summary_table <- function(
  dt,
  metric,
  by,
  rank_within,
  stat = c("median", "mean"),
  lower_is_better = TRUE
) {
  stopifnot(
    "rank_within must be non-empty character" = is.character(rank_within) &&
      length(rank_within) >= 1L,
    "rank_within must be subset of by" = all(rank_within %in% by)
  )

  agg <- aggregate_metric(
    dt = dt,
    metric = metric,
    by = by,
    stat = match.arg(stat),
    value_name = metric
  )
  if (lower_is_better) {
    agg[, rank := rank(get(metric), ties.method = "min"), by = rank_within]
  } else {
    agg[, rank := rank(-get(metric), ties.method = "min"), by = rank_within]
  }

  summary_cols <- setdiff(by, rank_within)
  summary_dt <- agg[,
    .(
      mean_rank = mean(rank),
      se_rank = if (.N > 1) sd(rank) / sqrt(.N) else 0,
      n_rank1 = sum(rank == 1),
      n_top2 = sum(rank <= 2),
      n = .N
    ),
    by = summary_cols
  ]

  list(cells = agg, summary = summary_dt)
}

mean_se_table <- function(dt, value_col, by) {
  stopifnot(
    "dt must be data.table" = data.table::is.data.table(dt),
    "value_col must be length-1 character" = is.character(value_col) &&
      length(value_col) == 1L,
    "by must be character" = is.character(by)
  )

  dt[
    !is.na(get(value_col)),
    .(
      mean = mean(get(value_col), na.rm = TRUE),
      se = if (.N > 1) sd(get(value_col), na.rm = TRUE) / sqrt(.N) else 0,
      median = median(get(value_col), na.rm = TRUE),
      n = .N
    ),
    by = by
  ]
}

fmt_mean_se <- function(mean, se, digits = 2, scale = 1, suffix = "") {
  sprintf(
    paste0("%.", digits, "f ± %.", digits, "f%s"),
    mean * scale,
    se * scale,
    suffix
  )
}

mc_summary <- function(dt, metric, by) {
  metric_chr <- as.character(substitute(metric))
  # Compute n_fail from unfiltered data, then summarize filtered
  fail_counts <- dt[,
    .(n_fail = sum(failure | is.na(get(metric_chr)))),
    by = by
  ]
  result <- dt[
    !is.na(get(metric_chr)) & failure == FALSE,
    .(
      mean = mean(get(metric_chr), na.rm = TRUE),
      median = median(get(metric_chr), na.rm = TRUE),
      sd = sd(get(metric_chr), na.rm = TRUE),
      mc_se = sd(get(metric_chr), na.rm = TRUE) / sqrt(.N),
      q25 = quantile(get(metric_chr), 0.25, na.rm = TRUE),
      q75 = quantile(get(metric_chr), 0.75, na.rm = TRUE),
      n = .N
    ),
    by = by
  ]
  merge(result, fail_counts, by = by)
}

# --- Ladders ------------------------------------------------------------------

amplitude_ladders <- function() {
  list(
    "harmonic + simple" = c(D01 = "none", D05 = "rank1", D09 = "highrank"),
    "harmonic + complex" = c(D02 = "none", D07 = "rank2", D10 = "highrank"),
    "wiggly + complex" = c(D03 = "none", D08 = "rank2", D12 = "highrank"),
    "wiggly + simple" = c(D14 = "none", D06 = "rank1", D11 = "highrank")
  )
}

warp_ladders <- function() {
  list(
    "harmonic + none" = c(D04 = "affine", D01 = "simple", D02 = "complex"),
    "harmonic + highrank" = c(D15 = "affine", D09 = "simple", D10 = "complex"),
    "wiggly + none" = c(D14 = "simple", D03 = "complex"),
    "wiggly + highrank" = c(D13 = "affine", D11 = "simple", D12 = "complex")
  )
}

# --- Degradation Ratios (Q6) -------------------------------------------------

degradation_ratios <- function(dt) {
  # Noise degradation: ratio of each noisy condition to clean (noise=0)
  noise_agg <- dt[
    failure == FALSE,
    .(warp_mise = median(warp_mise, na.rm = TRUE)),
    by = .(dgp, method, noise_sd, severity)
  ]
  noise_levels <- setdiff(
    as.character(sort(unique(noise_agg$noise_sd))),
    "0"
  )
  clean <- noise_agg[
    noise_sd == "0",
    .(dgp, method, severity, clean = warp_mise)
  ]
  noise_long <- noise_agg[noise_sd != "0"]
  noise_long <- merge(noise_long, clean, by = c("dgp", "method", "severity"))
  noise_long[, noise_ratio := warp_mise / clean]

  # Severity degradation
  sev_agg <- dt[
    failure == FALSE,
    .(warp_mise = median(warp_mise, na.rm = TRUE)),
    by = .(dgp, method, noise_sd, severity)
  ]
  sev_wide <- dcast(
    sev_agg,
    dgp + method + noise_sd ~ severity,
    value.var = "warp_mise"
  )
  setnames(sev_wide, c("0.5", "1"), c("low_sev", "high_sev"))
  sev_wide[, severity_ratio := high_sev / low_sev]

  list(noise = noise_long, severity = sev_wide)
}

# --- Failure Summary ----------------------------------------------------------

failure_summary <- function(dt) {
  dt[,
    .(
      n_total = .N,
      n_fail = sum(failure),
      rate = mean(failure)
    ),
    by = .(method, warp_type)
  ]
}

# --- Plotting Helpers ---------------------------------------------------------

theme_benchmark <- function(base_size = 11) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      strip.text = ggplot2::element_text(face = "bold"),
      legend.position = "bottom",
      panel.grid.minor = ggplot2::element_blank()
    )
}

method_colors <- function() {
  c(
    srvf = "#E41A1C",
    fda_default = "#377EB8",
    fda_crit1 = "#4DAF4A",
    affine_ss = "#984EA3",
    landmark_auto = "#FF7F00"
  )
}

method_labels <- function() {
  c(
    srvf = "SRVF",
    fda_default = "FDA (default)",
    fda_crit1 = "FDA (crit 1)",
    affine_ss = "Affine (S+S)",
    landmark_auto = "Landmark (auto)"
  )
}

#' Condition labeller for severity × noise facets
cond_labeller <- ggplot2::labeller(
  severity = function(x) paste0("sev=", x),
  noise_sd = function(x) paste0("noise=", x)
)

# --- Heatmap (rank-based coloring) -------------------------------------------

#' Heatmap: rank-based coloring per row.
#' Yellow (rank 1) = best in row, purple (rank 5) = worst.
#' Actual median values shown as white labels.
make_heatmap <- function(dt, metric, label = metric, lower_is_better = TRUE) {
  agg <- dt[
    failure == FALSE,
    .(value = median(get(metric), na.rm = TRUE)),
    by = .(dgp, method)
  ]
  agg[, label_text := fmt(value)]

  if (lower_is_better) {
    agg[, rank := rank(value, ties.method = "min"), by = dgp]
  } else {
    agg[, rank := rank(-value, ties.method = "min"), by = dgp]
  }

  dlabs <- dgp_short_labels()

  ggplot2::ggplot(agg, ggplot2::aes(method, dgp, fill = rank)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.5) +
    ggplot2::geom_label(
      ggplot2::aes(label = label_text),
      size = 2.3,
      fill = "white",
      alpha = 0.85,
      label.size = 0,
      label.padding = ggplot2::unit(1.5, "pt")
    ) +
    ggplot2::scale_x_discrete(labels = method_labels()) +
    ggplot2::scale_y_discrete(labels = dlabs) +
    ggplot2::scale_fill_viridis_c(
      option = "plasma",
      direction = -1,
      limits = c(1, 5),
      breaks = 1:5,
      labels = paste("Rank", 1:5),
      name = NULL
    ) +
    ggplot2::labs(x = NULL, y = NULL, subtitle = label) +
    theme_benchmark()
}

# --- Win-Rate Heatmap ---------------------------------------------------------

#' Win-rate heatmap: fraction of runs where each method has the best metric
#' value for a given DGP (across all severity/noise/rep conditions).
#' A "win" means rank 1 among the non-failed methods for that specific run.
make_winrate_heatmap <- function(
  dt,
  metric,
  label = metric,
  lower_is_better = TRUE
) {
  # For each (dgp, severity, noise_sd, rep), rank methods
  run_dt <- dt[
    failure == FALSE & !is.na(get(metric)),
    .(dgp, method, severity, noise_sd, rep, value = get(metric))
  ]
  if (lower_is_better) {
    run_dt[,
      rnk := rank(value, ties.method = "min"),
      by = .(dgp, severity, noise_sd, rep)
    ]
  } else {
    run_dt[,
      rnk := rank(-value, ties.method = "min"),
      by = .(dgp, severity, noise_sd, rep)
    ]
  }
  run_dt[, win := rnk == 1L]

  winrate <- run_dt[,
    .(win_pct = 100 * mean(win)),
    by = .(dgp, method)
  ]
  winrate[, label_text := paste0(round(win_pct), "%")]

  dlabs <- dgp_short_labels()

  ggplot2::ggplot(winrate, ggplot2::aes(method, dgp, fill = win_pct)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.5) +
    ggplot2::geom_label(
      ggplot2::aes(label = label_text),
      size = 2.3,
      fill = "white",
      alpha = 0.85,
      label.size = 0,
      label.padding = ggplot2::unit(1.5, "pt")
    ) +
    ggplot2::scale_x_discrete(labels = method_labels()) +
    ggplot2::scale_y_discrete(labels = dlabs) +
    ggplot2::scale_fill_gradient(
      low = "grey90",
      high = "#E41A1C",
      limits = c(0, 100),
      name = "Win %"
    ) +
    ggplot2::labs(x = NULL, y = NULL, subtitle = label) +
    theme_benchmark()
}

# --- Ladder Plot --------------------------------------------------------------

make_ladder_plot <- function(
  dt,
  ladder,
  ladder_name,
  metric = "warp_mise",
  metric_label = "Warp MISE"
) {
  dgps <- names(ladder)
  sub <- dt[dgp %in% dgps & failure == FALSE]
  sub[, ladder_label := ladder[as.character(dgp)]]
  sub[, ladder_label := factor(ladder_label, levels = ladder)]

  agg <- sub[,
    .(
      median = median(get(metric), na.rm = TRUE),
      q25 = quantile(get(metric), 0.25, na.rm = TRUE),
      q75 = quantile(get(metric), 0.75, na.rm = TRUE)
    ),
    by = .(method, ladder_label, severity, noise_sd)
  ]

  ggplot2::ggplot(
    agg,
    ggplot2::aes(
      x = ladder_label,
      y = median,
      color = method,
      group = method
    )
  ) +
    ggplot2::geom_line(linewidth = 0.7) +
    ggplot2::geom_point(size = 2) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = q25, ymax = q75),
      width = 0.1,
      alpha = 0.5
    ) +
    ggplot2::facet_grid(
      . ~ severity + noise_sd,
      labeller = cond_labeller
    ) +
    ggplot2::scale_color_manual(
      values = method_colors(),
      labels = method_labels()
    ) +
    ggplot2::scale_y_log10() +
    ggplot2::labs(
      title = ladder_name,
      x = NULL,
      y = metric_label,
      color = "Method"
    ) +
    theme_benchmark()
}

# --- Registration Ratio Plot (Q4) --------------------------------------------

make_reg_ratio_plot <- function(dt) {
  agg <- dt[
    failure == FALSE & !is.na(registration_ratio_median),
    .(
      median_ratio = median(registration_ratio_median, na.rm = TRUE),
      q25 = quantile(registration_ratio_median, 0.25, na.rm = TRUE),
      q75 = quantile(registration_ratio_median, 0.75, na.rm = TRUE)
    ),
    by = .(dgp, method, warp_type, severity, noise_sd)
  ]

  ggplot2::ggplot(
    agg,
    ggplot2::aes(y = method, x = median_ratio, color = method)
  ) +
    ggplot2::geom_vline(xintercept = 1, linetype = "dashed", color = "grey40") +
    ggplot2::geom_jitter(height = 0.15, size = 1.5, alpha = 0.5) +
    ggplot2::stat_summary(
      fun = median,
      geom = "crossbar",
      width = 0.4,
      color = "black",
      linewidth = 0.4
    ) +
    ggplot2::facet_grid(
      warp_type ~ severity + noise_sd,
      labeller = cond_labeller
    ) +
    ggplot2::scale_color_manual(
      values = method_colors(),
      labels = method_labels(),
      guide = "none"
    ) +
    ggplot2::scale_y_discrete(labels = method_labels()) +
    ggplot2::scale_x_continuous(
      trans = "log2",
      breaks = c(0.25, 0.5, 1, 2, 4),
      labels = c("1/4", "1/2", "1", "2", "4")
    ) +
    ggplot2::labs(
      x = expression("Registration ratio (" * log[2] * " scale)"),
      y = NULL
    ) +
    theme_benchmark()
}

# --- DGP Gallery --------------------------------------------------------------

#' Convert tf object to long data.frame for ggplot2
tf_to_df <- function(x) {
  m <- as.matrix(x)
  arg <- as.numeric(colnames(m))
  id <- rep(seq_len(nrow(m)), each = length(arg))
  data.frame(
    id = factor(id),
    arg = rep(arg, nrow(m)),
    value = as.vector(t(m))
  )
}

#' Generate example data plots for one DGP
#' Returns a patchwork of 3 panels: template, warps, observed data
make_dgp_example <- function(
  dgp_name,
  n = 20,
  severity = 0.5,
  noise_sd = 0,
  seed = 42
) {
  source(file.path(base_dir, "sim-dgp.R"))
  data <- generate_data(
    dgp_name,
    n = n,
    n_grid = 101,
    severity = severity,
    noise_sd = noise_sd,
    seed = seed
  )
  spec <- dgp_spec(dgp_name)
  desc <- dgp_desc_labels()[[dgp_name]]

  # Template
  tpl_df <- tf_to_df(data$template)
  p_template <- ggplot2::ggplot(tpl_df, ggplot2::aes(arg, value)) +
    ggplot2::geom_line() +
    ggplot2::labs(
      title = desc,
      x = "t",
      y = "template"
    ) +
    theme_benchmark(base_size = 9)

  # Warps
  warp_df <- tf_to_df(data$warps)
  p_warps <- ggplot2::ggplot(warp_df, ggplot2::aes(arg, value, group = id)) +
    ggplot2::geom_line(alpha = 0.4) +
    ggplot2::geom_abline(
      slope = 1,
      intercept = 0,
      linetype = "dashed",
      color = "grey40"
    ) +
    ggplot2::labs(
      title = sprintf("Warps (%s)", spec$phase),
      x = "t",
      y = "h(t)"
    ) +
    theme_benchmark(base_size = 9)

  # Observed data
  obs_df <- tf_to_df(data$x)
  p_data <- ggplot2::ggplot(obs_df, ggplot2::aes(arg, value, group = id)) +
    ggplot2::geom_line(alpha = 0.3) +
    ggplot2::labs(
      title = sprintf("Observed (amp: %s)", spec$amplitude),
      x = "t",
      y = "x(t)"
    ) +
    theme_benchmark(base_size = 9)

  p_template + p_warps + p_data + patchwork::plot_layout(nrow = 1)
}

#' Generate contaminated data example for one DGP × outlier type
#' Returns a patchwork of 3 panels: template, warps, observed data
#' Outlier curves shown in red (alpha = 0.5)
#' For phase outliers, the extreme warps (severity=3.0) are shown in red
#' in the warps panel (replacing the original warps for those curves).
make_contam_example <- function(
  dgp_name,
  outlier_type = "shape",
  contam_frac = 0.20,
  n = 20,
  severity = 0.5,
  noise_sd = 0.1,
  seed = 42
) {
  source(file.path(base_dir, "sim-dgp.R"))
  data <- generate_data(
    dgp_name,
    n = n,
    n_grid = 101,
    severity = severity,
    noise_sd = noise_sd,
    seed = seed
  )
  arg <- data$arg

  data <- contaminate_data(
    data,
    contam_frac = contam_frac,
    outlier_type = outlier_type,
    noise_sd = noise_sd,
    seed = seed
  )
  spec <- dgp_spec(dgp_name)
  desc <- dgp_desc_labels()[[dgp_name]]
  mask <- data$outlier_mask
  n_outliers <- sum(mask)

  # For phase outliers: regenerate the extreme warps for display
  # (contaminate_data replaces observed curves but doesn't store the warps)
  if (outlier_type == "phase" && n_outliers > 0) {
    old_seed <- if (exists(".Random.seed", envir = .GlobalEnv)) {
      .Random.seed
    } else {
      NULL
    }
    set.seed(seed + 10000L)
    # Advance RNG past sample.int (same call as contaminate_data)
    sample.int(n, n_outliers)
    # Now generate the same extreme warps
    extreme_warps <- generate_smooth_warps(
      n_outliers,
      arg,
      severity = 3.0,
      gp_type = "simple"
    )
    if (!is.null(old_seed)) .Random.seed <<- old_seed
  }

  # Template
  tpl_df <- tf_to_df(data$template)
  p_template <- ggplot2::ggplot(tpl_df, ggplot2::aes(arg, value)) +
    ggplot2::geom_line() +
    ggplot2::labs(
      title = sprintf("%s — %s outliers", desc, outlier_type),
      x = "t",
      y = "template"
    ) +
    theme_benchmark(base_size = 9)

  # Warps panel
  if (outlier_type == "phase" && n_outliers > 0) {
    # Show inlier warps in black + extreme outlier warps in red
    inlier_warps <- data$warps[!mask]
    warp_inlier_df <- tf_to_df(inlier_warps)
    warp_inlier_df$outlier <- FALSE

    warp_outlier_df <- tf_to_df(extreme_warps)
    warp_outlier_df$outlier <- TRUE
    # Re-index outlier IDs to avoid overlap
    warp_outlier_df$id <- factor(
      as.integer(warp_outlier_df$id) + length(inlier_warps)
    )

    warp_df <- rbind(warp_inlier_df, warp_outlier_df)
  } else {
    warp_df <- tf_to_df(data$warps)
    warp_df$outlier <- rep(mask, each = length(unique(warp_df$arg)))
  }

  p_warps <- ggplot2::ggplot(
    warp_df,
    ggplot2::aes(arg, value, group = id, color = outlier)
  ) +
    ggplot2::geom_line(
      data = warp_df[!warp_df$outlier, ],
      alpha = 0.4
    ) +
    ggplot2::geom_line(
      data = warp_df[warp_df$outlier, ],
      alpha = 0.5
    ) +
    ggplot2::geom_abline(
      slope = 1,
      intercept = 0,
      linetype = "dashed",
      color = "grey40"
    ) +
    ggplot2::scale_color_manual(
      values = c("FALSE" = "black", "TRUE" = "#E41A1C"),
      guide = "none"
    ) +
    ggplot2::labs(
      title = if (outlier_type == "phase") {
        "Warps (red = extreme, sev=3.0)"
      } else {
        sprintf("Warps (%s)", spec$phase)
      },
      x = "t",
      y = "h(t)"
    ) +
    theme_benchmark(base_size = 9)

  # Observed data with contamination
  obs_df <- tf_to_df(data$x)
  obs_df$outlier <- rep(mask, each = length(unique(obs_df$arg)))

  p_data <- ggplot2::ggplot(
    obs_df,
    ggplot2::aes(arg, value, group = id, color = outlier)
  ) +
    ggplot2::geom_line(
      data = obs_df[!obs_df$outlier, ],
      alpha = 0.3
    ) +
    ggplot2::geom_line(
      data = obs_df[obs_df$outlier, ],
      alpha = 0.5
    ) +
    ggplot2::scale_color_manual(
      values = c("FALSE" = "black", "TRUE" = "#E41A1C"),
      guide = "none"
    ) +
    ggplot2::labs(
      title = sprintf("Observed (%d%% contaminated)", round(100 * contam_frac)),
      x = "t",
      y = "x(t)"
    ) +
    theme_benchmark(base_size = 9)

  p_template + p_warps + p_data + patchwork::plot_layout(nrow = 1)
}

#' Case study: run all methods on clean vs contaminated data, return template
#' comparison plot
#'
#' @param dgp_name DGP identifier (e.g. "D02")
#' @param outlier_type "shape" or "phase"
#' @param contam_frac contamination fraction
#' @param noise_sd noise level
#' @param severity warp severity
#' @param n number of curves
#' @param seed random seed
#' @param methods character vector of method names to run
#' @return list with: plot (patchwork), metrics (data.frame)
make_contam_case_study <- function(
  dgp_name,
  outlier_type = "phase",
  contam_frac = 0.20,
  noise_sd = 0.1,
  severity = 0.5,
  n = 50,
  seed = 42,
  methods = c("srvf", "fda_default", "fda_crit1", "affine_ss", "landmark_auto")
) {
  source(file.path(base_dir, "sim-dgp.R"))
  source(file.path(base_dir, "sim-methods.R"))
  source(file.path(base_dir, "sim-metrics.R"))

  # Generate clean data
  data_clean <- generate_data(
    dgp_name,
    n = n,
    n_grid = 101,
    severity = severity,
    noise_sd = noise_sd,
    seed = seed
  )

  # Generate contaminated data (same base)
  data_contam <- generate_data(
    dgp_name,
    n = n,
    n_grid = 101,
    severity = severity,
    noise_sd = noise_sd,
    seed = seed
  )
  data_contam <- contaminate_data(
    data_contam,
    contam_frac = contam_frac,
    outlier_type = outlier_type,
    noise_sd = noise_sd,
    seed = seed
  )

  arg <- data_clean$arg
  template_true <- data_clean$template
  tpl_true_df <- tf_to_df(template_true)
  tpl_true_df$type <- "True template"

  mlabs <- method_labels()
  mcols <- method_colors()

  plots <- list()
  metrics_list <- list()

  for (m in methods) {
    # Run on clean data
    res_clean <- tryCatch(
      fit_method(data_clean, m),
      error = function(e)
        list(registration = NULL, time = NA, error = e$message)
    )
    # Run on contaminated data
    res_contam <- tryCatch(
      fit_method(data_contam, m),
      error = function(e)
        list(registration = NULL, time = NA, error = e$message)
    )

    # Extract metrics
    m_clean <- if (is.null(res_clean$error)) {
      extract_metrics(data_clean, res_clean)
    } else {
      failure_metrics(res_clean)
    }
    m_contam <- if (is.null(res_contam$error)) {
      extract_metrics(
        data_contam,
        res_contam,
        outlier_mask = data_contam$outlier_mask
      )
    } else {
      failure_metrics(res_contam)
    }
    metrics_list[[m]] <- list(clean = m_clean, contam = m_contam)

    # Build template comparison plot for this method
    for (condition in c("clean", "contam")) {
      res <- if (condition == "clean") res_clean else res_contam
      label <- if (condition == "clean") "Clean" else {
        sprintf("%d%% %s", round(100 * contam_frac), outlier_type)
      }

      if (!is.null(res$registration)) {
        tpl_est <- tf_template(res$registration)
        tpl_est_df <- tf_to_df(tpl_est)
        tpl_est_df$type <- "Estimated"

        aligned <- tf_aligned(res$registration)
        aligned_df <- tf_to_df(aligned)

        p <- ggplot2::ggplot() +
          ggplot2::geom_line(
            data = aligned_df,
            ggplot2::aes(arg, value, group = id),
            alpha = 0.15,
            color = "grey60"
          ) +
          ggplot2::geom_line(
            data = tpl_true_df,
            ggplot2::aes(arg, value),
            color = "black",
            linewidth = 1
          ) +
          ggplot2::geom_line(
            data = tpl_est_df,
            ggplot2::aes(arg, value),
            color = mcols[m],
            linewidth = 1,
            linetype = "dashed"
          ) +
          ggplot2::labs(
            title = if (condition == "clean") mlabs[m] else NULL,
            subtitle = label,
            x = "t",
            y = "x(t)"
          ) +
          theme_benchmark(base_size = 8)
      } else {
        p <- ggplot2::ggplot() +
          ggplot2::annotate("text", x = 0.5, y = 0.5, label = "FAILED") +
          ggplot2::labs(
            title = if (condition == "clean") mlabs[m] else NULL,
            subtitle = label
          ) +
          theme_benchmark(base_size = 8)
      }
      plots[[paste0(m, "_", condition)]] <- p
    }
  }

  # Arrange: methods as columns, clean/contam as rows
  plot_list <- list()
  for (m in methods) {
    plot_list <- c(
      plot_list,
      list(
        plots[[paste0(m, "_clean")]],
        plots[[paste0(m, "_contam")]]
      )
    )
  }

  combined <- patchwork::wrap_plots(
    plot_list,
    ncol = length(methods),
    byrow = FALSE
  )

  # Build metrics summary table
  metrics_df <- do.call(
    rbind,
    lapply(methods, function(m) {
      mc <- metrics_list[[m]]$clean
      mm <- metrics_list[[m]]$contam
      data.frame(
        method = mlabs[m],
        tmise_clean = mc$template_mise,
        tmise_contam = mm$template_mise,
        edist_clean = mc$template_elastic_dist,
        edist_contam = mm$template_elastic_dist,
        stringsAsFactors = FALSE
      )
    })
  )

  list(plot = combined, metrics = metrics_df)
}
