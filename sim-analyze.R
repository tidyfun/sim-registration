# sim-analyze.R -- Analysis Helpers for Registration Benchmark v2
#
# Provides:
#   load_study_a()        - Load and annotate Study A results
#   dgp_factors()         - DGP factor lookup table
#   dgp_desc_labels()     - Descriptive DGP labels for plots
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

# --- Aggregation with MC SEs -------------------------------------------------

mc_summary <- function(dt, metric, by) {
  metric_chr <- as.character(substitute(metric))
  dt[
    !is.na(get(metric_chr)) & failure == FALSE,
    .(
      mean = mean(get(metric_chr), na.rm = TRUE),
      median = median(get(metric_chr), na.rm = TRUE),
      sd = sd(get(metric_chr), na.rm = TRUE),
      mc_se = sd(get(metric_chr), na.rm = TRUE) / sqrt(.N),
      q25 = quantile(get(metric_chr), 0.25, na.rm = TRUE),
      q75 = quantile(get(metric_chr), 0.75, na.rm = TRUE),
      n = .N,
      n_fail = sum(is.na(get(metric_chr)))
    ),
    by = by
  ]
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
  dlabs <- dgp_short_labels()
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
    agg[noise_sd == "0" & severity == "0.5"],
    ggplot2::aes(y = method, x = median_ratio, color = method)
  ) +
    ggplot2::geom_vline(xintercept = 1, linetype = "dashed", color = "grey40") +
    ggplot2::geom_point(size = 2) +
    ggplot2::geom_errorbarh(
      ggplot2::aes(xmin = q25, xmax = q75),
      height = 0.2
    ) +
    ggplot2::facet_wrap(
      ~dgp,
      nrow = 2,
      labeller = ggplot2::as_labeller(dlabs)
    ) +
    ggplot2::scale_color_manual(
      values = method_colors(),
      labels = method_labels()
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
    theme_benchmark() +
    ggplot2::theme(legend.position = "none")
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
