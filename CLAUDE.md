# CLAUDE.md

Registration benchmark comparing 5 curve registration methods across 15 DGPs.
Companion study for the [tf](https://github.com/tidyfun/tf) R package.

## Project Structure

```
sim-dgp.R          DGP generation (templates, warps, amplitude, generate_data())
sim-methods.R       Method wrappers around tf_register()
sim-metrics.R       All metrics (warp MISE, alignment error, SRSF distance, etc.)
sim-config.R        Study designs (A–D), seed management
sim-run.R           Parallel runner with incremental saves
sim-validate.R      Known-answer tests (run before benchmarking)
sim-analyze.R       Analysis helpers: loading, aggregation, plotting
sim-report.qmd      Quarto report (Study A results)
slurm/              SLURM job scripts for LRZ cluster
attic/              Archived v1-era scripts and one-off analyses
results/v1/         Archived v1 results (13 DGPs, 2 noise levels)
BENCHMARK-PLAN.md   Full study design and analysis plan
```

## Running Things

```r
# Validate (from project root):
Rscript sim-validate.R

# Render report:
quarto render sim-report.qmd

# Run benchmark (local):
Rscript sim-run.R A 4

# Run benchmark (LRZ cluster):
sbatch slurm/study-a.slurm
```

## Key Dependencies

- `tf` package (dev version from `tidyfun/tf`)
- `fdasrvf` for SRVF registration
- `data.table`, `ggplot2`, `patchwork` for analysis/reporting
- `here` for path resolution (`.here` sentinel in project root)

## Conventions

- All scripts use `base_dir <- here::here()` for paths.
- Results are `.rds` files in `results/`, gitignored.
- DGPs are `D01`–`D15`; methods: `srvf`, `cc_default`, `cc_crit1`, `affine_ss`, `landmark_auto`.
- Severity (0.5, 1.0) and noise (0, 0.1, 0.3) are the experimental factors beyond DGP × method.
- Study A v2: 15 DGPs × 2 severity × 3 noise × 5 methods × 100 reps = 45,000 runs.
- Warp comparison uses forward warps (invert `tf_inv_warps()` via `tf_invert()`).
- SRSF elastic distance as secondary template quality metric.
