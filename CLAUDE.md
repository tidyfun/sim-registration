# CLAUDE.md

Registration benchmark comparing 5 curve registration methods across 13 DGPs.
Companion study for the [tf](https://github.com/tidyfun/tf) R package.

## Project Structure

```
sim-dgp.R          DGP generation (templates, warps, amplitude, generate_data())
sim-methods.R       Method wrappers around tf_register()
sim-metrics.R       All metrics (warp MISE, alignment error, registration ratio, etc.)
sim-config.R        Study designs (A–D), seed management
sim-run.R           Parallel runner with incremental saves
sim-validate.R      Known-answer tests (run before benchmarking)
sim-analyze.R       Analysis helpers: loading, aggregation, plotting
sim-report.qmd      Quarto report (Study A results)
run.sh              SLURM job script for LRZ cluster
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
sbatch run.sh
```

## Key Dependencies

- `tf` package (dev version from `tidyfun/tf`)
- `fdasrvf` for SRVF registration
- `data.table`, `ggplot2`, `patchwork` for analysis/reporting
- `here` for path resolution (`.here` sentinel in project root)

## Conventions

- All scripts use `base_dir <- here::here()` for paths.
- Results are `.rds` files in `results/`, gitignored.
- DGPs are `D01`–`D13`; methods: `srvf`, `fda_default`, `fda_crit1`, `affine_ss`, `landmark_auto`.
- Severity (0.5, 1.0) and noise (0, 0.1) are the two experimental factors beyond DGP × method.
