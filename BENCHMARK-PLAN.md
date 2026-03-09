# Registration Benchmark v2

Claude session: `claude --resume 2776a7a9-eda5-4ff6-9e56-a8b9b62a036c`

Monte Carlo benchmark comparing 5 registration methods across 15 DGPs with
3 warp types, 4 amplitude types, and new diagnostic measures. Redesign of v1
(archived in `attic/sim-registration-v1/`).

## 1. Templates (2 synthetic)

**Harmonic** (simple, existing):
```r
harmonic <- sin(2 * pi * t) + 0.5 * sin(4 * pi * t)
```
3 zero-crossings, symmetric, ~3-4 landmarks.

**Wiggly** (complex, new): sum of sinusoids at incommensurate frequencies:
```r
wiggly <- 0.6*sin(2*pi*t) + 0.5*sin(5*pi*t + 0.3) +
  0.45*cos(8*pi*t - 0.7) + 0.4*sin(11*pi*t + 1.2) +
  0.3*cos(14*pi*t)
```
~7 landmarks, pronounced wiggles throughout [0, 0.4] and [0.6, 1.0].
Both templates z-scored for comparability.

## 2. Phase Variation (3 types)

| Type | Construction | Domain-preserving? |
|------|-------------|-------------------|
| simple | GP(squareexp, scale=0.3) -> exp -> integrate -> normalize | Yes |
| complex | GP(matern, order=1.5, scale=0.05) -> exp -> integrate -> normalize | Yes |
| affine | h(t) = a*t + b, TruncNorm params, re-centered | No |

### Warp generation (smooth types)

GP-based warps via `tf_rgp()` -> positivity transform -> cumulative integral ->
endpoint normalization:

```r
gp <- tf_rgp(n, arg, cov = ..., scale = ..., nugget = 0.001)
positive <- exp(severity * (gp - mean(gp)))
warps <- tf_integrate(positive, definite = FALSE)
warps / as.numeric(warps[, max(arg)])   # normalize h(1) = 1
```

Simple warps use `squareexp` kernel (scale=0.3) for smooth, convex/concave
deviations from identity. Complex warps use `matern` kernel (order=1.5,
scale=0.05) for rapidly oscillating warps that snake around the identity.

### Affine warps

`h_i(t) = a_i * t + b_i` with `a_i ~ TruncNorm(1, sd_a)`,
`b_i ~ TruncNorm(0, sd_b)`. Re-centered for identifiability:
`a <- a - mean(a) + 1`, `b <- b - mean(b)`.

Parameters (doubled from initial values after calibration):
- `sd_a = 0.10 * severity` (was 0.05), bounds ±0.30*severity
- `sd_b = 0.06 * severity` (was 0.03), bounds ±0.15*severity
- At severity=0.5, this gives mean ISE from identity ~0.0017

Non-domain-preserving: values can extend outside [0,1]. Handled via
`approx(..., rule = 2)` (boundary extrapolation) when warping template.

## 3. Amplitude Variation (4 types)

| Type | Model | Parameters |
|------|-------|-----------|
| none | x_i(t) = m(h_i(t)) + noise | - |
| rank1_mult | x_i(t) = a_i * m(h_i(t)) | a_i ~ LogNormal(-s^2/2, s^2), s=0.3 |
| rank2 | x_i(t) = a_i * m(h_i(t)) + c_i | See calibrated params below |
| highrank | x_i(t) = m(h_i(t)) + sum_j c_ij * phi_j(t) | 3 FPCs from gait knee |

### Rank-2 amplitude calibration

Total target variance amp_sd^2 split equally between scale and shift so
total Var(a_i * m + c_i) ~ amp_sd^2, comparable to other types:
```r
scale_sd <- amp_sd / sqrt(2)   # ~0.212
shift_sd <- amp_sd / sqrt(2)   # ~0.212
a <- rlnorm(n, meanlog = -scale_sd^2/2, sdlog = scale_sd)
c_shift <- rnorm(n, 0, shift_sd)
```

### High-rank amplitude details

FPCs extracted from registered gait knee angles:
```r
reg <- tf_register(gait$knee_angle, method = "srvf")
fpc <- tfb_fpc(tf_aligned(reg), pve = 0.99)
```
Cached in `gait_knee_fpcs.rds` (20-point grid, 3 eigenfunctions).

Eigenfunctions interpolated from 20 to 101 points, then **re-orthonormalized
via QR decomposition** to ensure correct variance calibration:
```r
phi_interp <- approx(gait_arg, phi, xout = arg)$y
qr_phi <- qr(phi_interp * sqrt(dt))
phi_interp <- qr.Q(qr_phi) / sqrt(dt)
```

Eigenvalue decay flattened to 3:2:1 ratio, then scaled so total amplitude
variance matches `amp_sd^2` (default 0.3^2 = 0.09).

## 4. DGP Table (15 DGPs)

### Baseline (phase only)
| DGP | Template | Warp | Amplitude | Purpose |
|-----|----------|------|-----------|---------|
| D01 | harmonic | simple | none | Easy baseline |
| D02 | harmonic | complex | none | Warp complexity effect |
| D03 | wiggly | complex | none | Template + warp complexity |
| D04 | harmonic | affine | none | Affine warp baseline |
| D14 | wiggly | simple | none | Wiggly baseline (completes wiggly+simple ladder) |

### Low-rank amplitude
| DGP | Template | Warp | Amplitude | Purpose |
|-----|----------|------|-----------|---------|
| D05 | harmonic | simple | rank1_mult | Mild scenario |
| D06 | wiggly | simple | rank1_mult | Template complexity + amplitude |
| D07 | harmonic | complex | rank2 | Warp complexity + 2-param amplitude |
| D08 | wiggly | complex | rank2 | Everything moderately hard |

### High-rank amplitude (the novelty)
| DGP | Template | Warp | Amplitude | Purpose |
|-----|----------|------|-----------|---------|
| D09 | harmonic | simple | highrank | Easy phase, realistic amplitude |
| D10 | harmonic | complex | highrank | Hard phase + realistic amplitude |
| D11 | wiggly | simple | highrank | Hard template, easy phase, hard amplitude |
| D12 | wiggly | complex | highrank | Hardest overall |
| D13 | wiggly | affine | highrank | Non-domain-preserving + high-rank amplitude |
| D15 | harmonic | affine | highrank | Affine + highrank (completes harmonic+highrank ladder) |

Every template (2), warp (3), and amplitude (4) level appears at least twice.
D14 and D15 were added after council review to complete amplitude and warp
comparison ladders that previously had missing starting points.

## 5. Methods (5)

| Method | Call | Notes |
|--------|------|-------|
| srvf | `tf_register(x, method="srvf")` | SRVF elastic |
| fda_default | `tf_register(x, method="fda")` | FDA criterion 2 |
| fda_crit1 | `tf_register(x, method="fda", crit=1)` | FDA criterion 1 |
| affine_ss | `tf_register(x, method="affine", type="shift_scale")` | Shift+scale |
| landmark_auto | `tf_register(x, method="landmark", landmarks=...)` | Auto landmarks |

All methods estimate their own template (no oracle in main study).
Landmarks auto-detected via `tf_landmarks_extrema(tf_smooth(x), "both")`.
Template-based methods use `max_iter = 10` Procrustes iterations when estimating.

## 6. Performance Measures

### Primary
- **Warp MISE**: `mean(tf_integrate((h_est - h_true)^2))`
  Always compares **forward warps**. `tf_registration` stores inverse warps,
  so estimated inverse warps are inverted once via `tf_invert()` to recover
  the forward warp for comparison with the true forward warp.
- **Alignment error**: `mean(tf_integrate((aligned_i - true_prewarp_i)^2))`
- **Template MISE**: `tf_integrate((template_est - template_true)^2)` (L^2)
- **Template elastic distance** (SRSF/Fisher-Rao): `||SRSF(f_est) - SRSF(f_true)||_L2`
  where `SRSF q(t) = sign(f'(t)) * sqrt(|f'(t)|)`. Measures template shape
  distance invariant to reparameterization -- unlike L^2 MISE which conflates
  shape and phase error in the template estimate.

### Secondary
- **Alignment CC**: `mean(tf_crosscor(aligned, template_true))`
- **Time per curve**: `elapsed / n_curves`

### Diagnostics
- **Registration ratio**: `||h_est - id||^2 / ||h_true - id||^2` per curve.
  Median and IQR reported. < 1 = under-registration, > 1 = over-registration.
- **Amplitude variance ratio**: `mean_t(Var_i(aligned - template_est)) / mean_t(Var_i(true_prewarp - template_true))`.
  < 1 = over-registration absorbed amplitude into phase.
  > 1 = under-registration left phase variation in amplitude.
  NA for amplitude="none" DGPs (no true amplitude to compare).
- **Warp slope ratio**: `range(h'_est) / range(h'_true)` per curve.
  Detects overly smooth (compressed range) or overly wiggly (exaggerated range)
  estimated warps. NA for affine DGPs (constant derivative).

## 7. Study Structure & Computation

### Study A: Main (45,000 runs)
15 DGPs x 2 severity (0.5, 1.0) x 3 noise (0, 0.1, 0.3) x 5 methods = 450 cells x 100 reps

The third noise level (0.3) was added to test whether method rankings (especially
SRVF vs FDA) flip under high noise, since unregularized SRVF operates on SRSFs
(derivatives) which amplify high-frequency noise.

### Study B: Penalization (7,200 runs)
3 DGPs (data-driven from Study A) x 6 lambdas x 2 methods (srvf, fda_default)
x 2 noise x 2 severity = 144 cells x 50 reps

### Study C: Oracle paired (2,000 runs)
4 DGPs (data-driven from Study A) x 2 template modes x 5 methods
x noise=0.1 x severity=0.5 = 40 cells x 50 reps

### Study D: Outlier contamination (4,500 runs)
3 representative DGPs x 2 outlier types (shape, phase) x 3 fractions (5%, 10%, 20%)
x 5 methods x noise=0.1 x severity=0.5 = 90 cells x 50 reps

Shape outliers: random sinusoids scaled to ~2x template range.
Phase outliers: extreme warps (severity=3.0, 6x default).

**Total: ~58,700 runs**

## 8. Code Structure

```
sim-registration/            Standalone repo (~/fda/sim-registration/)
  BENCHMARK-PLAN.md          This file
  gait_knee_fpcs.rds         Cached FPCs from registered gait knee angles
  sim-dgp.R          ~500 lines   Templates, warps, amplitude, generate_data()
  sim-methods.R       ~110 lines   5 method wrappers returning tf_registration
  sim-metrics.R       ~300 lines   All metrics using tf API (incl. SRSF elastic dist)
  sim-config.R        ~220 lines   All study designs (A-D)
  sim-run.R           ~300 lines   Runner with incremental saves per DGP
  sim-validate.R      ~150 lines   Known-answer tests for v2
  sim-analyze.R       ~340 lines   Data loading, aggregation, plot helpers
  sim-report.qmd    ~1200 lines   Quarto report (Study A analysis)
  smoke-test.R         ~30 lines   LRZ smoke test (D01, 2 reps)
  study-a-run.R        ~20 lines   Study A production entry point
  run.sh               ~30 lines   Local launch script
```

### Key tf-syntax improvements (vs v1)
- All `trapez_int()` replaced with `tf_integrate()`
- `sapply(tf_evaluations(...), ...)` loops replaced with tf vector arithmetic
- `fit_method()` returns `tf_registration` object directly
- Metrics use `tf_integrate((a - b)^2)`, `tf_derive()`, `tf_crosscor()`
- Warp generation uses `tf_rgp()` -> `exp()` -> `tf_integrate(definite = FALSE)` -> normalize

## 9. Implementation Phases

### Phase 0: Scaffolding (done)
- Archived old code to `attic/sim-registration-v1/`
- Extracted and cached gait-knee FPCs to `gait_knee_fpcs.rds`
  (3 eigenfunctions, 20-point grid, score variances 437/171/104)

### Phase 1: DGP generation -- sim-dgp.R (done)
- Two templates (harmonic, wiggly), both z-scored
- Three warp types via `generate_warps()`:
  - Simple/complex: GP -> exp -> integrate -> endpoint normalization
  - Affine: TruncNorm with re-centering (doubled severity after calibration)
- Four amplitude types via `generate_amplitude()` + `apply_amplitude()`
  - rank2 calibrated: scale_sd = shift_sd = amp_sd/sqrt(2) for comparable total variance
- `true_prewarp_curves()` computes ground-truth aligned curves for metrics
- `contaminate_data()` for Study D (shape and phase outliers)
- `dgp_spec()` dispatches DGP name -> (template, phase, amplitude) tuple (D01-D15)

### Phase 2: Methods + Metrics (done)
- sim-methods.R: 5 thin wrappers around `tf_register()`, all return
  `{registration, time, error}` lists
- sim-metrics.R: all metrics via tf API, `safe_metric()` error handling
  - Warp comparison: always forward warps (invert estimated inverse via `tf_invert()`)
  - SRSF elastic distance added as secondary template quality metric

### Phase 3: Config + Runner (done)
- sim-config.R: all 4 study designs, `full_design()`, deterministic seeding
- sim-run.R: `create_tasks()` -> `run_one_task()` -> `run_batch()` with
  mclapply parallelization and incremental saves per DGP

### Council Review R1 (done)

Codex, Gemini, and Claude reviewed the full codebase. Key findings and fixes:

**Applied fixes:**
1. Added `na.rm = TRUE` to `median()`/`IQR()` in diagnostic ratios
2. Re-orthonormalize FPCs after interpolation via QR decomposition
3. Use endpoint evaluation `warps[, max(arg)]` instead of `tf_fmax()` for
   warp normalization (monotone warps have max at endpoint)
4. Max-iteration guard (10K) on `rtruncnorm` rejection sampler
5. `library(tf)` with `devtools::load_all()` fallback for HPC portability

**Acknowledged but deferred:**
- Affine recentering can technically push values outside truncation bounds
  (negligible with current severity range)
- `safe_metric` silently returns NA on failure (acceptable for benchmark)
- Warp slope ratio structurally NA for affine DGPs (by design)

### Council Review R2 (done)

Second council reviewed study design, parameter ranges, and report quality:

**Applied fixes (from R2):**
1. Added D14 (wiggly+simple+none) and D15 (harmonic+affine+highrank) to complete
   amplitude and warp comparison ladders
2. Doubled affine warp parameters (sd_a: 0.05 -> 0.10, sd_b: 0.03 -> 0.06) after
   calibration showed near-identity warps at severity=0.5
3. Calibrated rank2 amplitude variance: split amp_sd equally between scale and shift
   (scale_sd = shift_sd = amp_sd/sqrt(2)) so total variance matches other types
4. Added SRSF elastic distance metric for reparameterization-invariant template
   shape comparison
5. Added noise=0.3 as third noise level to test whether SRVF vs FDA rankings flip
   under high noise (derivative amplification)
6. Warp comparison always uses forward warps (was inconsistent: inverse for DP,
   forward for affine). Now always invert estimated inverse warps via `tf_invert()`.
7. All report findings rewritten as data-driven `results: asis` code chunks
   (6/7 original speculative findings were factually wrong)

### Phase 4: Validation (done)
- Mini pilot (D01 + D09, 2 reps): all 5 methods work, 0 failures
- D01 ranking: srvf (9.5e-6) > fda_crit1 (2.4e-5) > fda_default (2.1e-4) >
  landmark_auto (1.5e-3) > affine_ss (2.2e-3) -- as expected

### Phase 5: LRZ Environment (done)
- tf + all dependencies installed on CoolMUC-4 (R 4.3.3-gcc13-mkl)
- Smoke test passed: 0/10 failures, results match local pilot exactly
- Job 5125545 completed in ~4 minutes (serial_std partition, 1 core)

### Phase 6: Study A production run v1 (done, superseded)
- Ran with 13 DGPs x 2 noise levels (26K runs), all 65/65 batches complete
- Results used for initial analysis and report draft
- **Superseded** by updated design (15 DGPs, 3 noise levels, calibrated params)

### Phase 7: Study A production run v2 (pending)
- Updated design: 15 DGPs x 2 severity x 3 noise x 5 methods = 450 cells x 100 reps = 45K runs
- Changes from v1: +2 DGPs (D14, D15), +1 noise level (0.3), doubled affine severity,
  calibrated rank2 amplitude, SRSF elastic distance metric, forward warp comparison
- Batched by (DGP, method): 75 batches, each ~600 runs
- Results saved as `results_D{nn}_{method}_A.rds`

### Phase 8: Study A Analysis (in progress)
- sim-analyze.R: data loading, DGP factor annotation, aggregation helpers
- sim-report.qmd: Quarto report implementing analysis plan below
- All findings sections use data-driven `results: asis` code chunks

### Phases 9+: Sub-studies B, C, D (pending)
- DGPs for Studies B-D selected from Study A diagnostics
- Studies B/C/D runs, analysis, integrated report

## 10. Study A Analysis Plan

### Context

Study A (45K runs) compares 5 registration methods across
15 DGPs x 2 severity x 3 noise = 450 cells x 100 reps. The DGP design is NOT full
factorial (15 of 24 possible template x warp x amplitude combos), but it enables clean
**paired comparisons** that vary exactly one factor while holding the other two constant.

**Key principle**: Never average over heterogeneous settings. Each analysis conditions
on the relevant factors and answers a specific question.

### Available metrics (per run)

| Metric | Type | Interpretation |
|--------|------|----------------|
| `warp_mise` | Primary | ISE of estimated vs true forward warps (lower = better) |
| `alignment_error` | Primary | ISE of aligned curves vs true pre-warp curves |
| `template_mise` | Primary | L^2 ISE of estimated vs true template |
| `template_elastic_dist` | Primary | SRSF distance: shape distance invariant to reparameterization |
| `alignment_cc` | Secondary | Cross-correlation of aligned vs truth (higher = better) |
| `registration_ratio_median` | Diagnostic | `||h_est-id||^2 / ||h_true-id||^2` per curve, median. <1 = under-registered, >1 = over-registered |
| `registration_ratio_iqr` | Diagnostic | IQR of above -- warp recovery heterogeneity |
| `amp_variance_ratio` | Diagnostic | Residual amplitude variance ratio (meaningful only when amplitude != none) |
| `warp_slope_ratio_median` | Diagnostic | Estimated vs true warp derivative range ratio |
| `time_per_curve` | Practical | Seconds per curve |
| `failure` | Practical | Boolean: method errored |

### DGP factor structure

```
         template    warp      amplitude
D01      harmonic    simple    none
D02      harmonic    complex   none
D03      wiggly      complex   none
D04      harmonic    affine    none
D05      harmonic    simple    rank1
D06      wiggly      simple    rank1
D07      harmonic    complex   rank2
D08      wiggly      complex   rank2
D09      harmonic    simple    highrank
D10      harmonic    complex   highrank
D11      wiggly      simple    highrank
D12      wiggly      complex   highrank
D13      wiggly      affine    highrank
D14      wiggly      simple    none
D15      harmonic    affine    highrank
```

### Research questions -> metrics -> conditioning

**Q1: Which method best recovers true warps, conditional on warp type?**

- **Metric**: `warp_mise` (distribution, not just median)
- **Stratify by**: warp_type (simple / complex / affine)
- **Within each stratum**: boxplots/violins by method, faceted by severity
- **Aggregate within stratum**: OK to pool across templates and amplitude types IF
  the pattern is consistent (check with interaction plot first)
- **Visualization**: 3-panel figure (one per warp type), each showing method x severity

**Q2: Does amplitude variation degrade warp recovery? (Phase-amplitude separation)**

- **Amplitude ladders** (clean comparisons, each varies only amplitude):
  - harmonic + simple: D01(none) -> D05(rank1) -> D09(highrank)
  - harmonic + complex: D02(none) -> D07(rank2) -> D10(highrank)
  - wiggly + complex: D03(none) -> D08(rank2) -> D12(highrank)
  - wiggly + simple: D14(none) -> D06(rank1) -> D11(highrank)
- **Metrics**: warp_mise (primary) + amp_variance_ratio (secondary)
- **Conditioning**: one panel per ladder, one line/group per method
- **Visualization**: line plots showing warp_mise vs amplitude complexity

**Q3: How does warp complexity affect recovery?**

- **Warp ladders** (clean comparisons, each varies only warp type):
  - harmonic + none: D04(affine) -> D01(simple) -> D02(complex)
  - harmonic + highrank: D15(affine) -> D09(simple) -> D10(complex)
  - wiggly + none: D14(simple) -> D03(complex)
  - wiggly + highrank: D13(affine) -> D11(simple) -> D12(complex)
- **Metric**: warp_mise, also registration_ratio
- **Conditioning**: one panel per ladder, methods as groups

**Q4: Do methods over-register or under-register?**

- **Metric**: `registration_ratio_median`
- **Conditioning**: by method x DGP, grouped by warp_type
- **Threshold**: ratio ~ 1 is ideal, <1 = under, >1 = over
- **Visualization**: dot plot with reference line at 1, faceted by warp_type

**Q5: Template estimation quality**

- **Metrics**: `template_mise` (L^2 MISE) and `template_elastic_dist` (SRSF distance)
  - L^2 MISE conflates shape and phase error in the template estimate
  - SRSF elastic distance measures template shape distance invariant to
    reparameterization -- the "correct" metric but new/unfamiliar
  - Both reported; primary conclusions from elastic distance
- **Conditioning**: by method x severity x noise (NOT harmonic vs wiggly)
- **Key questions**: noise/severity sensitivity of template recovery

**Q6: Sensitivity to noise and severity**

- Report degradation ratios (noisy/clean, high/low severity) per method
- Three noise levels (0, 0.1, 0.3) enable testing whether SRVF degrades faster
  than FDA under high noise (SRVF operates on derivatives, amplifies HF noise)
- Degradation computed as ratio to baseline: noise=X vs noise=0, severity=1 vs severity=0.5
- Paired ratio dot plots -- one dot per DGP, grouped by method
- Absolute MISE values reported alongside ratios to contextualize

**Q7: Practical considerations**

- Runtime: `time_per_curve` by method
- Failure rate: by method x warp_type
- Error messages: tabulate failure modes

### Report structure (`sim-report.qmd`)

```
1. Setup
2. Overview: DGP x Method heatmaps (warp_mise, alignment_cc)
3. Warp Recovery (Q1-Q3)
   3.1 By warp type (Q1)
   3.2 Phase-amplitude interaction (Q2)
   3.3 Warp complexity effect (Q3)
4. Registration Diagnostics (Q4-Q5)
   4.1 Over/under-registration (Q4)
   4.2 Template estimation: L^2 MISE (Q5)
   4.3 Template estimation: SRSF elastic distance (Q5)
5. Sensitivity Summary (Q6)
6. Practical (Q7)
7. Summary table
```

All findings sections are **data-driven**: `results: asis` code chunks that
compute summary statistics and generate prose via `cat(sprintf(...))`. No
speculative findings (6/7 original speculative findings were factually wrong).

### Visualization choices

| Plot type | Where used | Why |
|-----------|-----------|-----|
| Heatmap (tile + text) | Section 2 overview | Compact DGP x method summary |
| Boxplot/violin | Q1 warp by type | Shows full distribution |
| Line plot (ladder) | Q2, Q3 paired comparisons | Reveals factor effects cleanly |
| Dot plot + ref line | Q4 registration ratio | Natural threshold at 1 |
| Paired ratio plot | Q6 sensitivity | Proportional effects |
| Bar chart (log scale) | Q7 runtime | Orders of magnitude difference |

## 11. Verification

1. **Unit**: Each DGP generates correct structure; metrics return 0/1 for identity/perfect alignment
2. **Pilot** (10 reps): failure rates < 5%, severity calibration, warp complexity distinguishable
3. **Sanity**: D01 (harmonic+simple+none) -- SRVF should dominate, MISE ranges comparable
4. **Seed reproducibility**: re-run 5 random cells, verify identical
5. **Diagnostic validity**: registration ratio ~ 1 for well-specified combos
6. **Affine sanity**: D04/D13/D15 -- affine_ss should outperform smooth methods
7. **Analysis verification**: ladder comparisons show monotone difficulty; MC SEs support rankings

## Open Questions (resolve during analysis)

- Which DGPs for Studies B, C, D (data-driven from Study A)
- Lambda grid for penalization arm (method-specific, possibly data-dependent)
- Whether `fda_crit1` differs meaningfully from `fda_default` in harder settings

## Reproducibility

- Deterministic seeding: `digest::digest2int(paste(dgp, rep))` per task
- Session info saved with results
- Incremental saves per DGP group (crash-safe)
- All code in `~/fda/sim-registration/`, data in `results/` subdirectory
- `.gitignore` in tf repo blocks `attic/sim-registration/` from being pushed
