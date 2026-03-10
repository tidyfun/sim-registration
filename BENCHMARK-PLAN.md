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

### Study B: Penalization (~19,000 runs)
4 DGPs x 8 lambdas x 3 methods (srvf, fda_default, fda_crit1)
x 3 noise (0, 0.1, 0.3) x 2 severity = 576 cells x ~33 reps

**Purpose**: Test whether regularization (lambda > 0) can rescue methods in regimes
where they struggle, particularly SRVF under high noise.

**DGP selection** (4 DGPs, covering the key method-ranking regimes):

| DGP | Why |
|-----|-----|
| D01 (harmonic+simple+none) | Easy baseline where SRVF dominates — does lambda *hurt* a winning method? |
| D02 (harmonic+complex+none) | Complex warps, no amplitude confound — over-registration regime, lambda might help |
| D09 (harmonic+simple+highrank) | Realistic amplitude — does lambda interact with amplitude variation? |
| D12 (wiggly+complex+highrank) | Hardest DGP — stress test for regularization |

**Lambda grid** (calibrated from pilot, job 5131464):
`c(0, 1e-4, 1e-3, 0.01, 0.05, 0.1, 1, 10)` — 8 values spanning 5 orders of magnitude.

Lambda pilot findings (2 DGPs × 10 lambdas × 3 methods × 3 noise × 10 reps = 1,800 runs):
- **Methods need very different lambda ranges.** fda_default optimal: 1e-4 (D02) to
  0.5–10 (D12). fda_crit1 optimal: 1e-3 (D02) to 0.05–0.1 (D12). SRVF optimal:
  1e-4 (clean) to 10 (noisy D02) — but lambda=0 is best on D12 regardless of noise.
- **SRVF lambda×noise interaction is enormous on D02**: λ*=1e-4 at noise=0, λ*=5e-3
  at noise=0.1, λ*=10 at noise=0.3. On D12 the optimal is λ=0 at all noise levels —
  regularization only hurts.
- **FDA benefits substantially on D12**: fda_crit1 improves 50–60%, fda_default
  improves 45–52%. Smaller improvement on D02 (10–25%).
- **fda_crit1 vs fda_default respond differently**: fda_crit1 optimal is consistently
  higher than fda_default on D02 (1e-3 vs 1e-4), but lower on D12 (0.05–0.1 vs
  0.5–10). The two criteria have qualitatively different lambda profiles.
- **A shared lambda grid is justified**: the 8-value grid covers all methods' optima.
  Method-specific grids are unnecessary since the shared grid includes each method's
  optimal region.

**Rationale for design choices**:
- **+fda_crit1**: Criterion 1 penalizes warp roughness directly, so its interaction
  with external lambda differs qualitatively from fda_default. Pilot confirmed:
  different lambda profiles. [Council R4: unanimous]
- **+noise=0.3**: The key regime where SRVF degrades (derivative amplification) and
  FDA rankings may flip. Pilot confirmed: SRVF optimal lambda jumps 5 orders of
  magnitude from noise=0 to noise=0.3.
- **+1 DGP** (D01): Need a case where the default (lambda=0) is already optimal to
  verify that the lambda grid includes "no regularization" as a viable choice.
- **4 DGPs instead of 3**: Covers easy/hard × with/without amplitude. The extra DGP
  adds modest cost but avoids confounding warp difficulty with amplitude.
- **~33 reps**: Paired design (same seeds across lambda values) makes within-lambda
  comparisons powerful even at 33 reps. The primary signal is the U-shape of the
  MISE-vs-lambda curve, not a point estimate. If ambiguous U-shapes are observed,
  run a targeted 50-rep extension for the 2–3 most interesting cells.

**Key analyses**:
- Lambda × noise interaction: does optimal lambda shift with noise level?
- Per-DGP optimal lambda profiles for SRVF vs FDA (and fda_crit1 vs fda_default)
- Warp MISE vs lambda curves (U-shape confirmed by pilot for FDA; monotone for SRVF on D12)
- fda_crit1 vs fda_default: criterion choice interacts with lambda (pilot-confirmed)
- Does lambda rescue SRVF at noise=0.3? (pilot suggests: only on D02, not D12)

### Study C: Grid resolution sensitivity (~11,250 runs)
5 DGPs x 3 grid sizes (51, 101, 201) x 3 noise (0, 0.1, 0.3)
x 5 methods x severity=1.0 = 225 cells x 50 reps

**Motivation**: Grid-size pilot (job 5131384) revealed that SRVF warp MISE degrades
1.82× at n_grid=201 vs 101, while FDA template MISE *improves* ~40%. This is because
SRVF operates on SRSFs (numerical derivatives) which amplify discretization artifacts
at finer grids. Grid resolution is a confound for method comparison that needs
systematic investigation.

**DGP selection** (5 DGPs, covering all warp types × both templates):

| DGP | Why |
|-----|-----|
| D01 (harmonic+simple+none) | Easy baseline — isolates pure grid effect without other confounds |
| D03 (wiggly+complex+none) | Hard warps + hard template — do complex Matérn warps need finer grids? |
| D04 (harmonic+affine+none) | Affine warps — no numerical integration in warp generation, so grid artifacts come purely from the method |
| D09 (harmonic+simple+highrank) | Realistic amplitude — does highrank amplitude interact with grid resolution? |
| D12 (wiggly+complex+highrank) | Hardest DGP — stress test, maximum grid sensitivity expected |

**Coverage**: Templates: harmonic (D01, D04, D09) + wiggly (D03, D12). Warps: simple
(D01, D09) + complex (D03, D12) + affine (D04). Amplitude: none (D01, D03, D04) +
highrank (D09, D12).

**Design choices**:
- **Grid levels**: 51, 101, 201 — span coarse-to-fine on a log scale (doubling).
  51 tests whether SRVF actually improves at coarser grids; 201 confirms pilot
  degradation with more reps and DGPs.
- **Severity=1.0 only**: Pilot used severity=1.0. A grid × severity interaction is
  plausible in principle (at low severity, small warps → discretization noise may
  dominate), but likely second-order vs grid × method and grid × noise. Adding
  severity would double runs for limited insight. [Council R4: noted as limitation]
- **All 3 noise levels**: Grid × noise interaction is the key question — does finer
  grid amplify noise more for derivative-based methods?
- **50 reps**: Sufficient for method × grid comparisons; paired design via common seeds.

**Key analyses**:
- Method × grid interaction plots (warp MISE, template MISE, elastic distance)
- Per-method degradation ratios (grid 201 / grid 101, grid 51 / grid 101)
- Grid × noise interaction: does noise=0.3 + grid=201 compound for SRVF?
- Wilcoxon signed-rank tests on paired within-DGP comparisons
- Recommendation for default grid resolution in tf_register()

**Grid × penalization cross-check** (~1,800 runs): If Study B finds that lambda
rescues SRVF at high noise, test whether it also compensates for grid artifacts.
Run D02 at grid={101, 201} × lambda={0, 0.01, 0.1} × noise={0, 0.3} for SRVF
only, 50 reps. Can be framed as a Study C extension or a small Study F.
[Council R4: raised by all 3 reviewers as a key cross-study interaction]

### Study D: Oracle template (~6,000 runs)
6 DGPs x 2 template modes (oracle, estimated) x 5 methods
x 2 conditions = 120 cells x 50 reps

**Purpose**: Decompose total registration error into template estimation error vs
warp estimation error. When given the true template, template MISE = 0 and any
remaining warp MISE is purely from the warp estimation algorithm.

**DGP selection** (6 DGPs, balanced template × difficulty):

| DGP | Why |
|-----|-----|
| D01 (harmonic+simple+none) | Easy baseline — template is easy to estimate, expect small oracle benefit |
| D03 (wiggly+complex+none) | Hard template — biggest expected oracle benefit (wiggly is hard to estimate) |
| D09 (harmonic+simple+highrank) | Amplitude confounds template estimation — oracle removes this confound |
| D10 (harmonic+complex+highrank) | Complex warps + amplitude — is template error or warp error dominant? |
| D13 (wiggly+affine+highrank) | Affine warps + wiggly template — different warp mechanism, non-domain-preserving |
| D14 (wiggly+simple+none) | Wiggly template without amplitude — isolates template difficulty effect |

**Two conditions** (easy anchor + hard stress test):

| Condition | noise | severity | Purpose |
|-----------|-------|----------|---------|
| Easy | 0.1 | 0.5 | Anchor: template estimation is relatively easy, oracle benefit should be small |
| Hard | 0.3 | 1.0 | Stress: high noise + large warps stress iterative template estimation → larger differential oracle benefit across methods |

**Rationale for changes vs original design**:
- **+2 DGPs** (D01, D13 or D14): Original selection (D03, D08, D10, D12) overweighted
  wiggly+complex. Need harmonic DGPs to test the method×template interaction (Study A
  showed SRVF loses template MISE specifically on harmonic). Need an easy case to anchor
  the oracle benefit scale.
- **D14 over D08**: D08 (wiggly+complex+rank2) is similar to D12 (wiggly+complex+highrank).
  D14 (wiggly+simple+none) provides a clean wiggly-only case without amplitude or complex
  warp confounds.
- **Two conditions instead of one**: A single easy condition (noise=0.1, severity=0.5)
  compresses oracle benefit differences — template estimation isn't stressed enough to
  reveal which methods are most template-sensitive. The hard condition (noise=0.3,
  severity=1.0) is where template estimation actually fails meaningfully, exposing the
  interesting differential oracle benefit across methods. [Council R4: 2/3 reviewers]

**Key analyses**:
- Oracle benefit = (warp MISE with estimated template) − (warp MISE with oracle template)
- Per-method oracle benefit: which methods are most sensitive to template quality?
- DGP-conditional oracle benefit: where does template estimation matter most?
- Condition-conditional: does oracle benefit grow from easy → hard as expected?
- Template MISE in estimated mode as sanity check (should be ~0 in oracle mode)

### Study E: Outlier contamination (~9,000 runs)
3 DGPs x 2 outlier types (shape, phase) x 3 fractions (10%, 20%, 30%)
x 5 methods x 2 noise (0.1, 0.3) x severity=0.5 = 180 cells x 50 reps

**Purpose**: Test method robustness when a fraction of curves are outliers (wrong shape
or extreme warping). Relevant for practical use where data cleaning is imperfect.

**DGP selection** (3 DGPs, spanning difficulty):

| DGP | Why |
|-----|-----|
| D02 (harmonic+complex+none) | Moderate difficulty, no amplitude — cleanest outlier signal. Complex warps mean phase outliers compete with legitimate warp variation. |
| D09 (harmonic+simple+highrank) | Realistic amplitude — do shape outliers get confused with legitimate amplitude variation? Key practical question. |
| D12 (wiggly+complex+highrank) | Hardest DGP — stress test. Methods already struggle here; how much worse do outliers make it? |

**Rationale for changes vs original design**:
- **D02 replaces D01**: D01 (harmonic+simple+none) is too easy — all methods handle it
  well, so outlier effects would be trivially small. D02 adds complex warps, making
  phase outliers harder to distinguish from legitimate variation.
- **D03 dropped**: D03 (wiggly+complex+none) overlaps with D12 on template+warp
  difficulty. D09 adds the amplitude dimension instead, which is the more interesting
  confound for outlier detection.
- **D12 replaces D03 as "hard" case**: Hardest DGP tests whether outliers break methods
  that are already near their limits.
- **+noise=0.3**: At high noise, legitimate observations already look noisy, making
  shape outliers harder to distinguish from noisy inliers. Breakdown curves at noise=0.1
  alone would overestimate robustness. [Council R4: Claude, Codex]
- **Contamination fractions raised to 10%/20%/30%**: At n=50 curves, 5% contamination
  means only 2–3 outliers — discrete rounding artifacts dominate. 10% (5 curves) is
  the practical minimum for stable results. [Council R4: Codex]

**Outlier types**:
- **Shape outliers**: Random sinusoids scaled to ~2x template range. These curves have
  the wrong shape entirely — a well-functioning method should downweight or exclude them.
- **Phase outliers**: Extreme warps (severity=3.0, 6× default). These curves have the
  right shape but extreme time warping — tests whether iterative methods (FDA, SRVF)
  get pulled toward extreme warps.

**Implementation note**: Verify that `contaminate_data()` uses deterministic outlier
selection given the same base seed, so that paired comparison across contamination
fractions is valid (same curves are selected as outliers at 10% and 20%, with 20%
adding more). [Council R4: Claude]

**Key analyses**:
- Breakdown curves: at what contamination fraction does each method fail?
- Shape vs phase outlier sensitivity: which type is more damaging per method?
- Noise × contamination interaction: is robustness worse at high noise?
- Template bias under contamination (outliers pull template estimate)
- Warp MISE degradation curves by contamination fraction
- Report failure rates explicitly — conditional-on-success comparisons are biased
  if failure probability changes with contamination [Council R4: Codex]

**Total: ~45,850 runs (Studies B–E + grid×lambda cross-check) + 45,000 (Study A) ≈ 91,000 runs**

## 8. Code Structure

```
sim-registration/            Standalone repo (~/fda/sim-registration/)
  BENCHMARK-PLAN.md          This file
  gait_knee_fpcs.rds         Cached FPCs from registered gait knee angles
  sim-dgp.R          ~500 lines   Templates, warps, amplitude, generate_data()
  sim-methods.R       ~110 lines   5 method wrappers returning tf_registration
  sim-metrics.R       ~300 lines   All metrics using tf API (incl. SRSF elastic dist)
  sim-config.R        ~220 lines   All study designs (A-E)
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
- `contaminate_data()` for Study E (shape and phase outliers)
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

### Council Review R3 (done)
Council of 4 AI agents (2× Claude, Codex, Gemini) reviewed Study A report. See
TODO section for full list. Key fixes applied: Type II ANOVA, elastic distance
metric (fdasrvf::elastic.distance), interaction plots, alignment error section.

### Council Review R4 (done)
Council of 3 AI agents (Claude, Codex, Gemini) reviewed Studies B–E design.
Applied changes:
1. Added fda_crit1 to Study B (unanimous) — different criterion, different lambda response
2. Added hard condition (noise=0.3, severity=1.0) to Study D (2/3 reviewers)
3. Added noise=0.3 to Study E + raised min contamination from 5% to 10% (Codex: discrete artifacts at n=50)
4. Added grid × lambda cross-check to Study C (~1,800 runs, all 3 reviewers)
5. Noted severity=1.0-only limitation in Study C, lambda-grid pilot for Study B
6. Added metric hierarchy question (template MISE vs elastic distance) to open questions
Acknowledged but deferred:
- 9 missing DGP cells (low priority)
- v1→v2 design differences may invalidate some B–E rationale (finalize after v2)
- Failure rate integration into B–E analyses (implement during analysis)

### Phase 7: Study A production run v2 (submitted)
- Updated design: 15 DGPs x 2 severity x 3 noise x 5 methods = 450 cells x 100 reps = 45K runs
- Changes from v1: +2 DGPs (D14, D15), +1 noise level (0.3), doubled affine severity,
  calibrated rank2 amplitude, Fisher-Rao elastic distance metric, forward warp comparison
- tf updated to 639d19d ("Fix CC warp monotonicity") — FDA results may differ from v1
- LRZ job 5131463 (2026-03-10), 32 cores, serial_std partition, 24h time limit
- Batched by (DGP, method): 75 batches, each ~600 runs
- Results saved as `results_D{nn}_{method}_A.rds`

### Phase 8: Study A Analysis (in progress)
- sim-analyze.R: data loading, DGP factor annotation, aggregation helpers
- sim-report.qmd: Quarto report implementing analysis plan below
- All findings sections use data-driven `results: asis` code chunks

### Phases 9+: Sub-studies B–E (pending)
- DGPs for Studies B–E selected from Study A diagnostics
- Studies B–E runs, analysis, integrated report

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

- ~~Lambda grid for Study B~~ **Resolved**: Pilot (job 5131464) calibrated grid to
  `c(0, 1e-4, 1e-3, 0.01, 0.05, 0.1, 1, 10)`. Shared grid covers all methods' optima.
  Methods have qualitatively different lambda profiles (documented in Section 7).
- Study B–E DGP selections are tentative (documented in Section 7). Finalize after
  Study A re-run based on actual method rankings, interaction patterns, and failure rates.
- Metric hierarchy: template MISE and elastic distance are strongly correlated (ρ=0.85).
  Clarify which is primary and which is confirmatory for Studies B–E to avoid redundant
  testing. [Council R4: Codex]

## TODO (from Council Review R3)

Council of 4 AI agents reviewed the Study A report and code. Key issues to address:

### Before Study A re-run
- [x] **Fix elastic distance metric**: `compute_template_elastic_dist()` in `sim-metrics.R`
  now uses `fdasrvf::elastic.distance()` for proper Fisher-Rao amplitude distance.
  Verified in grid-size pilot (job 5131384): values range 0.2–2.6, ρ=0.85 with
  template MISE, sensible DGP/method patterns.
- [ ] **Re-run Study A** to populate `template_elastic_dist` column (currently all NA).
  Template quality conclusions are provisional until this is done.
  Submitted as job 5131463 on LRZ (2026-03-10). Also includes tf "Fix CC warp
  monotonicity" (639d19d) — FDA method results may differ from v1.

### Report fixes (applied)
- [x] Switch ANOVA from Type I SS (`aov()`) to Type II SS (`car::Anova()`) — non-full-
  factorial design makes Type I η² order-dependent
- [x] Add noise-conditional heatmaps (noise=0 vs noise=0.3) to overview section
- [x] Replace MC SE of mean with median [Q25, Q75] in summary tables
- [x] Add alignment_error section (primary metric, was completely missing from report)
- [x] Add method × template and method × noise interaction plots
- [x] Qualify "best general-purpose method" claims with template/noise interaction caveats
- [x] Fix hardcoded "rankings largely stable" prose → data-driven rank concordance
- [x] Fix Q5 "invariant to reparameterization" claim → "SRSF L2 distance (upper bound)"
- [x] Explain FDA noise ratio < 1 (MISE decreasing with noise) — regularization effect

### Study B–E design (tentative, finalize after Study A re-run)
DGP selections, rationale, and council feedback documented in Section 7.
Pending implementation tasks:
- [ ] **Study B (penalization)**: Update `study_b_design()` in `sim-config.R` — add
  fda_crit1, noise=0.3, DGPs D01/D02/D09/D12, ~33 reps, 8 lambdas
  `c(0, 1e-4, 1e-3, 0.01, 0.05, 0.1, 1, 10)`. Lambda grid calibrated by pilot
  (job 5131464).
- [ ] **Study C (grid sensitivity)**: Implement `study_c_design()` in `sim-config.R` —
  DGPs D01/D03/D04/D09/D12, grid sizes 51/101/201, severity=1.0 only.
  Also implement grid × lambda cross-check (~1,800 runs).
- [ ] **Study D (oracle template)**: Update `study_c_design()` → `study_d_design()` in
  `sim-config.R` — DGPs D01/D03/D09/D10/D13/D14, two conditions (easy + hard).
- [ ] **Study E (outlier contamination)**: Update `study_d_design()` → `study_e_design()`
  in `sim-config.R` — DGPs D02/D09/D12, fractions 10%/20%/30%, 2 noise levels.
  Verify `contaminate_data()` seed determinism for paired fraction comparisons.

### Design considerations
- [ ] **9 missing DGP cells**: 15/24 is fine for ladders but limits global claims. Affine
  is particularly sparse (only none and highrank amplitude, no rank1/rank2 → no amplitude
  ladder for affine). Consider running the 9 missing cells at 1 condition × 50 reps
  (2,250 runs) as a completeness check.
- [ ] **Paired analysis**: Seeds are `(dgp, rep)`-based → same data across methods/noise/
  severity. Current ANOVA and summaries treat runs as independent. Paired within-replicate
  differences would be more powerful, especially for close races (FDA default vs crit 1).
- [ ] 100 reps sufficient for large separations but tight for close pooled calls.
  The top-3 methods are separated by ~1e-4 on harmonic DGPs — within MC noise without
  paired analysis. Consider 200 reps for the 3–6 most interesting DGPs.

## Reproducibility

- Deterministic seeding: `digest::digest2int(paste(dgp, rep))` per task
- Session info saved with results
- Incremental saves per DGP group (crash-safe)
- All code in `~/fda/sim-registration/`, data in `results/` subdirectory
- `.gitignore` in tf repo blocks `attic/sim-registration/` from being pushed
