# Test 15: Theory-Informed Direct q Model

## Purpose

Test 15 evaluates a new direct adaptive-REQ model for predicting the REQ quantile:

```matlab
q_pred = F(spectral_features, theory_q_candidates, optional_user_field_guess)
```

This is intentionally **not** a residual model. A residual version should be treated as a separate future test.

The key question is whether adding discrete REQ theory candidates improves the current data-driven hybrid model. The current hybrid model learns:

```matlab
spectral_features -> q_theory
```

Test 15 instead tests:

```matlab
spectral_features + q_dir2D + q_diffuse2D + q_projected3D + user_guess_prior -> q_theory
```

## Relationship to Test 12

Test 15 uses the already generated Test 12 dataset:

```matlab
test_12_cs_guess_window_sweep
```

It does not regenerate Test 12 data, and it does not modify Test 11, Test 12, or the previously trained LocalOnly/GlobalOnly/HybridLocalGlobal models.

Although the training dataset comes from Test 12, this analysis is conceptually a new test. For compatibility with previous scripts, the primary analysis also exists under the Test 12 run folder, but the outputs are mirrored into a standalone Test 15 location:

```matlab
outputs/test_15_theory_informed_direct_q/analysis/
```

## Files Added

Main analysis:

```matlab
experiments/analysis/analyze_test_15_theory_informed_direct_q.m
```

Feature builders:

```matlab
src/+adaptive_req/+analysis/build_theory_q_features.m
src/+adaptive_req/+analysis/build_user_field_guess_features.m
```

Model registry:

```matlab
src/+adaptive_req/+analysis/register_trained_model.m
outputs/model_registry/model_manifest.csv
```

Supporting fixes:

```matlab
src/+adaptive_req/+analysis/train_q_model_from_predictors.m
src/+adaptive_req/+analysis/Test12Analysis.m
```

The training helper was updated so that mixed numeric and categorical predictors can be encoded together. This was required for `user_field_guess`.

## Theory-q Candidate Features

For each unique REQ geometry, Test 15 computes discrete theoretical quantiles using:

```matlab
adaptive_req.theory.q_theory_REQ_discrete_shearUZ
```

The added theory features are:

```matlab
q_theory_dir2D
q_theory_diffuse2D
q_theory_projected3D
q_theory_mean_dir2D_projected3D
q_theory_mean_all
```

The theory calculation is cached by unique geometry so the same quantile candidates are not recomputed for every local window.

## User Field Guess Features

Test 15 also simulates a user-provided qualitative field prior:

```matlab
user_field_guess
```

with possible values:

```matlab
directional_like
partially_diffuse
diffuse_like
unknown
```

The corresponding prior quantile is:

```matlab
q_user_guess_prior
```

Mapping:

| user_field_guess | q_user_guess_prior |
|---|---|
| `directional_like` | `q_theory_dir2D` |
| `partially_diffuse` | `mean(q_theory_dir2D, q_theory_projected3D)` |
| `diffuse_like` | `q_theory_projected3D` |
| `unknown` | `mean(q_theory_dir2D, q_theory_diffuse2D, q_theory_projected3D)` |

The user guess is categorical. No continuous aperture variable is used as an operational predictor.

## Models Compared

### CurrentHybridBaseline

Existing deployed model:

```matlab
HybridLocalGlobal | NoCsGuess | bagged_trees
```

This model is loaded from the existing Test 12 deployment folder and is not retrained.

### TheoryCandidatesDirect

New operational direct model.

Predictors:

```matlab
Hybrid spectral features
q_theory_dir2D
q_theory_diffuse2D
q_theory_projected3D
q_theory_mean_dir2D_projected3D
q_theory_mean_all
```

### TheoryCandidatesPlusUserGuessDirect

New operational direct model.

Predictors:

```matlab
Hybrid spectral features
theory-q candidates
q_user_guess_prior
user_field_guess
```

### UserGuessPriorOnly

No-ML baseline:

```matlab
q_pred = q_user_guess_prior
```

This tests how much the qualitative user prior can do by itself.

### TheoryBestCandidateOracle

Diagnostic-only oracle. For each row, it picks the theory candidate closest to the true target `q_theory`.

This is not operational and should not be used as a deployable model.

## Anti-Leakage Rules

Operational models do not use:

```matlab
q_theory
q_true
q_global_theory
q_local_minus_global
abs_q_local_minus_global
M_eff_true_diag
cs_true
cs_pred
sws_error
abs_sws_error
sws_error_pct
abs_sws_error_pct
residual
abs_error
aperture_weight
solid_angle_weight
true_aperture_weight
```

`M_eff_true_diag` is kept only for diagnostic grouping after prediction.

## Training Setup

Dataset:

```matlab
test_12_cs_guess_window_sweep
```

Target:

```matlab
q_theory
```

Split:

```matlab
grouped condition split
```

Main model type:

```matlab
bagged_trees
```

The trained direct models are registered in:

```matlab
outputs/model_registry/
```

## Outputs

Standalone Test 15 output:

```matlab
outputs/test_15_theory_informed_direct_q/analysis/
```

Compatibility output:

```matlab
outputs/test_12_cs_guess_window_sweep/<run>/analysis/level_15_theory_informed_direct_q/
```

Tables:

```matlab
level15_predictions.csv
level15_q_metrics.csv
level15_sws_metrics.csv
level15_sws_metrics_by_user_guess.csv
level15_sws_metrics_by_M_eff.csv
level15_model_comparison.csv
level15_theory_q_cache.csv
```

Model file:

```matlab
level15_theory_informed_direct_q_models.mat
```

Figures:

```matlab
level15_mape_by_model.png
level15_high_error_by_model.png
level15_q_true_vs_pred.png
level15_delta_mape_vs_CurrentHybridBaseline.png
level15_mape_by_user_field_guess.png
level15_mape_by_M_eff_true_diag.png
level15_user_guess_prior_vs_ml_direct.png
```

## Main Results From Current Run

The current run produced:

| Model | MAPE (%) | High-error >20% (%) |
|---|---:|---:|
| CurrentHybridBaseline | 1.87 | 0.059 |
| TheoryCandidatesPlusUserGuessDirect | 2.82 | 0.204 |
| TheoryCandidatesDirect | 2.91 | 0.237 |
| UserGuessPriorOnly | 20.61 | 36.55 |

Interpretation:

1. The current hybrid baseline remains the best operational model in this run.
2. Adding theory-q candidates directly did not improve over the current hybrid baseline.
3. Adding `user_field_guess` improved the theory-informed direct model slightly relative to using theory candidates alone.
4. The user prior alone is not sufficient.
5. The theory-informed models still perform much better than the user-prior-only baseline, meaning the spectral features are doing real work.

## Scientific Interpretation

The theory candidates are physically meaningful, but this direct formulation does not yet beat the learned hybrid baseline. This suggests that simply adding candidate theoretical quantiles as extra predictors is not enough.

Possible explanations:

1. The hybrid spectral features already encode much of the information contained in the theory candidates.
2. The relationship between theoretical field classes and simulated finite-window spectra may be more naturally modeled as a correction or residual.
3. The user guess is too coarse to strongly improve performance when the spectral features are already informative.

This motivates a future Test 16:

```matlab
q_pred = q_theory_candidate_or_prior + residual_model(features)
```

That residual formulation may be a better way to use theory because it forces the model to learn corrections around a physically meaningful baseline.

## Current Takeaway

Test 15 is a useful negative/diagnostic result:

```text
Theory-informed direct q prediction is feasible, but it does not outperform the current HybridLocalGlobal baseline.
```

The next theory-informed step should be a residual model rather than another direct model.
