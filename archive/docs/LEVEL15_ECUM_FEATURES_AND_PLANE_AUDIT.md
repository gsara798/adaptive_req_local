# Level 15: Ecum Features and In-Plane Wave Audit

## Purpose

Level 15 tests two questions:

1. Does Fibonacci cone sampling include at least one wave in the imaging
   plane?
2. Do features derived from the REQ cumulative energy curve improve q
   prediction?
3. After predicting q, how does the error propagate into shear wave speed
   (`SWS`, `c_s`)?

## In-Plane Wave Audit

For an imaging plane `y = 0`, a wave is exactly in-plane when its direction
has `uy = 0`.

With the current `cone + fibonacci` sampling and `Nwaves = 2000`, only the
zero-aperture case has exact in-plane waves. For all nonzero apertures, the
nearest wave is very close to the plane but not exactly in it.

The new option `ForceInPlaneWave` replaces the closest sampled direction with
an exact in-plane direction. This is available in the simulator and defaults
to `false` to avoid silently changing old experiments.

Relevant files:

```text
src/+adaptive_req/+simulate/summarize_wave_direction_plane_coverage.m
src/+adaptive_req/+simulate/estimate_fibonacci_cone_plane_coverage.m
src/+adaptive_req/+simulate/simulate_rswe_plane.m
```

## Ecum Shape Features

REQ computes a cumulative radial energy function `Ecum(k)`. Level 15 extracts
shape descriptors from this curve without using `cs_true`, `q_true`, or a
fixed `cs_guess`.

The features include:

- normalized inverse-CDF locations: `k10`, `k25`, `k75`, `k90` relative to
  `k50`;
- cumulative-curve widths: `ecum_width_50_rel`, `ecum_width_80_rel`;
- tail asymmetry;
- transition asymmetry and lower/upper transition width ratio;
- normalized area under `Ecum(k)`;
- entropy and peak concentration of CDF increments;
- increment Gini concentration;
- maximum normalized slope, peak-to-mean slope, and slope spread.

Relevant file:

```text
src/+adaptive_req/+quantile/extract_ecum_shape_features.m
```

## Link to REQ Theory

The REQ method starts from the 2D spectrum, builds a radial mean spectrum, and
then computes the cumulative distribution

```text
Ecum(k) = integral_0^k Srad(k') dk' / integral_0^kmax Srad(k') dk'
```

The predicted quantile selects `kq` through `Ecum(kq) = q`, and the final speed
is

```text
c_s = 2*pi*f0/kq
```

The useful point for ML is that the shape of `Ecum(k)` changes with field type.
In the pasted REQ theory, directional fields behave close to a step around the
true wavenumber, while diffuse/projected fields spread energy over lower
projected radial wavenumbers and shift the theoretical quantile. This means
that cumulative-curve width, slope, entropy, and increment concentration are
physically meaningful predictors of q, not arbitrary statistical features.

## q Prediction Result

Using Test 08 and adding Ecum features from the stored lightweight
`req_mapping`, the best operational model improved:

| Model | Test RMSE_q |
|---|---:|
| Model C baseline, bagged trees | 0.05075 |
| Model H Ecum shape, bagged trees | 0.04191 |

The improvement also appears for the difficult `M = 2` case:

| Model | M=2 RMSE_q |
|---|---:|
| Model C baseline | 0.06023 |
| Model H Ecum shape | 0.05160 |

This supports using REQ-curve shape features as operational predictors.

## SWS Result

Level 15 now converts every `q_pred` into SWS using the local lightweight
`req_mapping`, so the same analysis contains both quantile and final-speed
performance.

For the test split:

| Model | SWS MAPE | SWS RMSE | Bias | R2 |
|---|---:|---:|---:|---:|
| Model H Ecum shape, bagged trees | 3.05% | 5.16% | 0.92% | 0.958 |
| Model I angular + Ecum, bagged trees | 3.10% | 5.18% | 0.91% | 0.958 |
| Model J Ecum by M, bagged trees | 3.32% | 5.73% | 1.30% | 0.955 |
| Model C baseline, bagged trees | 4.43% | 6.84% | 1.50% | 0.933 |

By `M`, the best Level 15 model is still weakest for `M = 2`:

| REQ M | SWS MAPE | SWS RMSE | p95 APE | Max APE |
|---:|---:|---:|---:|---:|
| 2 | 4.31% | 7.11% | 15.25% | 32.68% |
| 3 | 2.47% | 3.98% | 8.80% | 36.06% |
| 4 | 1.91% | 3.18% | 7.41% | 20.82% |

The large maximum errors are patch-level outliers. The typical error is much
lower: median APE is 1.58% for the best model.

The per-M model is now trained with the same global train/test split. It does
not outperform the global Ecum model. This suggests that `REQ_M` is already
being used effectively by the tree models, including nonlinear interactions
between `M` and the Ecum shape features.

## High-Error Rate

For the best model (`ModelH_ecum_shape`, bagged trees), 30 of 2400 test
predictions had absolute SWS error above 20%, i.e. 1.25%.

By `M`:

| REQ M | n | n > 20% | percent |
|---:|---:|---:|---:|
| 2 | 900 | 23 | 2.56% |
| 3 | 900 | 6 | 0.67% |
| 4 | 600 | 1 | 0.17% |

By true SWS:

| true c_s | n | n > 20% | percent |
|---:|---:|---:|---:|
| 2 m/s | 900 | 6 | 0.67% |
| 3 m/s | 900 | 18 | 2.00% |
| 4 m/s | 600 | 6 | 1.00% |

## Failure Interpretation

The dominant failure mechanism is q underestimation at high aperture and low
`M`. In the worst outliers, `q_pred` is too small relative to `q_true`; when
that smaller q is inverted through `Ecum(k)`, it selects a lower `kq`, and
because `c_s = 2*pi*f0/kq`, the predicted SWS becomes too high.

For the best model, the strongest associations with absolute SWS error were:

| Feature | Spearman rho with abs SWS error |
|---|---:|
| `q_abs_error` | 0.923 |
| `q_true` | 0.371 |
| `Omega_sr` | 0.362 |
| `ecum_width_80_rel` | 0.351 |
| `q_pred` | 0.330 |
| `ecum_width_50_rel` | 0.303 |

This suggests the main remaining risk is not random noise. The model still
struggles when the cumulative curve is broad, the aperture is large, and the
REQ window is small. Those are exactly the settings where projection/diffuse
broadening makes `Ecum(k)` less step-like and the q-to-SWS inversion more
sensitive.

## New Outputs

Level 15 now writes:

```text
level15_sws_predictions.csv
level15_sws_metrics_by_model.csv
level15_sws_error_by_M.csv
level15_sws_error_by_aperture.csv
level15_sws_error_by_condition.csv
level15_failure_correlations.csv
level15_failure_outliers.csv
level15_high_error_rate_by_model.csv
level15_high_error_rate_by_model_M.csv
level15_high_error_rate_by_model_cs_bg.csv
level15_sws_rmse_by_model.png
level15_sws_pred_vs_true_best_models.png
level15_sws_error_vs_true_speed_best_models.png
level15_sws_error_boxplots.png
level15_best_model_failure_diagnostics.png
level15_failure_correlation_summary.png
level15_high_error_rates.png
```

## Run

```matlab
run('experiments/analysis/analyze_test_08_level15_ecum_features_and_plane_audit.m')
```

Outputs are written to:

```text
outputs/test_08_advanced_angular_features/<run>/analysis/level_15_ecum_features_and_plane_audit
```
