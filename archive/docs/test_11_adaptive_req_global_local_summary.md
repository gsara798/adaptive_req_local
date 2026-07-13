# Test 11: Adaptive REQ With Local and Global Features

## 1. Objective

Test 11 evaluated **adaptive-q** models for REQ. The central idea was to learn the reference quantile `q_theory` from spectral information in the simulated wavefield, then convert each predicted quantile `q_pred` into SWS using the local `req_mapping` from each analysis window.

The dataset included:

- local features computed from each window;
- global features computed from the full-field power spectrum;
- the `LocalOnly`, `GlobalOnly`, and `HybridLocalGlobal` models;
- prediction of `q_theory`;
- downstream conversion `q_pred -> req_mapping -> k_pred -> cs_pred`.

## 2. Structure

Main files:

- `configs/test_11_global_req_features.m`
- `experiments/run_test_11_global_req_features.m`
- `experiments/analysis/analyze_test_11_level01_global_vs_local.m`
- `experiments/analysis/analyze_test_11_level04_grouped_generalization.m`
- `experiments/analysis/analyze_test_11_level05_outlier_diagnostics.m`

## 3. Evaluated Models

- `LocalOnly`: uses only local window features and predicts one `q` per window.
- `GlobalOnly`: uses only global features from the full field. It tests how much the global spectral context can explain, but it loses local spatial variation.
- `HybridLocalGlobal`: uses local + global features. It predicts one `q` per window while being informed by global context.
- Diagnostic/oracle models, when present, are not operational. In particular, any model using `q_global_theory`, `q_local_minus_global`, `M_eff_true_diag`, `cs_true`, residuals, or errors should be interpreted as diagnostic only.

## 4. Main Level 01 Conclusions

Level 01 compared local, global, and hybrid models using less strict splits than the later grouped generalization tests.

Conclusions:

- `GlobalOnly` is not sufficient as the main model for local quantile prediction. It can capture broad field-level trends, but it does not describe window-to-window variation well enough.
- `LocalOnly` works well and supports the hypothesis that local spectral features contain physical information about the REQ quantile.
- `HybridLocalGlobal` appears to be the most robust model overall because it combines the local information that defines the SWS map with global context from the full field.

## 5. Main Level 04 Conclusions

Level 04 evaluated clean generalization using grouped splits:

- `leave_one_frequency`: holds out an entire `SIM_f0` value.
- `leave_one_M`: holds out an entire `REQ_M` value.
- `leave_one_wave_model`: holds out either `planewave` or `spherical`.
- `leave_one_aperture`: holds out an entire `step_idx`.

These splits are stricter than a random window-level split because they avoid leakage between windows from the same condition or realization. This is especially important because global features are repeated across many local windows.

Conclusions:

- The model generalizes reasonably well to unseen frequencies.
- Generalization to unseen `REQ_M` is harder, especially for `M = 2`.
- It does not appear necessary to train separate models for `planewave` and `spherical`; the model seems to learn spectral geometry rather than details of the wave generator.
- `HybridLocalGlobal` is the most stable operational candidate.

In the Level 04 grouped-split metrics, `leave_one_M` was the hardest case. For `bagged_trees`, the approximate MAPE values in `leave_one_M` were:

- `LocalOnly`: 6.15%
- `GlobalOnly`: 5.78%
- `HybridLocalGlobal`: 4.82%

## 6. Main Level 05 Conclusions

Level 05 diagnosed Level 04 outliers and added occupancy/support analysis.

There are two useful ways to summarize the metrics:

1. Broad summary of predictions stored in Level 05:
   - `HybridLocalGlobal`: MAPE approximately 1.83%, high-error >20% approximately 0.31%.
   - `LocalOnly`: MAPE approximately 2.23%, high-error >20% approximately 0.87%.
   - `GlobalOnly`: MAPE approximately 3.99%, high-error >20% approximately 1.27%.

2. Support-aware summary using only `bagged_trees` test rows:
   - `HybridLocalGlobal`: high-error >20% approximately 1.15%.
   - `LocalOnly`: high-error >20% approximately 3.04%.
   - `GlobalOnly`: high-error >20% approximately 2.66%.

The difference comes from the subset being summarized. Exact values should be checked in the final Level 05 tables, especially:

- `level11_level05_outlier_rate_by_model.csv`
- `level11_level05_support_by_M_eff_true_diag.csv`

Main conclusions:

- Large outliers are globally rare.
- Outliers are not random: they concentrate in small effective windows.
- For `LocalOnly` and `HybridLocalGlobal`, the main problematic bin is `M_eff_true_diag_bin = [1.5, 2)`.
- That bin contains a small fraction of the test set, around 11.1%, but has a much higher high-error rate.
- This indicates a small but problematic bin, not only a large bin accumulating many errors because of occupancy.
- For `GlobalOnly`, outliers appear more related to occupancy plus a moderate error rate: a large fraction falls in `M_eff_true_diag_bin = [3, 3.5)`, which contains around 33.3% of the rows.

## 7. Physical Interpretation

`REQ_M` is the nominal window size in guessed wavelengths, but the quantity that matters physically is the effective window size relative to the true wavelength.

A useful approximation is:

```matlab
M_eff_true_diag ~= REQ_M * REQ_cs_guess / SIM_cs_bg
```

Interpretation:

- `M_eff_true_diag` is diagnostic/oracle because it uses the true simulated shear wave speed. It should not be used as an operational predictor for real data.
- `M_eff_guess` is closer to an operational/nominal quantity because it is tied to the window size the algorithm believes it is using.
- The main failure mode appears when the true effective window size is too small, especially around 1.5-2 wavelengths.
- This suggests insufficient effective spatial support for robust local spectrum and quantile estimation.

## 8. Current Limitations

- A fixed `cs_guess` artificially couples `REQ_M`, `cs_true`, and `M_eff_true_diag`.
- Test 11 does not fully decouple the effects of `REQ_M` and `cs_guess`.
- Noise, attenuation, and depth-dependent SNR should not be mixed into this stage yet. The effective-window problem should be isolated first.
- The current simulations are controlled. k-Wave, bilayer, and inclusion tests should be used later as more realistic validation, not as the first place to diagnose this mechanism.

## 9. Motivation for Test 12

Test 12 should decouple:

- `REQ_M`
- `REQ_cs_guess`
- `SIM_cs_bg`
- `SIM_f0`
- `M_eff_guess`
- `M_eff_true_diag`

The central question will be:

> Does the model fail because of nominal `REQ_M`, because of `cs_guess`, or because of the true effective window size?

To answer this, Test 12 explicitly sweeps `REQ.cs_guess` in addition to `REQ.M`, `SIM.cs_bg`, and `SIM.f0`, without adding noise or attenuation.

## 10. Next Steps

1. Run Test 12 without noise or attenuation.
2. Analyze generalization to unseen `REQ_cs_guess`.
3. Analyze error versus `M_eff_guess` and `M_eff_true_diag`.
4. Create Test 13 for SNR.
5. Create Test 14 for depth-dependent SNR and attenuation.
6. Then return to bilayer/k-Wave as more realistic validation.
