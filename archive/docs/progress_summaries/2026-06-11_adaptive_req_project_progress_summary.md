# Adaptive REQ Project Progress Summary

**Date:** 2026-06-11  
**Project:** `adaptive_req_local`  
**Purpose:** Track what has been implemented, tested, and learned so far while building a reproducible ML-assisted adaptive-quantile REQ workflow.

## 1. Big-Picture Goal

The project aims to make REQ more adaptive and reproducible by using machine learning to predict the REQ quantile `q` from spectral features of simulated shear-wave fields.

The intended operational workflow is:

```text
wavefield -> local/global spectral features -> q_pred -> REQ mapping -> k_pred -> SWS estimate
```

The main scientific question is:

> Can the REQ quantile be predicted automatically from measurable spectral features, so that SWS estimation works across directional, diffuse, and mixed wavefields?

## 2. Repository Organization Work

The project has been progressively organized into:

- `configs/`: experiment configuration files.
- `experiments/`: dataset generation scripts.
- `experiments/analysis/`: analysis levels for each test.
- `src/+adaptive_req/`: reusable project functions.
- `docs/`: project notes, summaries, interpretation, and progress tracking.
- `outputs/`: generated datasets, figures, tables, and trained models.

Important reproducibility choices:

- experiment configurations are separated from analysis scripts;
- generated outputs are timestamped under `outputs/`;
- heavy `req_curve` storage is avoided where possible;
- compact mappings such as `req_mapping` and `global_req_mapping` are kept instead;
- analysis scripts read saved datasets instead of regenerating data.

## 3. Early ML and SWS Pipeline Work

Earlier analysis levels established the base workflow:

- extraction of spectral features from the power spectrum;
- training models to predict `q_theory`;
- converting `q_pred` into SWS using the REQ mapping;
- evaluating error as SWS MAPE, RMSE, bias, p95 absolute error, max error, and high-error rates.

Additional features were added over time:

- angular entropy and angular concentration features;
- cumulative-energy curve features from `Ecum`;
- slope/width/asymmetry features from `Ecum`;
- Srad proxy features describing radial spectrum shape;
- residual-correction ideas;
- feature-audit and outlier-diagnostic plots.

Key intermediate result:

- Low effective window support, especially small `M`, repeatedly appeared as a major source of large SWS errors.

## 4. Level 13 to Level 17 Work

Several analysis layers were created before Test 11 to understand SWS performance more deeply.

Relevant docs:

- `docs/LEVEL13_SWS_METRICS.md`
- `docs/LEVEL15_ECUM_FEATURES_AND_PLANE_AUDIT.md`
- `docs/LEVEL16_ERROR_LEARNING_FEATURE_AUDIT.md`
- `docs/LEVEL17_MODEL_COMPARISON_AND_NEXT_TESTS.md`

Main takeaways:

- SWS error must be evaluated after converting `q_pred` through the local REQ mapping.
- Boxplots can reveal rare but important high-error outliers that are hidden by mean metrics.
- `Ecum` features are physically meaningful because the cumulative energy curve changes shape with spectral diffuseness.
- Some added features help, but feature usefulness must be verified by ablation and grouped generalization, not only random splits.
- A residual-corrected model was promising in earlier tests, but later work shifted toward cleaner physical/generalization diagnostics.

## 5. Test 11: Local vs Global REQ Features

### Goal

Test 11 was created to compare:

- local features from each REQ window;
- global features from the full-field power spectrum;
- hybrid local + global features.

The key models were:

- `LocalOnly`: local window features only.
- `GlobalOnly`: full-field global features only.
- `HybridLocalGlobal`: local + global features.

The prediction target remained:

```matlab
q_theory
```

The operational SWS conversion used:

```text
q_pred -> req_mapping -> k_pred -> cs_pred
```

### Main Files

- `configs/test_11_global_req_features.m`
- `experiments/run_test_11_global_req_features.m`
- `experiments/analysis/analyze_test_11_level01_global_vs_local.m`
- `experiments/analysis/analyze_test_11_level04_grouped_generalization.m`
- `experiments/analysis/analyze_test_11_level05_outlier_diagnostics.m`
- `docs/test_11_adaptive_req_global_local_summary.md`

### Level 01 Conclusions

Level 01 showed that:

- `GlobalOnly` is not enough as the main local quantile model;
- `LocalOnly` works well, supporting the idea that local spectral features carry physical information about `q`;
- `HybridLocalGlobal` is generally the strongest operational candidate because it keeps local information while adding global context.

### Level 04 Conclusions

Level 04 implemented strict grouped generalization:

- leave-one-frequency-out;
- leave-one-M-out;
- leave-one-wave-model-out;
- leave-one-aperture-out.

This was important because random window-level splits can leak information between windows from the same condition or realization, especially when global features are repeated across many rows.

Main conclusion:

- grouped generalization is reasonable overall;
- leave-one-`REQ_M` is harder, especially for low `M`;
- `HybridLocalGlobal` is the most stable operational candidate.

### Level 05 Outlier Diagnostics

Level 05 diagnosed high-error outliers and added support/occupancy analysis.

Main conclusion:

- large outliers are globally rare;
- outliers are not random;
- for `LocalOnly` and `HybridLocalGlobal`, the main problematic bin is:

```text
M_eff_true_diag_bin = [1.5, 2)
```

This bin is relatively small but has a much higher high-error rate, suggesting a true physical support problem rather than only a dataset-occupancy artifact.

Important distinction:

- `M_eff_true_diag` is diagnostic/oracle because it uses true simulated SWS;
- `M_eff_guess` is closer to operational because it represents the effective window size the algorithm believes it is using.

## 6. Bilayer / Heterogeneous Transfer Work

A bilayer transfer analysis was explored to test whether adaptive-q models trained on homogeneous simulations could be used in heterogeneous fields.

The intended comparison included:

- a theory-based diffuse-3D quantile baseline;
- a global-q model;
- a local model;
- a hybrid local-global model.

Several issues were identified and addressed:

- the correct REQ estimator should be reused rather than reimplementing a slow or inconsistent version;
- predicted SWS maps should use the same style and assumptions as the reference REQ workflow;
- plots need valid-region zooming, smaller titles, correct subscript formatting, and reusable figure helpers;
- ROI statistics should be extracted from soft and hard regions away from boundaries.

This work motivated later reusable plotting and configuration cleanup.

## 7. Theory Integration

REQ theory code was brought into the project instead of depending on an external `+REQ` folder.

The intended behavior is:

- if the estimator is asked for a theoretical discrete quantile, it should compute it internally from the requested field type;
- the code should not depend on external research folders;
- theory-based quantiles can serve as baselines against ML-predicted quantiles.

Theory-based quantiles remain useful as:

- sanity checks;
- baselines;
- diagnostic references.

They are not a replacement for validating adaptive ML models.

## 8. Test 12: cs_guess and Effective Window Sweep

### Motivation

Test 11 revealed that a fixed `cs_guess` coupled:

```text
REQ_M
REQ_cs_guess
SIM_cs_bg
M_eff_true_diag
```

Test 12 was created to decouple these variables.

Central question:

> Does the model fail because of nominal `REQ_M`, because of `REQ_cs_guess`, or because of the true effective window size?

### Config and Run Script

Created:

- `configs/test_12_cs_guess_window_sweep.m`
- `experiments/run_test_12_cs_guess_window_sweep.m`

Test 12 sweeps:

- `SIM.WaveModel`
- `SIM.f0`
- `SIM.cs_bg`
- `REQ.M`
- `REQ.cs_guess`

It intentionally does **not** include:

- SNR sweep;
- attenuation;
- depth-dependent SNR;
- bilayer;
- k-Wave.

This keeps the test focused on effective-window physics.

### Dataset Status

The loaded dataset contains:

- 112,500 rows;
- 178 columns;
- `REQ_cs_guess`;
- `M_eff_guess`;
- `M_eff_true_diag`;
- local `req_mapping`;
- `global_req_mapping`;
- local and global spectral/Ecum/Srad-proxy features.

## 9. Test 12 Analysis Scripts Prepared

Created:

- `experiments/analysis/analyze_test_12_level01_model_comparison.m`
- `experiments/analysis/analyze_test_12_level02_grouped_generalization.m`
- `experiments/analysis/analyze_test_12_level03_compare_test11.m`
- `src/+adaptive_req/+analysis/Test12Analysis.m`
- `docs/test_12_cs_guess_window_sweep_analysis_summary.md`

### Level 01: Model and Feature-Set Comparison

Purpose:

- compare `LocalOnly`, `GlobalOnly`, and `HybridLocalGlobal` within Test 12;
- test whether `REQ_cs_guess` and/or `M_eff_guess` improve the model.

Feature sets:

- `NoCsGuess`;
- `WithCsGuess`;
- `WithMeffGuess`;
- `WithCsGuessAndMeffGuess`;
- `DiagnosticWithMeffTrue`.

Important rule:

- `M_eff_true_diag` is diagnostic-only and is not used in operational models.

### Level 02: Grouped Generalization

Purpose:

- evaluate strict leave-one-group-out generalization inside Test 12.

Grouped splits:

- leave-one-frequency-out;
- leave-one-`REQ_M`-out;
- leave-one-`REQ_cs_guess`-out;
- leave-one-`SIM_cs_bg`-out;
- leave-one-aperture-out;
- leave-one-wave-model-out when more than one wave model exists.

Main question:

> Does adding `REQ_cs_guess`, `M_eff_guess`, or both improve generalization to unseen regimes?

### Level 03: Test 12 vs Test 11

Purpose:

- compare Test 12 against Test 11 globally;
- compare matched/subset conditions when possible.

The matched comparison is important because Test 12 changed the dataset design. It attempts to restrict comparisons to common values of:

- `SIM_f0`;
- `SIM_cs_bg`;
- `REQ_M`;
- `SIM_WaveModel`;
- `REQ_cs_guess == 3.0` in Test 12 when available.

### Validation Performed

The new Test 12 analysis code was syntax-checked with MATLAB `checkcode`.

Result:

```text
0 checkcode issues
```

A smoke test was also run on 500 Test 12 rows:

- loaded the dataset;
- built all 15 model/feature-set combinations;
- created a condition-level split;
- trained a minimal linear `LocalOnly / NoCsGuess` model;
- converted `q_pred` to SWS;
- confirmed finite `cs_pred`.

The full Test 12 analyses were not run yet because they train many models and can be time-consuming.

## 10. Current Scientific Interpretation

The current working interpretation is:

1. Adaptive REQ is promising because spectral features can predict useful quantiles.
2. `HybridLocalGlobal` is the best current operational direction.
3. `GlobalOnly` is useful diagnostically but not enough for local SWS mapping.
4. The major failure mode is not simply nominal `REQ_M`; it is likely effective spatial support.
5. Small effective windows, especially around 1.5-2 true wavelengths, are fragile.
6. Test 12 is needed to determine whether `REQ_cs_guess` and `M_eff_guess` can make the model aware of this support issue operationally.

## 11. Leakage Rules Established

Operational models must not use:

```text
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
```

Allowed operational additions in Test 12:

- `REQ_cs_guess`;
- `M_eff_guess`.

Diagnostic-only:

- `M_eff_true_diag`.

## 12. GitHub / Reproducibility Work

The project was prepared conceptually for GitHub version control:

- repository name discussed: `adaptive_req`;
- private repository preferred;
- goal is to version code, configs, docs, and lightweight metadata;
- generated outputs should be handled carefully because some result files are large.

Suggested eventual GitHub hygiene:

- track source code, configs, scripts, and docs;
- avoid committing large raw outputs unless intentionally curated;
- add or refine `.gitignore` for bulky outputs;
- keep dated progress summaries like this one.

## 13. What Is Ready To Run Next

Recommended order:

```matlab
experiments/analysis/analyze_test_12_level01_model_comparison.m
experiments/analysis/analyze_test_12_level02_grouped_generalization.m
experiments/analysis/analyze_test_12_level03_compare_test11.m
```

After running these, update:

- `docs/test_12_cs_guess_window_sweep_analysis_summary.md`;
- this progress summary or a new dated progress summary.

## 14. Next Experiments Not Yet Created

Do not mix these into Test 12. They should be separate tests:

1. SNR sweep.
2. Depth-dependent SNR.
3. Attenuation.
4. Bilayer transfer with improved REQ map generation.
5. Inclusion phantom transfer.
6. k-Wave validation.
7. Broader stress tests varying number of waves, source model, amplitude jitter, and aperture.

## 15. Current Project State

Implemented and ready:

- Test 11 dataset and analyses through Level 05.
- Test 12 config and run script.
- Test 12 analysis Levels 01-03.
- Shared Test 12 analysis helper utilities.
- English Test 11 closure document.
- Spanish Test 12 placeholder summary document.
- This dated English progress summary.

Still pending:

- run full Test 12 Level 01;
- run full Test 12 Level 02;
- run Test 12 Level 03 comparison after Level 01/02 outputs exist;
- fill in final Test 12 numerical conclusions;
- decide whether `REQ_cs_guess`, `M_eff_guess`, or both improve operational performance;
- decide whether Test 12 improves over Test 11 in the matched comparison.

## 16. Short Version

The project has moved from exploratory quantile prediction to a structured adaptive-REQ workflow. Test 11 established that local and hybrid spectral features can predict useful REQ quantiles, but also revealed that high-error outliers cluster in small effective-window regimes. Test 12 was designed and run to decouple nominal `REQ_M`, guessed speed, true speed, and effective window size. The analysis scripts are now ready to determine whether `REQ_cs_guess` and `M_eff_guess` genuinely improve the model and whether Test 12 improves over Test 11 under fair matched conditions.
