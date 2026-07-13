# Level 17: Model Comparison Plots and Next Tests

## Purpose

Level 17 is a plotting/reporting layer. It does not retrain models. It reads
Level 15 and Level 16 prediction CSVs and compares the operational candidates
with consistent labels, SWS metrics, q diagnostics, outlier plots, and
aperture/M breakdowns.

## Best Current Model

The best current model is:

```text
Level16_residual_corrected / bagged_trees
```

Test split summary:

| Model | MAPE | RMSE | p95 APE | Max APE |
|---|---:|---:|---:|---:|
| L16 residual corrected | 2.70% | 4.85% | 10.17% | 33.16% |
| H Ecum | 3.05% | 5.16% | 11.62% | 36.06% |
| I Angular + Ecum | 3.10% | 5.18% | 11.63% | 35.46% |
| L16 base Ecum + Srad | 3.06% | 5.24% | 11.60% | 34.54% |
| J Ecum by M | 3.32% | 5.73% | 12.78% | 46.22% |
| C baseline | 4.43% | 6.84% | 15.07% | 42.64% |

The residual corrector improves absolute SWS error in 64.83% of paired test
points. The median paired change is -0.17 percentage points, and the mean
paired change is -0.37 percentage points.

## Key Outputs

```text
level17_model_metric_dashboard.png
level17_sws_pred_vs_true_grid.png
level17_error_boxplots.png
level17_rmse_by_M_and_cs.png
level17_aperture_and_high_error.png
level17_q_diagnostics.png
level17_level16_paired_improvement.png
level17_worst_outliers.png
```

Tables:

```text
level17_model_comparison_metrics.csv
level17_model_comparison_by_M.csv
level17_model_comparison_by_cs_bg.csv
level17_model_comparison_by_aperture.csv
level17_high_error_rate_by_model.csv
level17_high_error_rate_by_model_M.csv
level17_level16_paired_error_delta.csv
level17_worst_outliers.csv
```

## How To Know If The Model Is Good

The current result is promising, but not enough by itself. It is good on the
current synthetic distribution. To claim robustness, the model should pass
tests where the data-generating process changes.

Recommended validation ladder:

1. Repeat train/test splits by condition with multiple random seeds.
2. Hold out full regimes, not just random conditions: one `M`, one SWS, one
   frequency, or one aperture range.
3. Add stress tests: SNR, number of waves, amplitude distributions, spatial
   sampling, window size, and out-of-plane content.
4. Test calibration: predicted error probability or uncertainty should increase
   where SWS error increases.
5. Compare against simple REQ fixed-q baselines and against q_true oracle SWS.

## Next Experiment Ideas

The cleanest next step is a new Test 09 dataset designed as a robustness
matrix rather than only an accuracy dataset.

Suggested factors:

```text
SNR_dB: Inf, 40, 30, 20, 10
Nwaves: 20, 50, 100, 500, 2000
WaveModel: planewave, spherical
REQ_M: 2, 3, 4
cs_bg: 2, 3, 4 m/s
f0: 400, 500, 600 Hz
aperture steps: current 10-step cone schedule
ForceInPlaneWave: false, true
```

To keep runtime controlled, do this in stages:

1. Small screening design: fewer realizations, all factors.
2. Identify fragile regimes: likely low `M`, low SNR, low Nwaves, high aperture.
3. Dense follow-up only in fragile regimes.
4. Train on broad conditions, then test on held-out stress regimes.

## Run

```matlab
run('experiments/analysis/analyze_test_08_level17_model_comparison_plots.m')
```

Outputs are written to:

```text
outputs/test_08_advanced_angular_features/<run>/analysis/level_17_model_comparison_plots
```
