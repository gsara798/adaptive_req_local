# Level 16: Error Learning and Feature Audit

## Purpose

Level 16 tests three ideas:

1. Can the model learn from its own q errors?
2. Can REQ theory about the radial spectrum `Srad(k)` provide additional
   operational features?
3. Are any predictors hurting the model?

## Theory-Derived Features

REQ defines `Ecum(k)` as the cumulative distribution of the radial mean
spectrum `Srad(k)`. Because `Srad(k)` is proportional to the local derivative
or increment of `Ecum(k)`, Level 16 extracts lightweight `Srad` proxies from
the stored `req_mapping`, without storing the heavy 2D spectrum.

New proxy features include:

```text
srad_proxy_centroid_k_norm
srad_proxy_std_k_norm
srad_proxy_skewness
srad_proxy_kurtosis
srad_proxy_peak_k_norm
srad_proxy_peak_to_centroid
srad_proxy_low_side_frac
srad_proxy_high_side_frac
```

These describe whether the radial distribution is sharp, broad, skewed toward
low projected wavenumbers, or peak-concentrated.

## Error Learning

The error-learning model is a two-stage model:

```text
q_base = F(features)
q_residual = G(features, q_base)
q_corrected = clip(q_base + q_residual, 0, 1)
```

The residual target is `q_true - q_base`, trained only on the training split and
applied to the test split.

## Result

| Model | SWS MAPE | SWS RMSE | p95 APE | Max APE | >20% |
|---|---:|---:|---:|---:|---:|
| Level16 residual corrected | 2.70% | 4.85% | 10.17% | 33.16% | 1.17% |
| Level16 base Ecum + Srad proxy | 3.06% | 5.24% | 11.60% | 34.54% | 1.29% |

The residual corrector improves the test split, so this is a promising path.
The high-error rate above 20% also drops from 31/2400 to 28/2400.

By `M`, the corrected model still concentrates large errors in `M = 2`:

| REQ M | n | n > 20% | percent |
|---:|---:|---:|---:|
| 2 | 900 | 25 | 2.78% |
| 3 | 900 | 3 | 0.33% |
| 4 | 600 | 0 | 0.00% |

## Feature Audit

The strongest feature association with q is `ecum_increment_gini`, followed by
angular/radial concentration features such as `circ_var`,
`ecum_increment_peak_frac`, `ang_entropy`, and `radial_entropy`.

Drop-one ablation should be interpreted cautiously because correlated features
can substitute for each other. Still, it gives a useful signal:

- Dropping `ecum_increment_gini` worsened SWS RMSE by about +0.51 percentage
  points, so it appears strongly useful.
- Dropping `REQ_M`, `radial_entropy`, `ecum_asymmetry_10_90`,
  `width_90_50_rel`, and `ecum_width_ratio_80_50` also worsened performance.
- Dropping `ecum_lower_upper_width_ratio`, `dom_dir_frac`,
  `REQ_Nbins_effective`, `ecum_auc_norm`, or `ecum_upper_tail_rel` slightly
  improved performance. These are candidates for pruning, but the effect is
  small and should be confirmed with repeated splits before removing them from
  the operational model.

## Run

```matlab
run('experiments/analysis/analyze_test_08_level16_error_learning_feature_audit.m')
```

Outputs are written to:

```text
outputs/test_08_advanced_angular_features/<run>/analysis/level_16_error_learning_feature_audit
```

Key files:

```text
level16_sws_metrics_by_model.csv
level16_high_error_rate_by_model.csv
level16_high_error_rate_by_model_M.csv
level16_high_error_rate_by_model_cs_bg.csv
level16_feature_q_error_associations.csv
level16_drop_one_ablation_sws_metrics.csv
level16_model_sws_rmse.png
level16_error_learning_diagnostics.png
level16_feature_association.png
level16_drop_one_ablation.png
```
