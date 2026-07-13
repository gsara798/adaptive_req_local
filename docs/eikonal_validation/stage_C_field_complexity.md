# Stage C: Field Complexity Transfer Validation

Stage C evaluates whether the frozen `baseline_minimal_v1` q/SWS models remain stable when clean Eikonal fields become less directional and more diffuse-like. No model is retrained, and no readout noise, risk mask, correction, or reliability layer is used.

## Design
- Geometries: homogeneous_cs3, homogeneous_cs4, inclusion_2_3_D24, inclusion_2_4_D24, bilayer_2_3.
- Frequencies: 300, 400, 500, 600 Hz.
- Field regimes: single_source_lateral, diffuse_like_4src_layoutA, diffuse_like_8src_layoutA, diffuse_like_16src_layoutA, diffuse_like_4src_layoutB, diffuse_like_8src_layoutB, diffuse_like_16src_layoutB.
- REQ: M=2, cs_guess=3.00 m/s, target step 1.00 mm, valid windows only.
- Clean-only fields are globally RMS-normalized over the central evaluation region so source-count effects are not confounded with global amplitude scaling. Local interference patterns are preserved.

## Models
- Primary: `q_spectrum_plus_composition`.
- Diagnostic: `q_spectrum_only`.

## Main Results
| model_name | N_valid_patches | MAPE_pct | signed_bias_pct | high_error_20_pct |
| --- | --- | --- | --- | --- |
| q_spectrum_only | 5.49e+04 | 2.02 | -1.73 | 2.13 |
| q_spectrum_plus_composition | 5.49e+04 | 2.02 | -1.75 | 2.13 |

### Primary Model By Geometry
| geometry | MAPE_pct | signed_bias_pct | high_error_20_pct |
| --- | --- | --- | --- |
| homogeneous_cs4 | 1.41 | -1.38 | 0.0146 |
| inclusion_2_4_D24 | 2.63 | -2.11 | 4.25 |

### Primary Model By Field Regime
| field_regime | N_sources | N_in_plane_sources | MAPE_pct | signed_bias_pct | high_error_20_pct |
| --- | --- | --- | --- | --- | --- |
| diffuse_like_16src_layoutA | 16 | 1 | 2.27 | -2.17 | 2.24 |
| diffuse_like_4src_layoutA | 4 | 1 | 2.39 | -2.28 | 2.31 |
| diffuse_like_8src_layoutA | 8 | 1 | 2.21 | -2.11 | 2.3 |
| single_source_lateral | 1 | 1 | 1.2 | -0.439 | 1.68 |

### Primary Model By Frequency
| f0 | MAPE_pct | signed_bias_pct | high_error_20_pct |
| --- | --- | --- | --- |
| 400 | 2.5 | -2.2 | 2.89 |
| 600 | 1.61 | -1.37 | 1.49 |

## Interpretation Guide
- If homogeneous 3 m/s stays accurate across all source counts, the feature extraction and frozen model are stable for clean Eikonal fields without hard-speed bias.
- If homogeneous 4 m/s becomes increasingly negative with source count, the issue is hard-speed plus field-regime shift, not interface mixing.
- If inclusion 2/4 degrades more than homogeneous 4 m/s, field complexity is interacting with the hard interface.
- If bilayer error is concentrated in interface ROIs, this supports a window-mixing interpretation.

## Worst Condition
Worst primary-model condition: `inclusion_2_4_D24__f400__clean__diffuse_like_4src_layoutA__M2`, inclusion 2/4, D=24 mm, 400 Hz, 4 sources, layout A, MAPE 3.76%, bias -3.54%.

## Figures
Key analysis figures are under `figures/`, including source-count trends, geometry-field heatmaps, ROI-field heatmaps, frequency trends, model comparisons, and representative maps from the runner.

## How To Run
```matlab
setenv('ADAPTIVE_REQ_EIKONAL_STAGE_C_MODE','validate_only'); run('experiments/runners/eikonal_validation/run_stage_c_field_complexity.m')
setenv('ADAPTIVE_REQ_EIKONAL_STAGE_C_MODE','quick'); run('experiments/runners/eikonal_validation/run_stage_c_field_complexity.m')
setenv('ADAPTIVE_REQ_EIKONAL_STAGE_C_MODE','full'); run('experiments/runners/eikonal_validation/run_stage_c_field_complexity.m')
run('experiments/analysis/eikonal_validation/analyze_stage_c_field_complexity.m')
```
