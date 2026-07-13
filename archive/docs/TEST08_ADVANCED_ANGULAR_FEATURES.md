# Test 08: Advanced Angular Features

## Scientific Goal

Test 08 evaluates whether REQ quantile prediction improves when the model is
given rotation-invariant features that describe angular spectral geometry.
The operational models must generalize across directional and diffuse fields
without receiving the simulated `WaveModel`, true SWS, or aperture.

## New Features

The feature extractor now stores:

- `ang_moment_1`: first circular Fourier moment.
- `ang_moment_2`: detects opposing or bidirectional lobes.
- `ang_moment_4`: detects higher-order angular structure.
- `ang_peak_count_rel`: number of separated lobes above a relative threshold.
- `ang_top1_window_frac`: energy in the strongest angular window.
- `ang_top2_window_frac`: energy in the two strongest separated windows.
- `ang_top3_window_frac`: energy in the three strongest separated windows.
- `ang_top2_to_top1`: relative strength of the second lobe.
- `ang_peak_separation_deg`: separation between the two strongest lobes.

The implementation lives in:

```text
src/+adaptive_req/+features/extract_angular_shape_features.m
```

## Why Simulations Must Be Rerun

Test 07 stored scalar features but did not store the complete feature
structures or angular spectra. Therefore, the new angular descriptors cannot
be reconstructed from the existing Test 07 output.

Test 08 preserves the same physical sweep as Test 07 and writes results to a
new output folder.

## Run the Dataset

From MATLAB:

```matlab
run('experiments/run_test_08_advanced_angular_features.m')
```

This runs all 54 conditions and stores the lightweight `req_mapping`.

## Compare Models

After Test 08 finishes:

```matlab
run('experiments/analysis/analyze_test_08_level14_compare_advanced_features.m')
```

Level 14 compares:

| Model | Role | Additional information |
|---|---|---|
| Model C | Operational baseline | Existing rich spectral features |
| Model E | Operational | Advanced angular-shape features |
| Model F | Operational | Advanced features plus `M_eff_guess` |
| Model G | Diagnostic aperture | Advanced features plus `Omega_sr` |
| Model D | Diagnostic true-SWS ceiling | Advanced features plus true `M_eff` |

`ModelG_diagnostic_with_aperture` must not be presented as the main
deployable model. It measures how much performance could improve if aperture
were explicitly supplied.

## Interpretation

The main scientific comparison is Model C versus Models E and F. Improvement
there supports the claim that the local spectrum contains enough information
to adapt REQ across directional and diffuse fields.

Models G and D provide diagnostic ceilings and help identify information that
the operational features still fail to recover.
