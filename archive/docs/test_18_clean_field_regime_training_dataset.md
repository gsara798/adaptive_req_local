# Test 18: Clean Field-Regime Training Dataset

## Purpose

Test 18 creates a homogeneous, noise-free dataset that explicitly represents
the physical field regimes used in later transfer tests. It addresses a gap in
the previous aperture-sweep training data: the older models may not have seen
enough canonical directional 2D, diffuse 2D, partially diffuse 3D, and diffuse
3D fields.

Test 18 generates data only. Model retraining and old-versus-new comparisons
are intentionally deferred to Test 19.

## Field regimes

- `directional_2D`: one spherical 2D wave.
- `diffuse_2D`: 32, 64, 128, or 256 in-plane waves.
- `partial_3D`: 4, 8, 16, or 32 Fibonacci-sphere sources, with at least one
  source forced into the xz plane.
- `diffuse_3D`: 128 or 256 Fibonacci-sphere sources.

## Physical grid

- SWS: 2, 3, and 4 m/s.
- Frequency: 400, 500, and 600 Hz.
- REQ M: 2, 3, and 4.
- Paired `dx = dz`: 0.1, 0.2, 0.3, and 0.5 mm.
- Three realizations and five patches per physical condition and M.
- SNR is infinite and no attenuation or amplitude jitter is included.

The full design contains 396 physical conditions and 17,820 dataset rows.
Each simulated field is reused for all three M values.

## Predictor policy

The following variables are diagnostic metadata and must not enter the main
operational models:

- `field_regime_label`
- `field_regime_variant`
- `SIM_Nwaves`
- `SIM_Is2D`
- `SIM_ForceInPlaneWave`

They may be used for grouping, validation, and a future user-guess analysis.
The primary model must infer the regime from operational local/global spectral
and cumulative-energy features.

## Outputs

The run is resumable and uses a stable output folder:

```text
outputs/test_18_clean_field_regime_training_dataset/dataset/
```

It stores per-condition checkpoints, the complete MAT dataset, a lightweight
CSV table, condition/status tables, predictor-policy metadata, and sanity-check
figures for representative fields, spectra, q distributions, and occupancy.

## Running

From the project root in MATLAB:

```matlab
run('experiments/run_test_18_clean_field_regime_training_dataset.m')
```

For a short validation:

```matlab
setenv('ADAPTIVE_REQ_TEST18_MODE', 'smoke');
run('experiments/run_test_18_clean_field_regime_training_dataset.m');
setenv('ADAPTIVE_REQ_TEST18_MODE', '');
```
