# Test 09 and Test 10 Plan

## Separation of Tests

These are intentionally separate:

- **Test 09 stress robustness**: homogeneous simulations with controlled stress
  factors such as SNR and number of waves.
- **Test 10 heterogeneous maps**: simple bilayer/inclusion phantoms for visual
  map inspection and ROI statistics.

This prevents stress robustness and heterogeneous-phantom behavior from being
mixed into one ambiguous result.

## Test 09: Stress Robustness

Config:

```text
configs/test_09_stress_robustness.m
```

Run:

```matlab
run('experiments/run_test_09_stress_robustness.m')
```

Analyze after it finishes:

```matlab
run('experiments/analysis/analyze_test_09_level01_stress_robustness.m')
```

Screening design:

```text
WaveModel: planewave, spherical
Nwaves: 50, 200, 2000
SNR_dB: Inf, 30, 20, 10
ForceInPlaneWave: false, true
REQ_M: 2, 3, 4
f0: 500 Hz
cs_bg: 3 m/s
steps: 10 aperture steps
realizations: 2
patches: 5
```

This produces 144 conditions and 14,400 patch rows. It is larger than the
heterogeneous test and may take a while, but it is still a screening dataset,
not the final robustness dataset.

## Test 10: Heterogeneous Maps

Run:

```matlab
run('experiments/run_test_10_heterogeneous_maps.m')
```

What it does:

1. Trains the current Level 16 residual-corrected model on Test 08.
2. Simulates two simple heterogeneous phantoms:
   - bilayer background/fast layer;
   - circular fast inclusion.
3. Extracts patch-level features on a 9 x 9 grid.
4. Predicts q, converts q to SWS, and compares to the median true SWS inside
   each patch.
5. Saves maps, ROI MAPE, ROI CoV, and patch-level error distributions.

Current run:

```text
outputs/test_10_heterogeneous_maps/test_10_heterogeneous_maps_2026-06-08_175652
```

Current ROI result:

| Phantom | ROI | MAPE | RMSE | Bias | CoV |
|---|---|---:|---:|---:|---:|
| bilayer | background/boundary | 14.72% | 26.03% | -11.89% | 26.64% |
| bilayer | fast | 22.98% | 36.43% | -22.31% | 39.61% |
| circular inclusion | background/boundary | 32.96% | 44.44% | -31.00% | 45.83% |
| circular inclusion | fast | 39.91% | 49.69% | -39.43% | 53.61% |

Interpretation: the homogeneous-trained model does not yet generalize well to
heterogeneous maps. It strongly underestimates fast regions, especially in the
inclusion case. This is useful evidence: the model is good on the current
homogeneous synthetic distribution, but heterogeneous phantoms need either
heterogeneous training examples, better map-level inference, or boundary-aware
features.

## Viability

The heterogeneous test is viable as a qualitative diagnostic now. It is not yet
a fair final validation because:

- patches near boundaries mix multiple true speeds;
- the training distribution is homogeneous;
- spherical propagation through the heterogeneous `k_map` is an approximation;
- the model has never seen inclusion or bilayer conditions.

The next fair step is to add simple heterogeneous phantoms into training or
fine-tuning, while holding out other heterogeneous phantoms for testing.

## Recommended Next Steps

1. Run Test 09 to quantify homogeneous robustness under noise and low wave
   count.
2. Use Test 09 to decide which stress factors must enter the main training set.
3. Improve Test 10 by adding:
   - more ROI definitions, especially eroded interior ROIs away from edges;
   - separate boundary ROI;
   - repeated seeds;
   - SNR and Nwaves variations.
4. Create Test 11 heterogeneous training/fine-tuning:
   - train on some bilayer/inclusion settings;
   - test on held-out geometry, contrast, radius, and boundary softness.
