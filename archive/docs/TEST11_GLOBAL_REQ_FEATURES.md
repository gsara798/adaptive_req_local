# Test 11: Local vs Global REQ Features

## Goal

Test 11 compares three ways of learning the REQ quantile:

- Local-only: each patch uses its own FFT, local Ecum/Srad-proxy features, and local `q_theory`.
- Global-only: each full simulated field gets one global FFT/REQ curve, and every patch row receives the same `global_*` descriptors.
- Hybrid: local patch descriptors plus global full-field descriptors.

It also stores `q_global_theory` and `q_local_minus_global` to diagnose whether a global q is close enough to the local q. This is especially important before using the model in heterogeneous fields, where a global q may erase spatial differences.

## Run

Generate the dataset:

```matlab
run('/Users/sara/Documents/adaptive_req_local/experiments/run_test_11_global_req_features.m')
```

Analyze local/global/hybrid models:

```matlab
run('/Users/sara/Documents/adaptive_req_local/experiments/analysis/analyze_test_11_level01_global_vs_local.m')
```

## Outputs

The dataset is saved under:

```text
/Users/sara/Documents/adaptive_req_local/outputs/test_11_global_req_features/
```

The analysis creates:

- `level11_sws_metrics.csv`
- `level11_sws_metrics_by_M.csv`
- `level11_sws_metrics_by_aperture.csv`
- `level11_q_local_global_gap_summary.csv`
- plots comparing local-only, global-only, hybrid, and diagnostic global-q models.

## Interpretation

Use `HybridLocalGlobal` as the main operational candidate if it improves over `LocalOnly`.

Use `GlobalQDirect` and `GlobalQDiagnostic` only as diagnostics. They use global q information that is not equivalent to a deployable local estimate, but they reveal whether global REQ contains useful information for the local target.
