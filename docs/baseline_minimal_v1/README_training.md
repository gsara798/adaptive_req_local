# Baseline Minimal v1 Training Recipe

This document describes the minimal clean baseline used to rebuild the
publishable q/SWS estimator after the exploratory Test 38/Test 53 line.

The goal is intentionally narrow:

- train a clean operational q estimator;
- keep the experiment reproducible through a JSON config;
- avoid dependencies on older T18/old/confidence/correction models;
- keep training separate from analysis.

## Entry Points

Training / experiment runner:

```matlab
experiments/runners/run_baseline_minimal_v1.m
```

Analysis-only companion:

```matlab
experiments/analysis/analyze_baseline_minimal_v1.m
```

Main config:

```matlab
configs/final_training/baseline_minimal_v1.json
```

Shared REQ defaults:

```matlab
configs/shared/req_defaults.json
configs/shared/feature_defaults.json
```

## Models Trained

Only two deployable q models are trained.

### 1. q_spectrum_only

`q_spectrum_only` is a bagged-tree regression model trained to predict
`q_oracle` from operational spectral/REQ features only.

It does not use:

- true SWS;
- true material label;
- true patch purity;
- distance to interface;
- prediction error;
- confidence/risk;
- previous T18/old model predictions;
- manually supplied oracle regions.

### 2. q_spectrum_plus_composition

`q_spectrum_plus_composition` uses the same operational spectral/REQ
features as `q_spectrum_only`, plus three internally predicted composition
features:

- `predicted_patch_purity`;
- `p_mixed`;
- `p_strong_mixed`.

These three auxiliary composition outputs are produced by models trained
inside the same experiment from the train split only. They are saved
separately so the composition block can be inspected or reused later.

Saved auxiliary artifact:

```matlab
outputs/baseline_minimal_v1[/quick|/medium]/models/composition_auxiliary_models.mat
```

Saved full model bundle:

```matlab
outputs/baseline_minimal_v1[/quick|/medium]/models/baseline_minimal_v1_q_models.mat
```

## Target

The q target is `q_oracle`.

For each REQ patch, the script computes:

```matlab
k_true = 2*pi*f0 / SWS_true(center)
```

Then `q_oracle` is found by inverting the local REQ cumulative mapping:
the selected q is the cumulative-energy quantile whose mapped wavenumber is
closest to `k_true`.

`SWS_true` is only used to define this supervised target and for evaluation.
It is not passed as an input feature.

## Training Data

The default full recipe uses clean synthetic wave fields with:

- `dx = dz = 0.2 mm`;
- `f0 = [200 300 400 500 600] Hz`;
- `M = [2 3]`;
- `cs_guess = 3.0 m/s`;
- `TargetStepM = 2.0 mm`.

The current config includes these geometry cases:

- homogeneous: `1.5`, `2.0`, `2.5`, `3.0`, `3.5`, `4.0 m/s`;
- bilayer: `2/3`, `2/4 m/s`;
- inclusion: `2/3`, `2/4 m/s`;
- ellipse `2/4 m/s`;
- off-center inclusion `2/4 m/s`;
- two inclusions `2/4 m/s`;
- three-material map `2/3/4 m/s`.

The current config includes these field regimes:

- `directional_2D_angle0`;
- `directional_2D_angle30`;
- `directional_2D_angle60`;
- `diffuse_2D_seed1`;
- `diffuse_3D_seed1`;
- `partial_3D_8src`;
- `partial_3D_16src`.

Every regime forces at least one in-plane source/wave when applicable.

## REQ Settings

REQ is extracted with:

```matlab
QuantileMode = "local_req"
Nbins = "auto"
Nbins_auto_oversample = 1
Nbins_min = 16
smooth_sigma = 1
Gamma = 1
PadFactor = 1
EdgeMode = "valid"
```

The spatial sampling of REQ windows is controlled in physical units:

```matlab
StepX = max(1, round(TargetStepM / dx))
StepZ = max(1, round(TargetStepM / dz))
```

For the default full config, `TargetStepM = 2.0 mm`, so the model is trained
on a sparse but physically consistent patch grid. Dense maps should be
generated in separate validation scripts, not by making this training run
unnecessarily dense.

## Train/Test Split

The runner uses a condition-level split, not a random patch split.

The default train fraction is:

```matlab
TrainFraction = 0.70
```

Condition grouping includes:

- geometry/case;
- field regime;
- frequency;
- M.

This is not the strongest possible validation. Stronger validation remains
the responsibility of the strong-split analysis line, but this baseline
runner keeps the train/test separation cleaner than patch-level random
splitting.

## Operational Predictors

The base predictor table is selected automatically from numeric feature
columns after applying a leakage guard.

Forbidden predictor name patterns include:

- `true`;
- `oracle`;
- `purity`;
- `mixed`;
- `confidence`;
- `error`;
- `pred`;
- `sws`;
- `cs_`;
- `k_true`;
- `q_local`;
- `q_pred`;
- `q_theory`;
- `req_mapping`;
- patch/map coordinates;
- interface distance;
- condition identifiers.

Allowed metadata-like operational predictors include:

- `REQ_M`;
- `M`;
- `SIM_f0`;
- `f0`;
- `dx`;
- `dz`;
- `REQ_Nbins_effective`.

## Metrics Saved

The runner saves patch-level predictions and summary tables under:

```matlab
outputs/baseline_minimal_v1[/quick|/medium]/tables/
```

Important tables:

- `baseline_minimal_v1_patch_level_predictions.csv`;
- `baseline_minimal_v1_summary_overall.csv`;
- `baseline_minimal_v1_summary_by_case.csv`;
- `baseline_minimal_v1_summary_by_family.csv`;
- `baseline_minimal_v1_summary_by_regime.csv`;
- `baseline_minimal_v1_summary_by_frequency.csv`;
- `baseline_minimal_v1_summary_by_M.csv`;
- `baseline_minimal_v1_summary_by_purity_bin.csv`;
- `baseline_minimal_v1_summary_by_geometry_frequency.csv`;
- `baseline_minimal_v1_summary_by_roi_frequency.csv`;
- `baseline_minimal_v1_summary_by_geometry_roi_frequency.csv`;
- `baseline_minimal_v1_worst_conditions.csv`.

ROI labels are evaluation-only and include:

- `homogeneous_center`;
- `soft_core`;
- `intermediate_core`;
- `hard_core`;
- `interface_0_1mm`;
- `interface_1_2mm`;
- `interface_2_4mm`.

The ROI labels use true material and distance information only after
prediction, for analysis. They are never used as model inputs.

## Why cs_guess Is Fixed Here

This baseline intentionally keeps:

```matlab
cs_guess = 3.0 m/s
```

Changing `cs_guess` changes the physical REQ window and the feature
distribution. That is a real modeling question, but it should be tested as a
separate controlled experiment. The clean baseline should remain simple and
auditable.

The next natural experiment is a two-pass or cs-guess sweep:

1. first pass with fixed `cs_guess = 3`;
2. estimate SWS;
3. recompute REQ with a predicted/adaptive `cs_guess`;
4. evaluate whether hard-core bias improves without damaging soft or
   homogeneous regions.

## How To Run

From the local project root:

```bash
cd /Users/sara/local/adaptive_req_local
```

Validation-only smoke test:

```bash
/Applications/MATLAB_R2025a.app/bin/matlab -batch "setenv('ADAPTIVE_REQ_BASELINE_MODE','quick'); setenv('ADAPTIVE_REQ_BASELINE_VALIDATE_ONLY','true'); run('experiments/runners/run_baseline_minimal_v1.m')"
```

Quick training run:

```bash
/Applications/MATLAB_R2025a.app/bin/matlab -batch "setenv('ADAPTIVE_REQ_BASELINE_MODE','quick'); setenv('ADAPTIVE_REQ_BASELINE_VALIDATE_ONLY','false'); run('experiments/runners/run_baseline_minimal_v1.m')"
```

Full training run:

```bash
/Applications/MATLAB_R2025a.app/bin/matlab -batch "setenv('ADAPTIVE_REQ_BASELINE_MODE','full'); setenv('ADAPTIVE_REQ_BASELINE_VALIDATE_ONLY','false'); setenv('ADAPTIVE_REQ_BASELINE_USE_PARALLEL_TRAINING','true'); run('experiments/runners/run_baseline_minimal_v1.m')"
```

Analysis after training:

```bash
/Applications/MATLAB_R2025a.app/bin/matlab -batch "setenv('ADAPTIVE_REQ_BASELINE_MODE','full'); run('experiments/analysis/analyze_baseline_minimal_v1.m')"
```

## Main Interpretation Rule

For paper use, prefer the simplest model that is stable across:

- frequency;
- geometry;
- field regime;
- ROI/core regions;
- patch-purity bins.

If `q_spectrum_plus_composition` improves mixed/interface regions but harms
homogeneous or hard-core regions, that tradeoff should be reported directly
rather than hidden by global MAPE.
