# Adaptive REQ Progress Summary: Confidence, Mixedness, q-Spectrum Models, and Transfer

**Date:** 2026-06-28  
**Project:** `adaptive_req_local`  
**Scope:** Progress from frozen confidence validation through Test 43.  

## 1. Current Status in One Sentence

The most stable current direction is to use a clean spectral/composition q model
as the primary estimator, especially `q_spectrum_plus_composition`, and treat
mixedness/composition as reliability information unless a future correction
model can decide both whether to correct and in which direction.

## 2. Big-Picture Evolution

The project started with adaptive REQ models that predicted the REQ quantile
`q` from local/global spectral features, then moved through confidence
detection, interface diagnostics, structure-aware corrections, mixedness-aware
correction, and finally clean q-spectrum models trained on richer synthetic
data.

The operational goal remains:

```text
wavefield -> spectral/composition features -> q_pred -> REQ mapping -> k_pred -> SWS map
```

Important policy:

- true SWS is used only for labels/evaluation;
- true patch purity, true material side, and distance-to-interface are
  diagnostic/evaluation variables only;
- frozen q/confidence/composition models are not silently retrained in
  validation or correction tests;
- any correction method must be judged separately on homogeneous, pure
  heterogeneous, mixed/interface, soft, and hard regions.

## 3. Confidence and Interface Diagnostics

### Test 22: Frozen Confidence External Validation

Test 22 evaluated previously trained confidence detectors without retraining.
It covered homogeneous maps, bilayer 2/3 m/s, circular inclusion 2/3 m/s,
multiple field regimes, frequencies, M values, and spatial resolutions.

Key result:

- confidence maps were useful for locating likely high-error regions;
- `ML_bagged_trees` became the main practical detector;
- confidence alone was not enough to correct SWS maps.

Runtime was later optimized with `pilot`, `medium`, and `full` modes, plus a
physical target step for REQ window spacing.

### Tests 23-25: Patch Contamination and Two-Radius Spectra

These tests diagnosed whether interface errors were caused by patches crossing
soft/hard interfaces.

Main findings:

- large errors increase when `patch_purity` decreases;
- hard-side underestimation is common near mixed patches;
- two-radius spectral mixture is visible near interfaces;
- pure heterogeneous patches can still have nontrivial error, so not all error
  is explained by simple soft/hard contamination.

Interpretation:

```text
mixed patches -> spectral contamination is real
pure heterogeneous patches -> additional model/domain/context bias exists
```

## 4. Correction Attempts Before the New q-Spectrum Family

### Test 26: Confidence-Gated Corrections

Test 26 compared physical donut filters, prior-donut corrections,
confidence-gated donut corrections, interpolation, and peak-candidate methods.

Finding:

- gated interpolation helped somewhat;
- interface errors remained too high;
- some corrections smoothed or distorted sharp interfaces.

### Test 27: Adaptive Window and Edge-Aware Correction

Test 27 compared Hybrid, Local, confidence blending, smaller windows, and
edge-aware interpolation.

Findings:

- `LocalOnly_T18` often beat Hybrid in pure/near-pure heterogeneous regions;
- M=2 was usually safer near interfaces than larger M;
- edge-aware interpolation improved some maps but did not fully solve mixed
  patches.

### Tests 28-29: Oracle/Half-Window/Graph-TV Diagnostics

These tests probed upper bounds and visual diagnostics.

Findings:

- same-material/oracle references showed that better side identification could
  improve interface estimates;
- half-window and graph/TV ideas were diagnostically useful but not yet robust
  as deployable estimators.

## 5. Theory-Based Structure and Simple Interpolation

### Test 30: Theory Structure + Local M2 Region Levels

Test 30 used `TheoryQDiscrete` to infer a simple structure, then estimated
region levels from Local M=2.

Representative full synthetic results:

- global MAPE: about `1.67%`;
- homogeneous MAPE: about `0.46%`;
- pure MAPE: about `0.72%`;
- near-pure MAPE: about `1.75%`;
- moderate mixed MAPE: about `2.46%`;
- strong mixed MAPE: about `9.67%`.

Interpretation:

- excellent for clean two-material synthetic maps;
- too rigid to trust as the final general estimator for realistic/unknown
  tissue;
- very useful as a structural reference.

### Test 31: Simple Confidence-Based q/SWS Interpolation

Test 31 compared copying/interpolating q or SWS from high-confidence regions.

Representative full synthetic results:

- Local baseline MAPE: about `4.69%`;
- SWS nearest high-confidence: about `3.71%`;
- q median global: about `4.10%`;
- Test30 region levels: about `1.75%`;
- edge-aware SWS: about `3.76%`;
- edge-aware q: about `4.23%`.

Finding:

- SWS interpolation helped more than q interpolation;
- simple interpolation did not match theory-structure region levels.

### Test 32: Structural-Confidence Hybrid

Test 32 blended flexible local/interpolated maps with the rigid Test30
structure prior.

Representative full synthetic results:

- Local MAPE: about `4.69%`;
- SWS nearest MAPE: about `3.71%`;
- Test30 MAPE: about `1.75%`;
- relaxed structural hybrids: about `2.7%`.

Finding:

- hybrids reduced some rigidity concerns but still did not beat Test30 on clean
  two-material synthetic maps.

## 6. Mixedness-Aware Correction and k-Wave Transfer

### Test 33: Mixedness-Aware q/log-k Correction

Test 33 trained an operational mixedness detector, patch-purity regressor,
log-k residual corrector, q-candidate selector, and posterior reliability
model.

Representative synthetic held-out results:

- Local MAPE: about `4.10%`;
- Test30 region levels: about `1.56%`;
- boundary-protected hybrid: about `2.49%`;
- mixedness log-k corrected: about `2.22%`;
- mixedness q selector: about `3.30%`.

Mixedness detection was strong:

- mixed `<0.95`: ROC AUC about `0.990`, PR AUC about `0.975`;
- strong mixed `<0.75`: ROC AUC about `0.989`, PR AUC about `0.948`;
- patch-purity MAE about `0.026`.

Interpretation:

- mixedness can be predicted operationally;
- mixedness-aware correction can help;
- choosing among existing q candidates is weaker than continuous log-k
  correction.

### Test 34: k-Wave Mixedness Transfer

Test 34 applied the mixedness-aware correction idea to k-Wave-like data.

Important outcome:

- Test34 mixedness-logk correction became one of the strongest k-Wave baselines;
- k-Wave transfer remained harder than clean synthetic external validation;
- hard-side underestimation remained a central failure mode.

## 7. New Independent q-Spectrum Model Family

### Test 35: Spectral Composition to q Model

Test 35 deliberately moved away from the older T18/Hybrid stack. It trained
new models from clean spectral features, with optional predicted composition
features.

Important conceptual point:

- true patch purity is used as a training label for composition;
- predicted composition can be used by downstream q models;
- true purity is not an inference input.

The model family included:

- `q_spectrum_only`;
- `q_spectrum_plus_composition`;
- `q_spectrum_plus_theory_composition`;
- `delta_q_theory_composition`;
- `delta_logk` variants.

### Test 37: OOD Validation of Test35 q-Spectrum Models

Test 37 added out-of-distribution simulations:

- unseen velocities;
- unseen frequency values;
- new geometries;
- new directional angles and diffuse seeds.

Finding:

- q-spectrum models looked promising;
- failures became easier to diagnose with denser diagnostic maps;
- validating only on the same synthetic family was not enough.

### Test 38: Velocity/Field-Diverse q Training

Test 38 expanded training diversity:

- more velocities;
- multiple directional angles;
- random diffuse 3D variants;
- partial 3D variants with different source counts;
- diffuse 2D included;
- medium/full modes to manage runtime.

Recommended bundle:

```text
outputs/model_registry/test38_velocity_field_diverse_q_training/
test38__velocity_field_diverse_q__medium_bundle.mat
```

Current practical model preference:

```text
q_spectrum_plus_composition
```

Reason:

- stable across Test39 OOD and k-Wave transfer;
- conceptually uses useful composition information;
- avoids the instability of aggressive correction stacks.

### Test 39: Frozen Test38 External Validation

Test 39 validated frozen Test38 models on external/OOD simulations.

Important plotting/runtime decision:

- validation tables can stay sparse and fast;
- dense diagnostic maps can be recomputed separately for visual inspection.

Finding:

- q-spectrum/composition models are substantially more robust than earlier
  estimators in many synthetic OOD conditions;
- summary figures now omit `theory_discrete` where needed so learned model
  differences are visible.

## 8. k-Wave Comparison and Correction Transfer

### Test 40: k-Wave Latest Model Comparison

Test 40 compared old Test34 strategies with new Test35/Test38 frozen q models
on the same k-Wave table.

Representative k-Wave global results:

- `T34_mixedness_logk_corrected`: MAPE about `7.36%`;
- `T34_theory_baseline`: MAPE about `8.11%`;
- best new Test38-style models: about `9%` MAPE;
- new models had lower bias than some old corrections but did not beat Test34
  globally.

Interpretation:

- the new q-spectrum family transfers reasonably to k-Wave;
- Test34 remains a very strong k-Wave-specific reference;
- domain shift is real.

### Test 41: Mixedness-Corrected Latest k-Wave

Test 41 trained a Test33-style correction on k-Wave/Test40 outputs, keeping
base q models frozen.

Representative result:

- k-Wave performance improved strongly when training and testing within the
  k-Wave family;
- this made Test41 look very good on k-Wave-like data.

Important caveat:

- this is not proof of general correction transfer;
- it is partly domain-specific.

### Test 42: k-Wave-Trained Corrector on Test39 OOD

Test 42 transferred the k-Wave-trained correction to Test39 OOD simulations.

Finding:

- the correction helped some strongly mixed hard/interface subsets;
- it degraded the already strong OOD baseline globally;
- the corrector was too domain-sensitive.

Conclusion:

```text
k-Wave-trained correction is diagnostic, not deployable as a global default.
```

### Test 43: Synthetic-Trained Composition-Aware Correction Transfer

Test 43 trained only residual correction layers on synthetic Test38 rows,
keeping base q and composition models frozen.

Medium-mode outcome:

#### Test39 external

Best models were still uncorrected baselines:

| Model | MAPE | Signed error | High-error >20% |
|---|---:|---:|---:|
| `q_spectrum_only` | about `5.44%` | about `-0.01%` | about `4.47%` |
| `q_spectrum_plus_composition` | about `5.46%` | about `-0.05%` | about `4.47%` |
| best T43 mixedness-weighted correction | about `5.91%` | about `-0.13%` | about `4.08%` |

Interpretation:

- correction reduced high-error rate slightly;
- it increased MAPE;
- baseline remains preferable for general OOD synthetic validation.

#### k-Wave transfer

Representative results:

| Model | MAPE | Signed error | High-error >20% |
|---|---:|---:|---:|
| `T34_mixedness_logk_corrected` | about `7.36%` | about `-4.21%` | about `14.01%` |
| `T34_theory_baseline` | about `8.11%` | about `-0.60%` | about `10.42%` |
| best T43/new correction | about `8.7-8.8%` | about `-1.6%` to `-2.2%` | about `11-12%` |
| `q_spectrum_plus_composition` baseline | about `8.99%` | about `-1.12%` | about `11.39%` |

Important regional finding:

- T43-style correction can improve hard/interface regions;
- the same correction can strongly degrade soft/background regions.

Interpretation:

```text
The correction signal is real, but the current gate does not know enough about
correction direction.
```

## 9. Current Best Model Choice

If one model must be chosen now, independent of simulation type, the safest
choice is:

```text
q_spectrum_plus_composition
```

Why:

- strong and stable on Test39 external simulations;
- reasonable on k-Wave transfer;
- uses predicted composition information without relying on a brittle
  post-hoc correction;
- avoids choosing different estimators depending on whether the data are
  synthetic or k-Wave.

Second-best conservative option:

```text
q_spectrum_only
```

It is slightly simpler and can be very competitive, but the composition-aware
version is more informative for reliability and future correction.

## 10. What We Know Scientifically

### Interface contamination is real

Mixed patches often contain spectral energy from both soft and hard materials.
This can shift the effective selected wavenumber and produce SWS bias.

### Contamination is not the whole story

Pure heterogeneous patches can still have elevated error compared with
homogeneous maps. This suggests domain/context effects beyond simple patch
mixture.

### Theory maps are powerful but rigid

`TheoryQDiscrete` and Test30-style region levels can be extremely accurate for
clean two-material maps, but they may fail if the real material structure is
not piecewise constant or if the theory-inferred structure is wrong.

### M=2 remains safer near interfaces

Larger windows can make inclusion/hard-region underestimation worse, likely
because they mix more interface information. M=2 often preserves regions better.

### Correction needs direction

A useful future corrector must decide:

```text
should correct?
if yes, should SWS go up or down?
how much?
how reliable is the corrected value?
```

The current correction families sometimes know that a patch is problematic,
but they do not always know the correct sign of the correction.

## 11. Recommended Next Step

The next clean experiment should be a direction-aware correction test.

Suggested Test 44:

```text
Test 44: direction-aware correction and reliability for q_spectrum_plus_composition
```

Core idea:

- use `q_spectrum_plus_composition` as the frozen base estimator;
- train a correction decision layer, not a replacement estimator;
- predict both correction usefulness and correction direction;
- apply correction only when expected gain is high and harm risk is low;
- output both corrected SWS and final reliability.

Candidate operational predictors:

- predicted patch purity;
- `p_mixed`, `p_strong_mixed`;
- q/SWS disagreement with theory;
- spectral peak asymmetry;
- high-k vs low-k energy ratios;
- spectral width;
- local SWS gradients;
- local consistency of neighboring predictions;
- M/frequency/regime metadata.

Variables that must remain evaluation-only:

- true SWS;
- true material side;
- true patch purity;
- distance-to-interface;
- signed error;
- high-error labels.

Success criteria:

- no degradation in homogeneous maps;
- no degradation in pure/near-pure heterogeneous regions;
- improvement in mixed/interface regions;
- improvement on hard-side underestimation without damaging soft/background;
- stable transfer across Test39 OOD and k-Wave.

## 12. Important Artifacts

### Documentation

- `docs/TESTS_INDEX.md`
- `docs/TESTS_22_34_CONFIDENCE_AND_CORRECTION_SUMMARY.md`
- `docs/LATEST_MODEL_REGISTRY.md`
- `docs/TEST34_KWAVE_MIXEDNESS_TRANSFER.md`

### Current recommended model registry

```text
outputs/model_registry/test38_velocity_field_diverse_q_training/
test38__velocity_field_diverse_q__medium_bundle.mat
```

### Recent validation outputs

```text
outputs/test_39_frozen_test38_external_validation/
outputs/test_40_kwave_latest_model_comparison/
outputs/test_41_mixedness_corrected_latest_kwave/
outputs/test_42_external_validation_mixedness_corrector/
outputs/test_43_synthetic_trained_composition_correction_transfer/
```

## 13. Practical Recommendation Right Now

For the next report/figure set, present:

1. `q_spectrum_plus_composition` as the main current estimator.
2. `q_spectrum_only` as the simpler ablation.
3. `T34_mixedness_logk_corrected` as the strong k-Wave-specific reference.
4. `TheoryQDiscrete` / Test30 as structural diagnostic baselines, not the
   deployable solution.
5. Mixedness/composition maps as reliability information.
6. Direction-aware correction as the next research step, not yet solved.

