# Tests 22-34: Confidence, Interface Contamination, and Correction Summary

This document summarizes the current diagnostic/correction sequence for
adaptive REQ q/SWS maps, especially around heterogeneous interfaces.

## High-Level Story

The project moved from frozen confidence validation to interface diagnostics,
then to several correction families. The core finding so far is:

- Homogeneous maps are generally low error.
- Heterogeneous pure patches are harder than homogeneous patches, but still
  much better than mixed/interface patches.
- The largest failures remain concentrated near interfaces and in hard regions,
  especially inclusion hard material.
- `TheoryQDiscrete` / Test 30 region-level maps are very strong for clean
  piecewise-constant 2/3 m/s synthetic maps, but they are structurally rigid.
- Mixedness-aware continuous log-k correction from Test 33 is the most promising
  flexible correction, especially for hard-side underestimation.
- Simple q-candidate selection is not enough; choosing among existing q maps
  does not solve the hard-region bias.
- k-Wave transfer is harder than synthetic external maps. The Test 33 log-k
  corrector improves global k-Wave MAPE, but some regimes/M values still show
  strong hard-side underestimation.

## Test 22: Frozen Confidence External Validation

Purpose:
Validate frozen confidence detectors trained earlier, without retraining, over
external synthetic conditions.

Key design:

- Homogeneous 2/3 m/s, bilayer 2/3, circular inclusion 2/3.
- Field regimes: directional 2D, diffuse 2D, partial 3D, diffuse 3D.
- Frequencies and M/dx sweeps with pilot/medium/full modes.
- Main detector: `ML_bagged_trees`, plus other rule/ML detectors.

Outcome:

- Confidence was informative but not a correction by itself.
- Confidence maps helped identify interface/mixed regions and large-error
  zones.
- Runtime was optimized with physical target step and mode controls.

## Tests 23-25: Interface Contamination and Spectral Diagnostics

Purpose:
Diagnose whether interface errors come from mixed patches containing two
wavenumber radii.

Key outputs:

- Signed distance to interface.
- Patch purity / mixedness.
- Error by soft/hard side.
- Spectral mixture diagnostics and representative radial spectra.

Findings:

- Large errors increase strongly as patch purity drops.
- Hard-side SWS underestimation is common near mixed patches.
- However, pure heterogeneous patches can still have nontrivial error, so the
  problem is not exclusively spectral two-radius contamination.
- In pure patches, errors likely involve model/domain/context bias, not simple
  patch contamination.

## Test 26: Confidence-Gated Corrections

Purpose:
Try simple correction strategies without retraining q models.

Strategies:

- Physical donut.
- Prior donut.
- Confidence-gated prior donut.
- Confidence-gated interpolation.
- Experimental peak-candidate selection.

Finding:

- Gated interpolation improved some metrics but did not solve mixed/interface
  errors enough.
- Some strategies smoothed or distorted interfaces.
- This motivated edge-aware and structure-aware strategies.

## Test 27: Adaptive Window and Edge-Aware Correction

Purpose:
Compare Local, Hybrid, confidence blending, smaller windows, and edge-aware
interpolation.

Findings:

- `LocalOnly_T18` often beat Hybrid in pure/near-pure heterogeneous regions.
- Hybrid was not automatically better in heterogeneous maps.
- Edge-aware interpolation helped somewhat but still left large mixed-interface
  errors.
- M=2 consistently looked safer than M=3 for heterogeneous interfaces.

## Test 28-29: Oracle / Half-Window / Graph-TV Ideas

Purpose:
Probe upper bounds and alternative operational corrections.

Findings:

- Same-material/oracle-style references showed there is room to improve if the
  correct side of the interface can be identified.
- Half-window and graph/TV ideas helped diagnostically but did not yet become
  a robust deployable solution.
- Additional representative maps and ROI error bars made clear that the
  inclusion hard region remained problematic.

## Test 30: Theory Structure + Local M2 Region Levels

Purpose:
Reverse the usual order: use TheoryQDiscrete to infer a simple structure, then
use Local M=2 levels within regions.

Result on full synthetic cached design:

- `userguess_region_levels` global MAPE: about `1.67%`.
- Homogeneous MAPE: about `0.46%`.
- Pure MAPE: about `0.72%`.
- Near-pure MAPE: about `1.75%`.
- Moderate mixed MAPE: about `2.46%`.
- Strong mixed MAPE: about `9.67%`.

Interpretation:

- Very strong for clean two-material synthetic maps.
- Rigid: likely risky for more realistic continuous or complex tissue.
- Excellent as a structural baseline / upper-ish reference, not a complete
  deployable solution by itself.

## Test 31: Simple Confidence-Based q/SWS Interpolation

Purpose:
Compare low-confidence q interpolation/copying vs SWS interpolation/copying,
and compare with Test 30.

Full synthetic results:

- Local MAPE: `4.69%`.
- SWS nearest high-confidence: `3.71%`.
- q median global: `4.10%`.
- Test30 region levels: `1.75%`.
- Edge-aware SWS: `3.76%`.
- Edge-aware q: `4.23%`.

Finding:

- Simple interpolation helps over Local but cannot match structure-based region
  levels.
- q interpolation was not clearly better than SWS interpolation.

## Test 32: Structural-Confidence Hybrid

Purpose:
Blend the flexible Local/SWS-nearest maps with the rigid Test30 region-level
prior using operational weights.

Full synthetic highlights:

- Local MAPE: `4.69%`.
- SWS nearest MAPE: `3.71%`.
- Test30 MAPE: `1.75%`.
- `hybrid_lowconf_region_else_sws_nearest`: `2.73%`.
- `hybrid_relaxed_region_blend`: `2.69%`.
- `hybrid_boundary_protected_region`: `2.73%`, with lowest >20% among several
  hybrids.

Geometry:

- Bilayer: Test30 remained best.
- Inclusion: Test30 remained strong, but hybrids reduced some rigidity concerns.

Finding:

- Hybrids improved over flexible baselines but did not beat Test30 on clean
  synthetic 2/3 maps.
- They provided useful ingredients for a less rigid correction stack.

## Test 33: Mixedness-Aware q/log-k Correction

Purpose:
Train an operational mixedness detector, patch-purity regressor, continuous
log-k residual corrector, q-candidate selector, and posterior reliability model.

Important:

- Frozen q/confidence models were not retrained.
- Mixedness and truth variables were labels/evaluation only.
- Inference used operational features: confidence, Local/SWS-nearest/Test30
  predicted maps, disagreement, estimated boundary distance, structure
  agreement, M/frequency/regime/geometry family.

Full synthetic held-out caveat:

- Full predictions include all regimes.
- Held-out split in `*_test.csv` contains directional 2D and partial 3D;
  diffuse 2D/3D are in train for this split and should be read as diagnostic
  when using patch-level all-row summaries.

Held-out results:

- Local MAPE: `4.10%`.
- Test30 region levels: `1.56%`.
- Boundary-protected hybrid: `2.49%`.
- Mixedness log-k corrected: `2.22%`.
- Mixedness q selector: `3.30%`.

Mixedness detector:

- Mixed `<0.95`: ROC AUC `0.990`, PR AUC `0.975`.
- Strong mixed `<0.75`: ROC AUC `0.989`, PR AUC `0.948`.
- Patch purity MAE: `0.026`.

Key interpretation:

- ML can learn mixed vs pure well from operational features.
- The continuous log-k correction strongly improves hard-side bias, especially
  inclusion hard regions.
- The q-candidate selector is not enough: choosing among existing q maps does
  not solve the hard-region bias.

Full held-out hard-side examples:

- Inclusion hard MAPE:
  - Local: `18.29%`.
  - Test30: `7.76%`.
  - Mixedness log-k: `3.99%`.
- Bilayer hard MAPE:
  - Local: `6.82%`.
  - Test30: `2.18%`.
  - Mixedness log-k: `1.71%`.

Main caveat:

- Log-k correction can overcorrect soft regions. A future candidate should be
  side-aware or structure-aware.

## Test 34: k-Wave Mixedness Transfer

Purpose:
Apply the frozen Test 33 stack to cached k-Wave inclusion simulations.

Inputs:

- Cached Test 12 Level 06 k-Wave `CASE_OUT`.
- Frozen Test 19 q models.
- Frozen Test 21 `ML_bagged_trees` confidence detector.
- Frozen Test 33 mixedness/log-k/reliability models.

Cases:

- Directional 2D.
- Diffuse 2D.
- Projected diffuse 3D.
- 3D-rev / partial 3D.
- M = 2, 3, 4.

Important implementation note:

- The cached k-Wave `q_theory_discrete` column was NaN, so Test34 recomputes
  TheoryQDiscrete using the analytic theory function, matching the Test22
  convention.
- Test33 was trained on M=2/3; M=4 k-Wave transfer is extrapolation.
- The q-candidate selector in Test34 is a proxy/fallback because Test31 q
  candidate maps are not available for k-Wave.

Overall k-Wave results:

- Local MAPE: `8.98%`, high-error >20: `15.59%`.
- Hybrid baseline MAPE: `10.02%`, high-error >20: `15.85%`.
- Theory baseline MAPE: `8.11%`, high-error >20: `10.42%`.
- Test30-style region levels MAPE: `9.00%`, high-error >20: `13.39%`.
- Boundary-protected hybrid MAPE: `9.08%`, high-error >20: `14.57%`.
- Mixedness log-k corrected MAPE: `7.36%`, high-error >20: `14.01%`.

By regime:

- Directional 2D:
  - Local `7.67%`.
  - Theory `5.97%`.
  - Mixedness log-k `5.13%`.
- Diffuse 2D:
  - Local `10.36%`.
  - Theory `5.27%`.
  - Mixedness log-k `8.04%`.
- Projected diffuse 3D:
  - Local `9.19%`.
  - Theory `10.21%`.
  - Mixedness log-k `8.41%`.
- 3D-rev / partial 3D:
  - Local `8.69%`.
  - Theory `10.97%`.
  - Mixedness log-k `7.85%`.

By M:

- M=2:
  - Local `6.90%`.
  - Theory `6.80%`.
  - Mixedness log-k `5.08%`.
- M=3:
  - Local `9.24%`.
  - Theory `8.11%`.
  - Mixedness log-k `7.51%`.
- M=4:
  - Local `13.01%`.
  - Theory `10.89%`.
  - Mixedness log-k `11.99%`.

Interpretation:

- k-Wave confirms M=2 is safest.
- Mixedness log-k transfer improves global MAPE and bias, especially for
  directional 2D and M=2.
- It does not fully solve hard-side underestimation in k-Wave, especially for
  projected diffuse/partial 3D and larger M.
- Theory baseline is surprisingly strong in some k-Wave regimes but fails in
  others, consistent with its rigidity/domain assumptions.

## Current Recommendation

For synthetic clean maps:

- Use Test30 region levels as the strongest piecewise-constant reference.
- Use Test33 mixedness/log-k as a flexible correction layer, especially for
  hard-side underestimation and inclusion cases.

For k-Wave transfer:

- Prefer M=2.
- Treat Test33 log-k correction as the best current learned transfer candidate,
  but not final.
- Do not rely on the q-candidate selector.
- Add a side-aware / structure-aware correction that protects soft regions and
  targets hard-side underestimation explicitly.
- Consider training a kWave-aware or spectrum-aware residual model using
  spectral features and mixedness labels, but keep a strict validation split by
  k-Wave case/regime.

