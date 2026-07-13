# Paper Closure Plan: Robust Learned q/SWS Estimator

## Core Claim

The v1 paper should focus on a robust learned REQ quantile estimator for SWS
mapping. The main model family is the clean spectral/composition branch from
Test 38, with Test 53 defining the final controlled retraining setup for the
paper.

The paper should not claim that confidence-gated correction fully solves
interfaces. Instead, interface degradation is a characterized limitation, and
confidence/SNR maps are used to decide where reporting SWS is reliable.

## Publication Model

- Main baseline: `q_spectrum_plus_composition`.
- Simple baseline: `q_spectrum_only`.
- Diagnostic reference: `q_spectrum_plus_theory_composition` and
  `theory_discrete`, but the paper should avoid making Theory a required
  operational dependency unless its advantage is decisive.
- Final clean retraining: Test 53.
- Realistic readout/SNR validation: Test 54.

## Final Training Design

Test 53 pins the paper-facing clean design:

- frequencies: `200, 300, 400, 500, 600 Hz`;
- REQ M: `2, 3`;
- primary resolution: `dx = dz = 0.2 mm`;
- clean synthetic fields only;
- no confidence maps, previous q models, true SWS, true patch purity, material
  side, or error labels as predictors;
- grouped/condition-level held-out evaluation inherited from the Test 38
  engine.

The 200 Hz cases are included for low-frequency/liver-like relevance. They
should be interpreted carefully, especially in heterogeneous maps where the
REQ window spans more physical wavelength cycles.

## Realistic Readout/SNR Validation

Test 54 uses `/Users/sara/Documents/wave_sim_project` as a validation source,
not a training source. It evaluates how much the clean model degrades when the
field includes ultrasound-like effects:

- clean: no readout noise and no acoustic shadow;
- moderate realistic: shear attenuation, spatially correlated readout noise,
  depth/SNR weighting;
- hard realistic: stronger readout noise, acoustic shadowing, backscatter
  contrast, and material-dependent attenuation.

The key metric is not only absolute MAPE, but the degradation relative to the
matched clean condition.

## Figures to Prioritize

1. Model ranking, excluding Theory from the main ranking.
2. MAPE vs frequency for M=2 and M=3.
3. MAPE vs M and geometry family.
4. Core/interface ROI metrics.
5. Dense representative SWS maps with absolute and signed error.
6. q vs omega/aperture sweep to show that the learned model follows expected
   physical trends.
7. Clean vs moderate/hard readout validation maps, including phase and
   amplitude QC.
8. Confidence/SNR estimability maps showing where SWS should or should not be
   reported.

## What Stays Out of the Main Claim

- Full correction/regularization of interface artifacts.
- FEM/breast-like validation unless it becomes clean quickly.
- Training directly on realistic readout/noisy simulations, unless the clean
  model clearly fails under moderate readout.
- Deployment-level liver claims without phantom or realistic validation.

## Recommended Narrative

1. Classical/theory q works in narrow regimes but struggles across field
   structure, geometry, and frequency.
2. A learned spectral q estimator trained on clean but diverse simulations
   improves generalization.
3. Composition-aware features help mixed/heterogeneous patches, but immediate
   interfaces remain the hardest case.
4. Realistic ultrasound readout introduces an estimability problem, not only an
   estimation problem.
5. The practical output should be SWS plus reliability/SNR mask.
