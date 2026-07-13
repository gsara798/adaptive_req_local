# Adaptive REQ Test Index

This file is the compact project-wide index of tests implemented or prepared
through Test 27. A script being listed here means that its code exists; it does
not imply that every full, computationally expensive analysis has been run.

## Software and workflow checks

| Test | Purpose | Main entry point |
|---|---|---|
| 00 | Default configuration construction | `tests/integration/test_00_config.m` |
| 01 | One-field simulation engine and patch geometry | `tests/integration/test_01_simulation_engine.m` |
| 02 | Aperture-sweep integration | `tests/integration/test_02_aperture_sweep.m` |
| 03 | Config-driven aperture sweep | `tests/integration/test_03_config_driven_aperture_sweep.m` |
| 04 | Monte Carlo configuration preview | `tests/integration/test_04_mc_config_preview.m` |
| 05 | Small Monte Carlo sweep | `tests/integration/test_05_mc_sweep.m` |

Unit coverage currently includes REQ mappings, angular-shape features,
cumulative-energy features, and wave-direction plane coverage under
`tests/unit/`.

## Scientific experiments and analyses

| Test | Scientific purpose | Main files / status |
|---|---|---|
| 06 | Baseline feature-to-q associations and controlled diagnostics | `configs/test_06_feature_q_baseline.m`; analysis Levels 01-11 exist |
| 07 | Feature-q dataset and first trained-q-to-SWS evaluation | `experiments/run_test_07_feature_q_dataset.m`; analysis Levels 12-13 exist |
| 08 | Advanced angular, Ecum, plane-coverage, error-learning, and model-comparison audit | `experiments/run_test_08_advanced_angular_features.m`; analysis Levels 14-17 exist |
| 09 | Homogeneous stress robustness | `experiments/run_test_09_stress_robustness.m`; Level 01 analysis exists |
| 10 | Heterogeneous bilayer/inclusion maps | `experiments/run_test_10_heterogeneous_maps.m` |
| 11 | Local-only, global-only, and hybrid local-global q models; bilayer transfer and grouped generalization | `experiments/run_test_11_global_req_features.m`; analysis Levels 01-05 exist |
| 12 | `cs_guess` / effective-window sweep, model comparison, generalization, maps, SNR, and k-Wave transfer | `experiments/run_test_12_cs_guess_window_sweep.m`; analysis Levels 01-07 exist |
| 13 | Reserved in the scientific sequence; no standalone Test 13 driver exists | SWS metric work appears as Test 07 Level 13 |
| 14 | Spatial-resolution sensitivity | `experiments/run_test_14_dx_dz_resolution_sweep.m`; Level 01 analysis exists |
| 15 | Theory-informed direct-q model | `experiments/analysis/analyze_test_15_theory_informed_direct_q.m` |
| 16 | Theory residual-q model | `experiments/analysis/analyze_test_16_theory_residual_q.m` |
| 17 | Synthetic heterogeneous external cases and model comparison | `experiments/analysis/analyze_test_17_synthetic_inclusion_kWave_like.m`; `analyze_test_17_model_comparison_heterogeneous_cases.m` |
| 18 | Clean, homogeneous, noise-free field-regime training dataset | `experiments/run_test_18_clean_field_regime_training_dataset.m` |
| 19 | Train the Test 18 operational q models and compare old/theory baselines | `experiments/analysis/analyze_test_19_train_clean_field_regime_models.m` |
| 20 | External heterogeneous validation and homogeneous aperture-q tracking; no training | `experiments/analysis/analyze_test_20_external_validation_and_aperture_q_tracking.m` |
| 21 | Train rule-based and ML confidence detectors for high SWS error | `experiments/analysis/analyze_test_21_interface_confidence_detector.m` |
| 22 | External validation of the frozen Test 21 confidence detectors over geometry, frequency, regime, M, and resolution; no training | `experiments/analysis/analyze_test_22_confidence_external_validation.m` |
| 23 | Signed-distance, patch-purity, and representative-spectrum diagnosis of interface contamination; no training or model inference | `experiments/analysis/analyze_test_23_interface_patch_contamination_step01.m` |
| 24 | Quantify interface spectral failure modes using k ratios, material-reference regions, purity/distance bins, and representative spectra | `experiments/analysis/analyze_test_24_interface_spectral_failure_modes.m` |
| 25 | Test the two-radius contamination hypothesis using matched hard/soft spectral bands in pure and mixed patches | `experiments/analysis/analyze_test_25_two_radius_spectral_mixture.m` |
| 26 | Compare frozen-q physical donuts, prior donuts, confidence-gated corrections, interpolation, and peak candidates without retraining | `experiments/analysis/analyze_test_26_confidence_gated_corrections.m` |
| 27 | Compare frozen Hybrid/Local switching, blending, adaptive small windows, and edge-aware low-confidence correction | `experiments/analysis/analyze_test_27_adaptive_window_edge_aware.m` |

Test 22 supports `ADAPTIVE_REQ_TEST22_SIZE=pilot|medium|full`. Its REQ map
stride is fixed in physical units within each profile rather than in pixels.
For quick pipeline debugging, `ADAPTIVE_REQ_TEST22_FAST_MODE=true` restricts
the pilot to the primary hybrid q model and three frozen confidence detectors.

Example MATLAB launches:

```matlab
setenv('ADAPTIVE_REQ_TEST22_SIZE','pilot');
setenv('ADAPTIVE_REQ_TEST22_FAST_MODE','true'); % optional debug mode
run('experiments/analysis/analyze_test_22_confidence_external_validation.m');
```

Use `medium` or `full` in the first line for the larger frozen validations.
Outputs are separated under `analysis/pilot_fast`, `analysis/pilot`,
`analysis/medium`, and the original `analysis` folder for `full`.

Test 23 Step 01 reads frozen Test 20 predictions, computes material purity
using the exact REQ window, and recomputes spectra only for ten representative
patches. Its diagnostic output is under
`outputs/test_23_interface_patch_contamination/analysis/level_23_interface_patch_contamination_step01/`.

Test 24 consumes the cached Test 23 joined table without retraining or rerunning
models. It writes patch-level k diagnostics, grouped summaries, and twelve
representative failure spectra under
`outputs/test_24_interface_spectral_failure_modes/`.

Test 25 extracts a stratified spectral sample from the saved Test 17 fields and
uses frozen Test 23 errors only for retrospective association. Its
`spectral_mixture_index` is the balanced energy in bands around `k_hard` and
`k_soft`; matched homogeneous fields with the same regime and M provide the
baseline needed to distinguish material contamination from ordinary spectral
broadening. Outputs are under
`outputs/test_25_two_radius_spectral_mixture/`.

Test 26 consumes the frozen Test 22 condition and field caches. It reuses the
two cached hybrid predictions directly; because Test 22 did not retain
LocalOnly/GlobalOnly predictions, it reconstructs their operational features
once per condition and checkpoints the resulting frozen-model predictions.
Radial donuts hold q fixed and change only the q-to-k inversion. Select the
run profile with `ADAPTIVE_REQ_TEST26_MODE=quick|full` (preferred); the legacy
`ADAPTIVE_REQ_TEST26_QUICK_MODE=true` remains supported. Validation-only mode uses
`ADAPTIVE_REQ_TEST26_VALIDATE_ONLY=true`. Correction maps use an approximately
physical stride of 1.0 mm in quick mode and 0.75 mm in full mode. Checkpoints
store single-precision prediction cubes; the patch CSV stores only six key
diagnostic strategies. `GlobalOnly_T18` and `ML_boosted_trees` are optional via
`ADAPTIVE_REQ_TEST26_INCLUDE_GLOBAL=true` and
`ADAPTIVE_REQ_TEST26_INCLUDE_BOOSTED=true`. By default, a complete eight-panel
map diagnostic is exported for every condition under
`figures/maps_by_condition_dx200um/<geometry>/<regime>/`; disable this only with
`ADAPTIVE_REQ_TEST26_SAVE_ALL_MAPS=false`. Outputs are under
`outputs/test_26_confidence_gated_corrections/` (or its `quick/` subfolder).
The default Test 26 design evaluates only `dx=0.2 mm` (`CORR.Dx=0.2e-3`), giving
16 quick conditions and 192 full conditions.

Test 27 consumes the compact Test 26 checkpoints and their spectral caches. It
builds operational Hybrid/Local switches and blends, aligns `M-1` as the
small-window candidate, and applies a two-cluster pseudo-segmentation for
edge-aware interpolation. Oracle material/purity information is used only by
the explicitly diagnostic `oracle_region_reference`. Select the profile with
`ADAPTIVE_REQ_TEST27_MODE=quick|full`, validate without analysis using
`ADAPTIVE_REQ_TEST27_VALIDATE_ONLY=true`, and control per-condition map export
with `ADAPTIVE_REQ_TEST27_SAVE_ALL_MAPS=true|false`. Outputs are under
`outputs/test_27_adaptive_window_edge_aware/` or its `quick/` subfolder.

Test 28 addresses mixed-interface patches before spatial interpolation. It
uses the measured radial spectrum to detect two separated components without
using expected material wavenumbers, assigns the low-k or high-k component
from an operational Local/Hybrid pseudo-segmentation, and optionally recomputes
REQ after a soft same-region mask is applied to the velocity patch. Both
corrections are confidence-gated and fall back to frozen LocalOnly predictions
when spectral evidence or numerical quality is insufficient. Select the
profile with `ADAPTIVE_REQ_TEST28_MODE=quick|full`, run isolated checks with
`ADAPTIVE_REQ_TEST28_VALIDATE_ONLY=true`, control field-level masked REQ with
`ADAPTIVE_REQ_TEST28_EDGE_MASKED=true|false`, and control map export with
`ADAPTIVE_REQ_TEST28_SAVE_ALL_MAPS=true|false`. Outputs are under
`outputs/test_28_two_component_edge_masked_req/` or its `quick/` subfolder.

Test 29 implements a staged go/no-go workflow for interface correction. Its
diagnostic-only oracle same-material mask measures the recoverable ceiling;
an operational bank of eight smoothly tapered half-windows searches for a
cleaner one-sided spectrum; and confidence-weighted bilateral graph/TV
reconstruction produces an edge-preserving continuous map. Run quick first
with `ADAPTIVE_REQ_TEST29_MODE=quick`. Full mode reads the saved quick gate and
is blocked unless an operational method improves mixed-patch MAPE by at least
0.5 points while keeping homogeneous and pure/near-pure degradation within
configured limits. `ADAPTIVE_REQ_TEST29_FORCE_FULL=true` explicitly overrides
that guard. Validation-only and map controls are
`ADAPTIVE_REQ_TEST29_VALIDATE_ONLY=true` and
`ADAPTIVE_REQ_TEST29_SAVE_ALL_MAPS=true|false`. Outputs are under
`outputs/test_29_oracle_halfwindow_graph_tv/` or its `quick/` subfolder.
The companion `analyze_test_29_extended_maps_and_rois.m` exports all 96
heterogeneous condition panels, a 24-case representative subset, and matched
6 x 6 mm ROI summaries for bilayer soft/hard and inclusion
background/core regions without rerunning REQ. `TheoryQDiscrete` is
reconstructed from the cached radial curves and the discrete theoretical q
for each frequency/regime/M condition, and is included as a diagnostic
baseline in Test 29 summary tables, maps, error panels, and ROI comparisons.

Test 30 reverses the usual adaptive-q ordering: `TheoryQDiscrete` at `M=2`
provides an initial structural segmentation, while frozen `LocalOnly_T18` at
`M=2` supplies robust quantitative region levels. It compares automatic
(`unknown`) and user-provided geometry families (`homogeneous`, `bilayer`, or
`inclusion`), region-level reconstruction, boundary-aware graph refinement,
and a diagnostic oracle reference. `M=3` and `M=4` are never used to replace
the output map; they only reduce reliability when their projected estimates
show multiscale instability or a monotonic SWS drop. The inclusion prior keeps
the largest connected component and fills holes, preventing isolated Theory
artifacts from becoming hard material. Run validation with
`ADAPTIVE_REQ_TEST30_VALIDATE_ONLY=true`, quick analysis with
`ADAPTIVE_REQ_TEST30_MODE=quick`, or the complete cached design with
`ADAPTIVE_REQ_TEST30_MODE=full`. Outputs, including every condition map, are
under `outputs/test_30_theory_structure_local_m2/`. The region-level methods
are explicitly intended for approximately piecewise-constant maps; their
performance is a topology/region-estimation result, not evidence that they
will preserve arbitrary continuous tissue gradients.

Test 31 is a simpler and more interpretable confidence-interpolation
diagnostic. It uses `dx=0.2 mm`, `M=2` and `M=3`, frozen `LocalOnly_T18` q/SWS,
frozen `ML_bagged_trees` confidence, cached REQ radial curves, and
`TheoryQDiscrete` to compare whether low-confidence pixels are better handled
by copying/interpolating q or by copying/interpolating SWS. The low-confidence
mask is fixed at `confidence < 0.8`. Strategies include nearest high-confidence
q/SWS copy, global and region-wise high-confidence q medians, local-structure
SWS interpolation, Theory-structure edge-aware q/SWS interpolation, and the
Test 30 region-level output as a structural baseline. No true SWS, material
label, patch purity, or interface distance is used for correction; these are
attached only for summaries, ROI plots, and interpretation. Run validation with
`ADAPTIVE_REQ_TEST31_VALIDATE_ONLY=true`, quick analysis with
`ADAPTIVE_REQ_TEST31_MODE=quick`, and the complete cached design with
`ADAPTIVE_REQ_TEST31_MODE=full`. Map export is controlled by
`ADAPTIVE_REQ_TEST31_SAVE_ALL_MAPS=true|false`. Outputs are under
`outputs/test_31_simple_confidence_interpolation/`, or its `quick/` subfolder,
including patch-level results, q-correction summaries, ROI SWS-vs-frequency
tables, strategy rankings, donor-distance maps, q correction maps, and one map
panel per condition when all-map export is enabled.

Test 32 builds a structural-confidence hybrid on top of Test 31 outputs. It
does not rerun REQ or train any model. Instead, it combines the flexible
`LocalOnly_T18` / high-confidence SWS-nearest maps with the stronger but more
rigid Test 30 Theory-structure region-level estimate. Hybrid weights are
operational only: confidence, estimated structure distance-to-boundary,
agreement between Local-derived and Theory/Test30-derived structures, and
Local-vs-region disagreement. Truth, patch purity, material side, and true
interface distance are used only for summaries. The main candidate relaxes
toward SWS-nearest near estimated edges and toward high-confidence Local when
the structural prior is less reliable, while allowing Test30-like levels in
low-confidence interiors. Run with `ADAPTIVE_REQ_TEST32_MODE=quick|full`,
validate with `ADAPTIVE_REQ_TEST32_VALIDATE_ONLY=true`, and control maps with
`ADAPTIVE_REQ_TEST32_SAVE_ALL_MAPS=true|false`. Outputs are under
`outputs/test_32_structural_confidence_hybrid/`.

Test 33 trains a mixedness-aware post-map layer on top of Test 32 compact
outputs. It does not retrain the frozen q estimators or confidence detectors.
Instead, it learns operational post-processing models from cached maps: a patch
mixedness detector (`patch_purity < 0.95` and `< 0.75` as labels), a patch-purity
regressor, a residual log-k corrector for the selected Test 32 base strategy,
and, when Test 31 q outputs are available, a mixedness-aware q-candidate
selector that chooses among previously computed q maps such as local q,
nearest high-confidence q, region median q, and edge-aware q. A posterior
reliability detector then estimates whether the corrected map is still likely
to have error >20%. Correction inputs are restricted to
operational quantities such as confidence, Local/SWS-nearest/Test30 predicted
maps, model disagreement, estimated boundary distance, structure agreement, M,
frequency, geometry family, and field regime. True SWS, patch purity, material
side, q_true, and true interface distance are used only as training labels or
evaluation variables. Run validation with `ADAPTIVE_REQ_TEST33_VALIDATE_ONLY=true`,
quick analysis with `ADAPTIVE_REQ_TEST33_MODE=quick`, and the complete cached
design with `ADAPTIVE_REQ_TEST33_MODE=full`. Outputs are under
`outputs/test_33_mixedness_aware_q_correction/`.

Test 34 transfers the frozen Test 33 mixedness/log-k correction stack to cached
k-Wave inclusion simulations from Test 12 Level 06 and compares the correction
layer against the cached Test 33 reference domain. It loads k-Wave REQ feature
tables, applies frozen Test 19 q models, frozen Test 21 `ML_bagged_trees`
confidence, and frozen Test 33 mixedness/posterior-reliability models. No model
is retrained. k-Wave truth is used only for ROI, soft/hard, distance, and map
evaluation. Outputs are under `outputs/test_34_kwave_mixedness_transfer/`,
including copied frozen models, patch-level predictions, ROI summaries,
soft/hard summaries, confidence/reliability summaries, and maps for all
k-Wave cases/M values.

Test 35 starts a clean modeling branch that does not use frozen Local/Hybrid
q models, confidence detectors, Test 30 structure maps, or previous correction
outputs as predictors. It rebuilds local REQ spectral/Ecum feature tables from
Test 22 field caches, learns a scalar dominant-material `patch_purity`
regressor plus mixed/strong-mixed classifiers directly from primitive spectral
features, then trains direct q/SWS models using only those primitive features,
physical parameters, optional analytic `TheoryQDiscrete`, and the predicted
composition outputs. The composition target is intentionally scalar rather
than soft/hard fractions, so it remains meaningful for windows containing
three or more materials. q labels are produced by inverting each patch's
`Ecum(k)` mapping at `k_true`; true SWS and patch purity are training/evaluation
labels only, not inference predictors. Run validation with
`ADAPTIVE_REQ_TEST35_VALIDATE_ONLY=true`, quick training with
`ADAPTIVE_REQ_TEST35_MODE=quick`, and full training with
`ADAPTIVE_REQ_TEST35_MODE=full`. k-Wave transfer is controlled by
`ADAPTIVE_REQ_TEST35_USE_KWAVE=true|false`, and representative map export by
`ADAPTIVE_REQ_TEST35_SAVE_MAPS=true|false`. Outputs are under
`outputs/test_35_spectral_composition_to_q_model/`, or its `quick/` subfolder,
including saved composition/q models, synthetic held-out summaries, purity-bin
summaries, k-Wave external summaries when enabled, and diagnostic figures.

Test 37 is an out-of-distribution validation of the frozen Test 35
`q_spectrum_only` model. It does not retrain. It generates new diagnostic
simulations outside the original 2/3 m/s training family, including unseen
homogeneous speeds, 2/4 and 1.5/3 contrasts, ellipse/off-center/two-inclusion
and oblique/thin/three-material geometries, unseen frequencies, new directional
angles, new diffuse seeds, and forced in-plane source alignment for every
regime. It applies the frozen Test 35 spectral q model, computes true patch
purity and boundary distance only for evaluation, saves every condition map,
and summarizes failures by case, family, regime, frequency, M, true purity bin,
and boundary-distance bin. Run quick diagnostics with
`ADAPTIVE_REQ_TEST37_MODE=quick`, full OOD validation with
`ADAPTIVE_REQ_TEST37_MODE=full`, and validation-only checks with
`ADAPTIVE_REQ_TEST37_VALIDATE_ONLY=true`. Outputs are under
`outputs/test_37_ood_q_spectrum_validation/`, or its `quick/` subfolder.

Test 38 retrains the clean spectral q family with broader velocity and field
diversity. It keeps the Test 35 philosophy: no frozen Local/Hybrid/confidence
models are predictors. Primitive REQ spectrum features plus physical metadata
are used to learn a scalar composition estimate (`predicted_patch_purity`,
`p_mixed`, `p_strong_mixed`), and the q/SWS models then compare spectrum-only,
theory-informed, composition-informed, theory+composition, `delta_q`, and
`delta_logk` variants. The training design includes multiple material speeds,
three directional angles in full mode, at least one diffuse 2D regime, three
random diffuse 3D seeds in full mode, and partial 3D source counts of 8/16/32
with at least one in-plane/aligned source. q labels are produced by inverting
the local `Ecum(k)` mapping at `k_true`; true SWS and true patch purity are
labels/evaluation variables only, never primitive predictors. Run lightweight
checks with `ADAPTIVE_REQ_TEST38_VALIDATE_ONLY=true`, quick training with
`ADAPTIVE_REQ_TEST38_MODE=quick`, medium training with
`ADAPTIVE_REQ_TEST38_MODE=medium` (all geometries/regimes, frequencies 350 and
550 Hz, M=2), and full training with `ADAPTIVE_REQ_TEST38_MODE=full`. The
current full design uses frequencies 350/450/550/650 Hz, M values
1.5/2/2.5/3/3.5/4, and spatial resolutions 0.15/0.2/0.25/0.3 mm by default.
Override resolutions with `ADAPTIVE_REQ_TEST38_DX_LIST_MM`, the physical REQ
center spacing with `ADAPTIVE_REQ_TEST38_TARGET_STEP_M`, and the post-extraction
row cap with `ADAPTIVE_REQ_TEST38_MAX_PATCHES_PER_CONDITION`; a cap of `0`
retains every valid StepX/StepZ grid center. Dense full runs can create several
million patches, so model fitting and held-out prediction tables are separately
subsampled by condition using `ADAPTIVE_REQ_TEST38_MAX_MODEL_TRAIN_ROWS` and
`ADAPTIVE_REQ_TEST38_MAX_MODEL_EVAL_ROWS` (`0` means no subsampling). Set
`ADAPTIVE_REQ_TEST38_SAVE_ALL_MAPS=true|false` to control per-condition map
export, and `ADAPTIVE_REQ_TEST38_MAP_INTERP_SCALE` to densify diagnostic
visualization.
Parallel extraction/training can be enabled with
`ADAPTIVE_REQ_TEST38_USE_PARFOR=true` and
`ADAPTIVE_REQ_TEST38_USE_PARALLEL_TRAINING=true`.
Outputs are under
`outputs/test_38_velocity_field_diverse_q_training/`, or its `quick/`
subfolder, including saved model bundles, patch-level predictions, grouped
summaries, composition diagnostics, model rankings, and per-condition maps.

Test 39 is the frozen external validation for the Test 38 model family. It
loads a saved Test 38 model bundle (`ADAPTIVE_REQ_TEST39_MODEL_SOURCE=medium`,
`full`, `quick`, or an explicit `.mat` path), never retrains, and applies all
frozen models (`q_spectrum_only`, theory/composition variants, `delta_q`,
`delta_logk`, and `theory_discrete`) to new OOD simulations. The OOD design
uses velocities, frequencies, directional angles, random diffuse seeds,
partial-3D source counts, and geometry parameters that differ from the Test 38
training runs. True SWS, true patch purity, and boundary distance remain
diagnostic/evaluation variables only. Run quick validation with
`ADAPTIVE_REQ_TEST39_MODE=quick`, medium validation with
`ADAPTIVE_REQ_TEST39_MODE=medium`, and full validation with
`ADAPTIVE_REQ_TEST39_MODE=full`. Optional controls include
`ADAPTIVE_REQ_TEST39_VALIDATE_ONLY=true`,
`ADAPTIVE_REQ_TEST39_SAVE_ALL_MAPS=true|false`,
`ADAPTIVE_REQ_TEST39_USE_PARFOR=true`, and
`ADAPTIVE_REQ_TEST39_MAP_STYLE=patch|interp`. Dense diagnostic map export is
controlled by `ADAPTIVE_REQ_TEST39_DENSE_MAPS=true|false` and
`ADAPTIVE_REQ_TEST39_DENSE_MAP_STEP_M`; this keeps the validation tables sparse
and fast while recomputing only figure maps on a denser, uncropped REQ grid.
The default `patch` style shows the REQ patch grid directly; `interp` is only
for smoother visual inspection and uses `ADAPTIVE_REQ_TEST39_MAP_INTERP_SCALE`.
Summary plots omit the diagnostic `theory_discrete` baseline so the learned
frozen models can be compared on a useful scale. Test 39 also writes
material-region, core/interface, frequency-region, and M-region summaries plus
region diagnostics for soft, hard, inclusion, background, and core ROIs.
Outputs are under
`outputs/test_39_frozen_test38_external_validation/`, or its `quick/` and
`medium/` subfolders.

Test 40 compares the k-Wave transfer results from Test 34 against the latest
frozen clean-q model bundles from Test 35 and Test 38. It does not train any
model. Instead, it loads the Test 34 k-Wave patch-level predictions, the Test
35 k-Wave REQ feature table, and the registered Test 35/Test 38 model bundles,
then evaluates all compatible learned q variants on the same k-Wave patches.
The resulting long table includes old strategies such as `local_baseline`,
`hybrid_baseline`, `test30_theory_region_levels`,
`mixedness_logk_corrected`, and `mixedness_q_candidate_selector`, plus the new
`q_spectrum_only`, theory-informed, composition-informed, theory+composition,
`delta_q`, and `delta_logk` variants from the registered bundles. No true SWS,
true composition, confidence, or correction output is used as a predictor by
the new frozen models. Test 40 writes overall, model-family, case/regime, M,
soft/hard, ROI, distance-to-interface, and analysis-region summaries, together
with ranking plots and per-condition maps. Run validation checks with
`ADAPTIVE_REQ_TEST40_VALIDATE_ONLY=true`. Control map export with
`ADAPTIVE_REQ_TEST40_SAVE_ALL_MAPS=true|false` and limit diagnostic maps with
`ADAPTIVE_REQ_TEST40_MAX_MAP_CONDITIONS`. Select the Test38 bundle with
`ADAPTIVE_REQ_TEST40_TEST38_SOURCE=quick|medium|full|/path`; the full bundle is
loaded from `outputs/test_38_velocity_field_diverse_q_training/models/` and can
take several minutes to load because it is a large frozen ensemble bundle. Set
`ADAPTIVE_REQ_TEST40_INCLUDE_TEST35=false` to compare only Test34/Test33-style
strategies against Test38, and restrict Test38 variants with
`ADAPTIVE_REQ_TEST40_TEST38_MODELS`, for example
`q_spectrum_only,q_spectrum_plus_composition,q_spectrum_plus_theory_composition`.
This script uses the cached Test35 k-Wave REQ feature table and converts
predicted q through the cached `req_mapping`; it does not re-run dense REQ on
the raw k-Wave fields. Outputs are written under
`outputs/test_40_kwave_latest_model_comparison/test38_<source>/` by default.

Test 41 applies a Test33-style mixedness-aware residual log-k correction to
the latest frozen clean-q models on the same k-Wave table produced by Test 40.
It keeps the base q models frozen and trains only a post-map corrector using
operational predictors: the frozen model SWS/q prediction, frozen confidence,
predicted composition/mixedness, M/f0/regime, and operational local/theory
reference predictions. True SWS is used only as the residual-correction label
and for evaluation; true material side, ROI, and distance-to-interface are
diagnostic grouping variables only. The script compares the corrected new
models against the old Test34 stack and writes Test33-style plots for
soft/hard errors, core ROI errors, M/regime dependence, distance-to-interface
behavior, and correction gain by region. Run quick mode with
`ADAPTIVE_REQ_TEST41_MODE=quick` and full mode with
`ADAPTIVE_REQ_TEST41_MODE=full`; use
`ADAPTIVE_REQ_TEST41_VALIDATE_ONLY=true`,
`ADAPTIVE_REQ_TEST41_SAVE_ALL_MAPS=true|false`, and
`ADAPTIVE_REQ_TEST41_MAX_MAP_CONDITIONS` to control checks and map export.
Outputs are under `outputs/test_41_mixedness_corrected_latest_kwave/`.

Test 42 externally validates the k-Wave-trained mixedness log-k corrector from
the Test 41 idea on the OOD simulation predictions produced by Test 39. The
base Test38 q models remain frozen. Test 42 trains only residual correctors on
the k-Wave/Test40 table, then applies those frozen correctors to Test39 quick,
medium, or full external tables (`ADAPTIVE_REQ_TEST42_TEST39_SOURCE`). To make
the validation honest, the corrector uses only predictors available in both
k-Wave and Test39 outputs: model SWS/q prediction, theory q/SWS reference,
predicted composition/mixedness, M/f0/regime, and model-theory disagreement.
It does not use true SWS, error, true patch purity, material side, ROI, or
distance-to-interface for inference. Test 42 writes overall, case/family,
regime, frequency, M, purity, material-region, region-zone, distance, gain,
and worst-condition tables plus Test33-style diagnostic plots and optional
maps. Parallel ensemble training can be enabled with
`ADAPTIVE_REQ_TEST42_USE_PARALLEL=true`. Quick results are intentionally a
stress test: the k-Wave-trained correction improves some strongly mixed
hard/interface subsets but degrades the already strong OOD baseline globally,
so it should be treated as a diagnostic correction rather than a deployable
default until a better gating/acceptance rule is developed. Outputs are under
`outputs/test_42_external_validation_mixedness_corrector/`, with `quick/` and
`medium/` subfolders when those modes are used.

Test 43 trains only composition-aware correction layers on synthetic Test38
training rows, while keeping both the Test38 composition/mixedness predictor
and the base q models frozen. It uses the already predicted composition
variables (`predicted_patch_purity`, `p_mixed`, `p_strong_mixed`) together
with operational q/SWS/theory disagreement features to learn residual
`delta_logk` corrections and an acceptance gate. The frozen correction stack
is evaluated on three domains: Test38 synthetic held-out rows, Test39 OOD
external simulations, and Test40 k-Wave transfer rows. No true SWS, true patch
purity, material labels, ROI labels, distance-to-interface, or error variables
are used at inference time. Runtime controls include
`ADAPTIVE_REQ_TEST43_MODE=quick|medium|full`,
`ADAPTIVE_REQ_TEST43_TRAIN_SOURCE`, `ADAPTIVE_REQ_TEST43_TEST39_SOURCE`,
`ADAPTIVE_REQ_TEST43_USE_PARALLEL=true|false`, and map export controls. Quick
results show the synthetic-trained correction improves synthetic held-out
MAPE, but the existing `q_spectrum_plus_composition` baseline remains best on
Test39 OOD and the old Test34 mixedness-logk correction remains best on the
current k-Wave transfer table. This suggests the learned residual correction
is still domain-sensitive and should not replace the frozen q baseline without
a stronger acceptance/gating policy. Outputs are under
`outputs/test_43_synthetic_trained_composition_correction_transfer/`, with
mode-specific subfolders.

Test 44 trains a conservative direction-aware correction and reliability layer
on top of the frozen `q_spectrum_plus_composition` estimator. The base q model
is not retrained. The script learns a direction classifier (`increase_sws`,
`decrease_sws`, `keep`), residual magnitude models, acceptance/harm gates, and
a posterior high-error reliability score using only operational predictors:
base q/SWS predictions, theory and q-spectrum disagreements, predicted
composition/mixedness, M/frequency/regime/geometry metadata, and local
variation proxies. True SWS, true purity, material side, ROI, and
distance-to-interface are labels or evaluation groups only and are explicitly
blocked from predictor tables. Strategies include the frozen baseline,
residual correction applied everywhere, direction-only diagnostics,
direction-gated correction, conservative direction-gated correction,
mixedness-only gate ablation, oracle apply-if-improves, and reliability-only.
Runtime controls include `ADAPTIVE_REQ_TEST44_MODE=quick|medium|full`,
`ADAPTIVE_REQ_TEST44_SOURCE=test39|kwave|both`,
`ADAPTIVE_REQ_TEST44_VALIDATE_ONLY=true|false`,
`ADAPTIVE_REQ_TEST44_SAVE_ALL_MAPS=true|false`, and
`ADAPTIVE_REQ_TEST44_USE_PARALLEL=true|false`. Outputs are under
`outputs/test_44_direction_aware_correction/`, with `quick/` and `medium/`
subfolders when those modes are used.

Test 45 evaluates a cleaner no-Theory correction layer on top of the frozen
`q_spectrum_plus_composition` estimator. It does not retrain the base q model
and does not use `TheoryQDiscrete`, true SWS, true material side, true purity,
or distance-to-interface as operational predictors. Instead it builds a small
set of deployable SWS candidates from the base model, `q_spectrum_only`, a
pseudo-region median, a same-region high-quality neighbor median, and a
base/region blend driven by predicted composition. A candidate selector,
conservative selector, mixedness-gated selector, oracle best-candidate upper
bound, and reliability-only baseline are compared across synthetic held-out,
Test39 OOD, and k-Wave/Test40 rows when available. The maps explicitly show
candidate SWS maps, pseudo-regions, correction masks, correction magnitude,
predicted mixedness, and final reliability so geometry preservation can be
inspected visually. It also writes a comparison table/figure against Test33
and Test34 reference summaries. Runtime controls include
`ADAPTIVE_REQ_TEST45_MODE=quick|medium|full`,
`ADAPTIVE_REQ_TEST45_SOURCE=test39|kwave|both`,
`ADAPTIVE_REQ_TEST45_VALIDATE_ONLY=true|false`,
`ADAPTIVE_REQ_TEST45_SAVE_ALL_MAPS=true|false`,
`ADAPTIVE_REQ_TEST45_USE_PARALLEL=true|false`,
`ADAPTIVE_REQ_TEST45_MAP_STYLE=patch|interp`, and
`ADAPTIVE_REQ_TEST45_MAX_MAP_CONDITIONS`. Outputs are under
`outputs/test_45_candidate_aware_no_theory_correction/`, with mode-specific
subfolders.

Test 46 returns to the smoother correction pattern that worked best in Test33.
It keeps the frozen `q_spectrum_plus_composition` estimator as the base map
and trains only a residual `delta log-k` correction layer using operational
features. The correction is blended by predicted mixedness, so pure-like and
homogeneous regions keep the base estimate while mixed/interface regions can
receive a continuous log-k adjustment. Unlike Test45, Test46 deliberately does
not construct pseudo-regions, does not use region medians, does not use nearest
neighbor interpolation, and does not use `TheoryQDiscrete` as an input. Maps
are saved in a Test33-like patch-grid layout to make geometry preservation
visible. Strategies include the frozen base, `q_spectrum_only` reference,
residual log-k applied everywhere, mixedness-gated log-k, conservative
mixedness-gated log-k, posterior reliability only, and an oracle delta-log-k
upper bound. Runtime controls include
`ADAPTIVE_REQ_TEST46_MODE=quick|medium|full`,
`ADAPTIVE_REQ_TEST46_SOURCE=test39|kwave|both`,
`ADAPTIVE_REQ_TEST46_VALIDATE_ONLY=true|false`,
`ADAPTIVE_REQ_TEST46_SAVE_ALL_MAPS=true|false`,
`ADAPTIVE_REQ_TEST46_USE_PARALLEL=true|false`,
`ADAPTIVE_REQ_TEST46_MAX_TRAIN_ROWS`, and
`ADAPTIVE_REQ_TEST46_MAX_MAP_CONDITIONS`. Outputs are under
`outputs/test_46_test33_style_logk_correction/`, with mode-specific
subfolders.

Test 47 is a post-hoc diagnostic analysis of the saved Test38 full prediction
table. It does not run REQ and does not train anything. The script reads
`test38_patch_level_predictions.csv`, keeps `theory_discrete` out of the main
rankings, and focuses on `q_spectrum_plus_composition`, `q_spectrum_only`, and
`q_spectrum_plus_theory_composition`. It writes grouped summaries by
frequency, M, dx, regime, case family, case id, purity bin,
distance-to-interface bin, and ROI/core/interface zone. ROI definitions are
diagnostic only: homogeneous center ROI is an 8 mm x 8 mm square, soft core is
soft material farther than 8 mm from an interface, hard core is hard material
farther than 4 mm from an interface, and interface bands are 0-0.5, 0.5-1,
1-2, 2-4, and >4 mm. Figures include model ranking without Theory, MAPE versus
frequency/M/dx, family/case failures, distance-to-interface curves, ROI error
summaries, and overlay diagrams showing where ROIs are located. Run with
`ADAPTIVE_REQ_TEST47_MODE=full` for the full saved CSV or
`ADAPTIVE_REQ_TEST47_MODE=quick` for a sampled check. Outputs are under
`outputs/test_47_test38_full_results_diagnostics/`, with a `quick/` subfolder
when quick mode is used.

Test 48 recalculates dense REQ maps for representative validation cases using
the frozen Test38 model bundle. It does not train or modify the saved models.
The script is intended for visual/ROI diagnosis rather than a combinatorial
benchmark: quick mode evaluates eight representative seen and unseen cases
with M=2 and M=3; full mode expands to the representative case list across
directional 15 degrees, new diffuse 2D/3D seeds, and partial 3D with 12
sources. Cases include homogeneous 2/3 m/s, bilayer 2/3, circular inclusion
2/3, thin layer 2/4, three material 2/3/4, unseen bilayer 2.25/3.75, unseen
inclusion 2.25/3.75, and an off-center ellipse 2/4. Maps show true SWS,
predicted SWS, absolute/signed error, predicted q, predicted patch purity,
predicted mixedness, true patch purity, distance-to-interface, and ROI
overlays. ROI tables report MAPE, signed error, mean/std SWS, N, and ROI size.
Controls include `ADAPTIVE_REQ_TEST48_MODE=quick|full`,
`ADAPTIVE_REQ_TEST48_MODEL_SOURCE=quick|medium|full|/path`,
`ADAPTIVE_REQ_TEST48_TARGET_STEP_M`, `ADAPTIVE_REQ_TEST48_SAVE_ALL_MAPS`, and
`ADAPTIVE_REQ_TEST48_VALIDATE_ONLY`. By default the dense REQ extraction now
uses `ADAPTIVE_REQ_TEST48_REQ_PROFILE=test38_training`, matching Test38
training settings (`Nbins='auto'`, `Nbins_auto_oversample=1`,
`Nbins_min=16`). Use `ADAPTIVE_REQ_TEST48_REQ_PROFILE=dense_default` to recover
the earlier denser-bin Test48 behavior, or override explicitly with
`ADAPTIVE_REQ_TEST48_NBINS_OVERSAMPLE` and `ADAPTIVE_REQ_TEST48_NBINS_MIN`.
Set `ADAPTIVE_REQ_TEST48_COMPARE_TEST33=true` to append cached Test33 reference
summaries for `local_baseline`, `mixedness_logk_corrected`, and
`mixedness_q_candidate_selector`; these are clearly marked as cached Test33
domain results because original Test33 strategies depend on Test31/Test32
intermediate maps and are not regenerated from the dense Test48 fields.
Quick/validate defaults to the medium model bundle to avoid repeatedly loading
the 5.3 GB full bundle; full mode defaults to the full bundle. Validation-only
runs write under `quick/validate/` or `validate/` so they do not overwrite
analysis summaries. Outputs are under
`outputs/test_48_test38_dense_req_validation_maps/`, with a `quick/` subfolder
when quick mode is used.

Test 49 performs the analogous dense validation on the raw k-Wave fields. It
loads `data/k-wave/<case>/data_500Hz.mat`, re-extracts REQ/features with the
Test38 training profile (`Nbins='auto'`, `Nbins_auto_oversample=1`,
`Nbins_min=16`, `smooth_sigma=1`, `gamma_win=1`, `pad_factor=1`, and
`EdgeMode='valid'`), and applies frozen Test38 q models through the newly
computed local `req_mapping`. No q model is retrained, and k-Wave truth is used
only for q-oracle/error/ROI diagnostics. Quick mode evaluates Directional 2D
and Projected diffuse 3D at M=2; medium mode evaluates three k-Wave regimes at
M=2/3; full mode evaluates Directional 2D, Diffuse 2D, Projected diffuse 3D,
and 3D-rev at M=2/3/4. Controls include
`ADAPTIVE_REQ_TEST49_MODE=quick|medium|full`,
`ADAPTIVE_REQ_TEST49_MODEL_SOURCE=quick|medium|full|/path`,
`ADAPTIVE_REQ_TEST49_TARGET_STEP_M` (default 0.5 mm),
`ADAPTIVE_REQ_TEST49_USE_PARFOR=true|false`,
`ADAPTIVE_REQ_TEST49_SAVE_ALL_MAPS=true|false`, and
`ADAPTIVE_REQ_TEST49_MODELS`. Feature caches are stored per condition under
`outputs/test_49_kwave_dense_req_test38_validation/<mode>/data/condition_cache`
for quick/medium and under the root Test49 output for full mode. Outputs
include patch-level predictions, summaries by case/M/ROI/soft-hard/distance,
and dense map panels with true SWS, theory SWS, predicted purity/mixedness,
model SWS maps, absolute error maps, distance-to-interface, and hard-region
mask.

Test 50 starts a clean, project-native k-Wave phase 1. Reusable helpers live
under `src/+adaptive_req/+kwave/`: toolbox discovery, controlled material maps,
2D elastic simulation, and harmonic phasor extraction. The runnable entry point
is `tests/integration/test_50_controlled_kwave_phase1.m`, not an experiment
script, so the simulation itself stays in the test layer. It does not move the
k-Wave installation; it validates and adds the toolbox root from
`ADAPTIVE_REQ_KWAVE_PATH` or the default
`/Users/sara/Documents/k-wave-toolbox-version-1.4.1`. `validate` mode checks
paths and map construction without running time stepping. `quick` mode runs a
small homogeneous and inclusion field with a single sine source and M=2. `full`
mode remains a controlled phase-1 set: homogeneous 2/3 m/s, bilayer 2/3,
inclusion 2/3, single/multisource/square-wave diagnostics, and M=2/3. Optional
REQ extraction uses the same project-native estimator and Test38-style profile
(`Nbins='auto'`, `Nbins_auto_oversample=1`, `Nbins_min=16`, `smooth_sigma=1`,
`gamma_win=1`, `pad_factor=1`). The selected k-Wave field component defaults
to `axial_shear`, i.e. the depth/axial shear velocity component returned by
k-Wave as `uy_split_s`/`uy_s`, then transposed to project convention `Uxz(z,x)`.
The source defaults to `SourceSide=left` and `SourcePolarization=axial`, so the
main propagation is lateral while the measured particle velocity is axial; this
better mimics ultrasound axial velocity observations of a transverse shear wave
than a bottom/top source pushing axially along its propagation direction.
The default analysis ROI is `exclude_source_buffer`: k-Wave is run on the full
grid, but the field passed to REQ is cropped after a background buffer on the
source side. This keeps the source inside the simulated material while moving
it outside the measured/estimated region, reducing near-field and boundary
contamination. The full harmonic field and absolute ROI origin are still saved
for diagnostics.
`lateral_shear` and `shear_magnitude` are available only as diagnostics.
Controls include
`ADAPTIVE_REQ_TEST50_MODE=validate|quick|full`,
`ADAPTIVE_REQ_TEST50_RUN_REQ=true|false`,
`ADAPTIVE_REQ_TEST50_TARGET_STEP_M`, `ADAPTIVE_REQ_TEST50_USE_PARFOR`,
`ADAPTIVE_REQ_TEST50_VELOCITY_COMPONENT`,
`ADAPTIVE_REQ_TEST50_SOURCE_SIDE`, `ADAPTIVE_REQ_TEST50_SOURCE_POLARIZATION`,
`ADAPTIVE_REQ_TEST50_ANALYSIS_ROI`,
`ADAPTIVE_REQ_TEST50_ANALYSIS_BUFFER_M`, `ADAPTIVE_REQ_TEST50_ANALYSIS_MARGIN_M`,
and
`ADAPTIVE_REQ_TEST50_SAVE_TIME_SERIES`. Outputs are under
`outputs/test_50_controlled_kwave_phase1/<mode>/`.

Test 51 applies frozen Test38 q models to the controlled k-Wave fields produced
by Test 50. It lives in `tests/integration/test_51_controlled_kwave_test38_validation.m`
and does not train any model. Because CSV cannot preserve the per-patch
`req_mapping` needed for q-to-SWS conversion, Test 51 reloads Test50 field
`.mat` files, re-extracts REQ/features with the Test38 profile, and caches the
full feature table as `.mat` under its own output folder. It then applies
`q_spectrum_only`, `q_spectrum_plus_composition`, and
`q_spectrum_plus_theory_composition` by default, with optional
`theory_discrete` as a diagnostic baseline. Evaluation uses true SWS,
material side, patch purity, ROI, and distance to interface only after
inference. Controls include `ADAPTIVE_REQ_TEST51_MODE=validate|quick|full`,
`ADAPTIVE_REQ_TEST51_TEST50_SOURCE=quick|full|/path`,
`ADAPTIVE_REQ_TEST51_MODEL_SOURCE=quick|medium|full|/path`,
`ADAPTIVE_REQ_TEST51_TARGET_STEP_M`, `ADAPTIVE_REQ_TEST51_USE_PARFOR`,
`ADAPTIVE_REQ_TEST51_SAVE_ALL_MAPS`, and `ADAPTIVE_REQ_TEST51_MODELS`. Outputs
are under `outputs/test_51_controlled_kwave_test38_validation/<mode>/`.

Test 52 directly compares the frozen Test35 and Test38 clean-q model bundles
on identical dense synthetic REQ maps. It lives in
`experiments/analysis/analyze_test_52_test35_vs_test38_dense_comparison.m` and
does not train any model. The script reuses the representative Test48-style
simulation design and extracts one dense REQ feature table per condition, then
applies both frozen bundles to the same rows. Model names are prefixed with
`T35_` or `T38_` so rankings, ROI summaries, frequency/M plots, and maps are
honest same-map comparisons. This test answers whether Test35's strong
patch-level results survive dense-map evaluation or whether the degradation
seen in later maps is mostly due to dense sampling and harder conditions.
Controls include `ADAPTIVE_REQ_TEST52_MODE=quick|full`,
`ADAPTIVE_REQ_TEST52_VALIDATE_ONLY=true|false`,
`ADAPTIVE_REQ_TEST52_TEST35_MODEL_SOURCE=quick|full|/path`,
`ADAPTIVE_REQ_TEST52_TEST38_MODEL_SOURCE=quick|medium|full|/path`,
`ADAPTIVE_REQ_TEST52_TARGET_STEP_M`, `ADAPTIVE_REQ_TEST52_REQ_PROFILE`,
`ADAPTIVE_REQ_TEST52_USE_PARFOR`, and
`ADAPTIVE_REQ_TEST52_SAVE_ALL_MAPS`. Outputs are under
`outputs/test_52_test35_vs_test38_dense_comparison/`.

Test 53 is the paper-facing clean retraining entry point. It delegates to the
Test 38 training engine but pins the controlled publication design:
frequencies `200, 300, 400, 500, 600 Hz`, `M=[2 3]`, `dx=0.2 mm`, clean
synthetic fields, and outputs under
`outputs/test_53_paper_final_clean_q_training/`. It exists so the paper model
can be rerun without editing Test 38. Use
`ADAPTIVE_REQ_TEST53_MODE=quick|medium|full`,
`ADAPTIVE_REQ_TEST53_USE_PARALLEL_TRAINING=true`, and optional row caps via
`ADAPTIVE_REQ_TEST53_MAX_MODEL_TRAIN_ROWS` /
`ADAPTIVE_REQ_TEST53_MAX_MODEL_EVAL_ROWS`.

Test 53 strong-splits is the paper-facing leakage audit and OOD validation for
the Test 53 baseline. It lives in
`experiments/analysis/analyze_test_53_strong_splits_q_training.m` and retrains
the q baselines inside each fold without reusing composition predictions from
the full dataset. It evaluates grouped condition splits, leave-one-frequency
out, leave-one-geometry-family out, field-regime OOD, and optional leave-one-M
out. The fold-level composition predictors for `q_spectrum_plus_composition`
(`predicted_patch_purity`, `p_mixed`, `p_strong_mixed`) are trained only on the
training rows of that fold. True SWS, oracle q, true purity, coordinates,
distance-to-interface, errors, and confidence variables are excluded from
predictors by an explicit leakage guard. Outputs are under
`outputs/test53_strong_splits/<quick|full>/`, including split summaries,
frequency/M/geometry/regime/purity/ROI metrics, representative maps, runtime
logs, and `README_results.md`. The script defaults to estimate-only mode; run
with `ADAPTIVE_REQ_TEST53_STRONG_ESTIMATE_ONLY=false` only when ready to launch
the selected quick or full fold set.

Test 54 validates a frozen Test 38/Test 53 q model under realistic
ultrasound-like readout using `/Users/sara/Documents/wave_sim_project`.
It does not train. It generates clean, moderate-realistic, and hard-realistic
fields, re-extracts REQ with the Test 38 profile
(`Nbins='auto'`, oversample 1, minimum 16 bins, `smooth_sigma=1`,
`gamma_win=1`, `pad_factor=1`, valid edges), applies frozen q models, and
reports clean-to-realistic degradation. The realistic levels include shear
attenuation, spatially correlated readout noise, depth/readout SNR, acoustic
shadowing, backscatter contrast, and phase/amplitude QC maps. Run with
`ADAPTIVE_REQ_TEST54_MODE=validate|quick|full` and choose the frozen bundle via
`ADAPTIVE_REQ_TEST54_MODEL_SOURCE=full|test53|/path/to/bundle.mat`.

Test 55 is a controlled field-frequency ROI sweep for frozen Test 38/Test 53
models. It lives in
`experiments/analysis/analyze_test_55_controlled_field_frequency_roi_sweep.m`
and does not train. It explicitly evaluates selected homogeneous, bilayer, and
inclusion geometries over a balanced frequency x field-regime grid so that
directional 2D, diffuse 2D, partial 3D, and diffuse 3D are present at every
requested frequency. It re-extracts REQ with the matched Test 38 profile by
default (`Nbins='auto'`, oversample 1, minimum 16 bins, `smooth_sigma=1`,
`gamma_win=1`, `pad_factor=1`, valid edges), applies frozen
`q_spectrum_only` and `q_spectrum_plus_composition`, and saves ROI/core/interface
summaries plus SWS-vs-frequency curves separated by field type. Use
`ADAPTIVE_REQ_TEST55_MODE=quick|full`, `ADAPTIVE_REQ_TEST55_MODEL_SOURCE`,
`ADAPTIVE_REQ_TEST55_M_LIST`, `ADAPTIVE_REQ_TEST55_FREQUENCIES_HZ`,
`ADAPTIVE_REQ_TEST55_TARGET_STEP_M`, and
`ADAPTIVE_REQ_TEST55_SAVE_ALL_MAPS`.

Test 56 trains an operational M-effective residual correction on top of the
frozen `q_spectrum_plus_composition` prediction, using Test55-style patch
results as input. It does not retrain the base q model and does not use
`M_eff_true`, true SWS, true purity, ROI, material side, or distance as
predictors. The second-pass predictors are operational first-pass quantities:
predicted SWS/q, predicted `M_eff`, predicted wavelength/window length,
predicted composition/mixedness, M, frequency, dx/dz, and encoded regime/case
metadata. The target is `delta_logk_true = log(k_true)-log(k_firstpass)`, used
only as a supervised label. Strategies include baseline, residual applied
everywhere, residual only in high predicted purity, conservative residual, and
diagnostic oracle apply-if-improves. Run with
`ADAPTIVE_REQ_TEST56_MODE=quick|full`,
`ADAPTIVE_REQ_TEST56_SOURCE_CSV=/path/to/test55_patch_level_results.csv`, and
`ADAPTIVE_REQ_TEST56_VALIDATE_ONLY=true|false`. Outputs are under
`outputs/test_56_meff_operational_correction/`.

Test 57 evaluates a multi-window ablation without training. It expects a
Test55-style patch-level CSV generated with multiple M values, such as
`ADAPTIVE_REQ_TEST55_M_LIST='1.5 2 2.5 3'`. It compares fixed-M performance and,
when rows align across windows, operational multi-window strategies: mean SWS,
median SWS, lowest predicted mixedness, highest predicted purity, and predicted
M-effective value closest to a target. If the source CSV contains only one M,
it writes fixed-M summaries and prints the rerun command needed for the full
ablation. Run with `ADAPTIVE_REQ_TEST57_MODE=quick|full`,
`ADAPTIVE_REQ_TEST57_SOURCE_CSV`, and
`ADAPTIVE_REQ_TEST57_TARGET_MEFF`. Outputs are under
`outputs/test_57_multiwindow_ablation/`.

Test 58 trains an estimability/high-error risk mask. It does not correct SWS.
The clean Test55 predictions are used to train classifiers for high-error
>10% and >20% using only operational predictors: predicted SWS/q, predicted
M-effective features, predicted composition/mixedness, M, frequency, dx/dz,
and encoded regime/case metadata. Error labels are supervised targets only.
If a Test54 realistic-readout patch table exists, the same frozen risk model is
also evaluated there to test transfer under noise/readout artifacts. Outputs
include ROC/PR AUC, threshold summaries, high-confidence versus low-confidence
MAPE, risk-bin summaries, and calibration plots. Run with
`ADAPTIVE_REQ_TEST58_MODE=quick|full`,
`ADAPTIVE_REQ_TEST58_CLEAN_CSV`, `ADAPTIVE_REQ_TEST58_REALISTIC_CSV`, and
`ADAPTIVE_REQ_TEST58_VALIDATE_ONLY=true|false`. Outputs are under
`outputs/test_58_estimability_risk_mask/`.

Test 59 combines the Test56 M-effective residual correction, the Test57
M=2/M=3 multi-window candidate, and a Test58-style high-error risk gate. It
uses the frozen `q_spectrum_plus_composition` predictions from a Test55
multi-M patch table and trains only second-pass operational layers: a
`delta_logk` residual, a high-error risk model, and a candidate selector that
chooses keep M=2, switch to M=3, or apply the residual. True SWS/k/error labels
are used only for supervised labels and evaluation; true purity, material
region, ROI, and distance are diagnostic only. The primary success target is
reducing hard-core subestimation without degrading homogeneous, soft-core, or
interface regions. Run with `ADAPTIVE_REQ_TEST59_MODE=quick|full`,
`ADAPTIVE_REQ_TEST59_SOURCE_CSV`, `ADAPTIVE_REQ_TEST59_VALIDATE_ONLY`, and
`ADAPTIVE_REQ_TEST59_SAVE_ALL_MAPS`. When map saving is enabled, Test59 now
saves one panel for every evaluated physical condition; use
`ADAPTIVE_REQ_TEST59_MAX_MAP_CONDITIONS=0` for no limit, or set a positive
integer for debug runs. M=2/M=3 rows are aligned by nearest physical center
using `ADAPTIVE_REQ_TEST59_ALIGN_TOLERANCE_M` because exact StepX/StepZ centers
can differ across M/frequency. Outputs are under
`outputs/test_59_hard_aware_meff_multiwindow_gate/`.

Test 60 is the conservative follow-up to Test59. It keeps the same frozen
`q_spectrum_plus_composition` base predictions and the same M=2/M=3 alignment,
but adds two deployable selector variants: an interface-penalized selector and
a spatial interface-penalized selector. The interface penalty increases the
required training gain for candidate corrections in pixels that look risky,
mixed, or low-purity using operational predictors only. The spatial variant
then removes small isolated correction islands per condition using a local
support rule on the patch grid; it does not use true SWS, material masks,
distance-to-interface, or ROI labels at inference. Runtime controls include
`ADAPTIVE_REQ_TEST60_MODE=quick|full`,
`ADAPTIVE_REQ_TEST60_VALIDATE_ONLY=true|false`,
`ADAPTIVE_REQ_TEST60_SAVE_ALL_MAPS=true|false`,
`ADAPTIVE_REQ_TEST60_MAX_MAP_CONDITIONS`,
`ADAPTIVE_REQ_TEST60_INTERFACE_EXTRA_GAIN_MARGIN_PCT`,
`ADAPTIVE_REQ_TEST60_PENALIZED_SELECTOR_RISK_THRESHOLD`,
`ADAPTIVE_REQ_TEST60_SPATIAL_GATE_RADIUS_PX`, and
`ADAPTIVE_REQ_TEST60_SPATIAL_GATE_MIN_SUPPORT`. Outputs are under
`outputs/test_60_interface_penalized_spatial_gate/`.

Test 61 replaces the hand-designed acceptance rule with a learned
benefit/harm correction gate. It still keeps the frozen
`q_spectrum_plus_composition` estimator fixed and reuses the Test59/Test60 M=2
versus M=3 alignment. For each pixel, the script builds candidate rows for M=3
and the Test56-style residual correction, then trains a regression model for
expected MAPE gain and a classifier for probability of harm using only
operational features: candidate q/SWS/log-k/M-effective values, predicted
composition/mixedness, risk scores, M/frequency/dx/regime/case metadata, and
candidate/base disagreement. True SWS and errors are labels only. Strategies
include fixed M=2, fixed M=3, Test56 residual, Test59 learned selector, Test60
interface-penalized selectors, learned benefit/harm gate, conservative
benefit/harm gate, and oracle best-candidate. Runtime controls include
`ADAPTIVE_REQ_TEST61_MODE=quick|full`,
`ADAPTIVE_REQ_TEST61_VALIDATE_ONLY=true|false`,
`ADAPTIVE_REQ_TEST61_SAVE_ALL_MAPS=true|false`,
`ADAPTIVE_REQ_TEST61_BENEFIT_MIN_GAIN_PCT`,
`ADAPTIVE_REQ_TEST61_BENEFIT_MAX_HARM_PROB`,
`ADAPTIVE_REQ_TEST61_BENEFIT_HARM_PENALTY_PCT`,
`ADAPTIVE_REQ_TEST61_CONSERVATIVE_MIN_GAIN_PCT`,
`ADAPTIVE_REQ_TEST61_CONSERVATIVE_MAX_HARM_PROB`, and
`ADAPTIVE_REQ_TEST61_CONSERVATIVE_HARM_PENALTY_PCT`. Outputs are under
`outputs/test_61_learned_benefit_harm_gate/`.

Test 62 validates the frozen Test61 correction gate on dense representative
REQ maps. It recalculates dense REQ features using the Test38-matched settings
(`Nbins='auto'`, `Nbins_auto_oversample=1`, `Nbins_min=16`,
`smooth_sigma=1`, `gamma_win=1`, `pad_factor=1`) for M=2 and M=3, applies the
frozen `q_spectrum_plus_composition` estimator from Test38, aligns M2/M3
centers by physical position, and then applies the already-trained Test61
strategies without retraining. The purpose is to check whether the learned
benefit/harm gate still behaves well on visually dense maps rather than sparse
training/evaluation rows. Tables include dense patch-level strategy results,
base M2/M3 predictions, summaries by case/regime/frequency/ROI, and correction
harm. Figures include per-condition dense maps with true SWS, fixed M2/M3,
interface-penalized selector, learned benefit/harm gate, conservative
benefit/harm gate, signed errors, risk maps, and correction masks. Runtime
controls include `ADAPTIVE_REQ_TEST62_MODE=quick|full`,
`ADAPTIVE_REQ_TEST62_VALIDATE_ONLY=true|false`,
`ADAPTIVE_REQ_TEST62_MODEL_SOURCE=medium|full|<path>`,
`ADAPTIVE_REQ_TEST62_TEST61_SOURCE=full|<path>`,
`ADAPTIVE_REQ_TEST62_TARGET_STEP_M`,
`ADAPTIVE_REQ_TEST62_SAVE_ALL_MAPS=true|false`, and
`ADAPTIVE_REQ_TEST62_MAX_MAP_CONDITIONS`. Outputs are under
`outputs/test_62_dense_test61_benefit_harm_validation/`.

Test 63 is the controlled dense matrix version of Test62. It uses the same
frozen Test38 `q_spectrum_plus_composition` model and frozen Test61
correction/reliability gate, but runs a factorial design so frequency,
geometry, field regime, and M effects can be compared directly. The full
controlled grid is: geometries `homogeneous_cs2`, `homogeneous_cs3`,
`homogeneous_cs4`, `bilayer_2_4`, `circular_inclusion_2_4`, and
`three_material_2_3_4`; frequencies `200, 300, 400, 500 Hz`; regimes
`directional_2D_angle0`, `directional_2D_angle30`, `diffuse_2D_seed1`,
`partial_3D_12src`, and `diffuse_3D_seed1`; and M values `[2 3]`. Quick mode
uses a smaller factorial subset but keeps the same design logic. Outputs add
controlled summaries by case/frequency/regime, case/frequency,
case/regime, ROI/frequency, ROI/regime, and an alignment coverage table
showing how many M2 and M3 dense centers were available and aligned per
condition. This coverage table is important at low frequencies because M=3
can have very few valid centers when the physical REQ window approaches the
field size. Runtime controls include `ADAPTIVE_REQ_TEST63_MODE=quick|full`,
`ADAPTIVE_REQ_TEST63_VALIDATE_ONLY=true|false`,
`ADAPTIVE_REQ_TEST63_MODEL_SOURCE=medium|full|<path>`,
`ADAPTIVE_REQ_TEST63_TEST61_SOURCE=full|<path>`,
`ADAPTIVE_REQ_TEST63_TARGET_STEP_M`,
`ADAPTIVE_REQ_TEST63_SAVE_ALL_MAPS=true|false`, and
`ADAPTIVE_REQ_TEST63_MAX_MAP_CONDITIONS`. Outputs are under
`outputs/test_63_controlled_dense_frequency_regime_matrix/`.

Test 64 audits `/Users/sara/Documents/wave_sim_project` before using it as a
realistic validation source. It does not run REQ and does not train or apply q
models. The primary path is `simcore.simulateEikonalVectorReadout2p5D`, using
the axial complex phasor `out.U` as the future REQ input. It checks clean,
shear-attenuation, readout-depth-noise, and readout-depth-plus-acoustic-shadow
layers; homogeneous phase-derived SWS; inclusion and bilayer masks;
signed-distance conventions; ROI placement/purity; source/polarization axial
measurement weights; and compact storage. The primary vector 2.5D medium
builder now reuses `simcore.medium.masks.alpha2D`, so it supports homogeneous,
circle/inclusion, bilayer, s-curve, and custom masks in the eikonal/readout
path. Materials are resolved through
`simcore.medium.resolveMaterialAtFrequency(material,f0)`, which applies
frequency-dependent attenuation and Kelvin-Voigt properties when configured.
The acoustic readout path now separates total ultrasound loss, homogeneous
background loss, heterogeneity-only extra shadow, echo amplitude, tracking SNR
from echo, Boukraa depth weight, echo/SNR noise weight, and total Boukraa noise
weight. The acoustic attenuation map follows
`alpha_us(x,z,f_us)=alpha0_material(x,z)*f_us_MHz^power_material`, and the ray
integral can use a single plane wave or plane-wave compounding. In the
readout-depth-plus-acoustic-shadow layer, acoustic shadow is now represented as
lower local echo amplitude / lower tracking SNR; Boukraa readout noise is then
scaled from that SNR map instead of treating shadow as an arbitrary standalone
noise multiplier. The companion `wave_sim_project` Test 01 documents and checks
this material-resolution and readout-weight path with explicit diagnostic
figures. Runtime controls include
`ADAPTIVE_REQ_TEST64_MODE=validate|quick|full`,
`ADAPTIVE_REQ_TEST64_WAVE_SIM_PATH`,
`ADAPTIVE_REQ_TEST64_SAVE_ALL_MAPS=true|false`, and
`ADAPTIVE_REQ_TEST64_SAVE_TIMESERIES_DIAG=true|false`. Outputs are under
`outputs/test_64_realistic_simulation_audit/`.

Test 65 audits finite lateral source size and readout geometry before running
adaptive REQ. It lives in `adaptive_req_local` and uses `wave_sim_project` only
as a simulation backend. It does not extract REQ features and does not train or
apply q models. It compares a point-like source with finite Hann-aperture
sources, using an external lateral source with axial displacement readout, and
checks clean versus readout-depth-plus-acoustic-shadow conditions. The plots
separate lateral profiles from axial/depth profiles because the shear wave
propagates primarily from the lateral source into the field; shear displacement
amplitude loss is therefore expected mainly along x, while ultrasound echo
amplitude and tracking SNR vary with z and with acoustic shadow. It also records
finite-source weights, clean/readout amplitude and phase, phase error, shear
attenuation at `f0`, ultrasound attenuation at `fUS`, echo amplitude, tracking
SNR from echo, total acoustic loss, heterogeneity-only extra shadow, Boukraa
depth weight, echo/SNR noise weight, and total Boukraa noise weight. Runtime
controls include `ADAPTIVE_REQ_TEST65_MODE=validate|quick|full`,
`ADAPTIVE_REQ_TEST65_WAVE_SIM_PATH`, and
`ADAPTIVE_REQ_TEST65_SAVE_ALL_MAPS=true|false`. Outputs are under
`outputs/test_65_finite_source_readout_audit/<mode>/`.

Test 66 auxiliary shadow audit isolates the double-looking shear-amplitude
shadow seen in Test 65. It does not run REQ and does not use any q or confidence
model. It compares a matched 2x2 design for an inclusion field: point versus
4 mm finite Hann source, and `straight_ray_integral` versus
`eikonal_ray_integral` shear attenuation. The script is
`experiments/analysis/analyze_test_66_shear_shadow_source_integral_comparison.m`
and outputs are under
`outputs/test_66_shear_shadow_source_integral_comparison/<mode>/`.

Test 66 zero-shot Eikonal realistic transfer validation is the paper-facing
transfer test for frozen Test53/Test38 q models. The script
`experiments/analysis/analyze_test_66_eikonal_realistic_transfer_validation.m`
loads `q_spectrum_only`, `q_spectrum_plus_composition`, and the frozen
composition/purity auxiliaries without retraining anything. It evaluates clean,
shear-attenuated, and readout-noisy 2.5D Eikonal fields generated by
`wave_sim_project_clean`, including homogeneous 2/3/4 m/s, bilayer, inclusion,
and a three-material bilayer+inclusion geometry with labels 0/1/2. It extracts
REQ with `Nbins='auto'`, `Nbins_min=16`, `smooth_sigma=1`, `cs_guess=3`, and
`TargetStepM=0.5 mm`, then reports MAPE, bias, high-error rates, q error,
patch purity, SNR/amplitude/depth bins, distance-to-interface bins, and ROI
summaries. Runtime controls include
`ADAPTIVE_REQ_TEST66_MODE=validate|quick|full_a|full_b`,
`ADAPTIVE_REQ_TEST66_MODEL_SOURCE=test53|full|medium|quick|<bundle.mat>`,
`ADAPTIVE_REQ_TEST66_WAVE_SIM_PATH`,
`ADAPTIVE_REQ_TEST66_TARGET_STEP_M`, and
`ADAPTIVE_REQ_TEST66_SAVE_ALL_MAPS=true|false`. Outputs are under
`outputs/test_66_eikonal_realistic_transfer_validation/<mode>/`.

Test 68 trains an operational estimability/risk mask on existing Test66
patch-level outputs without retraining or modifying any q, composition, or
correction model. The script
`experiments/analysis/analyze_test_68_test66_estimability_mask.m` uses only
inference-available predictors such as frozen q/SWS predictions, frozen
composition probabilities, q/theory disagreement, local q/SWS neighborhood
variation, amplitude/SNR readout diagnostics, source distance, frequency, M,
and field regime. True SWS, q-oracle, true patch purity, distance to interface,
ROI labels, and error maps are excluded from predictors and used only as
training labels or diagnostics. It reports ROC/PR AUC, calibration,
threshold tradeoffs, accuracy-coverage curves, and grouped summaries by
geometry, frequency, realism, field regime, ROI, true purity, and
distance/window-radius. Runtime controls include
`ADAPTIVE_REQ_TEST68_MODE=quick|full`,
`ADAPTIVE_REQ_TEST68_TEST66_SOURCE=full_a|quick|<run_dir>`,
`ADAPTIVE_REQ_TEST68_MAX_TRAIN_ROWS`,
`ADAPTIVE_REQ_TEST68_MAX_EVAL_ROWS`, and
`ADAPTIVE_REQ_TEST68_SAVE_ALL_MAPS=true|false`. Outputs are under
`outputs/test_68_test66_estimability_mask/`.

Test 69 applies the frozen Test 68 estimability mask back onto dense Test66
map outputs. The script
`experiments/analysis/analyze_test_69_dense_estimability_maps.m` does not
retrain any q, composition, correction, or risk model. It loads the frozen
Test68 high-error >20% bagged-tree detectors, reconstructs the same operational
features on selected dense Test66 patch maps, averages risk across the frozen
folds, and writes per-condition panels with true SWS, predicted SWS, masked
SWS, signed/absolute error, risk, reliability, SNR/amplitude diagnostics,
predicted mixedness, true patch purity, and distance/window-radius diagnostics.
Runtime controls include `ADAPTIVE_REQ_TEST69_MODE=quick|full`,
`ADAPTIVE_REQ_TEST69_TEST66_SOURCE=full_a|quick|<run_dir>`,
`ADAPTIVE_REQ_TEST69_TEST68_BUNDLE=<compact.mat>`,
`ADAPTIVE_REQ_TEST69_MODEL=q_spectrum_plus_composition`,
`ADAPTIVE_REQ_TEST69_M=2`, `ADAPTIVE_REQ_TEST69_RISK_THRESHOLD=0.5`,
`ADAPTIVE_REQ_TEST69_MAX_CONDITIONS`, and
`ADAPTIVE_REQ_TEST69_SAVE_ALL_MAPS=true|false`. Outputs are under
`outputs/test_69_dense_estimability_maps/<mode>/`.

Test 70 evaluates a conservative high-reliability M-effective M2/M3 selector
on dense Test66 Eikonal maps. The script
`experiments/analysis/analyze_test_70_eikonal_meff_m2_m3_selector.m` does not
retrain q, composition, risk, or correction models. It loads existing Test66
M=2 and M=3 dense patch predictions, aligns patch centers by physical
coordinates, applies the frozen Test68 high-error risk model to M=2, and
switches from M=2 to M=3 only when the patch is high-reliability, predicted
pure/non-mixed, hard-like, and the M=3 estimate increases SWS smoothly. It
compares fixed M=2, fixed M=3, the operational reliable-hard M2->M3 switch,
and a diagnostic oracle best-of-M2/M3. True SWS, ROI labels, true patch purity,
and interface distance are used only for evaluation. Test66 now supports
filtered, labeled runs via `ADAPTIVE_REQ_TEST66_RUN_LABEL` plus geometry,
frequency, realism, field-regime, and seed filters, and includes
`bilayer_2_4` for this validation. Runtime controls include
`ADAPTIVE_REQ_TEST70_MODE=quick|full`,
`ADAPTIVE_REQ_TEST70_TEST66_SOURCE=full_a|<run_dir>`,
`ADAPTIVE_REQ_TEST70_RISK_BUNDLE=<compact.mat>`,
`ADAPTIVE_REQ_TEST70_RISK_THRESHOLD`,
`ADAPTIVE_REQ_TEST70_SAVE_ALL_MAPS=true|false`, and optional filters for
geometries, frequencies, realism levels, and field regimes. Outputs are under
`outputs/test_70_eikonal_meff_m2_m3_selector/<mode>/`.

## Current operational chain

The current intended sequence is:

```text
Test 18 clean dataset
  -> Test 19 q-model training
  -> Test 20 external q/SWS validation
  -> Test 21 confidence-detector training
  -> Test 22 frozen confidence-detector external validation
```

Test 22 deliberately keeps operational predictors, diagnostic variables, and
error targets separate. True SWS, interface distance, SWS error, and high-error
labels are never detector inputs.

## Detailed documentation

- `docs/TEST08_ADVANCED_ANGULAR_FEATURES.md`
- `docs/TEST09_STRESS_AND_TEST10_HETEROGENEOUS_PLAN.md`
- `docs/TEST11_GLOBAL_REQ_FEATURES.md`
- `docs/test_11_adaptive_req_global_local_summary.md`
- `docs/test_12_cs_guess_window_sweep_analysis_summary.md`
- `docs/TEST15_THEORY_INFORMED_DIRECT_Q.md`
- `docs/test_18_clean_field_regime_training_dataset.md`
- `docs/progress_summaries/2026-06-11_adaptive_req_project_progress_summary.md`
