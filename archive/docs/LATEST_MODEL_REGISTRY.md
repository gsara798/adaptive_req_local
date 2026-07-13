# Latest Model Registry

Created: 2026-06-27

Local registry root:

`outputs/model_registry/`

The latest model registry uses hard links where possible, so the registry paths
remain valid even if the original output folders are later deleted, without
duplicating the large `.mat` payloads on disk.

## Recommended Current Bundle

Use this as the current frozen clean q/SWS model family:

`outputs/model_registry/test38_velocity_field_diverse_q_training/test38__velocity_field_diverse_q__medium_bundle.mat`

This is the Test 38 medium velocity/field-diverse clean spectral q bundle used
by Test 39 external validation. Test 39 itself is validation-only and does not
train a new model.

## Registered Bundles

| Model ID | Source Test | Role | Registry File | Summary |
|---|---|---|---|---|
| `test33__mixedness_aware_q_correction__full_bundle` | Test 33 | mixedness-aware q correction and reliability | `outputs/model_registry/test33_mixedness_aware_q_correction/test33__mixedness_aware_q_correction__full_bundle.mat` | research correction bundle |
| `test35__spectral_composition_q__quick_bundle` | Test 35 quick | clean spectral-composition q family | `outputs/model_registry/test35_spectral_composition_to_q_model/test35__spectral_composition_q__quick_bundle.mat` | best quick held-out MAPE about 3.43% |
| `test35__spectral_composition_q__full_bundle` | Test 35 full | clean spectral-composition q family | `outputs/model_registry/test35_spectral_composition_to_q_model/test35__spectral_composition_q__full_bundle.mat` | best synthetic held-out MAPE about 2.44% |
| `test38__velocity_field_diverse_q__quick_bundle` | Test 38 quick | velocity/field-diverse clean q family | `outputs/model_registry/test38_velocity_field_diverse_q_training/test38__velocity_field_diverse_q__quick_bundle.mat` | best quick held-out MAPE about 2.56% |
| `test38__velocity_field_diverse_q__medium_bundle` | Test 38 medium | recommended velocity/field-diverse clean q family | `outputs/model_registry/test38_velocity_field_diverse_q_training/test38__velocity_field_diverse_q__medium_bundle.mat` | best medium held-out MAPE about 4.21% |

## Local Manifests

Detailed local manifests:

- `outputs/model_registry/model_manifest.csv`
- `outputs/model_registry/latest_models_manifest.csv`
- `outputs/model_registry/LATEST_MODELS_README.md`

## Predictor Policy

For the clean Test 35 and Test 38 q bundles, primitive predictors are spectral
REQ/Ecum features, physical metadata such as `M`, `f0`, `dx`, `dz`, optional
theory prior terms, and predicted composition outputs when used by the specific
model variant. True SWS, true patch purity, error labels, and confidence maps
are not primitive inference predictors.
