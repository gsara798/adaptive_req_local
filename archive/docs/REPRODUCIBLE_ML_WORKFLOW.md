# Reproducible ML workflow

## Scientific objective

The production pipeline is:

```text
local wavefield patch
-> spectral features
-> trained ML model
-> predicted quantile q
-> local REQ mapping
-> predicted wavenumber k
-> predicted SWS
```

The ML target is `q_theory`. The final scientific metric is the error of the
SWS reconstructed from `q_pred`.

## Dataset outputs

For training and evaluation, store:

- Stable row identifiers: condition, aperture step, realization, and patch.
- Scalar spectral features.
- `q_theory`.
- Simulation and REQ metadata needed to define groups and reproduce runs.
- `req_mapping`, the lightweight 1D relation between cumulative energy and
  radial wavenumber.

Do not store `req_curve` for the complete ML dataset. It includes large 2D
diagnostic arrays and should only be enabled for small validation runs.

Recommended output configuration:

```matlab
CFG.OUTPUT.store_req_curve = false;
CFG.OUTPUT.store_req_mapping = true;
CFG.OUTPUT.store_req_metadata = true;
CFG.OUTPUT.store_feature_struct = false;
```

## Why `req_mapping` is sufficient

Converting a predicted quantile to SWS uses:

```matlab
k_pred = adaptive_req.quantile.quantile_to_k(req_mapping, q_pred);
cs_pred = 2*pi*f0/k_pred;
```

Only `req_mapping.Ecum` and `req_mapping.k_cent` are required. The mapping is
stored in single precision to reduce file size.

## Reproducibility rules

1. Keep each experiment definition in one configuration profile.
2. Record the profile name, resolved configuration, seeds, and sweep matrix.
3. Split ML train and test sets by simulation condition, not by random rows.
4. Keep full REQ curves only for a small diagnostic subset.
5. Treat generated outputs as disposable; the configuration and code should
   be sufficient to regenerate them.
6. Report both q-prediction metrics and final SWS metrics.

Run the fast automated checks with:

```matlab
results = run_unit_tests();
```

## Recommended experiment layers

```text
configs/
  Defines reproducible experiment inputs.

src/+adaptive_req/
  Reusable simulation, feature, quantile, ML, and evaluation code.

experiments/
  Thin scripts that call reusable code for one scientific question.

tests/
  Automated correctness and regression tests.

outputs/
  Generated results only.
```
