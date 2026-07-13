# adaptive_req_local

Local MATLAB workspace for the regenerated adaptive REQ baseline and Eikonal validation line.

## What is tracked

This repository tracks the maintainable project code:

- `src/`: reusable adaptive REQ MATLAB functions.
- `experiments/`: experiment runners and analysis scripts.
- `configs/`: JSON configuration files for reproducible runs.
- `docs/`: documentation and result summaries.
- `archive/`: selected archived scripts/docs from earlier project stages.
- `setup_adaptive_req.m`: MATLAB setup helper.

Generated outputs are intentionally ignored by git because they can be very large.

## Main current lines

- `baseline_minimal_v1`: rebuilt minimal baseline training and diagnostics.
- `eikonal_validation`: staged transfer validation of frozen baseline models on clean/realistic Eikonal simulations.

## Output policy

The `outputs/` folder is local-only. It contains models, cached simulations, figures, CSV tables, and other generated artifacts. Keep important interpretation summaries in `docs/` so they can be versioned without uploading large data.

## Setup

In MATLAB:

```matlab
cd('/Users/sara/local/adaptive_req_local')
run('setup_adaptive_req.m')
```
