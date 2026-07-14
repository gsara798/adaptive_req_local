# Adaptive REQ

Research-oriented MATLAB toolbox for **adaptive radial energy-quantile (REQ) estimation** and local shear-wave-speed mapping in harmonic shear-wave elastography.

This repository contains reusable analysis functions, experiment runners, configuration files, and technical documentation developed to study how the REQ quantile should adapt to the local wavefield. The current research focuses on directional, diffuse, and projected wavefields, as well as transfer validation using controlled Eikonal simulations.

> **Research status:** this is active research software. Interfaces, model definitions, and experiment organization may change before the first stable release. The code is not intended for clinical use or medical decision-making.

## Method overview

REQ estimates a representative local wavenumber from the cumulative radial energy of a windowed wavefield patch. Starting from the local radial spectrum, the normalized cumulative energy is

```text
E_cum(k) = integral_0^k S_rad(k') dk' / integral_0^kmax S_rad(k') dk'.
```

A quantile `q` selects the wavenumber `k_q` through

```text
E_cum(k_q) = q,
```

and the corresponding shear-wave speed is

```text
c_s = 2*pi*f0/k_q.
```

The adaptive workflow investigates whether spectral, cumulative-energy, acquisition, and local-composition features can select an operational value of `q` without using the true shear-wave speed at inference time.

## Scientific scope

The repository currently supports research workflows for:

- local radial-spectrum and cumulative-energy analysis;
- extraction of physically motivated REQ and `E_cum` features;
- training and evaluation of adaptive quantile models;
- conversion of predicted quantiles into local shear-wave-speed estimates;
- accuracy, bias, precision, robustness, and high-error analysis;
- validation across frequency, spatial sampling, wavefield dimensionality, and geometry;
- transfer evaluation on clean and perturbed Eikonal simulations.

The main active development lines are:

- **Minimal adaptive baseline:** compact operational models for quantile and SWS estimation.
- **Eikonal transfer validation:** evaluation of frozen adaptive models under controlled propagation, attenuation, and readout conditions.

## Repository structure

```text
adaptive_req_local/
├── src/                 Reusable MATLAB functions
├── experiments/         Training, validation, and analysis runners
├── configs/             Versioned JSON experiment configurations
├── docs/                Method notes and result summaries
├── archive/             Selected historical material; not part of the stable API
├── setup_adaptive_req.m MATLAB path setup helper
└── outputs/             Local generated artifacts; ignored by git
```

Code under `src/` is the intended reusable layer. Scripts under `experiments/` may depend on specific configurations, cached simulations, or trained model artifacts.

## Getting started

Clone the repository and open MATLAB in the repository root:

```bash
git clone https://github.com/gsara798/adaptive_req_local.git
cd adaptive_req_local
```

Then initialize the MATLAB path:

```matlab
root_dir = setup_adaptive_req();
```

The setup function adds the repository's `src/` directory to the MATLAB path and returns the resolved project root. No machine-specific absolute path is required.

## Running research workflows

Experiment entry points are organized under `experiments/`, with corresponding parameters under `configs/` and interpretation notes under `docs/`.

Before running an experiment:

1. Read the associated documentation and configuration.
2. Confirm that any required simulation cache or trained model bundle is available locally.
3. Run `setup_adaptive_req()` from the repository root.
4. Execute the selected experiment runner from MATLAB.
5. Record the configuration, MATLAB release, and Git commit used for the run.

Some long-running workflows generate large simulation caches and model files that are intentionally not stored in GitHub. A small self-contained public demonstration is part of the planned first release.

## Output and data policy

Generated content is written under `outputs/` and excluded from version control. This can include:

- simulation caches;
- trained `.mat` model bundles;
- intermediate feature tables;
- figures and result tables;
- logs and temporary analysis files.

This policy keeps the repository lightweight while preserving source code, configurations, and interpretation summaries. Results intended for long-term reference should be summarized in `docs/` together with the exact configuration and commit used.

## Evaluation principles

Performance claims should be based on held-out or explicitly external validation data. Training-set results are useful for diagnostics but should not be interpreted as unbiased estimates of generalization.

The project distinguishes among:

- **accuracy:** MAPE, median absolute percentage error, RMSE, and bias;
- **intramap precision:** coefficient of variation and robust spatial variability;
- **inter-realization robustness:** variability across seeds or repeated simulated realizations;
- **failure behavior:** high-error rate, tail error, and invalid-prediction fraction.

These quantities should be considered jointly because a low global mean error does not guarantee stable local performance near interfaces, at low SNR, or under distribution shift.

## Companion simulation repository

Controlled Eikonal wavefields used in transfer studies are developed separately in [`wave_sim_project_clean`](https://github.com/gsara798/wave_sim_project_clean). That repository models propagation and readout effects; this repository focuses on adaptive REQ feature extraction, model development, and performance evaluation.

## Limitations

- This repository is a research codebase, not a validated clinical product.
- Eikonal simulations are controlled phenomenological approximations and are not full elastodynamic solvers.
- Estimated SWS can depend on local window size, frequency, spatial sampling, field dimensionality, geometry, attenuation, and readout SNR.
- Near-interface and low-SNR results require separate reliability analysis.
- Large generated datasets and trained model bundles are not currently distributed with the repository.

## Roadmap

Planned public-release improvements include:

- a portable end-to-end example with a small synthetic wavefield;
- documented MATLAB and toolbox requirements;
- automated unit and regression tests;
- continuous integration for core numerical functions;
- tagged releases, citation metadata, and a stable public API;
- clearer separation between research experiments and reusable library functions.

## Contributing

Issues and pull requests are welcome, particularly for:

- reproducibility problems;
- numerical validation and edge cases;
- documentation improvements;
- portable examples and tests;
- MATLAB compatibility issues.

For a reproducible bug report, include the MATLAB release, operating system, configuration file, Git commit, and the smallest example that demonstrates the problem.

## License

This project is licensed under the Apache License 2.0. See the
[LICENSE](LICENSE) file for details.

## Maintainer

**Sara Gómez** — [@gsara798](https://github.com/gsara798)
