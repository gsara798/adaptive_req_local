# Changelog

All notable public changes to Adaptive REQ will be documented in this file.

The format follows the principles of [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and the project intends to use semantic versioning after the first stable public
API is defined.

## [Unreleased]

### Added

- Public-facing project documentation.
- Apache License 2.0.
- Citation metadata for scientific software.
- Contribution and community guidelines.

### Planned

- Portable end-to-end REQ example using a small synthetic wavefield.
- Documented MATLAB and toolbox requirements.
- Unit tests for radial-spectrum and cumulative-energy calculations.
- Regression tests for representative adaptive-REQ workflows.
- Continuous integration for core numerical functions.
- Versioned public model and configuration metadata.

## [0.1.0] - YYYY-MM-DD

### Added

- Initial public research release.
- Reusable MATLAB functions under `src/`.
- Experiment runners and analysis scripts under `experiments/`.
- Versioned experiment configurations under `configs/`.
- Documentation for adaptive REQ baselines and Eikonal transfer validation.

### Notes

- This release should be created only after confirming that the public example
  runs from a clean clone and that its required dependencies are documented.
- Large generated outputs, simulation caches, and trained model artifacts are
  not included in the repository.
