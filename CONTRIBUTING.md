# Contributing to Adaptive REQ

Thank you for considering a contribution. Adaptive REQ is active research
software for wavefield analysis and shear-wave-speed estimation. Contributions
that improve reproducibility, numerical reliability, documentation, testing,
and portability are especially welcome.

## Ways to contribute

You can contribute by:

- reporting reproducible bugs;
- proposing or implementing numerical tests;
- improving documentation and examples;
- identifying MATLAB compatibility problems;
- validating edge cases or physical assumptions;
- suggesting improvements to the public API.

For major methodological changes, please open an issue before preparing a pull
request so the scientific scope and validation plan can be discussed first.

## Before opening an issue

Please check whether a similar issue already exists. For bug reports, include:

- MATLAB release;
- operating system;
- relevant MATLAB toolboxes;
- repository commit or release;
- configuration file or parameter values;
- the smallest script or dataset that reproduces the problem;
- expected and observed behavior;
- complete error message and stack trace, when applicable.

Do not upload confidential, clinical, proprietary, or personally identifiable
data.

## Development workflow

1. Fork the repository.
2. Create a descriptive branch:

   ```bash
   git checkout -b fix/radial-spectrum-normalization
   ```

3. Run the repository setup in MATLAB:

   ```matlab
   setup_adaptive_req();
   ```

4. Make a focused change.
5. Add or update tests and documentation when appropriate.
6. Run the relevant tests or minimal validation scripts.
7. Commit with a concise description.
8. Open a pull request against `main`.

## Pull-request expectations

A pull request should explain:

- what changed;
- why the change is needed;
- which files or workflows are affected;
- how the change was validated;
- any known limitations or compatibility implications.

Please avoid combining unrelated changes in one pull request. Generated outputs,
large `.mat` files, caches, and local paths should not be committed.

## MATLAB style

- Use descriptive function and variable names.
- Prefer functions over scripts for reusable numerical operations.
- Document units in names, comments, or function help.
- Validate important inputs and fail with informative error messages.
- Avoid machine-specific absolute paths.
- Preserve deterministic random seeds in reproducibility tests.
- Separate physical assumptions, numerical approximations, and visualization.
- Use `fullfile` for portable paths.

## Scientific validation

Methodological changes should be supported by an appropriate validation case.
Depending on the change, this may include:

- analytical or limiting-case behavior;
- homogeneous-medium simulations;
- controlled directional and diffuse wavefields;
- sensitivity to frequency, sampling, or window size;
- comparison against a frozen reference output;
- held-out or external validation.

Training-set performance alone should not be presented as evidence of
generalization.

## Research and clinical scope

This repository is research software and is not intended for clinical diagnosis,
treatment, or medical decision-making. Contributions must not represent
experimental presets as universal clinical constants.

## Licensing

By submitting a contribution, you agree that it may be distributed under the
repository's Apache License 2.0.
