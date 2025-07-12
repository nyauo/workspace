# Development Guidelines

This repository contains MATLAB/Octave scripts for diode model fitting.
Follow these steps when modifying code or adding new features.

## Manual Testing
1. Run `main.m` to manually verify changes.
   This loads configuration, performs fitting, and interacts with the user.
2. Save result artifacts when prompted. They include `.mat`, `.csv`, `.png`
   and `.txt` files produced by `saveResults` or `saveAdjustedParameters`.
   Commit these under a new directory called `artifacts/` so that results
   from tests can be inspected later.

## Validation
Before committing any changes, validate the code using MATLAB or GNU
Octave:
1. Start MATLAB/Octave in this repository root.
2. Run `main` and check that it finishes without errors.
3. Examine the generated artifacts under `artifacts/` to ensure the
   outputs make sense.

If Octave is not installed locally, install it first (e.g. `apt-get
install octave`).

## Style and Comments
- Use four spaces for indentation and avoid tabs.
- Keep line width under 80 characters where possible.
- Write clear comments explaining complex operations. Use English comments
  unless the surrounding code uses another language consistently.
- Document any new functions with a short description, input/output
  details, and examples if helpful.
