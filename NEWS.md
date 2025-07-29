# TwoPhaseCorR 1.1.1

_Resubmission_: This update addresses a namespace masking issue that could cause  
`ginv` from **MASS** to be unavailable in examples under certain environments.  
The fix is applied by explicitly importing `ginv` and reordering package Imports.

## New Features
- Added tolerance-based classification for Efficiency Factor (E ≈ 1).
- Function `TwoPhaseDesign()` now prints Efficiency Factor (E) at the specified `rho` directly to the console.
- Enhanced `eff_plot` with improved axis labels, title, and visual clarity using `ggplot2`.

## Enhancements
- Improved output structure and modular layout for better interpretability.
- Streamlined discovery of efficient designs over a range of `v` and `rho`.
- Enhanced numerical stability in efficiency computations.

## Bug Fixes
- Fixed namespace masking issue for `ginv` by explicitly importing from **MASS**.
- Fixed global variable declaration issues for `Class`, `Rho`, and `Efficiency` to meet CRAN compliance (`utils::globalVariables()` used).
- Removed problematic Unicode characters (e.g., `≥`, `≈`) from `.Rd` files that previously caused LaTeX PDF manual generation errors.

## Documentation
- Added detailed argument and return value descriptions in function documentation.
- Cleaned and standardized `.Rd` files for LaTeX compatibility.
- Resolved all CRAN LaTeX PDF warnings by ensuring proper rendering under `pdflatex`.

## CRAN Readiness
- Successfully passes `R CMD check` with **no errors, warnings, or notes**.
- Verified compatibility with R 4.5.1 and R-devel across platforms.
