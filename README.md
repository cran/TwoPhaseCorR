# TwoPhaseCorR

**Two-Phase Experimental Designs under Correlated Observations**

`TwoPhaseCorR` is an R package for constructing and evaluating two-phase experimental designs, particularly in the presence of correlated observations between phases. It supports the computation of information matrices and Canonical Efficiency Factors (CEF), helping researchers identify efficient designs under various intra-block correlation scenarios.

---

## Features

- Constructs cyclic two-phase designs for a specified number of treatments (`v`).
- Models intra-block correlation via a correlation parameter (`rho`).
- Computes information matrices for treatment effects and interactions.
- Estimates Canonical Efficiency Factors (CEF) for assessing design efficiency.
- Visualizes the relationship between CEF and `rho` to guide optimal design choices.

---

## Installation

### From CRAN (after approval)

```r
install.packages("TwoPhaseCorR")
