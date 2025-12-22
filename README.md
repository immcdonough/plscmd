# plscmd: Partial Least Squares Analysis for R

An R package for Partial Least Squares (PLS) analysis with permutation testing and bootstrap resampling. This is an R implementation of the MATLAB PLS command line tools originally developed at the Rotman Research Institute.

## Features

- **Six PLS Methods:**
  1. Mean-Centering Task PLS
  2. Non-Rotated Task PLS
  3. Regular Behavior PLS
  4. Multiblock PLS
  5. Non-Rotated Behavior PLS
  6. Non-Rotated Multiblock PLS

- **Statistical Inference:**
  - Permutation testing for significance of latent variables
  - Bootstrap resampling for confidence intervals
  - Split-half validation for reliability assessment

- **Flexible Options:**
  - Four mean-centering types
  - Four correlation modes (Pearson, covariance, cosine, dot product)
  - Stratified and non-stratified bootstrap

## Installation

```r
# Install from local source
install.packages("/path/to/PLS-R", repos = NULL, type = "source")

# Or using devtools
devtools::install_local("/path/to/PLS-R")
```

## Quick Start

### Task PLS Example

```r
library(plscmd)

# Simulate data: 2 groups, 3 conditions, 10 subjects each, 100 brain regions
set.seed(42)
datamat1 <- matrix(rnorm(10 * 3 * 100), nrow = 30, ncol = 100)
datamat2 <- matrix(rnorm(10 * 3 * 100), nrow = 30, ncol = 100)

# Run analysis
result <- pls_analysis(
  datamat_lst = list(datamat1, datamat2),
  num_subj_lst = c(10, 10),
  num_cond = 3,
  option = list(
    method = 1,          # Mean-Centering Task PLS
    num_perm = 500,      # Permutation test
    num_boot = 500       # Bootstrap
  )
)

# View results
print(result)
result$perm_result$sprob  # p-values for LVs
```

### Behavior PLS Example

```r
# Single group with behavioral measures
datamat <- matrix(rnorm(15 * 2 * 50), nrow = 30, ncol = 50)
behavdata <- matrix(rnorm(15 * 2 * 3), nrow = 30, ncol = 3)  # 3 behavior measures

result <- pls_analysis(
  datamat_lst = list(datamat),
  num_subj_lst = c(15),
  num_cond = 2,
  option = list(
    method = 3,          # Behavior PLS
    stacked_behavdata = behavdata,
    num_perm = 500,
    num_boot = 500
  )
)

# Brain-behavior correlations
result$lvcorrs
```

## Data Organization

Data matrices should be organized with subjects in rows (stacked by condition) and variables (voxels/electrodes) in columns:

```
Group 1:
  Condition 1: Subject 1, Subject 2, ..., Subject n
  Condition 2: Subject 1, Subject 2, ..., Subject n
  Condition 3: Subject 1, Subject 2, ..., Subject n
Group 2:
  Condition 1: Subject 1, Subject 2, ..., Subject m
  ...
```

## Key Functions

| Function | Description |
|----------|-------------|
| `pls_analysis()` | Main PLS analysis function |
| `normalize()` | Normalize vectors to unit length |
| `rri_boot_order()` | Generate bootstrap resampling orders |
| `rri_perm_order()` | Generate permutation orders |
| `rri_distrib()` | Calculate bootstrap confidence intervals |

## Result Structure

The `pls_analysis()` function returns a list with:

- `u`: Brain latent variables (saliences)
- `s`: Singular values
- `v`: Design/behavior latent variables
- `usc`: Brain scores
- `vsc`: Design/behavior scores
- `perm_result`: Permutation test results
  - `sprob`: p-values for each LV
- `boot_result`: Bootstrap results
  - `compare_u`: Bootstrap ratios
  - `u_se`: Standard errors

## PLS Methods

### Method 1: Mean-Centering Task PLS
Standard approach that identifies brain patterns differentiating experimental conditions.

### Method 2: Non-Rotated Task PLS
Uses user-specified design contrasts instead of data-driven rotation.

### Method 3: Behavior PLS
Identifies brain patterns that correlate with behavioral/cognitive measures.

### Method 4: Multiblock PLS
Combines task and behavior blocks in a single analysis.

### Methods 5 & 6: Non-Rotated Versions
Behavior and Multiblock PLS with user-specified contrasts.

## References

McIntosh, A.R., & Lobaugh, N.J. (2004). Partial least squares analysis of neuroimaging data: applications and advances. NeuroImage, 23, S250-S263.

Krishnan, A., Williams, L.J., McIntosh, A.R., & Abdi, H. (2011). Partial Least Squares (PLS) methods for neuroimaging: A tutorial and review. NeuroImage, 56(2), 455-475.

## License

GPL-3

## Original MATLAB Code

This R package is a translation of the MATLAB PLS command line tools developed by Jimmy Shen at the Rotman Research Institute, Toronto, Canada.
