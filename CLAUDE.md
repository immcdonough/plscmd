# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Package Overview

**plscmd** is an R package for Partial Least Squares (PLS) analysis of neuroimaging data. It's a translation of the MATLAB PLS command line tools from the Rotman Research Institute.

## Build and Development Commands

```r
# Load package for development (run after any code changes)
devtools::load_all()

# Rebuild Roxygen documentation
devtools::document()

# Run all tests
devtools::test()

# Run specific test file
devtools::test(filter = "pls_analysis")

# Check package integrity (runs R CMD check)
devtools::check()

# Build package
devtools::build()

# Install locally
devtools::install()

# Lint package
lintr::lint_package()
```

## Architecture

### Core Pipeline

The package follows a pipeline architecture centered on `pls_analysis()`:

```
pls_analysis() [R/pls_analysis.R]
    ↓
stacking_datamat() → Combines group matrices
    ↓
rri_get_covcor() [R/covcor.R] → Mean-centering + correlation calculation
    ├── rri_task_mean() → Average across subjects per condition
    └── rri_xcor() → Core correlation function
        └── robust_cor() [R/robust.R] → Optional robust methods
    ↓
SVD computation → Latent variables
    ↓
Score calculation [R/scores.R]
    ↓
[Optional] Permutation testing [R/permutation.R]
    ↓
[Optional] Bootstrap resampling [R/bootstrap.R]
    ↓
pls_result S3 object
```

### Key Files

- **R/pls_analysis.R** - Main entry point, orchestrates all analysis steps
- **R/covcor.R** - `rri_get_covcor()` handles mean-centering and correlation modes
- **R/bootstrap.R** - Bootstrap resampling with `rri_boot_order()`, `rri_distrib()`
- **R/permutation.R** - Permutation testing with `rri_perm_order()`
- **R/scores.R** - LV score calculation with `rri_get_behavscores()`
- **R/robust.R** - Robust correlation methods (biweight, winsorized, percentage bend)
- **R/visualization.R** - ggplot2-based plotting functions
- **R/utils.R** - General utilities (`normalize()`, `percentile()`)
- **R/missing.R** - NA-aware statistics (`missmean()`, `misssvd()`)
- **R/imputation.R** - Missing data handling (`pls_check_missing()`, `pls_impute()`, `pls_analysis_mi()`, Rubin's rules pooling)

### Six PLS Methods

1. Mean-Centering Task PLS (default) - Data-driven condition effects
2. Non-Rotated Task PLS - User-specified design contrasts
3. Regular Behavior PLS - Brain-behavior correlations
4. Multiblock PLS - Combined task + behavior blocks
5. Non-Rotated Behavior PLS - Behavior PLS with contrasts
6. Non-Rotated Multiblock PLS - Multiblock with contrasts

### Two Subject Structure Formats

- **Vector** (`num_subj_lst = c(10, 10)`): Equal subjects per condition per group
- **List of vectors** (`num_subj_lst = list(c(10, 8, 12), c(9, 11, 10))`): Variable subjects per condition (triggers Split-Subject-Block/SSB code paths with `ssb_*` function variants)

### Mean-Centering Types (option$meancentering_type)

- **0**: Remove group-condition means (boosts condition differences)
- **1**: Remove grand condition means (boosts group differences)
- **2**: Remove grand mean
- **3**: Remove all main effects (interaction only)

### Correlation Modes (option$cormode)

- **0**: Pearson correlation (default, robust methods apply only here)
- **1**: Covariance
- **2**: Cosine angle
- **3**: Dot product

## Testing

Tests use testthat (edition 3). Test files are in `tests/testthat/`:
- `test-pls_analysis.R` - Core PLS function tests
- `test-utils.R` - Utility function tests

## Data Format

Data matrices: subjects in rows (stacked by condition), variables in columns.

```
Group 1:
  Condition 1: Subject 1-n rows
  Condition 2: Subject 1-n rows
  ...
Group 2:
  Condition 1: Subject 1-m rows
  ...
```

## Missing Data Handling

### Architecture (R/imputation.R)

```
pls_check_missing() → Diagnose missing data patterns
    ↓
pls_impute() → Single or multiple imputation
    ├── impute_mean() → Column mean imputation
    └── missRanger::missRanger() → Random forest imputation
    ↓
[If MI] pls_analysis_mi() → Run PLS on each imputed dataset
    ↓
pool_pls_results() → Combine results
    ├── align_pls_results() → Procrustes alignment across imputations
    ├── pool_estimates() → Rubin's rules for point estimates
    └── pool_pvalues() → Median p-value method
    ↓
pls_mi_result S3 class
```

### Key Functions

- `pls_check_missing()` - Returns `pls_missing_report` with per-variable/subject stats
- `pls_impute()` - Single (m=1) or multiple (m>1) imputation via mean or RF
- `pls_analysis_mi()` - Full MI workflow with automatic pooling
- `pool_estimates()` - Rubin's rules: θ̄ = mean(θᵢ), T = W̄ + (1+1/m)B
- `pool_pvalues()` - Median of p-values across imputations
- `pool_pls_results()` - Aligns LVs, pools all PLS outputs

### Warning Thresholds

- **>5%**: Informational message
- **>20%**: Warning
- **>50%**: Strong warning (critical)

### Integration with pls_analysis()

- `option$impute_method`: "none" (default, error on NA), "mean", or "rf"
- `option$check_missing`: Print diagnostics (default TRUE)
- If NA detected without impute_method set → informative error

## Key Implementation Notes

- Bootstrap uses Procrustes rotation (`rri_bootprocrust()`) to align with original solution
- Permutation uses two-step approach: within-subject condition shuffle, then between-group subject shuffle
- Robust correlations only work with `cormode = 0` (Pearson mode)
- Stratified bootstrap (default) respects group structure
- MI uses same Procrustes rotation to align LVs across imputations before pooling
