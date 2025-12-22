# plscmd: Partial Least Squares Analysis for R

An R package for Partial Least Squares (PLS) analysis with permutation testing and bootstrap resampling. This is an R implementation of the MATLAB PLS command line tools originally developed at the Rotman Research Institute.

## Requirements

- **R version:** >= 4.0.0
- **Core dependencies:** stats, utils (included with R)
- **Visualization dependencies:** ggplot2 (>= 3.4.0), patchwork (>= 1.1.0)

Tested on R 4.3.x and R 4.4.x.

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

- **Built-in Visualization:**
  - Publication-quality plots with ggplot2
  - Customizable colors and font sizes
  - Summary tables and multi-panel figures

## Installation

```r
# Install from GitHub (recommended)
install.packages("devtools")
devtools::install_github("immcdonough/plscmd")

# Install visualization dependencies
install.packages(c("ggplot2", "patchwork"))
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

### Using Robust Correlations

For data with outliers or influential points, use robust correlation methods:

```r
# Behavior PLS with biweight midcorrelation (robust to outliers)
result <- pls_analysis(
  datamat_lst = list(brain_data),
  num_subj_lst = c(n_subjects),
  num_cond = n_conditions,
  option = list(
    method = 3,
    stacked_behavdata = behavior_data,
    num_perm = 500,
    num_boot = 500,
    robust_method = "biweight"  # Use robust correlation
  )
)
```

**Available Robust Methods:**

| Method | Description | Best For |
|--------|-------------|----------|
| `"none"` | Standard Pearson correlation (default) | Clean data |
| `"spearman"` | Spearman rank correlation | Non-linear relationships, ordinal data |
| `"winsorized"` | Trims extreme values before correlation | Moderate outliers |
| `"biweight"` | Biweight midcorrelation | Strong outliers, high breakdown point |
| `"percentage_bend"` | Percentage bend correlation | Balance of robustness and efficiency |

**Additional Options:**
- `robust_trim`: Trim proportion for winsorized method (default: 0.1 = 10% from each tail)
- `robust_beta`: Bend constant for percentage bend (default: 0.2)

```r
# Winsorized correlation with 15% trimming
result <- pls_analysis(
  ...,
  option = list(
    method = 3,
    robust_method = "winsorized",
    robust_trim = 0.15
  )
)
```

**Note:** Robust methods only apply when `cormode = 0` (Pearson correlation mode). They are ignored for covariance, cosine, or dot product modes.

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
| `pls_summary_table()` | Create summary table of LV statistics |
| `plot_lv_summary()` | Create multi-panel summary figure |
| `plot_all_significant_lvs()` | Generate figures for all significant LVs |
| `plot_behavior_correlations()` | Bar plot of brain-behavior correlations |
| `plot_brain_behavior_scatter()` | Scatter plot of brain vs behavior scores |
| `plot_bootstrap_ratios()` | Bar plot of bootstrap ratios |
| `plot_lv_correlations()` | Bar plot of LV correlations |
| `robust_cor()` | Compute robust correlations |
| `biweight_midcor()` | Biweight midcorrelation |

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

## Visualizing Results

### Summary Table

```r
# Create and print summary table
summary_table <- pls_summary_table(result)
print(summary_table)

#   LV Singular_Value Variance_Explained_Pct Cumulative_Pct P_Value Significant
# 1  1         15.234                  45.23          45.23  0.0020           *
# 2  2          8.127                  24.56          69.79  0.0340           *
# 3  3          5.891                  18.12          87.91  0.2150
```

### Quick Visualization of All Significant LVs

```r
library(ggplot2)
library(patchwork)

# Define variable names
behav_names <- c("Age", "IQ", "Memory", "Processing Speed")
region_names <- paste0("ROI_", 1:ncol(brain_data))

# Generate and save figures for all significant LVs
figures <- plot_all_significant_lvs(
  result,
  behav_names = behav_names,
  region_names = region_names,
  output_dir = "figures"
)
```

### Individual Plots

#### Brain-Behavior Correlations with Bootstrap CIs

```r
p1 <- plot_behavior_correlations(
  result,
  lv = 1,
  behav_names = c("Age", "IQ", "Memory", "Speed")
)
print(p1)
```

#### Brain vs Behavior Scores Scatter

```r
p2 <- plot_brain_behavior_scatter(result, lv = 1)
print(p2)
```

#### Bootstrap Ratios for Brain Regions

```r
# Only shows regions where CI lower bound exceeds threshold (conservative)
p3 <- plot_bootstrap_ratios(
  result,
  lv = 1,
  region_names = paste0("ROI_", 1:100),
  threshold = 1.96  # p < .05
)
print(p3)
```

#### LV Correlations from Bootstrap

```r
p4 <- plot_lv_correlations(
  result,
  lv = 1,
  behav_names = c("Age", "IQ", "Memory", "Speed")
)
print(p4)
```

### Combined Summary Figure

```r
# Create multi-panel figure for a single LV
fig <- plot_lv_summary(
  result,
  lv = 1,
  behav_names = c("Age", "IQ", "Memory", "Speed"),
  region_names = paste0("ROI_", 1:100),
  save_path = "LV1_summary.png",
  width = 12,
  height = 10,
  dpi = 300
)
```

### Customizing Plot Appearance

#### Custom Colors

```r
# Get default colors and modify
my_colors <- pls_colors()
my_colors$positive <- "#E41A1C"      # Red for positive BSR
my_colors$negative <- "#377EB8"      # Blue for negative BSR
my_colors$significant <- "#4DAF4A"   # Green for significant
my_colors$nonsignificant <- "#E0E0E0" # Light gray for non-significant
my_colors$points <- "#984EA3"        # Purple for scatter points

# Use custom colors in plots
p <- plot_bootstrap_ratios(result, lv = 1, colors = my_colors)
```

#### Custom Font Sizes

```r
# Define font sizes
my_fonts <- list(
  base_size = 14,
  title_size = 18,
  axis_title_size = 14,
  axis_text_size = 12
)

# Use custom fonts in plots
p <- plot_behavior_correlations(
  result,
  lv = 1,
  font_sizes = my_fonts
)

# Or use the theme directly for further customization
p + theme_pls(base_size = 16, title_size = 20)
```

#### Full Customization Example

```r
library(ggplot2)
library(patchwork)

# Custom colors - colorblind-friendly palette
my_colors <- pls_colors()
my_colors$positive <- "#D55E00"   # Orange
my_colors$negative <- "#0072B2"   # Blue
my_colors$significant <- "#009E73" # Green
my_colors$points <- "#CC79A7"     # Pink

# Larger fonts for presentations
presentation_fonts <- list(
  base_size = 16,
  title_size = 20,
  axis_title_size = 16,
  axis_text_size = 14
)

# Generate all figures with custom styling
figures <- plot_all_significant_lvs(
  result,
  behav_names = c("Reaction Time", "Accuracy", "Working Memory"),
  region_names = roi_labels,
  colors = my_colors,
  font_sizes = presentation_fonts,
  output_dir = "presentation_figures",
  width = 14,
  height = 12,
  dpi = 300
)
```

### Available Color Options

The `pls_colors()` function returns a list with the following customizable colors:

| Color Name | Default | Used For |
|------------|---------|----------|
| `positive` | #d73027 (red) | Positive bootstrap ratios |
| `negative` | #4575b4 (blue) | Negative bootstrap ratios |
| `significant` | gray40 | Significant correlations |
| `nonsignificant` | gray80 | Non-significant correlations |
| `points` | #4575b4 (blue) | Scatter plot points |
| `line` | black | Regression lines |
| `ci_fill` | gray80 | Confidence interval bands |
| `bar_fill` | gray50 | General bar fill |
| `bar_outline` | black | Bar outlines |

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
