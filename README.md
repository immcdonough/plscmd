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
# Install from GitHub (recommended)
install.packages("devtools")
devtools::install_github("immcdonough/plscmd")

# Or install from local source
install.packages("/path/to/PLS-R", repos = NULL, type = "source")
```

### Dependencies for Visualization

```r
# For publication-quality plots
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

## Visualizing Results

### Summary Table of LV Statistics

```r
# Create summary table of p-values and variance explained
pls_summary_table <- function(result) {
  s <- result$s
  prop_explained <- (s^2) / sum(s^2) * 100
  pvals <- result$perm_result$sprob

  summary_df <- data.frame(
    LV = seq_along(s),
    `Singular Value` = round(s, 3),
    `Variance Explained (%)` = round(prop_explained, 2),
    `Cumulative (%)` = round(cumsum(prop_explained), 2),
    `p-value` = round(pvals, 4),
    Significant = ifelse(pvals < 0.05, "*", "")
  )

  return(summary_df)
}

# Print the summary
summary_table <- pls_summary_table(result)
print(summary_table)

# For a nicer formatted table (requires knitr)
# knitr::kable(summary_table, align = "c")
```

### Publication-Quality Plots with ggplot2

```r
library(ggplot2)
library(patchwork)  # for combining plots

# Set a clean theme for all plots
theme_pls <- function() {
  theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
      axis.title = element_text(face = "bold", size = 12),
      axis.text = element_text(size = 10),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      legend.position = "bottom"
    )
}
```

### Plot 1: Brain-Behavior Correlations with Confidence Intervals

```r
plot_behavior_correlations <- function(result, lv = 1, behav_names = NULL) {
  # Get brain scores and behavior data
  usc <- result$usc[, lv]
  behav <- result$stacked_behavdata
  n_behav <- ncol(behav)

  # Set default behavior names
  if (is.null(behav_names)) {
    behav_names <- paste0("Behavior ", 1:n_behav)
  }

  # Calculate correlations with bootstrap CIs
  correlations <- sapply(1:n_behav, function(b) cor(usc, behav[, b], use = "complete.obs"))

  # Bootstrap for CIs
  n_boot <- 1000
  boot_corrs <- matrix(NA, n_boot, n_behav)
  n <- length(usc)

  for (i in 1:n_boot) {
    idx <- sample(1:n, n, replace = TRUE)
    for (b in 1:n_behav) {
      boot_corrs[i, b] <- cor(usc[idx], behav[idx, b], use = "complete.obs")
    }
  }

  ci_low <- apply(boot_corrs, 2, quantile, probs = 0.025, na.rm = TRUE)
  ci_high <- apply(boot_corrs, 2, quantile, probs = 0.975, na.rm = TRUE)

  # Create data frame for plotting
  plot_df <- data.frame(
    Behavior = factor(behav_names, levels = behav_names),
    Correlation = correlations,
    CI_low = ci_low,
    CI_high = ci_high,
    Significant = sign(ci_low) == sign(ci_high)
  )

  # Create plot
  ggplot(plot_df, aes(x = Behavior, y = Correlation, fill = Significant)) +
    geom_bar(stat = "identity", color = "black", width = 0.7) +
    geom_errorbar(aes(ymin = CI_low, ymax = CI_high),
                  width = 0.2, linewidth = 0.8) +
    geom_hline(yintercept = 0, linetype = "solid", linewidth = 0.5) +
    scale_fill_manual(values = c("TRUE" = "gray40", "FALSE" = "gray80"),
                      labels = c("TRUE" = "p < .05", "FALSE" = "n.s.")) +
    labs(title = sprintf("Brain-Behavior Correlations (LV%d)", lv),
         x = "", y = "Correlation (r)", fill = "") +
    theme_pls() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

# Usage:
# p1 <- plot_behavior_correlations(result, lv = 1,
#         behav_names = c("RT", "Accuracy", "Memory"))
```

### Plot 2: Brain vs Behavior Scores Scatter

```r
plot_brain_behavior_scatter <- function(result, lv = 1, group_labels = NULL) {
  usc <- result$usc[, lv]  # Brain scores
  vsc <- result$vsc[, lv]  # Behavior scores

  # Remove NAs
  valid <- !is.na(usc) & !is.na(vsc)
  usc <- usc[valid]
  vsc <- vsc[valid]

  # Calculate correlation
  r_val <- cor(vsc, usc)
  p_val <- cor.test(vsc, usc)$p.value

  # Create data frame
  plot_df <- data.frame(
    Behavior_Score = vsc,
    Brain_Score = usc
  )

  # Create plot
  ggplot(plot_df, aes(x = Behavior_Score, y = Brain_Score)) +
    geom_point(size = 3, alpha = 0.7, color = "#4575b4", shape = 16) +
    geom_smooth(method = "lm", se = TRUE, color = "black",
                fill = "gray80", linewidth = 1) +
    annotate("text", x = min(vsc) + 0.05 * diff(range(vsc)),
             y = max(usc) - 0.05 * diff(range(usc)),
             label = sprintf("r = %.2f, p = %.3g", r_val, p_val),
             hjust = 0, fontface = "bold", size = 4) +
    labs(title = sprintf("LV%d: Brain vs Behavior Scores", lv),
         x = "Behavior Score", y = "Brain Score") +
    theme_pls()
}

# Usage:
# p2 <- plot_brain_behavior_scatter(result, lv = 1)
```

### Plot 3: Bootstrap Ratios for Brain Regions

```r
plot_bootstrap_ratios <- function(result, lv = 1, region_names = NULL,
                                   threshold = 1.96, top_n = 30) {
  bsr <- result$boot_result$compare_u[, lv]
  se <- result$boot_result$u_se[, lv]

  n_regions <- length(bsr)
  if (is.null(region_names)) {
    region_names <- paste0("Region ", 1:n_regions)
  }

  # Create data frame
  plot_df <- data.frame(
    Region = region_names,
    BSR = bsr,
    SE = se,
    Significant = abs(bsr) > threshold
  )

  # Filter to significant regions or top N by absolute BSR
  if (sum(plot_df$Significant) > 0) {
    plot_df <- plot_df[plot_df$Significant, ]
  } else {
    plot_df <- plot_df[order(abs(plot_df$BSR), decreasing = TRUE)[1:min(top_n, n_regions)], ]
  }

  # Sort by BSR value
  plot_df <- plot_df[order(plot_df$BSR, decreasing = TRUE), ]
  plot_df$Region <- factor(plot_df$Region, levels = plot_df$Region)

  # Create plot
  ggplot(plot_df, aes(x = Region, y = BSR, fill = BSR > 0)) +
    geom_bar(stat = "identity", color = "black", width = 0.7) +
    geom_errorbar(aes(ymin = BSR - SE, ymax = BSR + SE),
                  width = 0.2, linewidth = 0.6) +
    geom_hline(yintercept = c(-threshold, threshold),
               linetype = "dashed", color = "gray40", linewidth = 0.8) +
    geom_hline(yintercept = c(-3.3, 3.3),
               linetype = "dotted", color = "gray40", linewidth = 0.8) +
    scale_fill_manual(values = c("TRUE" = "#d73027", "FALSE" = "#4575b4"),
                      guide = "none") +
    labs(title = sprintf("Bootstrap Ratios (LV%d)", lv),
         subtitle = "Dashed: p < .05, Dotted: p < .001",
         x = "", y = "Bootstrap Ratio") +
    theme_pls() +
    coord_flip()  # Horizontal bars for readability
}

# Usage:
# p3 <- plot_bootstrap_ratios(result, lv = 1,
#         region_names = c("ROI1", "ROI2", ...))
```

### Plot 4: Behavior Correlations with Bootstrap CIs

```r
plot_lv_correlations <- function(result, lv = 1, behav_names = NULL) {
  # Original correlations
  orig_corr <- result$boot_result$orig_corr[, lv]
  ul_corr <- result$boot_result$ulcorr[, lv]
  ll_corr <- result$boot_result$llcorr[, lv]

  n_behav <- length(orig_corr)
  if (is.null(behav_names)) {
    behav_names <- paste0("Behavior ", 1:n_behav)
  }

  # Determine significance (CI doesn't cross zero)
  significant <- sign(ul_corr) == sign(ll_corr)

  # Create data frame
  plot_df <- data.frame(
    Behavior = factor(behav_names, levels = behav_names),
    Correlation = orig_corr,
    CI_low = ll_corr,
    CI_high = ul_corr,
    Significant = significant
  )

  # Sort by correlation
  plot_df <- plot_df[order(plot_df$Correlation, decreasing = TRUE), ]
  plot_df$Behavior <- factor(plot_df$Behavior, levels = plot_df$Behavior)

  # Create plot
  ggplot(plot_df, aes(x = Behavior, y = Correlation, fill = Significant)) +
    geom_bar(stat = "identity", color = "black", width = 0.7) +
    geom_errorbar(aes(ymin = CI_low, ymax = CI_high),
                  width = 0.2, linewidth = 0.6) +
    geom_hline(yintercept = 0, linewidth = 0.5) +
    scale_fill_manual(values = c("TRUE" = "gray40", "FALSE" = "white"),
                      labels = c("TRUE" = "p < .05", "FALSE" = "n.s."),
                      name = "") +
    labs(title = sprintf("LV%d Behavior Correlations", lv),
         x = "", y = "Correlation (r)") +
    theme_pls() +
    coord_flip()
}

# Usage:
# p4 <- plot_lv_correlations(result, lv = 1,
#         behav_names = c("Age", "IQ", "Memory", "Speed"))
```

### Combining Plots into a Figure

```r
# Create combined figure for significant LVs
create_lv_summary_figure <- function(result, lv = 1, behav_names = NULL,
                                      region_names = NULL, save_path = NULL) {
  # Check if LV is significant
  if (result$perm_result$sprob[lv] >= 0.05) {
    message(sprintf("LV%d is not significant (p = %.3f)", lv, result$perm_result$sprob[lv]))
    return(NULL)
  }

  # Create individual plots
  p1 <- plot_behavior_correlations(result, lv, behav_names)
  p2 <- plot_brain_behavior_scatter(result, lv)
  p3 <- plot_bootstrap_ratios(result, lv, region_names)

  # Combine with patchwork
  combined <- (p1 | p2) / p3 +
    plot_annotation(
      title = sprintf("Latent Variable %d Summary (p = %.4f)",
                      lv, result$perm_result$sprob[lv]),
      theme = theme(plot.title = element_text(face = "bold", size = 16, hjust = 0.5))
    )

  # Save if path provided
  if (!is.null(save_path)) {
    ggsave(save_path, combined, width = 12, height = 10, dpi = 300)
    message(sprintf("Figure saved to: %s", save_path))
  }

  return(combined)
}

# Usage:
# fig <- create_lv_summary_figure(result, lv = 1,
#          behav_names = c("RT", "Accuracy", "Memory"),
#          region_names = paste0("ROI_", 1:100),
#          save_path = "LV1_summary.png")
```

### Loop Through All Significant LVs

```r
# Process all significant LVs
plot_all_significant_lvs <- function(result, behav_names = NULL,
                                      region_names = NULL, output_dir = ".") {
  sig_lvs <- which(result$perm_result$sprob < 0.05)

  if (length(sig_lvs) == 0) {
    message("No significant LVs found")
    return(NULL)
  }

  message(sprintf("Found %d significant LV(s): %s",
                  length(sig_lvs), paste(sig_lvs, collapse = ", ")))

  figures <- list()
  for (lv in sig_lvs) {
    save_path <- file.path(output_dir, sprintf("LV%d_summary.png", lv))
    figures[[lv]] <- create_lv_summary_figure(
      result, lv, behav_names, region_names, save_path
    )
  }

  return(figures)
}

# Usage:
# figures <- plot_all_significant_lvs(result,
#              behav_names = my_behav_names,
#              region_names = my_region_names,
#              output_dir = "pls_figures")
```

### Complete Example Workflow

```r
library(plscmd)
library(ggplot2)
library(patchwork)

# Run analysis
result <- pls_analysis(
  datamat_lst = list(brain_data),
  num_subj_lst = c(n_subjects),
  num_cond = n_conditions,
  option = list(
    method = 3,
    stacked_behavdata = behavior_data,
    num_perm = 1000,
    num_boot = 1000
  )
)

# Print summary table
summary_table <- pls_summary_table(result)
print(summary_table)

# Generate figures for all significant LVs
behav_names <- c("Age", "IQ", "Memory", "Processing Speed")
region_names <- paste0("ROI_", 1:ncol(brain_data))

figures <- plot_all_significant_lvs(
  result,
  behav_names = behav_names,
  region_names = region_names,
  output_dir = "output/figures"
)
```

## References

McIntosh, A.R., & Lobaugh, N.J. (2004). Partial least squares analysis of neuroimaging data: applications and advances. NeuroImage, 23, S250-S263.

Krishnan, A., Williams, L.J., McIntosh, A.R., & Abdi, H. (2011). Partial Least Squares (PLS) methods for neuroimaging: A tutorial and review. NeuroImage, 56(2), 455-475.

## License

GPL-3

## Original MATLAB Code

This R package is a translation of the MATLAB PLS command line tools developed by Jimmy Shen at the Rotman Research Institute, Toronto, Canada.
