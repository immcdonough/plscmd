# Example: Task PLS Analysis
#
# This example demonstrates how to run a simple task PLS analysis
# with permutation testing and bootstrap resampling.

library(plscmd)

# Simulate data
# --------------
# Imagine an fMRI study with:
# - 2 groups (e.g., patients vs controls)
# - 3 conditions (e.g., different task demands)
# - 10 subjects per group
# - 100 brain regions (voxels/ROIs)

set.seed(42)

n_subj <- 10
n_cond <- 3
n_vars <- 100  # Number of brain regions

# Create simulated data with a condition effect
# Group 1: Control group
datamat1 <- matrix(rnorm(n_subj * n_cond * n_vars),
                   nrow = n_subj * n_cond,
                   ncol = n_vars)

# Add a condition effect to first 20 regions
for (c in 1:n_cond) {
  rows <- ((c - 1) * n_subj + 1):(c * n_subj)
  datamat1[rows, 1:20] <- datamat1[rows, 1:20] + (c - 2) * 0.5  # Linear effect
}

# Group 2: Patient group - stronger effect
datamat2 <- matrix(rnorm(n_subj * n_cond * n_vars),
                   nrow = n_subj * n_cond,
                   ncol = n_vars)

for (c in 1:n_cond) {
  rows <- ((c - 1) * n_subj + 1):(c * n_subj)
  datamat2[rows, 1:20] <- datamat2[rows, 1:20] + (c - 2) * 1.0  # Stronger effect
}

# Run PLS Analysis
# ----------------
result <- pls_analysis(
  datamat_lst = list(datamat1, datamat2),
  num_subj_lst = c(n_subj, n_subj),
  num_cond = n_cond,
  option = list(
    method = 1,           # Mean-Centering Task PLS
    num_perm = 500,       # Number of permutations
    num_boot = 500,       # Number of bootstrap samples
    clim = 95,            # 95% confidence interval
    meancentering_type = 0,  # Remove group condition means
    verbose = TRUE
  )
)

# View Results
# ------------
print(result)

# Singular values (effect sizes)
cat("\nSingular Values:\n")
print(result$s)

# Permutation p-values
cat("\nPermutation p-values (significance of each LV):\n")
print(result$perm_result$sprob)

# Significant LVs (p < 0.05)
sig_lvs <- which(result$perm_result$sprob < 0.05)
cat("\nSignificant LVs (p < 0.05):", sig_lvs, "\n")

# Bootstrap Ratios (salience / SE)
# Values > 2 or < -2 are typically considered reliable
cat("\nBootstrap Ratios for LV 1 (first 10 brain regions):\n")
print(head(result$boot_result$compare_u[, 1], 10))

# Brain Saliences (which regions contribute to each LV)
cat("\nBrain Saliences (u) for LV 1 (first 10 regions):\n")
print(head(result$u[, 1], 10))

# Design Saliences (condition patterns)
cat("\nDesign Saliences (v) for LV 1:\n")
print(result$v[, 1])
# Note: This shows how conditions load on LV 1
# Rows are organized as: Group1-Cond1, G1-C2, G1-C3, G2-C1, G2-C2, G2-C3

# Brain Scores (individual subject expression of LV)
cat("\nBrain Scores dimensions:", dim(result$usc), "\n")
cat("(rows = subjects*conditions*groups, cols = LVs)\n")

# Confidence Intervals for Brain Scores
cat("\nConfidence intervals for brain scores (LV 1):\n")
cat("Upper limit:\n")
print(result$boot_result$ulusc[, 1])
cat("Lower limit:\n")
print(result$boot_result$llusc[, 1])

# Example: Visualize results (requires ggplot2)
# ---------------------------------------------
# if (require(ggplot2)) {
#   # Plot bootstrap ratios for first LV
#   bsr <- result$boot_result$compare_u[, 1]
#   df <- data.frame(Region = 1:length(bsr), BSR = bsr)
#   df$Significant <- abs(df$BSR) > 2
#
#   ggplot(df, aes(x = Region, y = BSR, fill = Significant)) +
#     geom_bar(stat = "identity") +
#     geom_hline(yintercept = c(-2, 2), linetype = "dashed") +
#     labs(title = "Bootstrap Ratios - LV 1",
#          x = "Brain Region", y = "Bootstrap Ratio") +
#     theme_minimal()
# }
