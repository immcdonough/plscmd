# Example: Behavior PLS Analysis
#
# This example demonstrates how to run a behavior PLS analysis
# to identify brain patterns that correlate with behavioral measures.

library(plscmd)

# Simulate data
# --------------
# Study design:
# - 1 group of subjects
# - 3 conditions (e.g., easy, medium, hard task)
# - 15 subjects
# - 100 brain regions
# - 2 behavioral measures (e.g., reaction time, accuracy)

set.seed(123)

n_subj <- 15
n_cond <- 3
n_vars <- 100
n_behav <- 2

# Create brain data
datamat <- matrix(rnorm(n_subj * n_cond * n_vars),
                  nrow = n_subj * n_cond,
                  ncol = n_vars)

# Create behavioral data with brain-behavior relationship
behavdata <- matrix(rnorm(n_subj * n_cond * n_behav),
                    nrow = n_subj * n_cond,
                    ncol = n_behav)

# Add correlation between brain activity and behavior
# Regions 1-15 correlate with behavior measure 1
# Regions 16-30 correlate with behavior measure 2
for (i in 1:(n_subj * n_cond)) {
  datamat[i, 1:15] <- datamat[i, 1:15] + 0.5 * behavdata[i, 1]
  datamat[i, 16:30] <- datamat[i, 16:30] + 0.5 * behavdata[i, 2]
}

# Run Behavior PLS Analysis
# -------------------------
result <- pls_analysis(
  datamat_lst = list(datamat),
  num_subj_lst = c(n_subj),
  num_cond = n_cond,
  option = list(
    method = 3,           # Regular Behavior PLS
    stacked_behavdata = behavdata,
    num_perm = 500,
    num_boot = 500,
    cormode = 0,          # Pearson correlation
    verbose = TRUE
  )
)

# View Results
# ------------
print(result)

# Correlation maps (brain-behavior correlations)
cat("\nDatamat correlations list (one per group):\n")
cat("Shape:", dim(result$datamatcorrs_lst[[1]]), "\n")
cat("(rows = n_behav * n_cond, cols = n_vars)\n")

# LV correlations (correlation of behavior with brain scores)
cat("\nLV Correlations shape:", dim(result$lvcorrs), "\n")
cat("These show how behavioral measures correlate with brain LV scores\n")
cat("across conditions.\n\n")

# Significant LVs
sig_lvs <- which(result$perm_result$sprob < 0.05)
cat("Significant LVs (p < 0.05):",
    if (length(sig_lvs) > 0) sig_lvs else "None", "\n")

# Bootstrap results
if (!is.null(result$boot_result)) {
  cat("\nBootstrap Ratios for first LV (brain saliences):\n")
  bsr <- result$boot_result$compare_u[, 1]

  # Regions with reliable saliences (|BSR| > 2)
  reliable <- which(abs(bsr) > 2)
  cat("Regions with reliable saliences (|BSR| > 2):",
      head(reliable, 10), "...\n")

  # Confidence intervals for LV correlations
  cat("\nConfidence intervals for LV correlations (first LV):\n")
  cat("Original correlations:\n")
  print(result$boot_result$orig_corr[, 1])
  cat("Upper 95% CI:\n")
  print(result$boot_result$ulcorr[, 1])
  cat("Lower 95% CI:\n")
  print(result$boot_result$llcorr[, 1])
}

# Interpreting Behavior PLS Results
# ---------------------------------
# 1. Brain saliences (u): Show which regions contribute to the LV
#    - Positive saliences: Increase activity correlates positively with
#      the behavioral pattern
#    - Negative saliences: Decrease activity correlates with behavioral pattern
#
# 2. Behavior saliences (v): Show how behavioral measures load on each LV
#    - The pattern of behavior measure loadings indicates what combination
#      of behaviors the brain pattern relates to
#
# 3. LV correlations: Correlation between behavior measures and brain scores
#    - Shows the strength of brain-behavior relationship for each condition
#
# 4. Bootstrap ratios: Reliability of brain saliences
#    - |BSR| > 2 indicates reliable contribution to the LV
