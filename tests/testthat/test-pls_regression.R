# Tests for PLS regression functions

# ===========================================================================
# Basic Algorithm Tests
# ===========================================================================

test_that("pls_simpls produces correct dimensions", {
  set.seed(42)
  n <- 50
  p <- 20
  q <- 3
  ncomp <- 5

  X <- matrix(rnorm(n * p), n, p)
  Y <- matrix(rnorm(n * q), n, q)

  result <- pls_regression(X, Y, ncomp = ncomp, algorithm = "simpls", verbose = FALSE)

  expect_s3_class(result, "pls_regression_result")
  expect_equal(dim(result$X_loadings), c(p, ncomp))
  expect_equal(dim(result$Y_loadings), c(q, ncomp))
  expect_equal(dim(result$X_scores), c(n, ncomp))
  expect_equal(dim(result$Y_scores), c(n, ncomp))
  expect_equal(dim(result$X_weights), c(p, ncomp))
  expect_equal(dim(result$coefficients), c(p, q, ncomp))
  expect_equal(dim(result$vip), c(p, ncomp))
  expect_equal(result$ncomp, ncomp)
  expect_equal(result$n, n)
  expect_equal(result$p, p)
  expect_equal(result$q, q)
})


test_that("pls_nipals produces correct dimensions", {
  set.seed(42)
  n <- 50
  p <- 20
  q <- 3
  ncomp <- 5

  X <- matrix(rnorm(n * p), n, p)
  Y <- matrix(rnorm(n * q), n, q)

  result <- pls_regression(X, Y, ncomp = ncomp, algorithm = "nipals", verbose = FALSE)

  expect_s3_class(result, "pls_regression_result")
  expect_equal(dim(result$X_loadings), c(p, ncomp))
  expect_equal(dim(result$Y_loadings), c(q, ncomp))
  expect_equal(result$algorithm, "nipals")
})


test_that("SIMPLS and NIPALS produce similar R2", {
  set.seed(42)
  n <- 100
  p <- 10

  # Create data with known structure
  X <- matrix(rnorm(n * p), n, p)
  beta <- c(rep(1, 3), rep(0, p - 3))
  Y <- X %*% beta + rnorm(n, sd = 0.5)

  result_simpls <- pls_regression(X, Y, ncomp = 3, algorithm = "simpls", verbose = FALSE)
  result_nipals <- pls_regression(X, Y, ncomp = 3, algorithm = "nipals", verbose = FALSE)

  # R2 should be similar (within 10% relative difference)
  expect_equal(result_simpls$R2Y_cum[3], result_nipals$R2Y_cum[3], tolerance = 0.1)
})


test_that("R2 values are between 0 and 1", {
  set.seed(42)
  n <- 50
  p <- 10

  X <- matrix(rnorm(n * p), n, p)
  Y <- X[, 1:3] %*% c(1, 1, 1) + rnorm(n, sd = 0.5)

  result <- pls_regression(X, Y, ncomp = 5, verbose = FALSE)

  expect_true(all(result$R2X >= 0 & result$R2X <= 1))
  expect_true(all(result$R2Y >= 0 & result$R2Y <= 1))
  expect_true(all(result$R2Y_cum >= 0 & result$R2Y_cum <= 1))
  # R2Y_cum should be non-decreasing
  expect_true(all(diff(result$R2Y_cum) >= -1e-10))
})


# ===========================================================================
# Prediction Tests
# ===========================================================================

test_that("predict.pls_regression_result works correctly", {
  set.seed(42)
  n <- 100
  p <- 10

  X <- matrix(rnorm(n * p), n, p)
  beta <- c(rep(1, 3), rep(0, p - 3))
  Y <- X %*% beta + rnorm(n, sd = 0.5)

  # Train on first 80, test on last 20
  train_idx <- 1:80
  test_idx <- 81:100

  fit <- pls_regression(X[train_idx, ], Y[train_idx], ncomp = 3, verbose = FALSE)
  Y_pred <- predict(fit, X[test_idx, ])

  expect_equal(length(Y_pred), 20)

  # Predictions should be correlated with true values
  cor_test <- cor(Y_pred, Y[test_idx])
  expect_true(cor_test > 0.5)
})


test_that("predict works with matrix Y", {
  set.seed(42)
  n <- 100
  p <- 10
  q <- 3

  X <- matrix(rnorm(n * p), n, p)
  Y <- matrix(rnorm(n * q), n, q)
  Y[, 1] <- Y[, 1] + X[, 1]

  fit <- pls_regression(X, Y, ncomp = 3, verbose = FALSE)
  Y_pred <- predict(fit, X[1:10, ])

  expect_equal(dim(Y_pred), c(10, q))
})


test_that("predict respects ncomp argument", {
  set.seed(42)
  n <- 50
  p <- 10

  X <- matrix(rnorm(n * p), n, p)
  Y <- X[, 1:3] %*% c(1, 1, 1) + rnorm(n)

  fit <- pls_regression(X, Y, ncomp = 5, verbose = FALSE)

  pred1 <- predict(fit, X[1:5, ], ncomp = 1)
  pred5 <- predict(fit, X[1:5, ], ncomp = 5)

  # Different ncomp should give different predictions
  expect_false(all(pred1 == pred5))
})


# ===========================================================================
# VIP Tests
# ===========================================================================

test_that("VIP scores identify important variables", {
  set.seed(42)
  n <- 100
  p <- 10

  # Create data where first 3 variables are truly predictive
  X <- matrix(rnorm(n * p), n, p)
  beta <- c(rep(2, 3), rep(0, p - 3))
  Y <- X %*% beta + rnorm(n, sd = 0.5)

  result <- pls_regression(X, Y, ncomp = 3, verbose = FALSE)
  vip <- result$vip[, 3]

  # First 3 variables should have highest VIP
  top3 <- order(vip, decreasing = TRUE)[1:3]
  expect_true(all(top3 %in% 1:3))

  # VIP of important variables should be > 1
  expect_true(all(vip[1:3] > 1))
})


test_that("VIP values are non-negative", {
  set.seed(42)
  n <- 50
  p <- 10

  X <- matrix(rnorm(n * p), n, p)
  Y <- rnorm(n)

  result <- pls_regression(X, Y, ncomp = 3, verbose = FALSE)

  expect_true(all(result$vip >= 0))
})


# ===========================================================================
# Cross-Validation Tests
# ===========================================================================

test_that("K-fold CV runs correctly", {
  set.seed(42)
  n <- 100
  p <- 10

  X <- matrix(rnorm(n * p), n, p)
  Y <- X[, 1:3] %*% c(1, 1, 1) + rnorm(n, sd = 0.5)

  cv_result <- pls_cv(X, Y, ncomp_max = 5, method = "kfold", k = 5, verbose = FALSE)

  expect_s3_class(cv_result, "pls_cv_result")
  expect_equal(length(cv_result$RMSEP), 5)
  expect_equal(length(cv_result$Q2), 5)
  expect_true(cv_result$optimal_ncomp >= 1)
  expect_true(cv_result$optimal_ncomp <= 5)
  expect_equal(cv_result$k, 5)
})


test_that("LOO CV runs correctly", {
  set.seed(42)
  n <- 30  # Smaller for LOO speed
  p <- 5

  X <- matrix(rnorm(n * p), n, p)
  Y <- X[, 1:2] %*% c(1, 1) + rnorm(n, sd = 0.5)

  cv_result <- pls_cv(X, Y, ncomp_max = 3, method = "loo", verbose = FALSE)

  expect_equal(cv_result$k, n)  # LOO has k = n
  expect_true(all(cv_result$R2_cv >= -0.5))  # Can be slightly negative
  expect_true(all(cv_result$Q2 >= -0.5))
})


test_that("CV metrics are consistent", {
  set.seed(42)
  n <- 50
  p <- 10

  X <- matrix(rnorm(n * p), n, p)
  Y <- X[, 1:3] %*% c(1, 1, 1) + rnorm(n, sd = 0.3)

  cv_result <- pls_cv(X, Y, ncomp_max = 5, k = 5, verbose = FALSE)

  # Q2 and R2_cv should be equal
  expect_equal(cv_result$Q2, cv_result$R2_cv)

  # RMSEP should be positive
  expect_true(all(cv_result$RMSEP > 0))
})


# ===========================================================================
# Bootstrap Tests
# ===========================================================================

test_that("Bootstrap produces CIs", {
  set.seed(42)
  n <- 50
  p <- 5

  X <- matrix(rnorm(n * p), n, p)
  Y <- X[, 1:2] %*% c(1, 1) + rnorm(n, sd = 0.5)

  result <- pls_regression_boot(X, Y, ncomp = 2, num_boot = 50, verbose = FALSE)

  expect_true(!is.null(result$boot_result))
  expect_equal(result$boot_result$num_boot, 50)
  expect_true(all(result$boot_result$B_ci_low <= result$boot_result$B_ci_high))
  expect_true(all(result$boot_result$B_se >= 0))
  expect_true(all(result$boot_result$vip_ci_low <= result$boot_result$vip_ci_high))
})


test_that("Bootstrap ratios identify significant coefficients", {
  set.seed(42)
  n <- 100
  p <- 5

  X <- matrix(rnorm(n * p), n, p)
  Y <- X[, 1] * 3 + rnorm(n, sd = 0.5)  # Strong effect on X1

  result <- pls_regression_boot(X, Y, ncomp = 2, num_boot = 100, verbose = FALSE)

  # First coefficient should have large bootstrap ratio
  compare_B <- result$boot_result$compare_B
  expect_true(abs(compare_B[1, 1]) > 2)  # Should be significant
})


# ===========================================================================
# Permutation Tests
# ===========================================================================

test_that("Permutation test produces p-value", {
  set.seed(42)
  n <- 50
  p <- 5

  # Strong signal
  X <- matrix(rnorm(n * p), n, p)
  Y <- X[, 1] * 2 + rnorm(n, sd = 0.5)

  result <- pls_regression_perm(X, Y, ncomp = 2, num_perm = 50, verbose = FALSE)

  expect_true(!is.null(result$perm_result))
  expect_true(result$perm_result$p_value >= 0)
  expect_true(result$perm_result$p_value <= 1)

  # Strong signal should give low p-value
  expect_true(result$perm_result$p_value < 0.1)
})


test_that("Permutation with no signal gives high p-value", {
  set.seed(42)
  n <- 50
  p <- 5

  # No signal - X and Y are independent
  X <- matrix(rnorm(n * p), n, p)
  Y <- rnorm(n)

  result <- pls_regression_perm(X, Y, ncomp = 1, num_perm = 50, verbose = FALSE)

  # No signal should give high p-value (>0.05 usually, but due to randomness allow >0.01)
  expect_true(result$perm_result$p_value > 0.01)
})


# ===========================================================================
# Coefficient Extraction Tests
# ===========================================================================

test_that("coef.pls_regression_result works correctly", {
  set.seed(42)
  n <- 50
  p <- 10

  X <- matrix(rnorm(n * p), n, p)
  Y <- X[, 1:3] %*% c(1, 1, 1) + rnorm(n)

  fit <- pls_regression(X, Y, ncomp = 3, verbose = FALSE)
  B <- coef(fit)

  expect_equal(length(B), p)  # Single Y returns vector

  # Different ncomp should give different coefficients
  B1 <- coef(fit, ncomp = 1)
  B3 <- coef(fit, ncomp = 3)
  expect_false(all(B1 == B3))
})


test_that("coef with matrix Y returns matrix", {
  set.seed(42)
  n <- 50
  p <- 10
  q <- 3

  X <- matrix(rnorm(n * p), n, p)
  Y <- matrix(rnorm(n * q), n, q)

  fit <- pls_regression(X, Y, ncomp = 3, verbose = FALSE)
  B <- coef(fit)

  expect_equal(dim(B), c(p, q))
})


# ===========================================================================
# Edge Cases
# ===========================================================================

test_that("Single Y variable works", {
  set.seed(42)
  n <- 50
  p <- 10

  X <- matrix(rnorm(n * p), n, p)
  Y <- rnorm(n)  # Vector, not matrix

  result <- pls_regression(X, Y, ncomp = 3, verbose = FALSE)

  expect_equal(result$q, 1)
  expect_equal(nrow(result$Y_loadings), 1)
})


test_that("More components than rank handled", {
  set.seed(42)
  n <- 20
  p <- 50  # More predictors than observations

  X <- matrix(rnorm(n * p), n, p)
  Y <- rnorm(n)

  # Should handle gracefully
  result <- pls_regression(X, Y, ncomp = 15, verbose = FALSE)

  expect_true(result$ncomp <= min(n - 1, p))
})


test_that("Scaling and centering options work", {
  set.seed(42)
  n <- 50
  p <- 10

  X <- matrix(rnorm(n * p), n, p)
  Y <- rnorm(n)

  # With scaling
  result_scaled <- pls_regression(X, Y, ncomp = 3, scale = TRUE, center = TRUE, verbose = FALSE)
  expect_true(!is.null(result_scaled$X_scale))
  expect_true(!is.null(result_scaled$Y_scale))

  # Without scaling
  result_unscaled <- pls_regression(X, Y, ncomp = 3, scale = FALSE, center = TRUE, verbose = FALSE)
  expect_null(result_unscaled$X_scale)
  expect_null(result_unscaled$Y_scale)
})


# ===========================================================================
# Robust Methods Tests
# ===========================================================================

test_that("Robust methods work in pls_regression", {
  set.seed(42)
  n <- 50
  p <- 5

  X <- matrix(rnorm(n * p), n, p)
  Y <- X[, 1] + rnorm(n)

  # Add outliers
  X[1:3, ] <- X[1:3, ] + 10

  # Spearman should complete without error
  result_spearman <- pls_regression(X, Y, ncomp = 2, robust_method = "spearman",
                                     verbose = FALSE)
  expect_s3_class(result_spearman, "pls_regression_result")

  # Winsorized should complete without error
  result_winsorized <- pls_regression(X, Y, ncomp = 2, robust_method = "winsorized",
                                       verbose = FALSE)
  expect_s3_class(result_winsorized, "pls_regression_result")
})


test_that("Robust methods only work with SIMPLS", {
  set.seed(42)
  n <- 30
  p <- 5

  X <- matrix(rnorm(n * p), n, p)
  Y <- rnorm(n)

  # NIPALS with robust method should warn
  expect_warning(
    pls_regression(X, Y, ncomp = 2, algorithm = "nipals",
                    robust_method = "spearman", verbose = FALSE)
  )
})


# ===========================================================================
# Print and Summary Tests
# ===========================================================================

test_that("print.pls_regression_result works", {
  set.seed(42)
  X <- matrix(rnorm(50 * 10), 50, 10)
  Y <- rnorm(50)

  result <- pls_regression(X, Y, ncomp = 3, verbose = FALSE)

  # Should not error
  expect_output(print(result), "PLS Regression Result")
  expect_output(print(result), "Components: 3")
})


test_that("summary.pls_regression_result works", {
  set.seed(42)
  X <- matrix(rnorm(50 * 10), 50, 10)
  Y <- rnorm(50)

  result <- pls_regression(X, Y, ncomp = 3, verbose = FALSE)

  # Should not error
  expect_output(summary(result), "PLS Regression Summary")
  expect_output(summary(result), "VIP Scores")
})


# ===========================================================================
# MI Tests
# ===========================================================================

test_that("pls_regression_mi works with missing data", {
  skip_if_not_installed("missRanger")

  set.seed(42)
  n <- 50
  p <- 10

  X <- matrix(rnorm(n * p), n, p)
  Y <- X[, 1:3] %*% c(1, 1, 1) + rnorm(n)

  # Add 10% missing
  X[sample(length(X), 0.1 * length(X))] <- NA

  result <- pls_regression_mi(X, Y, ncomp = 3, m = 3,
                               impute_method = "mean", verbose = FALSE)

  expect_s3_class(result, "pls_mi_regression_result")
  expect_s3_class(result, "pls_regression_result")
  expect_equal(result$m, 3)
  expect_true(!is.null(result$coefficients_pooled))
})


test_that("pls_regression_mi with no missing data runs standard PLS", {
  set.seed(42)
  n <- 30
  p <- 5

  X <- matrix(rnorm(n * p), n, p)
  Y <- rnorm(n)

  result <- pls_regression_mi(X, Y, ncomp = 2, verbose = FALSE)

  # Should return standard pls_regression_result, not MI result
  expect_s3_class(result, "pls_regression_result")
  expect_false(inherits(result, "pls_mi_regression_result"))
})


# ===========================================================================
# CV Plotting Test
# ===========================================================================

test_that("plot.pls_cv_result returns ggplot object", {
  skip_if_not_installed("ggplot2")

  set.seed(42)
  n <- 50
  p <- 10

  X <- matrix(rnorm(n * p), n, p)
  Y <- X[, 1:3] %*% c(1, 1, 1) + rnorm(n)

  cv_result <- pls_cv(X, Y, ncomp_max = 5, k = 5, verbose = FALSE)

  p <- plot(cv_result)
  expect_s3_class(p, "ggplot")
})
