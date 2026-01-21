# Test utility functions

test_that("normalize works correctly", {
  # Test column normalization (default)
  mat <- matrix(c(3, 4, 1, 0), nrow = 2, ncol = 2)
  result <- normalize(mat)

  # First column should have norm 1
  expect_equal(sqrt(sum(result[, 1]^2)), 1, tolerance = 1e-10)

  # Second column (was [1, 0]) should be [1, 0] after normalization
  expect_equal(result[1, 2], 1, tolerance = 1e-10)
  expect_equal(result[2, 2], 0, tolerance = 1e-10)

  # Test row normalization
  result_row <- normalize(mat, dim = 2)
  expect_equal(sqrt(sum(result_row[1, ]^2)), 1, tolerance = 1e-10)

  # Test zero vector handling
  mat_zero <- matrix(c(0, 0, 1, 2), nrow = 2, ncol = 2)
  result_zero <- normalize(mat_zero)
  expect_equal(result_zero[, 1], c(0, 0))
})


test_that("percentile works correctly", {
  x <- 1:100
  expect_equal(percentile(x, 50), 50.5, tolerance = 0.1)
  expect_equal(percentile(x, 0), 1, tolerance = 0.1)
  expect_equal(percentile(x, 100), 100, tolerance = 0.1)

  # Multiple percentiles
  result <- percentile(x, c(25, 50, 75))
  expect_equal(length(result), 3)
  expect_true(result[1] < result[2])
  expect_true(result[2] < result[3])
})


test_that("cumulative_gaussian and inverse work correctly", {
  # Test standard normal
  expect_equal(cumulative_gaussian(0), 0.5, tolerance = 1e-10)
  expect_equal(cumulative_gaussian(1.96), 0.975, tolerance = 0.001)
  expect_equal(cumulative_gaussian(-1.96), 0.025, tolerance = 0.001)

  # Test inverse
  expect_equal(cumulative_gaussian_inv(0.5), 0, tolerance = 1e-10)
  expect_equal(cumulative_gaussian_inv(0.975), 1.96, tolerance = 0.01)

  # Round-trip
  x <- c(-2, -1, 0, 1, 2)
  expect_equal(cumulative_gaussian_inv(cumulative_gaussian(x)), x, tolerance = 1e-10)

  # Boundary cases
  expect_equal(cumulative_gaussian_inv(0), -Inf)
  expect_equal(cumulative_gaussian_inv(1), Inf)
})


test_that("rri_xcor computes correlations correctly", {
  # Simple test
  x <- matrix(c(1, 2, 3, 4, 5, 6), nrow = 3, ncol = 2)
  y <- matrix(c(1, 2, 3, 2, 4, 6), nrow = 3, ncol = 2)

  result <- rri_xcor(x, y, cormode = 0)
  expect_equal(dim(result), c(2, 2))

  # Perfect correlation case
  result_perfect <- rri_xcor(x, x, cormode = 0)
  expect_equal(diag(result_perfect), c(1, 1), tolerance = 1e-10)

  # Covariance mode
  result_cov <- rri_xcor(x, y, cormode = 2)
  expect_equal(dim(result_cov), c(2, 2))
})


test_that("rri_task_mean computes condition means correctly", {
  # 2 conditions, 3 subjects each
  # Column 1: condition 1 = [1,2,3], condition 2 = [4,5,6]
  # Column 2: condition 1 = [7,8,9], condition 2 = [10,11,12]
  datamat <- matrix(c(
    1, 2, 3, 4, 5, 6,      # Column 1
    7, 8, 9, 10, 11, 12    # Column 2
  ), nrow = 6, ncol = 2)

  result <- rri_task_mean(datamat, n = 3)

  expect_equal(nrow(result), 2)  # 2 conditions
  expect_equal(result[1, 1], 2)  # Mean of 1, 2, 3
  expect_equal(result[2, 1], 5)  # Mean of 4, 5, 6
  expect_equal(result[1, 2], 8)  # Mean of 7, 8, 9
  expect_equal(result[2, 2], 11) # Mean of 10, 11, 12
})
