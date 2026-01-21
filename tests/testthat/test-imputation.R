# Tests for missing data handling and imputation functions

test_that("pls_check_missing detects missing data correctly", {
  # Create test data with known missing pattern
  set.seed(42)
  data1 <- matrix(rnorm(300), 30, 10)  # 30 rows, 10 cols
  data1[1:3, 1] <- NA  # 10% missing in column 1
  data1[10, 5] <- NA   # Single missing value

  report <- pls_check_missing(list(data1), verbose = FALSE)

  expect_true(report$has_missing)
  expect_equal(report$total_missing, 4)
  expect_equal(report$total_cells, 300)
  expect_equal(report$overall_pct, 4/300)

  # Check variable-level stats
  expect_equal(report$by_variable$n_missing[1], 3)
  expect_equal(report$by_variable$n_missing[5], 1)
  expect_equal(sum(report$by_variable$n_missing), 4)
})

test_that("pls_check_missing handles no missing data", {
  data1 <- matrix(rnorm(100), 10, 10)

  report <- pls_check_missing(list(data1), verbose = FALSE)

  expect_false(report$has_missing)
  expect_equal(report$total_missing, 0)
})

test_that("pls_check_missing handles multiple groups", {
  set.seed(42)
  data1 <- matrix(rnorm(150), 15, 10)
  data2 <- matrix(rnorm(150), 15, 10)
  data1[1, 1] <- NA
  data2[1:3, 2] <- NA

  report <- pls_check_missing(list(data1, data2), verbose = FALSE)

  expect_true(report$has_missing)
  expect_equal(report$total_missing, 4)
  expect_equal(length(report$by_group), 2)
  expect_equal(report$by_group[[1]]$n_missing, 1)
  expect_equal(report$by_group[[2]]$n_missing, 3)
})

test_that("pls_check_missing generates warnings at thresholds", {
  # Create data with high missingness
  data1 <- matrix(rnorm(100), 10, 10)
  data1[1:6, 1] <- NA  # 60% missing in column 1

  # Should produce warnings about critical threshold
  expect_warning(
    report <- pls_check_missing(list(data1), critical_threshold = 0.5, verbose = TRUE)
  )
  # Verify warnings were recorded

  expect_true(length(report$warnings) > 0)
  expect_true(any(grepl("critical|>50%", report$warnings)))
})

test_that("pls_check_missing handles behavioral data", {
  data1 <- matrix(rnorm(100), 10, 10)
  behav <- matrix(rnorm(30), 10, 3)
  behav[1, 1] <- NA

  report <- pls_check_missing(list(data1), behavdata = behav, verbose = FALSE)

  expect_false(is.null(report$behavdata_stats))
  expect_equal(report$behavdata_stats$n_missing, 1)
})

test_that("impute_mean replaces NA with column means", {
  mat <- matrix(c(1, 2, NA, 4, 5, 6), nrow = 3, ncol = 2)
  # Column 1: 1, 2, NA -> mean = 1.5
  # Column 2: 4, 5, 6 -> no NA

  # impute_mean is an internal function, call it directly when sourced
  # or via the namespace when installed
  if (exists("impute_mean", mode = "function")) {
    result <- impute_mean(mat)
  } else {
    result <- plscmd:::impute_mean(mat)
  }

  expect_equal(result[3, 1], 1.5)
  expect_equal(result[, 2], c(4, 5, 6))
  expect_false(any(is.na(result)))
})

test_that("pls_impute with mean method removes all NAs", {
  set.seed(42)
  data1 <- matrix(rnorm(300), 30, 10)
  data1[sample(300, 15)] <- NA

  imputed <- pls_impute(list(data1), method = "mean", m = 1, verbose = FALSE)

  expect_false(any(is.na(imputed$datamat_lst[[1]])))
  expect_equal(dim(imputed$datamat_lst[[1]]), c(30, 10))
  expect_true(imputed$imputed)
})

test_that("pls_impute returns unchanged data when no missing", {
  data1 <- matrix(rnorm(100), 10, 10)

  imputed <- pls_impute(list(data1), method = "mean", m = 1, verbose = FALSE)

  expect_equal(imputed$datamat_lst[[1]], data1)
  expect_false(imputed$imputed)
})

test_that("pls_impute multiple returns m datasets", {
  set.seed(42)
  data1 <- matrix(rnorm(300), 30, 10)
  data1[sample(300, 15)] <- NA

  mi_data <- pls_impute(list(data1), method = "mean", m = 5, verbose = FALSE)

  expect_s3_class(mi_data, "pls_mi_data")
  expect_equal(mi_data$m, 5)
  expect_equal(length(mi_data$imputations), 5)

  # All imputations should have no missing data
  for (i in 1:5) {
    expect_false(any(is.na(mi_data$imputations[[i]]$datamat_lst[[1]])))
  }
})

test_that("pool_estimates correctly applies Rubin's rules", {
  # Create simple test case
  # 3 parameters, 2 imputations
  estimates <- matrix(c(1, 2, 3, 1.2, 2.1, 3.3), nrow = 3, ncol = 2)
  variances <- matrix(c(0.1, 0.2, 0.3, 0.12, 0.18, 0.32), nrow = 3, ncol = 2)

  result <- pool_estimates(estimates, variances)

  # Pooled point estimate should be mean
  expect_equal(result$pooled[1], mean(c(1, 1.2)))
  expect_equal(result$pooled[2], mean(c(2, 2.1)))
  expect_equal(result$pooled[3], mean(c(3, 3.3)))

  # Within variance should be mean of variances
  expect_equal(result$within_var[1], mean(c(0.1, 0.12)))

  # Check structure
  expect_equal(result$m, 2)
  expect_true(!is.null(result$between_var))
  expect_true(!is.null(result$total_var))
  expect_true(!is.null(result$fmi))
})

test_that("pool_estimates handles arrays", {
  # 2x2 matrix, 3 imputations
  estimates <- array(1:12, dim = c(2, 2, 3))

  result <- pool_estimates(estimates)

  expect_equal(dim(result$pooled), c(2, 2))
  expect_equal(result$m, 3)
})

test_that("pool_pvalues returns median", {
  pvals <- matrix(c(0.01, 0.02, 0.03,
                    0.05, 0.06, 0.04), nrow = 2, byrow = TRUE)

  result <- pool_pvalues(pvals)

  expect_equal(result[1], median(c(0.01, 0.02, 0.03)))
  expect_equal(result[2], median(c(0.05, 0.06, 0.04)))
})

test_that("pls_analysis errors on NA without impute_method", {
  set.seed(42)
  data1 <- matrix(rnorm(300), 30, 10)
  data1[1, 1] <- NA

  expect_error(
    pls_analysis(list(data1), c(10), 3, option = list(verbose = FALSE)),
    "Missing data"
  )
})

test_that("pls_analysis with impute_method='mean' handles missing data", {
  set.seed(42)
  data1 <- matrix(rnorm(300), 30, 10)
  data1[1, 1] <- NA

  # Should not error with impute_method set
  result <- pls_analysis(
    list(data1),
    c(10),
    3,
    option = list(
      impute_method = "mean",
      verbose = FALSE,
      check_missing = FALSE
    )
  )

  expect_s3_class(result, "pls_result")
  expect_true(!is.null(result$imputation_info))
  expect_equal(result$imputation_info$method, "mean")
})

# Skip RF tests if missRanger not available
test_that("pls_impute with RF method works", {
  skip_if_not_installed("missRanger")

  set.seed(42)
  data1 <- matrix(rnorm(300), 30, 10)
  data1[sample(300, 15)] <- NA

  imputed <- pls_impute(
    list(data1),
    method = "rf",
    m = 1,
    maxiter = 3,
    seed = 42,
    verbose = FALSE
  )

  expect_false(any(is.na(imputed$datamat_lst[[1]])))
  expect_true(imputed$imputed)
})

test_that("pls_analysis_mi runs with missing data", {
  skip_if_not_installed("missRanger")

  set.seed(42)
  data1 <- matrix(rnorm(300), 30, 10)
  data1[sample(300, 10)] <- NA  # ~3% missing

  # Use fewer permutations/bootstraps for speed
  result <- suppressMessages(pls_analysis_mi(
    datamat_lst = list(data1),
    num_subj_lst = c(10),
    num_cond = 3,
    option = list(
      method = 1,
      num_perm = 10,
      num_boot = 0,
      verbose = FALSE
    ),
    mi_option = list(
      method = "rf",
      m = 3,
      maxiter = 3,
      seed = 42
    )
  ))

  expect_s3_class(result, "pls_mi_result")
  expect_equal(result$m, 3)
  expect_true(!is.null(result$perm_result))
  expect_equal(length(result$perm_result$sprob), length(result$s))
})

test_that("print methods work for imputation objects", {
  # Test print.pls_missing_report
  data1 <- matrix(rnorm(100), 10, 10)
  data1[1, 1] <- NA

  report <- pls_check_missing(list(data1), verbose = FALSE)
  expect_output(print(report), "Missing Data Report")

  # Test print.pls_mi_data
  mi_data <- pls_impute(list(data1), method = "mean", m = 3, verbose = FALSE)
  expect_output(print(mi_data), "Multiple Imputation Data")
})
