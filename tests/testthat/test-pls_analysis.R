# Test main PLS analysis function

test_that("pls_analysis runs task PLS correctly", {
  set.seed(42)

  # Create simulated data: 2 groups, 3 conditions, 5 subjects each
  n_subj <- 5
  n_cond <- 3
  n_vars <- 50

  # Group 1: Add condition effect
  datamat1 <- matrix(rnorm(n_subj * n_cond * n_vars), nrow = n_subj * n_cond, ncol = n_vars)
  for (c in 1:n_cond) {
    rows <- ((c - 1) * n_subj + 1):(c * n_subj)
    datamat1[rows, 1:10] <- datamat1[rows, 1:10] + c  # Add condition effect
  }

  # Group 2: Similar structure
  datamat2 <- matrix(rnorm(n_subj * n_cond * n_vars), nrow = n_subj * n_cond, ncol = n_vars)
  for (c in 1:n_cond) {
    rows <- ((c - 1) * n_subj + 1):(c * n_subj)
    datamat2[rows, 1:10] <- datamat2[rows, 1:10] + c
  }

  # Run analysis
  result <- pls_analysis(
    datamat_lst = list(datamat1, datamat2),
    num_subj_lst = c(n_subj, n_subj),
    num_cond = n_cond,
    option = list(
      method = 1,
      verbose = FALSE
    )
  )

  # Check structure
  expect_s3_class(result, "pls_result")
  expect_equal(result$method, 1)
  expect_equal(length(result$s), min(n_cond * 2, n_vars))  # Number of LVs
  expect_true(all(result$s >= 0))  # Singular values should be non-negative
  expect_true(is.matrix(result$u))
  expect_true(is.matrix(result$v))
})


test_that("pls_analysis runs with permutation test", {
  set.seed(42)

  n_subj <- 4
  n_cond <- 2
  n_vars <- 20

  datamat1 <- matrix(rnorm(n_subj * n_cond * n_vars), nrow = n_subj * n_cond, ncol = n_vars)
  datamat2 <- matrix(rnorm(n_subj * n_cond * n_vars), nrow = n_subj * n_cond, ncol = n_vars)

  result <- pls_analysis(
    datamat_lst = list(datamat1, datamat2),
    num_subj_lst = c(n_subj, n_subj),
    num_cond = n_cond,
    option = list(
      method = 1,
      num_perm = 10,  # Small number for testing
      verbose = FALSE
    )
  )

  expect_true(!is.null(result$perm_result))
  expect_equal(result$perm_result$num_perm, 10)
  expect_true(all(result$perm_result$sprob >= 0))
  expect_true(all(result$perm_result$sprob <= 1))
})


test_that("pls_analysis runs with bootstrap", {
  set.seed(42)

  n_subj <- 4
  n_cond <- 2
  n_vars <- 20

  datamat1 <- matrix(rnorm(n_subj * n_cond * n_vars), nrow = n_subj * n_cond, ncol = n_vars)
  datamat2 <- matrix(rnorm(n_subj * n_cond * n_vars), nrow = n_subj * n_cond, ncol = n_vars)

  result <- pls_analysis(
    datamat_lst = list(datamat1, datamat2),
    num_subj_lst = c(n_subj, n_subj),
    num_cond = n_cond,
    option = list(
      method = 1,
      num_boot = 10,  # Small number for testing
      verbose = FALSE
    )
  )

  expect_true(!is.null(result$boot_result))
  expect_true(result$boot_result$num_boot >= 1)
  expect_true(is.matrix(result$boot_result$compare_u))
  expect_true(is.matrix(result$boot_result$u_se))
})


test_that("behavior PLS runs correctly", {
  set.seed(42)

  n_subj <- 5
  n_cond <- 2
  n_vars <- 30
  n_behav <- 2

  # Create brain data
  datamat1 <- matrix(rnorm(n_subj * n_cond * n_vars), nrow = n_subj * n_cond, ncol = n_vars)

  # Create behavior data correlated with brain
  behavdata <- matrix(rnorm(n_subj * n_cond * n_behav), nrow = n_subj * n_cond, ncol = n_behav)

  result <- pls_analysis(
    datamat_lst = list(datamat1),
    num_subj_lst = c(n_subj),
    num_cond = n_cond,
    option = list(
      method = 3,
      stacked_behavdata = behavdata,
      verbose = FALSE
    )
  )

  expect_equal(result$method, 3)
  expect_true(!is.null(result$lvcorrs))
  expect_true(!is.null(result$datamatcorrs_lst))
})


test_that("non-rotated PLS runs correctly", {
  set.seed(42)

  n_subj <- 4
  n_cond <- 3
  n_vars <- 20
  n_groups <- 2

  datamat1 <- matrix(rnorm(n_subj * n_cond * n_vars), nrow = n_subj * n_cond, ncol = n_vars)
  datamat2 <- matrix(rnorm(n_subj * n_cond * n_vars), nrow = n_subj * n_cond, ncol = n_vars)

  # Create design contrast (must have num_groups * num_cond rows)
  designdata <- matrix(c(
    1, -1,  # Group 1, Cond 1
    0,  1,  # Group 1, Cond 2
    -1, 0,  # Group 1, Cond 3
    1, -1,  # Group 2, Cond 1
    0,  1,  # Group 2, Cond 2
    -1, 0   # Group 2, Cond 3
  ), nrow = n_groups * n_cond, ncol = 2, byrow = TRUE)

  result <- pls_analysis(
    datamat_lst = list(datamat1, datamat2),
    num_subj_lst = c(n_subj, n_subj),
    num_cond = n_cond,
    option = list(
      method = 2,
      stacked_designdata = designdata,
      verbose = FALSE
    )
  )

  expect_equal(result$method, 2)
  expect_true(!is.null(result$lvintercorrs))
})
