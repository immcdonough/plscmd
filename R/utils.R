#' Normalize Euclidean Distance of Vectors to Unit 1
#'
#' Normalize vectors in a matrix to unit Euclidean length.
#'
#' @param origin A numeric matrix where each column (or row) is a vector to normalize
#' @param dim Integer specifying the direction of vectors. 1 = column vectors (default),
#'   2 = row vectors
#'
#' @return A matrix with vectors normalized to unit length
#'
#' @examples
#' mat <- matrix(c(3, 4, 1, 2), nrow = 2)
#' normalize(mat)
#'
#' @export
normalize <- function(origin, dim = 1) {
  if (!is.matrix(origin)) {
    origin <- as.matrix(origin)
  }

  if (dim == 1) {
    # Column vectors - compute norm for each column
    normal_base <- sqrt(colSums(origin^2))
    normal_base <- matrix(rep(normal_base, each = nrow(origin)),
                          nrow = nrow(origin), ncol = ncol(origin))
  } else if (dim == 2) {
    # Row vectors - compute norm for each row
    normal_base <- sqrt(rowSums(origin^2))
    normal_base <- matrix(rep(normal_base, times = ncol(origin)),
                          nrow = nrow(origin), ncol = ncol(origin))
  } else {
    stop("dim must be 1 or 2")
  }

  # Handle zero vectors
  zero_items <- which(normal_base == 0)
  normal_base[zero_items] <- 1

  normal <- origin / normal_base
  normal[zero_items] <- 0

  return(normal)
}


#' Calculate Percentile
#'
#' The kth percentile Pk is that value of X which corresponds to a
#' cumulative frequency of Nk/100.
#'
#' @param X A numeric vector or matrix. If matrix, percentiles are computed
#'   for each column.
#' @param Nk Percentile(s) to compute (0-100). Can be a scalar or vector.
#'
#' @return Percentile value(s). If X is a matrix and Nk is a scalar, returns
#'   a vector of percentiles for each column.
#'
#' @examples
#' x <- 1:100
#' percentile(x, 50)  # Should return ~50.5
#' percentile(x, c(25, 50, 75))
#'
#' @export
percentile <- function(X, Nk) {
  X <- drop(X)  # Equivalent to squeeze in MATLAB

  if (is.matrix(Nk) || length(dim(Nk)) > 1) {
    stop("Nk must either be a scalar or a vector")
  }

  # If X is a 1D row vector, transpose it
  if (is.vector(X)) {
    X <- matrix(X, ncol = 1)
  } else if (is.matrix(X) && nrow(X) == 1) {
    X <- t(X)
  }

  r <- nrow(X)

  # Create x coordinates for interpolation
  x <- c(0, (seq(0.5, r - 0.5, by = 1) / r) * 100, 100)

  # For each column, create y values and interpolate
  result <- matrix(NA, nrow = length(Nk), ncol = ncol(X))

  for (j in seq_len(ncol(X))) {
    col_data <- X[, j]
    y <- c(min(col_data), sort(col_data), max(col_data))
    result[, j] <- approx(x, y, Nk, rule = 2)$y
  }

  if (ncol(X) == 1) {
    return(drop(result))
  }
  return(result)
}


#' Cumulative Gaussian Distribution Function
#'
#' Computes the cumulative distribution function (CDF) of the Gaussian
#' (normal) distribution. This gives the probability that a variate will
#' assume a value <= x.
#'
#' @param x Numeric value(s) at which to evaluate the CDF
#' @param mu Mean of the distribution (default = 0)
#' @param sigma Standard deviation of the distribution (default = 1)
#'
#' @return Probability value(s) between 0 and 1
#'
#' @details
#' The CDF is computed as: D(x) = (1 + erf((x-mu)/(sigma*sqrt(2)))) / 2
#'
#' @examples
#' cumulative_gaussian(0)  # Should return 0.5
#' cumulative_gaussian(1.96)  # Should return ~0.975
#'
#' @export
cumulative_gaussian <- function(x, mu = 0, sigma = 1) {
  # Expand mu and sigma if they are scalars
  if (length(mu) == 1) {
    mu <- rep(mu, length(x))
  }
  if (length(sigma) == 1) {
    sigma <- rep(sigma, length(x))
  }

  D <- rep(0, length(x))

  # Handle bad data (sigma <= 0)
  bad_idx <- which(sigma <= 0)
  if (length(bad_idx) > 0) {
    x[bad_idx] <- NA
  }

  # Compute CDF using error function
  # D = (1 + erf((x - mu) / (sigma * sqrt(2)))) / 2
  # In R, we can use pnorm directly
  good_idx <- seq_along(x)
  D[good_idx] <- pnorm(x[good_idx], mean = mu[good_idx], sd = sigma[good_idx])

  # Clamp values to [0, 1]
  D[D < 0] <- 0
  D[D > 1] <- 1

  return(D)
}


#' Inverse Cumulative Gaussian Distribution Function
#'
#' Computes the inverse of the cumulative distribution function (quantile
#' function) of the Gaussian (normal) distribution.
#'
#' @param D Probability value(s) between 0 and 1
#' @param mu Mean of the distribution (default = 0)
#' @param sigma Standard deviation of the distribution (default = 1)
#'
#' @return Quantile value(s)
#'
#' @details
#' Given D = (1 + erf((x-mu)/(sigma*sqrt(2))))/2, solves for x:
#' x = sigma*sqrt(2)*erfinv(2*D-1) + mu
#'
#' @examples
#' cumulative_gaussian_inv(0.5)  # Should return 0
#' cumulative_gaussian_inv(0.975)  # Should return ~1.96
#'
#' @export
cumulative_gaussian_inv <- function(D, mu = 0, sigma = 1) {
  # Expand mu and sigma if they are scalars
  if (length(mu) == 1) {
    mu <- rep(mu, length(D))
  }
  if (length(sigma) == 1) {
    sigma <- rep(sigma, length(D))
  }

  x <- rep(0, length(D))

  # Handle bad data
  bad_idx <- which(sigma <= 0 | D < 0 | D > 1)
  if (length(bad_idx) > 0) {
    x[bad_idx] <- NA
  }

  # Compute inverse CDF
  good_idx <- which(sigma > 0 & D > 0 & D < 1)
  if (length(good_idx) > 0) {
    x[good_idx] <- qnorm(D[good_idx], mean = mu[good_idx], sd = sigma[good_idx])
  }

  # Handle boundary cases
  x[D == 0] <- -Inf
  x[D == 1] <- Inf

  return(x)
}


#' Cross-Correlation of Two Matrices
#'
#' Computes cross-correlation or related measures between two matrices.
#'
#' @param design Design matrix
#' @param datamat Data matrix (must have same number of rows as design)
#' @param cormode Correlation mode:
#'   \itemize{
#'     \item 0 = Pearson correlation (default)
#'     \item 2 = Covariance
#'     \item 4 = Cosine angle
#'     \item 6 = Dot product
#'   }
#' @param robust_method Robust correlation method (only applies when cormode = 0):
#'   \itemize{
#'     \item "none" or "pearson" = Standard Pearson correlation (default)
#'     \item "spearman" = Spearman rank correlation
#'     \item "winsorized" = Winsorized correlation
#'     \item "biweight" = Biweight midcorrelation
#'     \item "percentage_bend" = Percentage bend correlation
#'   }
#' @param trim Trim proportion for winsorized correlation (default: 0.1)
#' @param beta Bend constant for percentage bend correlation (default: 0.2)
#'
#' @return Cross-correlation matrix
#'
#' @export
rri_xcor <- function(design, datamat, cormode = 0, robust_method = "none",
                     trim = 0.1, beta = 0.2) {
  if (!is.matrix(design)) design <- as.matrix(design)
  if (!is.matrix(datamat)) datamat <- as.matrix(datamat)

  r <- nrow(datamat)
  dr <- nrow(design)

  if (r != dr) {
    stop("Error in rri_xcor: input matrices must have same number of rows")
  }

  # Normalize robust_method
  robust_method <- tolower(robust_method)
  if (robust_method == "none") robust_method <- "pearson"

  # Use robust correlation for cormode = 0 with non-Pearson method
  if (cormode == 0 && robust_method != "pearson") {
    return(rri_xcor_robust(design, datamat, cormode = 0,
                           robust_method = robust_method,
                           trim = trim, beta = beta))
  }

  if (cormode == 0) {
    # Pearson correlation
    avg <- colMeans(datamat)
    stdev <- apply(datamat, 2, sd)

    # Handle zero standard deviation
    checknan <- which(stdev == 0)
    if (length(checknan) > 0) {
      datamat[, checknan] <- 0
      avg[checknan] <- 0
      stdev[checknan] <- 1
    }

    # Standardize data
    datamat <- sweep(datamat, 2, avg, "-")
    datamat <- sweep(datamat, 2, stdev, "/")

    davg <- colMeans(design)
    dstdev <- apply(design, 2, sd)

    checknan <- which(dstdev == 0)
    if (length(checknan) > 0) {
      design[, checknan] <- 0
      davg[checknan] <- 0
      dstdev[checknan] <- 1
    }

    design <- sweep(design, 2, davg, "-")
    design <- sweep(design, 2, dstdev, "/")

    xprod <- t(design) %*% datamat
    outmat <- xprod / (r - 1)

  } else if (cormode == 2) {
    # Covariance
    avg <- colMeans(datamat)
    davg <- colMeans(design)

    datamat <- sweep(datamat, 2, avg, "-")
    design <- sweep(design, 2, davg, "-")

    xprod <- t(design) %*% datamat
    outmat <- xprod / (r - 1)

  } else if (cormode == 4) {
    # Cosine angle
    stdev <- apply(datamat, 2, sd)
    checknan <- which(stdev == 0)
    if (length(checknan) > 0) {
      datamat[, checknan] <- 0
      stdev[checknan] <- 1
    }

    dstdev <- apply(design, 2, sd)
    checknan <- which(dstdev == 0)
    if (length(checknan) > 0) {
      design[, checknan] <- 0
      dstdev[checknan] <- 1
    }

    datamat <- sweep(datamat, 2, stdev, "/")
    design <- sweep(design, 2, dstdev, "/")

    outmat <- t(design) %*% datamat

  } else if (cormode == 6) {
    # Dot product
    outmat <- t(design) %*% datamat

  } else {
    stop("cormode must be 0, 2, 4, or 6")
  }

  return(outmat)
}


#' Calculate Task Means
#'
#' Returns a matrix of task means for an array with data for each task
#' stacked on top of one another.
#'
#' @param inmat Input matrix with tasks stacked vertically
#' @param n Number of subjects per condition
#'
#' @return Matrix of task means (one row per condition)
#'
#' @export
rri_task_mean <- function(inmat, n) {
  if (!is.matrix(inmat)) inmat <- as.matrix(inmat)

  m1 <- nrow(inmat)
  m <- ncol(inmat)
  k <- m1 / n  # Number of conditions

  if (k != floor(k)) {
    stop("Number of rows is not evenly divisible by n")
  }

  meanmat <- matrix(0, nrow = k, ncol = m)

  for (i in seq_len(k)) {
    start_row <- 1 + (n * (i - 1))
    end_row <- n * i
    temp <- inmat[start_row:end_row, , drop = FALSE]
    meanmat[i, ] <- colMeans(temp)
  }

  return(meanmat)
}


#' Calculate Task Means for Split-Subject-Blocks
#'
#' Version of rri_task_mean for variable sample sizes per condition.
#'
#' @param inmat Input matrix with tasks stacked vertically
#' @param n Vector of number of subjects per condition
#'
#' @return Matrix of task means (one row per condition)
#'
#' @export
ssb_rri_task_mean <- function(inmat, n) {
  if (!is.matrix(inmat)) inmat <- as.matrix(inmat)

  k <- length(n)
  m <- ncol(inmat)

  meanmat <- matrix(0, nrow = k, ncol = m)

  accum <- 0
  for (i in seq_len(k)) {
    temp <- inmat[(accum + 1):(accum + n[i]), , drop = FALSE]
    meanmat[i, ] <- colMeans(temp)
    accum <- accum + n[i]
  }

  return(meanmat)
}


#' Create Correlation Maps
#'
#' Creates image-wide correlation map for k conditions with behavioral data.
#'
#' @param behav Behavioral data matrix
#' @param datamat Data matrix
#' @param n Number of subjects
#' @param k Number of conditions
#' @param cormode Correlation mode (see rri_xcor)
#' @param robust_method Robust correlation method (see rri_xcor)
#' @param trim Trim proportion for winsorized correlation
#' @param beta Bend constant for percentage bend correlation
#'
#' @return Stacked correlation maps
#'
#' @export
rri_corr_maps <- function(behav, datamat, n, k, cormode = 0,
                          robust_method = "none", trim = 0.1, beta = 0.2) {
  if (!is.matrix(behav)) behav <- as.matrix(behav)
  if (!is.matrix(datamat)) datamat <- as.matrix(datamat)

  maps <- NULL

  for (i in seq_len(k)) {
    start_row <- 1 + (n * (i - 1))
    end_row <- n * i

    behav_subset <- behav[start_row:end_row, , drop = FALSE]
    data_subset <- datamat[start_row:end_row, , drop = FALSE]

    temp <- rri_xcor(behav_subset, data_subset, cormode,
                     robust_method = robust_method, trim = trim, beta = beta)
    maps <- rbind(maps, temp)
  }

  return(maps)
}


#' Stack Data Matrices
#'
#' Stacks data matrices from multiple groups into a single matrix.
#'
#' @param datamat_lst List of data matrices, one per group
#' @param single_cond_lst Optional single condition list (for special cases)
#' @param verbose Logical, whether to print progress
#'
#' @return Stacked data matrix
#'
#' @export
stacking_datamat <- function(datamat_lst, single_cond_lst = NULL, verbose = TRUE) {
  num_groups <- length(datamat_lst)
  stacked_datamat <- NULL

  if (verbose) {
    cat("Stacking datamat from group:")
  }

  for (g in seq_len(num_groups)) {
    if (verbose) {
      cat(" ", g)
    }

    if (is.null(single_cond_lst)) {
      datamat <- datamat_lst[[g]]
      stacked_datamat <- rbind(stacked_datamat, datamat)
    } else if (g == 1) {
      datamat <- single_cond_lst[[1]]
      stacked_datamat <- rbind(stacked_datamat, datamat)
    }
  }

  if (verbose) {
    cat("\n")
  }

  return(stacked_datamat)
}


#' Procrustes Rotation for Bootstrap
#'
#' Rotates a matrix to align with a reference matrix using Procrustes analysis.
#' Used in bootstrap to align bootstrapped loadings with original loadings.
#'
#' @param ref Reference matrix
#' @param target Target matrix to rotate
#'
#' @return Rotation matrix
#'
#' @keywords internal
rri_bootprocrust <- function(ref, target) {
  # Compute rotation matrix using SVD
  # target_rotated = target %*% rotatemat should be aligned with ref

  cross <- t(ref) %*% target
  svd_result <- svd(cross)

  rotatemat <- svd_result$v %*% t(svd_result$u)

  return(rotatemat)
}


#' Create Mask for Behavior Scan Selection
#'
#' Creates a mask for selecting specific conditions in multiblock PLS.
#'
#' @param num_subj_lst Vector of number of subjects per group
#' @param num_cond Number of conditions
#' @param bscan Vector of conditions to include (default: all)
#'
#' @return Vector of row indices
#'
#' @keywords internal
mask4bscan <- function(num_subj_lst, num_cond, bscan = NULL) {
  if (is.null(bscan)) {
    bscan <- seq_len(num_cond)
  }

  mask <- NULL

  if (!is.list(num_subj_lst)) {
    start_subj <- 0

    for (g in seq_along(num_subj_lst)) {
      n <- num_subj_lst[g]
      mat <- matrix(0, nrow = n, ncol = num_cond)
      mat[, bscan] <- 1
      mask <- c(mask, start_subj + which(as.vector(mat) == 1))
      start_subj <- start_subj + length(mat)
    }
  } else {
    for (g in seq_along(num_subj_lst)) {
      n <- num_subj_lst[[g]]

      for (i in seq_along(n)) {
        if (i %in% bscan) {
          mask <- c(mask, rep(1, n[i]))
        } else {
          mask <- c(mask, rep(0, n[i]))
        }
      }
    }

    mask <- which(mask == 1)
  }

  return(mask)
}
