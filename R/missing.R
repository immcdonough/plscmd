#' Mean Ignoring NA Values
#'
#' Calculate mean while ignoring NA values.
#'
#' @param x Numeric vector or matrix
#' @param dim Dimension to operate on (1 = columns, 2 = rows). If NULL,
#'   operates on entire vector/matrix.
#'
#' @return Mean value(s)
#'
#' @export
missmean <- function(x, dim = NULL) {
  if (is.null(dim)) {
    return(mean(x, na.rm = TRUE))
  }

  if (!is.matrix(x)) x <- as.matrix(x)

  if (dim == 1) {
    return(colMeans(x, na.rm = TRUE))
  } else if (dim == 2) {
    return(rowMeans(x, na.rm = TRUE))
  } else {
    stop("dim must be 1 or 2")
  }
}


#' Standard Deviation Ignoring NA Values
#'
#' Calculate standard deviation while ignoring NA values.
#'
#' @param x Numeric vector or matrix
#' @param dim Dimension to operate on (1 = columns, 2 = rows). If NULL,
#'   operates on entire vector/matrix.
#'
#' @return Standard deviation value(s)
#'
#' @export
missstd <- function(x, dim = NULL) {
  if (is.null(dim)) {
    return(sd(x, na.rm = TRUE))
  }

  if (!is.matrix(x)) x <- as.matrix(x)

  if (dim == 1) {
    return(apply(x, 2, sd, na.rm = TRUE))
  } else if (dim == 2) {
    return(apply(x, 1, sd, na.rm = TRUE))
  } else {
    stop("dim must be 1 or 2")
  }
}


#' Variance Ignoring NA Values
#'
#' Calculate variance while ignoring NA values.
#'
#' @param x Numeric vector or matrix
#' @param dim Dimension to operate on (1 = columns, 2 = rows). If NULL,
#'   operates on entire vector/matrix.
#'
#' @return Variance value(s)
#'
#' @export
missvar <- function(x, dim = NULL) {
  if (is.null(dim)) {
    return(var(as.vector(x), na.rm = TRUE))
  }

  if (!is.matrix(x)) x <- as.matrix(x)

  if (dim == 1) {
    return(apply(x, 2, var, na.rm = TRUE))
  } else if (dim == 2) {
    return(apply(x, 1, var, na.rm = TRUE))
  } else {
    stop("dim must be 1 or 2")
  }
}


#' Weighted Sum Ignoring NA Values
#'
#' Calculate weighted sum while ignoring NA values.
#'
#' @param x Numeric vector or matrix
#' @param w Optional weights
#' @param dim Dimension to operate on (1 = columns, 2 = rows). If NULL,
#'   operates on entire vector/matrix.
#'
#' @return Sum value(s)
#'
#' @export
misssum <- function(x, w = NULL, dim = NULL) {
  if (is.null(w)) {
    w <- matrix(1, nrow = nrow(as.matrix(x)), ncol = ncol(as.matrix(x)))
  }

  # Replace NA with 0 for summation, but track positions
  x_clean <- x
  x_clean[is.na(x)] <- 0
  w[is.na(x)] <- 0

  if (is.null(dim)) {
    return(sum(x_clean * w))
  }

  if (!is.matrix(x_clean)) x_clean <- as.matrix(x_clean)
  if (!is.matrix(w)) w <- as.matrix(w)

  if (dim == 1) {
    return(colSums(x_clean * w))
  } else if (dim == 2) {
    return(rowSums(x_clean * w))
  } else {
    stop("dim must be 1 or 2")
  }
}


#' SVD with Missing Data Handling
#'
#' Performs Singular Value Decomposition while handling missing (NA) values
#' through iterative imputation.
#'
#' @param x Input matrix that may contain NA values
#' @param ncomp Number of components to return (0 = all)
#' @param max_iter Maximum iterations for convergence
#' @param tol Convergence tolerance
#'
#' @return If nu and nv are not specified, returns a list with components:
#'   \itemize{
#'     \item u: Left singular vectors
#'     \item d: Singular values (as diagonal matrix)
#'     \item v: Right singular vectors
#'   }
#'
#' @export
misssvd <- function(x, ncomp = 0, max_iter = 100, tol = 1e-12) {
  if (!is.matrix(x)) x <- as.matrix(x)

  # Find missing values
  na_idx <- which(is.na(x))

  if (length(na_idx) == 0) {
    # No missing values, use standard SVD
    svd_result <- svd(x)
    if (ncomp > 0 && ncomp < length(svd_result$d)) {
      return(list(
        u = svd_result$u[, 1:ncomp, drop = FALSE],
        d = diag(svd_result$d[1:ncomp], nrow = ncomp),
        v = svd_result$v[, 1:ncomp, drop = FALSE]
      ))
    }
    return(list(
      u = svd_result$u,
      d = diag(svd_result$d, nrow = length(svd_result$d)),
      v = svd_result$v
    ))
  }

  # Initial imputation: replace NA with column means
  x_filled <- x
  col_means <- colMeans(x, na.rm = TRUE)
  for (j in seq_len(ncol(x))) {
    x_filled[is.na(x_filled[, j]), j] <- col_means[j]
  }

  # Handle case where column mean is still NA
  x_filled[is.na(x_filled)] <- 0

  # Iterative SVD
  old_x <- x_filled

  for (iter in seq_len(max_iter)) {
    svd_result <- svd(x_filled)

    # Reconstruct matrix from SVD
    if (ncomp > 0 && ncomp < length(svd_result$d)) {
      reconstructed <- svd_result$u[, 1:ncomp, drop = FALSE] %*%
        diag(svd_result$d[1:ncomp], nrow = ncomp) %*%
        t(svd_result$v[, 1:ncomp, drop = FALSE])
    } else {
      reconstructed <- svd_result$u %*% diag(svd_result$d) %*% t(svd_result$v)
    }

    # Replace missing values with reconstructed values
    x_filled[na_idx] <- reconstructed[na_idx]

    # Check convergence
    diff <- max(abs(x_filled - old_x))
    if (diff < tol) {
      break
    }

    old_x <- x_filled
  }

  # Final SVD
  svd_result <- svd(x_filled)

  if (ncomp > 0 && ncomp < length(svd_result$d)) {
    return(list(
      u = svd_result$u[, 1:ncomp, drop = FALSE],
      d = diag(svd_result$d[1:ncomp], nrow = ncomp),
      v = svd_result$v[, 1:ncomp, drop = FALSE]
    ))
  }

  return(list(
    u = svd_result$u,
    d = diag(svd_result$d, nrow = length(svd_result$d)),
    v = svd_result$v
  ))
}


#' Cross-Correlation with Missing Data Handling
#'
#' Computes cross-correlation between two matrices, handling missing values.
#'
#' @param design Design matrix
#' @param datamat Data matrix
#' @param cormode Correlation mode (see rri_xcor)
#'
#' @return Cross-correlation matrix
#'
#' @keywords internal
missnk_rri_xcor <- function(design, datamat, cormode = 0) {
  # This version handles NA values by using pairwise complete observations
  if (!is.matrix(design)) design <- as.matrix(design)
  if (!is.matrix(datamat)) datamat <- as.matrix(datamat)

  r <- nrow(datamat)
  dr <- nrow(design)

  if (r != dr) {
    stop("Input matrices must have same number of rows")
  }

  # For now, use listwise deletion
  complete_rows <- complete.cases(cbind(design, datamat))

  if (sum(complete_rows) < 3) {
    warning("Too few complete cases for correlation")
    return(matrix(NA, nrow = ncol(design), ncol = ncol(datamat)))
  }

  return(rri_xcor(design[complete_rows, , drop = FALSE],
                  datamat[complete_rows, , drop = FALSE],
                  cormode))
}


#' Check for Low Variability in Resampled Data
#'
#' Checks if resampled data has low variability (less than 50% unique values).
#'
#' @param resampled Resampled data vector
#' @param original Original data vector
#'
#' @return Logical indicating if variability is low
#'
#' @keywords internal
rri_islowvariability <- function(resampled, original) {
  # Count unique values in resampled relative to original
  n_unique_resampled <- length(unique(resampled[!is.na(resampled)]))
  n_unique_original <- length(unique(original[!is.na(original)]))

  # Low variability if less than 50% of original unique values
  return(n_unique_resampled < (n_unique_original * 0.5))
}
