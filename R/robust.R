#' Robust Correlation Functions
#'
#' Functions for computing robust correlations that downweight influential points.
#'
#' @name robust
#' @keywords internal
NULL

#' Compute Robust Correlation
#'
#' Computes correlation using various robust methods that are resistant to outliers.
#'
#' @param x Numeric vector or matrix
#' @param y Numeric vector or matrix (if NULL, computes correlation matrix of x)
#' @param method Robust method to use:
#'   \itemize{
#'     \item "pearson" = Standard Pearson correlation (default)
#'     \item "spearman" = Spearman rank correlation
#'     \item "winsorized" = Winsorized correlation (trims extreme values)
#'     \item "biweight" = Biweight midcorrelation (downweights outliers)
#'     \item "percentage_bend" = Percentage bend correlation
#'   }
#' @param trim Trim proportion for winsorized correlation (default: 0.1)
#' @param beta Bend constant for percentage bend correlation (default: 0.2)
#'
#' @return Correlation coefficient(s)
#'
#' @details
#' \strong{Spearman:} Uses ranks instead of raw values, resistant to outliers
#' and monotonic transformations.
#'
#' \strong{Winsorized:} Replaces extreme values with less extreme values
#' (at the trim percentile) before computing Pearson correlation.
#'
#' \strong{Biweight midcorrelation:} Uses a weighted correlation where weights
#' decrease for observations far from the median. Very robust to outliers.
#'
#' \strong{Percentage bend:} Similar to Winsorized but uses a different
#' bending function. Good balance of robustness and efficiency.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' x <- c(1, 2, 3, 4, 100)  # Contains outlier
#' y <- c(1, 2, 3, 4, 5)
#' cor(x, y)  # Pearson: affected by outlier
#' robust_cor(x, y, method = "spearman")
#' robust_cor(x, y, method = "winsorized")
#' robust_cor(x, y, method = "biweight")
#' }
robust_cor <- function(x, y = NULL, method = "pearson", trim = 0.1, beta = 0.2) {
  method <- tolower(method)

  if (method == "pearson") {
    if (is.null(y)) {
      return(stats::cor(x))
    } else {
      return(stats::cor(x, y))
    }
  } else if (method == "spearman") {
    if (is.null(y)) {
      return(stats::cor(x, method = "spearman"))
    } else {
      return(stats::cor(x, y, method = "spearman"))
    }
  } else if (method == "winsorized") {
    return(winsorized_cor(x, y, trim = trim))
  } else if (method == "biweight") {
    return(biweight_midcor(x, y))
  } else if (method == "percentage_bend") {
    return(percentage_bend_cor(x, y, beta = beta))
  } else {
    stop("Unknown robust method: ", method,
         ". Use 'pearson', 'spearman', 'winsorized', 'biweight', or 'percentage_bend'")
  }
}


#' Winsorized Correlation
#'
#' Computes correlation after winsorizing (trimming) extreme values.
#'
#' @param x Numeric vector
#' @param y Numeric vector (optional)
#' @param trim Proportion to trim from each tail (default: 0.1 = 10%)
#'
#' @return Winsorized correlation coefficient
#'
#' @export
#'
#' @examples
#' \dontrun{
#' x <- c(1, 2, 3, 4, 100)
#' y <- c(1, 2, 3, 4, 5)
#' winsorized_cor(x, y)
#' }
winsorized_cor <- function(x, y = NULL, trim = 0.1) {
  winsorize <- function(v, trim) {
    v <- as.numeric(v)
    n <- length(v)
    lo <- stats::quantile(v, trim, na.rm = TRUE)
    hi <- stats::quantile(v, 1 - trim, na.rm = TRUE)
    v[v < lo] <- lo
    v[v > hi] <- hi
    return(v)
  }

  if (is.null(y)) {
    # Correlation matrix of x
    if (!is.matrix(x)) x <- as.matrix(x)
    x_wins <- apply(x, 2, winsorize, trim = trim)
    return(stats::cor(x_wins))
  } else {
    x_wins <- winsorize(x, trim)
    y_wins <- winsorize(y, trim)
    return(stats::cor(x_wins, y_wins))
  }
}


#' Biweight Midcorrelation
#'
#' Computes the biweight midcorrelation, a robust correlation measure that
#' downweights observations far from the median.
#'
#' @param x Numeric vector
#' @param y Numeric vector (optional)
#' @param c Tuning constant (default: 9 for 95% efficiency)
#'
#' @return Biweight midcorrelation coefficient
#'
#' @details
#' The biweight midcorrelation uses the biweight midvariance and midcovariance,
#' which are robust measures based on a weighting function that downweights
#' observations based on their distance from the median.
#'
#' @references
#' Wilcox, R.R. (2012). Introduction to Robust Estimation and Hypothesis Testing.
#' Academic Press.
#'
#' @export
biweight_midcor <- function(x, y = NULL, c = 9) {
  # Biweight midvariance
  biweight_midvar <- function(v, c = 9) {
    v <- as.numeric(v)
    v <- v[!is.na(v)]
    n <- length(v)
    med <- stats::median(v)
    mad_val <- stats::mad(v, constant = 1)

    if (mad_val == 0) {
      return(stats::var(v))
    }

    u <- (v - med) / (c * mad_val)
    a <- as.numeric(abs(u) < 1)

    num <- sum(a * (v - med)^2 * (1 - u^2)^4)
    denom <- abs(sum(a * (1 - u^2) * (1 - 5 * u^2)))

    if (denom == 0) {
      return(stats::var(v))
    }

    return(n * num / denom^2)
  }

  # Biweight midcovariance
  biweight_midcov <- function(x, y, c = 9) {
    x <- as.numeric(x)
    y <- as.numeric(y)

    # Remove NAs pairwise
    valid <- !is.na(x) & !is.na(y)
    x <- x[valid]
    y <- y[valid]
    n <- length(x)

    med_x <- stats::median(x)
    med_y <- stats::median(y)
    mad_x <- stats::mad(x, constant = 1)
    mad_y <- stats::mad(y, constant = 1)

    if (mad_x == 0 || mad_y == 0) {
      return(stats::cov(x, y))
    }

    u <- (x - med_x) / (c * mad_x)
    v <- (y - med_y) / (c * mad_y)

    a_x <- as.numeric(abs(u) < 1)
    a_y <- as.numeric(abs(v) < 1)

    num <- sum(a_x * a_y * (x - med_x) * (y - med_y) * (1 - u^2)^2 * (1 - v^2)^2)
    denom_x <- sum(a_x * (1 - u^2) * (1 - 5 * u^2))
    denom_y <- sum(a_y * (1 - v^2) * (1 - 5 * v^2))

    if (denom_x == 0 || denom_y == 0) {
      return(stats::cov(x, y))
    }

    return(n * num / (denom_x * denom_y))
  }

  if (is.null(y)) {
    # Correlation matrix of x
    if (!is.matrix(x)) x <- as.matrix(x)
    p <- ncol(x)
    cor_mat <- matrix(1, nrow = p, ncol = p)

    for (i in 1:(p-1)) {
      for (j in (i+1):p) {
        cov_ij <- biweight_midcov(x[, i], x[, j], c)
        var_i <- biweight_midvar(x[, i], c)
        var_j <- biweight_midvar(x[, j], c)

        if (var_i > 0 && var_j > 0) {
          r <- cov_ij / sqrt(var_i * var_j)
          r <- max(-1, min(1, r))  # Clamp to [-1, 1]
          cor_mat[i, j] <- r
          cor_mat[j, i] <- r
        }
      }
    }
    return(cor_mat)
  } else {
    cov_xy <- biweight_midcov(x, y, c)
    var_x <- biweight_midvar(x, c)
    var_y <- biweight_midvar(y, c)

    if (var_x <= 0 || var_y <= 0) {
      return(stats::cor(x, y))
    }

    r <- cov_xy / sqrt(var_x * var_y)
    return(max(-1, min(1, r)))
  }
}


#' Percentage Bend Correlation
#'
#' Computes the percentage bend correlation, a robust correlation measure.
#'
#' @param x Numeric vector
#' @param y Numeric vector (optional)
#' @param beta Bending constant (default: 0.2)
#'
#' @return Percentage bend correlation coefficient
#'
#' @details
#' The percentage bend correlation downweights a proportion (beta) of the
#' observations in each tail. It provides a good balance between robustness
#' and statistical efficiency.
#'
#' @references
#' Wilcox, R.R. (1994). The percentage bend correlation coefficient.
#' Psychometrika, 59, 601-616.
#'
#' @export
percentage_bend_cor <- function(x, y = NULL, beta = 0.2) {
  # Compute percentage bend for a single vector
  pbend_transform <- function(v, beta) {
    v <- as.numeric(v)
    v <- v[!is.na(v)]
    n <- length(v)

    med <- stats::median(v)
    sorted_abs_dev <- sort(abs(v - med))
    m <- floor((1 - beta) * n + 0.5)

    if (m < 1) m <- 1
    if (m > n) m <- n

    omega <- sorted_abs_dev[m]

    if (omega == 0) {
      return(list(transformed = v - med, omega = 1))
    }

    # Bend values outside omega
    psi <- (v - med) / omega
    psi[psi < -1] <- -1
    psi[psi > 1] <- 1

    return(list(transformed = psi, omega = omega))
  }

  if (is.null(y)) {
    # Correlation matrix of x
    if (!is.matrix(x)) x <- as.matrix(x)
    p <- ncol(x)
    cor_mat <- matrix(1, nrow = p, ncol = p)

    for (i in 1:(p-1)) {
      for (j in (i+1):p) {
        valid <- !is.na(x[, i]) & !is.na(x[, j])
        xi <- x[valid, i]
        xj <- x[valid, j]

        pb_i <- pbend_transform(xi, beta)
        pb_j <- pbend_transform(xj, beta)

        r <- stats::cor(pb_i$transformed, pb_j$transformed)
        cor_mat[i, j] <- r
        cor_mat[j, i] <- r
      }
    }
    return(cor_mat)
  } else {
    valid <- !is.na(x) & !is.na(y)
    x <- x[valid]
    y <- y[valid]

    pb_x <- pbend_transform(x, beta)
    pb_y <- pbend_transform(y, beta)

    return(stats::cor(pb_x$transformed, pb_y$transformed))
  }
}


#' Robust Cross-Correlation of Two Matrices
#'
#' Computes robust cross-correlation between two matrices using various
#' robust methods.
#'
#' @param design Design matrix
#' @param datamat Data matrix (must have same number of rows as design)
#' @param cormode Correlation mode (0 = correlation, 2 = covariance, etc.)
#' @param robust_method Robust method: "pearson", "spearman", "winsorized",
#'   "biweight", or "percentage_bend"
#' @param trim Trim proportion for winsorized correlation
#' @param beta Bend constant for percentage bend correlation
#'
#' @return Cross-correlation matrix
#'
#' @keywords internal
rri_xcor_robust <- function(design, datamat, cormode = 0,
                            robust_method = "pearson",
                            trim = 0.1, beta = 0.2) {
  if (!is.matrix(design)) design <- as.matrix(design)
  if (!is.matrix(datamat)) datamat <- as.matrix(datamat)

  r <- nrow(datamat)
  dr <- nrow(design)

  if (r != dr) {
    stop("Error in rri_xcor_robust: input matrices must have same number of rows")
  }

  robust_method <- tolower(robust_method)

  # For non-correlation modes (covariance, cosine, dot product), use standard method
  if (cormode != 0 || robust_method == "pearson") {
    # Use standard rri_xcor logic
    if (cormode == 0) {
      # Pearson correlation
      avg <- colMeans(datamat)
      stdev <- apply(datamat, 2, sd)

      checknan <- which(stdev == 0)
      if (length(checknan) > 0) {
        datamat[, checknan] <- 0
        avg[checknan] <- 0
        stdev[checknan] <- 1
      }

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
      avg <- colMeans(datamat)
      davg <- colMeans(design)

      datamat <- sweep(datamat, 2, avg, "-")
      design <- sweep(design, 2, davg, "-")

      xprod <- t(design) %*% datamat
      outmat <- xprod / (r - 1)

    } else if (cormode == 4) {
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
      outmat <- t(design) %*% datamat

    } else {
      stop("cormode must be 0, 2, 4, or 6")
    }

    return(outmat)
  }

  # Robust correlation methods (only for cormode = 0)
  n_design <- ncol(design)
  n_data <- ncol(datamat)

  outmat <- matrix(0, nrow = n_design, ncol = n_data)

  if (robust_method == "spearman") {
    # Spearman: rank-transform then compute Pearson
    design_ranked <- apply(design, 2, rank, na.last = "keep")
    datamat_ranked <- apply(datamat, 2, rank, na.last = "keep")

    # Standardize
    design_ranked <- scale(design_ranked)
    datamat_ranked <- scale(datamat_ranked)

    # Handle zero variance columns
    design_ranked[is.nan(design_ranked)] <- 0
    datamat_ranked[is.nan(datamat_ranked)] <- 0

    outmat <- t(design_ranked) %*% datamat_ranked / (r - 1)

  } else if (robust_method == "winsorized") {
    # Winsorize each column
    winsorize <- function(v, trim) {
      lo <- stats::quantile(v, trim, na.rm = TRUE)
      hi <- stats::quantile(v, 1 - trim, na.rm = TRUE)
      v[v < lo] <- lo
      v[v > hi] <- hi
      return(v)
    }

    design_wins <- apply(design, 2, winsorize, trim = trim)
    datamat_wins <- apply(datamat, 2, winsorize, trim = trim)

    # Standardize
    design_wins <- scale(design_wins)
    datamat_wins <- scale(datamat_wins)

    # Handle zero variance columns
    design_wins[is.nan(design_wins)] <- 0
    datamat_wins[is.nan(datamat_wins)] <- 0

    outmat <- t(design_wins) %*% datamat_wins / (r - 1)

  } else if (robust_method == "biweight") {
    # Compute column by column (slower but robust)
    for (i in seq_len(n_design)) {
      for (j in seq_len(n_data)) {
        outmat[i, j] <- biweight_midcor(design[, i], datamat[, j])
      }
    }

  } else if (robust_method == "percentage_bend") {
    # Compute column by column
    for (i in seq_len(n_design)) {
      for (j in seq_len(n_data)) {
        outmat[i, j] <- percentage_bend_cor(design[, i], datamat[, j], beta = beta)
      }
    }

  } else {
    stop("Unknown robust method: ", robust_method)
  }

  return(outmat)
}
