#' PLS Regression Analysis
#'
#' Performs Partial Least Squares (PLS) regression, a method for relating
#' predictor variables to response variables by finding latent components
#' that maximize covariance.
#'
#' @name pls_regression
#' @keywords regression
NULL


#' PLS Regression
#'
#' Fit a PLS regression model using SIMPLS or NIPALS algorithm.
#'
#' @param X Predictor matrix (n x p) where n is number of observations and
#'   p is number of predictor variables
#' @param Y Response matrix (n x q) or vector. If vector, treated as single response.
#' @param ncomp Number of components to extract. If NULL (default), uses
#'   min(n-1, p, 10)
#' @param algorithm Algorithm to use: "simpls" (default) or "nipals"
#' @param scale Logical, scale X and Y to unit variance (default: TRUE)
#' @param center Logical, mean-center X and Y (default: TRUE)
#' @param max_iter Maximum iterations for NIPALS convergence (default: 500)
#' @param tol Convergence tolerance for NIPALS (default: 1e-6)
#' @param robust_method Robust correlation method (only for SIMPLS):
#'   \itemize{
#'     \item "none" = Standard method (default)
#'     \item "spearman" = Spearman rank correlation
#'     \item "winsorized" = Winsorized correlation
#'     \item "biweight" = Biweight midcorrelation
#'     \item "percentage_bend" = Percentage bend correlation
#'   }
#' @param trim Trim proportion for winsorized correlation (default: 0.1)
#' @param beta Bend constant for percentage bend correlation (default: 0.2)
#' @param verbose Print progress messages (default: TRUE)
#'
#' @return A list with class "pls_regression_result" containing:
#'   \describe{
#'     \item{X_loadings}{X loadings matrix (p x ncomp)}
#'     \item{Y_loadings}{Y loadings matrix (q x ncomp)}
#'     \item{X_scores}{X scores matrix (n x ncomp)}
#'     \item{Y_scores}{Y scores matrix (n x ncomp)}
#'     \item{X_weights}{X weights matrix (p x ncomp)}
#'     \item{coefficients}{Regression coefficients array (p x q x ncomp)}
#'     \item{Y_fitted}{Fitted Y values for final ncomp}
#'     \item{Y_residuals}{Residuals for final ncomp}
#'     \item{R2X}{Variance explained in X per component}
#'     \item{R2Y}{Variance explained in Y per component}
#'     \item{R2Y_cum}{Cumulative R-squared for Y}
#'     \item{vip}{VIP scores (p x ncomp)}
#'     \item{algorithm}{Algorithm used}
#'     \item{ncomp}{Number of components}
#'     \item{X_center}{X centering values}
#'     \item{X_scale}{X scaling values (NULL if not scaled)}
#'     \item{Y_center}{Y centering values}
#'     \item{Y_scale}{Y scaling values (NULL if not scaled)}
#'     \item{n}{Number of observations}
#'     \item{p}{Number of predictors}
#'     \item{q}{Number of responses}
#'   }
#'
#' @details
#' PLS regression finds latent components that maximize the covariance between

#' X and Y while also explaining variance in both matrices. It is particularly
#' useful when predictors are highly collinear or when p > n.
#'
#' \strong{SIMPLS} (default): Faster algorithm that produces orthogonal X-scores.
#' Uses sequential orthogonalization of the covariance matrix.
#'
#' \strong{NIPALS}: Classical iterative algorithm. Can handle missing data
#' natively but is slower than SIMPLS.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Simulated data
#' set.seed(42)
#' n <- 100
#' p <- 20
#' X <- matrix(rnorm(n * p), n, p)
#' Y <- X[, 1:3] %*% c(1, 1, 1) + rnorm(n, sd = 0.5)
#'
#' # Fit PLS regression
#' fit <- pls_regression(X, Y, ncomp = 5)
#' print(fit)
#'
#' # Predictions
#' pred <- predict(fit, X[1:10, ])
#' }
pls_regression <- function(X, Y, ncomp = NULL, algorithm = "simpls",
                           scale = TRUE, center = TRUE,
                           max_iter = 500, tol = 1e-6,
                           robust_method = "none",
                           trim = 0.1, beta = 0.2,
                           verbose = TRUE) {

  # Input validation
  if (!is.matrix(X)) X <- as.matrix(X)
  if (!is.matrix(Y)) Y <- as.matrix(Y)

  n <- nrow(X)
  p <- ncol(X)
  q <- ncol(Y)

  if (nrow(Y) != n) {
    stop("X and Y must have the same number of rows")
  }

  # Handle missing data check
  if (any(is.na(X)) || any(is.na(Y))) {
    if (algorithm == "nipals") {
      if (verbose) message("Note: NIPALS can handle missing data")
    } else {
      stop("Missing values detected. Use algorithm='nipals' or impute first.")
    }
  }

  # Determine number of components
  if (is.null(ncomp)) {
    ncomp <- min(n - 1, p, 10)
  }
  ncomp <- min(ncomp, n - 1, p)

  if (ncomp < 1) {
    stop("ncomp must be at least 1")
  }

  # Validate algorithm
 algorithm <- tolower(algorithm)
  if (!algorithm %in% c("simpls", "nipals")) {
    stop("algorithm must be 'simpls' or 'nipals'")
  }

  # Validate robust method
  robust_method <- tolower(robust_method)
  if (!robust_method %in% c("none", "pearson", "spearman", "winsorized",
                             "biweight", "percentage_bend")) {
    stop("Unknown robust_method: ", robust_method)
  }

  if (robust_method != "none" && algorithm == "nipals") {
    warning("Robust methods only apply to SIMPLS. Ignoring robust_method for NIPALS.")
    robust_method <- "none"
  }

  if (verbose) {
    message(sprintf("Fitting PLS regression: %d obs, %d predictors, %d responses, %d components",
                    n, p, q, ncomp))
    message(sprintf("Algorithm: %s", toupper(algorithm)))
  }

  # Centering and scaling
  X_center <- colMeans(X, na.rm = TRUE)
  Y_center <- colMeans(Y, na.rm = TRUE)

  X0 <- sweep(X, 2, X_center, "-")
  Y0 <- sweep(Y, 2, Y_center, "-")

  if (scale) {
    X_scale <- apply(X, 2, sd, na.rm = TRUE)
    Y_scale <- apply(Y, 2, sd, na.rm = TRUE)

    # Handle zero variance
    X_scale[X_scale == 0] <- 1
    Y_scale[Y_scale == 0] <- 1

    X0 <- sweep(X0, 2, X_scale, "/")
    Y0 <- sweep(Y0, 2, Y_scale, "/")
  } else {
    X_scale <- NULL
    Y_scale <- NULL
  }

  # Run algorithm
  if (algorithm == "simpls") {
    pls_result <- pls_simpls(X0, Y0, ncomp,
                              robust_method = robust_method,
                              trim = trim, beta = beta)
  } else {
    pls_result <- pls_nipals(X0, Y0, ncomp,
                              max_iter = max_iter, tol = tol)
  }

  # Extract components
  T_mat <- pls_result$T  # X scores
  U_mat <- pls_result$U  # Y scores
  P_mat <- pls_result$P  # X loadings
  Q_mat <- pls_result$Q  # Y loadings
  W_mat <- pls_result$W  # X weights

  # Calculate regression coefficients for each number of components
  coefficients <- calculate_coefficients(W_mat, P_mat, Q_mat, ncomp)

  # Calculate fitted values and residuals (using all components)
  Y_fitted <- X0 %*% coefficients[, , ncomp, drop = FALSE]
  dim(Y_fitted) <- c(n, q)
  Y_residuals <- Y0 - Y_fitted

  # Unscale fitted values for output
  if (!is.null(Y_scale)) {
    Y_fitted_unscaled <- sweep(Y_fitted, 2, Y_scale, "*")
  } else {
    Y_fitted_unscaled <- Y_fitted
  }
  Y_fitted_unscaled <- sweep(Y_fitted_unscaled, 2, Y_center, "+")

  # Calculate variance explained
  SSX_total <- sum(X0^2, na.rm = TRUE)
  SSY_total <- sum(Y0^2, na.rm = TRUE)

  R2X <- numeric(ncomp)
  R2Y <- numeric(ncomp)
  R2Y_cum <- numeric(ncomp)

  for (a in seq_len(ncomp)) {
    # X variance explained by component a
    t_a <- T_mat[, a, drop = FALSE]
    p_a <- P_mat[, a, drop = FALSE]
    X_explained_a <- t_a %*% t(p_a)
    R2X[a] <- sum(X_explained_a^2, na.rm = TRUE) / SSX_total

    # Y variance explained (cumulative)
    Y_pred_a <- X0 %*% coefficients[, , a, drop = FALSE]
    dim(Y_pred_a) <- c(n, q)
    SS_res_a <- sum((Y0 - Y_pred_a)^2, na.rm = TRUE)
    R2Y_cum[a] <- 1 - SS_res_a / SSY_total

    if (a == 1) {
      R2Y[a] <- R2Y_cum[a]
    } else {
      R2Y[a] <- R2Y_cum[a] - R2Y_cum[a - 1]
    }
  }

  # Calculate VIP scores
  vip <- calculate_vip(W_mat, Q_mat, T_mat, ncomp)

  # Build result object
  result <- list(
    X_loadings = P_mat,
    Y_loadings = Q_mat,
    X_scores = T_mat,
    Y_scores = U_mat,
    X_weights = W_mat,
    coefficients = coefficients,
    Y_fitted = Y_fitted_unscaled,
    Y_residuals = Y_residuals,
    R2X = R2X,
    R2Y = R2Y,
    R2Y_cum = R2Y_cum,
    vip = vip,
    algorithm = algorithm,
    ncomp = ncomp,
    X_center = X_center,
    X_scale = X_scale,
    Y_center = Y_center,
    Y_scale = Y_scale,
    n = n,
    p = p,
    q = q
  )

  class(result) <- "pls_regression_result"

  if (verbose) {
    message(sprintf("Done. Final R2Y = %.4f", R2Y_cum[ncomp]))
  }

  return(result)
}


#' SIMPLS Algorithm for PLS Regression
#'
#' Implements the SIMPLS algorithm (de Jong, 1993) for PLS regression.
#'
#' @param X Centered/scaled predictor matrix (n x p)
#' @param Y Centered/scaled response matrix (n x q)
#' @param ncomp Number of components
#' @param robust_method Robust correlation method
#' @param trim Trim parameter for winsorized
#' @param beta Beta parameter for percentage bend
#'
#' @return List with T, U, P, Q, W matrices
#'
#' @references
#' de Jong, S. (1993). SIMPLS: An alternative approach to partial least squares
#' regression. Chemometrics and Intelligent Laboratory Systems, 18, 251-263.
#'
#' @keywords internal
pls_simpls <- function(X, Y, ncomp, robust_method = "none",
                        trim = 0.1, beta = 0.2) {
  n <- nrow(X)
  p <- ncol(X)
  q <- ncol(Y)

  # Storage matrices
  T_mat <- matrix(0, n, ncomp)      # X scores
  U_mat <- matrix(0, n, ncomp)      # Y scores
  P_mat <- matrix(0, p, ncomp)      # X loadings
  Q_mat <- matrix(0, q, ncomp)      # Y loadings
  W_mat <- matrix(0, p, ncomp)      # X weights
  V_mat <- matrix(0, p, ncomp)      # Orthonormal basis

  # Initial covariance matrix
  if (robust_method == "none" || robust_method == "pearson") {
    S <- t(X) %*% Y
  } else {
    # Use robust cross-correlation scaled to covariance
    S <- rri_xcor_robust(X, Y, cormode = 0, robust_method = robust_method,
                          trim = trim, beta = beta)
    S <- S * (n - 1)
  }

  for (a in seq_len(ncomp)) {
    # SVD of covariance matrix to get weight and loading
    if (q == 1) {
      # Single response: weight is proportional to S
      w <- S / sqrt(sum(S^2))
    } else {
      # Multiple responses: use first left singular vector
      svd_result <- svd(S, nu = 1, nv = 1)
      w <- svd_result$u[, 1]
    }

    # X score
    t_vec <- X %*% w
    t_norm <- sqrt(sum(t_vec^2))
    if (t_norm > 0) {
      t_vec <- t_vec / t_norm
    }

    # X loading
    p_vec <- t(X) %*% t_vec

    # Y loading
    q_vec <- t(Y) %*% t_vec

    # Y score (not strictly needed for regression but useful)
    u_vec <- Y %*% q_vec
    u_norm <- sqrt(sum(u_vec^2))
    if (u_norm > 0) {
      u_vec <- u_vec / u_norm
    }

    # Update orthonormal basis
    if (a == 1) {
      v <- p_vec
    } else {
      # Orthogonalize p against previous v's
      v <- p_vec - V_mat[, 1:(a-1), drop = FALSE] %*%
           (t(V_mat[, 1:(a-1), drop = FALSE]) %*% p_vec)
    }
    v_norm <- sqrt(sum(v^2))
    if (v_norm > 0) {
      v <- v / v_norm
    }

    # Deflate covariance matrix
    S <- S - v %*% (t(v) %*% S)

    # Store results
    T_mat[, a] <- t_vec * t_norm  # Undo normalization for scores
    U_mat[, a] <- u_vec * u_norm
    P_mat[, a] <- p_vec
    Q_mat[, a] <- q_vec
    W_mat[, a] <- w
    V_mat[, a] <- v
  }

  list(T = T_mat, U = U_mat, P = P_mat, Q = Q_mat, W = W_mat)
}


#' NIPALS Algorithm for PLS Regression
#'
#' Implements the NIPALS algorithm for PLS regression.
#'
#' @param X Centered/scaled predictor matrix (n x p)
#' @param Y Centered/scaled response matrix (n x q)
#' @param ncomp Number of components
#' @param max_iter Maximum iterations per component
#' @param tol Convergence tolerance
#'
#' @return List with T, U, P, Q, W matrices
#'
#' @keywords internal
pls_nipals <- function(X, Y, ncomp, max_iter = 500, tol = 1e-6) {
  n <- nrow(X)
  p <- ncol(X)
  q <- ncol(Y)

  # Storage matrices
  T_mat <- matrix(0, n, ncomp)      # X scores
  U_mat <- matrix(0, n, ncomp)      # Y scores
  P_mat <- matrix(0, p, ncomp)      # X loadings
  Q_mat <- matrix(0, q, ncomp)      # Y loadings
  W_mat <- matrix(0, p, ncomp)      # X weights

  # Working copies for deflation
  E <- X
  F <- Y

  for (a in seq_len(ncomp)) {
    # Initialize u as column of F with max variance
    col_var <- colSums(F^2, na.rm = TRUE)
    u <- F[, which.max(col_var), drop = FALSE]

    converged <- FALSE

    for (iter in seq_len(max_iter)) {
      u_old <- u

      # X weight: w = X'u / ||X'u||
      w <- t(E) %*% u
      w_norm <- sqrt(sum(w^2, na.rm = TRUE))
      if (w_norm > 0) {
        w <- w / w_norm
      }

      # X score: t = Xw
      t_vec <- E %*% w

      # Y loading: c = F't / t't
      t_ss <- sum(t_vec^2, na.rm = TRUE)
      if (t_ss > 0) {
        c_vec <- t(F) %*% t_vec / t_ss
      } else {
        c_vec <- t(F) %*% t_vec
      }

      # Y score: u = Fc
      u <- F %*% c_vec

      # Check convergence
      u_diff <- sqrt(sum((u - u_old)^2, na.rm = TRUE))
      if (u_diff < tol) {
        converged <- TRUE
        break
      }
    }

    if (!converged) {
      warning(sprintf("Component %d did not converge in %d iterations", a, max_iter))
    }

    # X loading: p = E't / t't
    t_ss <- sum(t_vec^2, na.rm = TRUE)
    if (t_ss > 0) {
      p_vec <- t(E) %*% t_vec / t_ss
    } else {
      p_vec <- t(E) %*% t_vec
    }

    # Deflate X and F
    E <- E - t_vec %*% t(p_vec)
    F <- F - t_vec %*% t(c_vec)

    # Store results
    T_mat[, a] <- t_vec
    U_mat[, a] <- u
    P_mat[, a] <- p_vec
    Q_mat[, a] <- c_vec
    W_mat[, a] <- w
  }

  list(T = T_mat, U = U_mat, P = P_mat, Q = Q_mat, W = W_mat)
}


#' Calculate PLS Regression Coefficients
#'
#' Computes regression coefficients from PLS weights and loadings.
#'
#' @param W X weights matrix (p x ncomp)
#' @param P X loadings matrix (p x ncomp)
#' @param Q Y loadings matrix (q x ncomp)
#' @param ncomp Number of components
#'
#' @return Array of coefficients (p x q x ncomp)
#'
#' @details
#' Coefficients are computed as: B = W(P'W)^{-1}Q'
#' This gives coefficients for centered/scaled data.
#'
#' @keywords internal
calculate_coefficients <- function(W, P, Q, ncomp) {
  p <- nrow(W)
  q <- nrow(Q)

  # Coefficients for each number of components
  B <- array(0, dim = c(p, q, ncomp))

  for (a in seq_len(ncomp)) {
    W_a <- W[, 1:a, drop = FALSE]
    P_a <- P[, 1:a, drop = FALSE]
    Q_a <- Q[, 1:a, drop = FALSE]

    # B = W * (P'W)^-1 * Q'
    PW <- t(P_a) %*% W_a

    # Use pseudo-inverse for numerical stability
    if (a == 1) {
      PW_inv <- 1 / PW[1, 1]
      B[, , a] <- W_a %*% PW_inv %*% t(Q_a)
    } else {
      PW_inv <- tryCatch(
        solve(PW),
        error = function(e) {
          MASS::ginv(PW)
        }
      )
      B[, , a] <- W_a %*% PW_inv %*% t(Q_a)
    }
  }

  B
}


#' Calculate Variable Importance in Projection (VIP) Scores
#'
#' Computes VIP scores indicating the importance of each X variable.
#'
#' @param W X weights matrix (p x ncomp)
#' @param Q Y loadings matrix (q x ncomp)
#' @param T X scores matrix (n x ncomp)
#' @param ncomp Number of components
#'
#' @return Matrix of VIP scores (p x ncomp)
#'
#' @details
#' VIP scores summarize the importance of each X variable in the PLS model.
#' Variables with VIP > 1 are typically considered important.
#'
#' Formula: VIP_j = sqrt(p * sum_a((SSY_a/SSY_total) * (w_ja/||w_a||)^2))
#'
#' @keywords internal
calculate_vip <- function(W, Q, T, ncomp) {
  p <- nrow(W)

  # Sum of squares explained by Y for each component
  # SSY_a = sum(q_a^2) * sum(t_a^2)
  SSY <- colSums(Q^2) * colSums(T^2)
  SSY_cum <- cumsum(SSY)

  # VIP for each number of components
  VIP <- matrix(0, p, ncomp)

  for (a in seq_len(ncomp)) {
    # Normalized weights
    W_a <- W[, 1:a, drop = FALSE]
    W_norm <- sweep(W_a, 2, sqrt(colSums(W_a^2)), "/")

    # Handle zero norm columns
    W_norm[is.nan(W_norm)] <- 0

    # VIP calculation
    # VIP_j^2 = p * sum_a((SSY_a / SSY_cum) * w_ja^2)
    weights_sq <- W_norm^2
    ss_weights <- SSY[1:a] / SSY_cum[a]

    VIP[, a] <- sqrt(p * rowSums(sweep(weights_sq, 2, ss_weights, "*")))
  }

  VIP
}


#' Predict Method for PLS Regression
#'
#' Predict response values for new predictor data.
#'
#' @param object A pls_regression_result object
#' @param newdata New X data matrix (m x p)
#' @param ncomp Number of components to use. If NULL (default), uses all.
#' @param ... Additional arguments (ignored)
#'
#' @return Predicted Y values (m x q matrix or vector if q=1)
#'
#' @export
predict.pls_regression_result <- function(object, newdata, ncomp = NULL, ...) {
  if (is.null(ncomp)) {
    ncomp <- object$ncomp
  }

  if (ncomp < 1 || ncomp > object$ncomp) {
    stop(sprintf("ncomp must be between 1 and %d", object$ncomp))
  }

  if (!is.matrix(newdata)) {
    newdata <- as.matrix(newdata)
  }

  if (ncol(newdata) != object$p) {
    stop(sprintf("newdata must have %d columns (predictors)", object$p))
  }

  # Center and scale using training parameters
  newdata <- sweep(newdata, 2, object$X_center, "-")
  if (!is.null(object$X_scale)) {
    newdata <- sweep(newdata, 2, object$X_scale, "/")
  }

  # Predict
  Y_pred <- newdata %*% object$coefficients[, , ncomp, drop = FALSE]
  dim(Y_pred) <- c(nrow(newdata), object$q)

  # Unscale Y
  if (!is.null(object$Y_scale)) {
    Y_pred <- sweep(Y_pred, 2, object$Y_scale, "*")
  }
  Y_pred <- sweep(Y_pred, 2, object$Y_center, "+")

  if (object$q == 1) {
    return(as.vector(Y_pred))
  }
  return(Y_pred)
}


#' Print Method for PLS Regression Result
#'
#' @param x A pls_regression_result object
#' @param ... Additional arguments (ignored)
#'
#' @export
print.pls_regression_result <- function(x, ...) {
  cat("PLS Regression Result\n")
  cat("=====================\n\n")
  cat("Algorithm:", toupper(x$algorithm), "\n")
  cat("Components:", x$ncomp, "\n")
  cat("Observations:", x$n, "\n")
  cat("Predictors (p):", x$p, "\n")
  cat("Responses (q):", x$q, "\n\n")

  cat("Variance Explained:\n")
  var_df <- data.frame(
    Component = seq_len(x$ncomp),
    R2X = round(x$R2X * 100, 2),
    R2Y = round(x$R2Y * 100, 2),
    R2Y_cumulative = round(x$R2Y_cum * 100, 2)
  )
  print(var_df, row.names = FALSE)

  cat("\n")

  # VIP summary
  vip_final <- x$vip[, x$ncomp]
  n_important <- sum(vip_final > 1)
  cat(sprintf("Variables with VIP > 1: %d of %d (%.1f%%)\n",
              n_important, x$p, 100 * n_important / x$p))

  invisible(x)
}


#' Summary Method for PLS Regression Result
#'
#' @param object A pls_regression_result object
#' @param top_n Number of top VIP variables to show (default: 10)
#' @param ... Additional arguments (ignored)
#'
#' @export
summary.pls_regression_result <- function(object, top_n = 10, ...) {
  cat("PLS Regression Summary\n")
  cat("======================\n\n")

  # Basic info
  print.pls_regression_result(object)

  # Top VIP variables
  cat("\nTop VIP Scores (Component", object$ncomp, "):\n")
  vip_final <- object$vip[, object$ncomp]
  top_idx <- order(vip_final, decreasing = TRUE)[1:min(top_n, length(vip_final))]

  vip_df <- data.frame(
    Variable = top_idx,
    VIP = round(vip_final[top_idx], 3)
  )
  print(vip_df, row.names = FALSE)

  # Model fit statistics
  cat("\nModel Fit:\n")
  cat(sprintf("  Total R2Y: %.4f\n", object$R2Y_cum[object$ncomp]))
  cat(sprintf("  Total R2X: %.4f\n", sum(object$R2X)))

  # Residual statistics
  cat("\nResiduals:\n")
  res_summary <- summary(as.vector(object$Y_residuals))
  cat(sprintf("  Min: %.4f, Median: %.4f, Max: %.4f\n",
              res_summary["Min."], res_summary["Median"], res_summary["Max."]))
  cat(sprintf("  RMSE: %.4f\n", sqrt(mean(object$Y_residuals^2))))

  invisible(object)
}


#' Extract Coefficients from PLS Regression
#'
#' @param object A pls_regression_result object
#' @param ncomp Number of components to use. If NULL (default), uses all.
#' @param ... Additional arguments (ignored)
#'
#' @return Coefficient matrix (p x q)
#'
#' @export
coef.pls_regression_result <- function(object, ncomp = NULL, ...) {
  if (is.null(ncomp)) {
    ncomp <- object$ncomp
  }

  if (ncomp < 1 || ncomp > object$ncomp) {
    stop(sprintf("ncomp must be between 1 and %d", object$ncomp))
  }

  B <- object$coefficients[, , ncomp]

  # Adjust for scaling
  if (!is.null(object$X_scale) && !is.null(object$Y_scale)) {
    B <- sweep(B, 1, object$X_scale, "/")
    B <- sweep(B, 2, object$Y_scale, "*")
  } else if (!is.null(object$X_scale)) {
    B <- sweep(B, 1, object$X_scale, "/")
  } else if (!is.null(object$Y_scale)) {
    B <- sweep(B, 2, object$Y_scale, "*")
  }

  if (object$q == 1) {
    return(as.vector(B))
  }
  return(B)
}


# ===========================================================================
# Bootstrap and Permutation Inference
# ===========================================================================

#' Bootstrap PLS Regression
#'
#' Perform bootstrap resampling to obtain confidence intervals for
#' PLS regression coefficients and VIP scores.
#'
#' @param X Predictor matrix (n x p)
#' @param Y Response matrix (n x q) or vector
#' @param ncomp Number of components
#' @param num_boot Number of bootstrap samples (default: 1000)
#' @param clim Confidence level as percentage (default: 95)
#' @param algorithm PLS algorithm: "simpls" or "nipals"
#' @param scale Scale X and Y (default: TRUE)
#' @param center Center X and Y (default: TRUE)
#' @param seed Random seed for reproducibility
#' @param verbose Print progress (default: TRUE)
#'
#' @return A pls_regression_result object with additional boot_result field:
#'   \describe{
#'     \item{num_boot}{Number of bootstrap samples}
#'     \item{clim}{Confidence level}
#'     \item{B_se}{Standard errors of coefficients}
#'     \item{B_ci_low}{Lower confidence limits for coefficients}
#'     \item{B_ci_high}{Upper confidence limits for coefficients}
#'     \item{compare_B}{Bootstrap ratios (coefficient / SE)}
#'     \item{vip_se}{Standard errors of VIP scores}
#'     \item{vip_ci_low}{Lower confidence limits for VIP}
#'     \item{vip_ci_high}{Upper confidence limits for VIP}
#'   }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' set.seed(42)
#' X <- matrix(rnorm(100 * 20), 100, 20)
#' Y <- X[, 1:3] %*% c(1, 1, 1) + rnorm(100, sd = 0.5)
#'
#' fit <- pls_regression_boot(X, Y, ncomp = 3, num_boot = 500)
#' # Check which coefficients are significant
#' significant <- abs(fit$boot_result$compare_B) > 2
#' }
pls_regression_boot <- function(X, Y, ncomp, num_boot = 1000, clim = 95,
                                 algorithm = "simpls", scale = TRUE, center = TRUE,
                                 seed = NULL, verbose = TRUE) {

  if (!is.matrix(X)) X <- as.matrix(X)
  if (!is.matrix(Y)) Y <- as.matrix(Y)

  n <- nrow(X)
  p <- ncol(X)
  q <- ncol(Y)

  if (!is.null(seed)) set.seed(seed)

  if (verbose) message("Fitting original model...")

  # Fit original model
  fit <- pls_regression(X, Y, ncomp = ncomp, algorithm = algorithm,
                         scale = scale, center = center, verbose = FALSE)

  if (verbose) message(sprintf("Running %d bootstrap samples...", num_boot))

  # Bootstrap storage
  B_boot <- array(NA, dim = c(p, q, num_boot))
  vip_boot <- matrix(NA, p, num_boot)

  for (b in seq_len(num_boot)) {
    if (verbose && b %% 100 == 0) {
      message(sprintf("  Bootstrap %d of %d", b, num_boot))
    }

    # Resample with replacement
    idx <- sample(n, n, replace = TRUE)

    # Fit on bootstrap sample
    fit_b <- tryCatch({
      pls_regression(X[idx, , drop = FALSE], Y[idx, , drop = FALSE],
                      ncomp = ncomp, algorithm = algorithm,
                      scale = scale, center = center, verbose = FALSE)
    }, error = function(e) NULL)

    if (is.null(fit_b)) next

    # Store coefficients and VIP
    B_boot[, , b] <- fit_b$coefficients[, , ncomp]
    vip_boot[, b] <- fit_b$vip[, ncomp]
  }

  # Remove failed bootstraps
  valid <- !is.na(B_boot[1, 1, ])
  num_valid <- sum(valid)

  if (num_valid < num_boot * 0.9) {
    warning(sprintf("Only %d of %d bootstrap samples succeeded", num_valid, num_boot))
  }

  if (num_valid < 10) {
    stop("Too few valid bootstrap samples")
  }

  # Calculate CIs
  ll <- (100 - clim) / 2
  ul <- 100 - ll

  # Coefficient CIs and SEs
  B_ci_low <- apply(B_boot[, , valid, drop = FALSE], c(1, 2),
                     function(x) quantile(x, ll / 100, na.rm = TRUE))
  B_ci_high <- apply(B_boot[, , valid, drop = FALSE], c(1, 2),
                      function(x) quantile(x, ul / 100, na.rm = TRUE))
  B_se <- apply(B_boot[, , valid, drop = FALSE], c(1, 2), sd, na.rm = TRUE)

  # Handle zero SE
  B_se[B_se == 0] <- .Machine$double.eps

  # Bootstrap ratios
  compare_B <- fit$coefficients[, , ncomp] / B_se

  # VIP CIs and SEs
  vip_ci_low <- apply(vip_boot[, valid, drop = FALSE], 1,
                       function(x) quantile(x, ll / 100, na.rm = TRUE))
  vip_ci_high <- apply(vip_boot[, valid, drop = FALSE], 1,
                        function(x) quantile(x, ul / 100, na.rm = TRUE))
  vip_se <- apply(vip_boot[, valid, drop = FALSE], 1, sd, na.rm = TRUE)

  # Add boot_result to fit
  fit$boot_result <- list(
    num_boot = num_boot,
    num_valid = num_valid,
    clim = clim,
    B_se = B_se,
    B_ci_low = B_ci_low,
    B_ci_high = B_ci_high,
    compare_B = compare_B,
    vip_se = vip_se,
    vip_ci_low = vip_ci_low,
    vip_ci_high = vip_ci_high
  )

  if (verbose) {
    n_sig <- sum(abs(compare_B) > 2, na.rm = TRUE)
    message(sprintf("Done. %d coefficients with |BSR| > 2", n_sig))
  }

  fit
}


#' Permutation Test for PLS Regression
#'
#' Perform permutation testing to assess statistical significance of
#' the PLS regression model.
#'
#' @param X Predictor matrix (n x p)
#' @param Y Response matrix (n x q) or vector
#' @param ncomp Number of components
#' @param num_perm Number of permutations (default: 1000)
#' @param algorithm PLS algorithm: "simpls" or "nipals"
#' @param scale Scale X and Y (default: TRUE)
#' @param center Center X and Y (default: TRUE)
#' @param seed Random seed for reproducibility
#' @param verbose Print progress (default: TRUE)
#'
#' @return A pls_regression_result object with additional perm_result field:
#'   \describe{
#'     \item{num_perm}{Number of permutations}
#'     \item{R2_observed}{Observed R-squared}
#'     \item{R2_perm}{Distribution of permuted R-squared values}
#'     \item{p_value}{Permutation p-value}
#'     \item{per_component}{List with p-values for each component}
#'   }
#'
#' @details
#' Permutation testing breaks the association between X and Y by randomly
#' shuffling Y. The proportion of permuted models with R-squared >= observed
#' R-squared gives the p-value.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' set.seed(42)
#' X <- matrix(rnorm(100 * 20), 100, 20)
#' Y <- X[, 1:3] %*% c(1, 1, 1) + rnorm(100, sd = 0.5)
#'
#' fit <- pls_regression_perm(X, Y, ncomp = 3, num_perm = 500)
#' print(fit$perm_result$p_value)
#' }
pls_regression_perm <- function(X, Y, ncomp, num_perm = 1000,
                                 algorithm = "simpls", scale = TRUE, center = TRUE,
                                 seed = NULL, verbose = TRUE) {

  if (!is.matrix(X)) X <- as.matrix(X)
  if (!is.matrix(Y)) Y <- as.matrix(Y)

  n <- nrow(X)

  if (!is.null(seed)) set.seed(seed)

  if (verbose) message("Fitting original model...")

  # Fit original model
  fit <- pls_regression(X, Y, ncomp = ncomp, algorithm = algorithm,
                         scale = scale, center = center, verbose = FALSE)

  obs_R2 <- fit$R2Y_cum[ncomp]
  obs_R2_per_comp <- fit$R2Y_cum

  if (verbose) message(sprintf("Running %d permutations...", num_perm))

  # Permutation storage
  R2_perm <- numeric(num_perm)
  R2_perm_per_comp <- matrix(NA, num_perm, ncomp)

  for (p_idx in seq_len(num_perm)) {
    if (verbose && p_idx %% 100 == 0) {
      message(sprintf("  Permutation %d of %d", p_idx, num_perm))
    }

    # Permute Y rows
    perm_idx <- sample(n)
    Y_perm <- Y[perm_idx, , drop = FALSE]

    # Fit on permuted data
    fit_p <- tryCatch({
      pls_regression(X, Y_perm, ncomp = ncomp, algorithm = algorithm,
                      scale = scale, center = center, verbose = FALSE)
    }, error = function(e) NULL)

    if (is.null(fit_p)) {
      R2_perm[p_idx] <- NA
      next
    }

    R2_perm[p_idx] <- fit_p$R2Y_cum[ncomp]
    R2_perm_per_comp[p_idx, ] <- fit_p$R2Y_cum
  }

  # Remove failed permutations
  valid <- !is.na(R2_perm)
  num_valid <- sum(valid)

  if (num_valid < num_perm * 0.9) {
    warning(sprintf("Only %d of %d permutations succeeded", num_valid, num_perm))
  }

  # Calculate p-value: proportion of permuted R2 >= observed R2
  # Add 1 to numerator and denominator for unbiased estimate
  p_value <- (sum(R2_perm[valid] >= obs_R2) + 1) / (num_valid + 1)

  # P-values per component
  p_per_comp <- numeric(ncomp)
  for (a in seq_len(ncomp)) {
    p_per_comp[a] <- (sum(R2_perm_per_comp[valid, a] >= obs_R2_per_comp[a]) + 1) /
                      (num_valid + 1)
  }

  # Add perm_result to fit
  fit$perm_result <- list(
    num_perm = num_perm,
    num_valid = num_valid,
    R2_observed = obs_R2,
    R2_perm = R2_perm[valid],
    p_value = p_value,
    per_component = list(
      R2_observed = obs_R2_per_comp,
      p_values = p_per_comp
    )
  )

  if (verbose) {
    message(sprintf("Done. P-value = %.4f", p_value))
  }

  fit
}


# ===========================================================================
# Multiple Imputation Support
# ===========================================================================

#' PLS Regression with Multiple Imputation
#'
#' Perform PLS regression on data with missing values using multiple imputation
#' and pool results using Rubin's rules.
#'
#' @param X Predictor matrix (n x p), may contain NA values
#' @param Y Response matrix (n x q) or vector, may contain NA values
#' @param ncomp Number of components
#' @param m Number of imputations (default: 5)
#' @param impute_method Imputation method: "mean" or "rf" (random forest, default)
#' @param algorithm PLS algorithm: "simpls" (default) or "nipals"
#' @param scale Scale X and Y (default: TRUE)
#' @param center Center X and Y (default: TRUE)
#' @param seed Random seed for reproducibility
#' @param verbose Print progress (default: TRUE)
#' @param ... Additional arguments passed to missRanger (for rf method)
#'
#' @return A list with class "pls_mi_regression_result" containing:
#'   \describe{
#'     \item{All fields from pls_regression_result}{Pooled estimates}
#'     \item{m}{Number of imputations}
#'     \item{impute_method}{Imputation method used}
#'     \item{coefficients_between_var}{Between-imputation variance for coefficients}
#'     \item{vip_between_var}{Between-imputation variance for VIP}
#'     \item{individual_results}{List of individual imputation results}
#'   }
#'
#' @details
#' Multiple imputation creates m complete datasets by imputing missing values,
#' runs PLS regression on each, and pools results using Rubin's rules.
#'
#' \strong{Imputation methods}:
#' \itemize{
#'   \item "mean": Simple mean imputation (fast but may underestimate variance)
#'   \item "rf": Random forest imputation via missRanger (recommended)
#' }
#'
#' \strong{Pooling}: Point estimates are averaged across imputations.
#' Between-imputation variance is computed for uncertainty quantification.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' set.seed(42)
#' X <- matrix(rnorm(100 * 20), 100, 20)
#' Y <- X[, 1:3] %*% c(1, 1, 1) + rnorm(100, sd = 0.5)
#'
#' # Add 10% missing data
#' X[sample(length(X), 0.1 * length(X))] <- NA
#'
#' # Fit with MI
#' fit <- pls_regression_mi(X, Y, ncomp = 3, m = 5)
#' print(fit)
#' }
pls_regression_mi <- function(X, Y, ncomp, m = 5,
                               impute_method = "rf",
                               algorithm = "simpls",
                               scale = TRUE, center = TRUE,
                               seed = NULL, verbose = TRUE, ...) {

  if (!is.matrix(X)) X <- as.matrix(X)
  if (!is.matrix(Y)) Y <- as.matrix(Y)

  n <- nrow(X)
  p <- ncol(X)
  q <- ncol(Y)

  if (!is.null(seed)) set.seed(seed)

  # Check for missing data
  has_na_X <- any(is.na(X))
  has_na_Y <- any(is.na(Y))

  if (!has_na_X && !has_na_Y) {
    if (verbose) message("No missing data detected. Running standard PLS regression.")
    return(pls_regression(X, Y, ncomp = ncomp, algorithm = algorithm,
                           scale = scale, center = center, verbose = verbose))
  }

  # Validate imputation method
  impute_method <- tolower(impute_method)
  if (!impute_method %in% c("mean", "rf")) {
    stop("impute_method must be 'mean' or 'rf'")
  }

  if (impute_method == "rf" && !requireNamespace("missRanger", quietly = TRUE)) {
    warning("missRanger not available, falling back to mean imputation")
    impute_method <- "mean"
  }

  if (verbose) {
    pct_miss_X <- 100 * sum(is.na(X)) / length(X)
    pct_miss_Y <- 100 * sum(is.na(Y)) / length(Y)
    message(sprintf("Missing data: X = %.1f%%, Y = %.1f%%", pct_miss_X, pct_miss_Y))
    message(sprintf("Running %d imputations using %s method...", m, impute_method))
  }

  # Storage for imputation results
  results_list <- vector("list", m)
  coef_array <- array(NA, dim = c(p, q, m))
  vip_mat <- matrix(NA, p, m)
  R2Y_cum_mat <- matrix(NA, ncomp, m)

  for (i in seq_len(m)) {
    if (verbose) message(sprintf("  Imputation %d of %d", i, m))

    # Impute X and Y
    if (impute_method == "mean") {
      X_imp <- impute_mean_matrix(X)
      Y_imp <- impute_mean_matrix(Y)
    } else {
      # Random forest imputation
      XY <- cbind(X, Y)
      XY_df <- as.data.frame(XY)

      XY_imp <- missRanger::missRanger(XY_df, verbose = 0, seed = seed + i, ...)
      XY_imp <- as.matrix(XY_imp)

      X_imp <- XY_imp[, 1:p]
      Y_imp <- XY_imp[, (p + 1):(p + q), drop = FALSE]
    }

    # Fit PLS on imputed data
    fit_i <- tryCatch({
      pls_regression(X_imp, Y_imp, ncomp = ncomp, algorithm = algorithm,
                      scale = scale, center = center, verbose = FALSE)
    }, error = function(e) {
      warning(sprintf("Imputation %d failed: %s", i, e$message))
      return(NULL)
    })

    if (is.null(fit_i)) next

    results_list[[i]] <- fit_i
    coef_array[, , i] <- fit_i$coefficients[, , ncomp]
    vip_mat[, i] <- fit_i$vip[, ncomp]
    R2Y_cum_mat[, i] <- fit_i$R2Y_cum
  }

  # Check how many succeeded
  valid <- !sapply(results_list, is.null)
  num_valid <- sum(valid)

  if (num_valid < m * 0.5) {
    stop(sprintf("Too few imputations succeeded: %d of %d", num_valid, m))
  }

  if (num_valid < m) {
    warning(sprintf("Only %d of %d imputations succeeded", num_valid, m))
  }

  # Pool results using Rubin's rules
  # Use pool_estimates from imputation.R for coefficients
  coef_pooled <- pool_estimates(coef_array[, , valid, drop = FALSE])
  vip_pooled <- pool_estimates(t(vip_mat[, valid, drop = FALSE]))

  # Pool other components by averaging
  ref_fit <- results_list[[which(valid)[1]]]

  # Average loadings, scores, etc.
  X_loadings_avg <- Reduce(`+`, lapply(results_list[valid], `[[`, "X_loadings")) / num_valid
  Y_loadings_avg <- Reduce(`+`, lapply(results_list[valid], `[[`, "Y_loadings")) / num_valid
  X_scores_avg <- Reduce(`+`, lapply(results_list[valid], `[[`, "X_scores")) / num_valid
  Y_scores_avg <- Reduce(`+`, lapply(results_list[valid], `[[`, "Y_scores")) / num_valid
  X_weights_avg <- Reduce(`+`, lapply(results_list[valid], `[[`, "X_weights")) / num_valid
  R2X_avg <- Reduce(`+`, lapply(results_list[valid], `[[`, "R2X")) / num_valid
  R2Y_avg <- Reduce(`+`, lapply(results_list[valid], `[[`, "R2Y")) / num_valid
  R2Y_cum_avg <- Reduce(`+`, lapply(results_list[valid], `[[`, "R2Y_cum")) / num_valid

  # Build pooled result object
  result <- list(
    X_loadings = X_loadings_avg,
    Y_loadings = Y_loadings_avg,
    X_scores = X_scores_avg,
    Y_scores = Y_scores_avg,
    X_weights = X_weights_avg,
    coefficients = ref_fit$coefficients,  # Keep full array structure
    Y_fitted = ref_fit$Y_fitted,
    Y_residuals = ref_fit$Y_residuals,
    R2X = R2X_avg,
    R2Y = R2Y_avg,
    R2Y_cum = R2Y_cum_avg,
    vip = ref_fit$vip,
    algorithm = algorithm,
    ncomp = ncomp,
    X_center = ref_fit$X_center,
    X_scale = ref_fit$X_scale,
    Y_center = ref_fit$Y_center,
    Y_scale = ref_fit$Y_scale,
    n = n,
    p = p,
    q = q,
    # MI-specific fields
    m = m,
    num_valid = num_valid,
    impute_method = impute_method,
    coefficients_pooled = coef_pooled$pooled,
    coefficients_between_var = coef_pooled$between_var,
    vip_pooled = vip_pooled$pooled,
    vip_between_var = vip_pooled$between_var,
    individual_results = results_list[valid]
  )

  # Update coefficients array with pooled values
  result$coefficients[, , ncomp] <- coef_pooled$pooled
  result$vip[, ncomp] <- vip_pooled$pooled

  class(result) <- c("pls_mi_regression_result", "pls_regression_result")

  if (verbose) {
    message(sprintf("Done. Pooled R2Y = %.4f (from %d imputations)", R2Y_cum_avg[ncomp], num_valid))
  }

  result
}


#' Mean Imputation for Matrix
#'
#' Replace NA values with column means.
#'
#' @param X Matrix with potential NA values
#'
#' @return Matrix with NA values replaced by column means
#'
#' @keywords internal
impute_mean_matrix <- function(X) {
  for (j in seq_len(ncol(X))) {
    na_idx <- is.na(X[, j])
    if (any(na_idx)) {
      X[na_idx, j] <- mean(X[, j], na.rm = TRUE)
    }
  }
  X
}


#' Print Method for MI PLS Regression Result
#'
#' @param x A pls_mi_regression_result object
#' @param ... Additional arguments (ignored)
#'
#' @export
print.pls_mi_regression_result <- function(x, ...) {
  cat("PLS Regression Result (Multiple Imputation)\n")
  cat("============================================\n\n")
  cat("Algorithm:", toupper(x$algorithm), "\n")
  cat("Components:", x$ncomp, "\n")
  cat("Observations:", x$n, "\n")
  cat("Predictors (p):", x$p, "\n")
  cat("Responses (q):", x$q, "\n\n")

  cat("Multiple Imputation:\n")
  cat("  Imputations:", x$num_valid, "of", x$m, "\n")
  cat("  Method:", x$impute_method, "\n\n")

  cat("Variance Explained (Pooled):\n")
  var_df <- data.frame(
    Component = seq_len(x$ncomp),
    R2X = round(x$R2X * 100, 2),
    R2Y = round(x$R2Y * 100, 2),
    R2Y_cumulative = round(x$R2Y_cum * 100, 2)
  )
  print(var_df, row.names = FALSE)

  cat("\n")

  # VIP summary
  vip_final <- x$vip[, x$ncomp]
  n_important <- sum(vip_final > 1)
  cat(sprintf("Variables with VIP > 1: %d of %d (%.1f%%)\n",
              n_important, x$p, 100 * n_important / x$p))

  invisible(x)
}
