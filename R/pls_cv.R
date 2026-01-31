#' Cross-Validation for PLS Regression
#'
#' Functions for cross-validation of PLS regression models to determine
#' optimal number of components.
#'
#' @name pls_cv
#' @keywords regression
NULL


#' Cross-Validation for PLS Regression
#'
#' Perform cross-validation to evaluate PLS regression model performance
#' and determine optimal number of components.
#'
#' @param X Predictor matrix (n x p)
#' @param Y Response matrix (n x q) or vector
#' @param ncomp_max Maximum number of components to test. If NULL (default),
#'   uses min(n-1, p, 10)
#' @param method CV method: "kfold" (default) or "loo" (leave-one-out)
#' @param k Number of folds for k-fold CV (default: 10). Ignored if method="loo"
#' @param algorithm PLS algorithm: "simpls" (default) or "nipals"
#' @param scale Scale X and Y (default: TRUE)
#' @param center Center X and Y (default: TRUE)
#' @param seed Random seed for fold assignment (default: NULL)
#' @param verbose Print progress (default: TRUE)
#'
#' @return A list with class "pls_cv_result" containing:
#'   \describe{
#'     \item{PRESS}{Predicted residual sum of squares per component}
#'     \item{RMSEP}{Root mean squared error of prediction per component}
#'     \item{R2_cv}{Cross-validated R-squared per component}
#'     \item{Q2}{Q-squared (predictive R-squared) per component}
#'     \item{optimal_ncomp}{Optimal number of components (minimizes RMSEP)}
#'     \item{method}{CV method used}
#'     \item{k}{Number of folds}
#'     \item{ncomp_max}{Maximum components tested}
#'     \item{fold_ids}{Fold assignments for each observation}
#'   }
#'
#' @details
#' Cross-validation helps determine the optimal number of components by
#' evaluating prediction performance on held-out data.
#'
#' \strong{K-fold CV}: Data is split into k folds. Each fold is used once
#' as test data while the remaining k-1 folds are used for training.
#'
#' \strong{Leave-one-out (LOO)}: Special case of k-fold where k equals n.
#' Each observation is left out once. More thorough but slower.
#'
#' \strong{Metrics}:
#' \itemize{
#'   \item PRESS: Predicted Residual Error Sum of Squares
#'   \item RMSEP: Root Mean Squared Error of Prediction
#'   \item Q2: 1 - PRESS/SS_total (predictive R-squared)
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' set.seed(42)
#' X <- matrix(rnorm(100 * 20), 100, 20)
#' Y <- X[, 1:3] %*% c(1, 1, 1) + rnorm(100, sd = 0.5)
#'
#' # K-fold CV
#' cv <- pls_cv(X, Y, ncomp_max = 10, method = "kfold", k = 5)
#' print(cv)
#' plot(cv)
#'
#' # LOO CV
#' cv_loo <- pls_cv(X, Y, ncomp_max = 5, method = "loo")
#' }
pls_cv <- function(X, Y, ncomp_max = NULL, method = "kfold", k = 10,
                   algorithm = "simpls", scale = TRUE, center = TRUE,
                   seed = NULL, verbose = TRUE) {

  # Input validation
  if (!is.matrix(X)) X <- as.matrix(X)
  if (!is.matrix(Y)) Y <- as.matrix(Y)

  n <- nrow(X)
  p <- ncol(X)
  q <- ncol(Y)

  if (nrow(Y) != n) {
    stop("X and Y must have the same number of rows")
  }

  # Determine max components
  if (is.null(ncomp_max)) {
    ncomp_max <- min(n - 1, p, 10)
  }

  # Validate method
  method <- tolower(method)
  if (!method %in% c("kfold", "loo")) {
    stop("method must be 'kfold' or 'loo'")
  }

  # For LOO, k = n
  if (method == "loo") {
    k <- n
    if (verbose) message(sprintf("Leave-one-out CV with %d observations", n))
  } else {
    if (k < 2 || k > n) {
      stop(sprintf("k must be between 2 and %d", n))
    }
    if (verbose) message(sprintf("%d-fold CV", k))
  }

  # Adjust ncomp_max based on fold sizes
  min_fold_size <- floor(n / k)
  ncomp_max <- min(ncomp_max, min_fold_size - 1, p)

  if (ncomp_max < 1) {
    stop("Not enough observations for cross-validation with this k")
  }

  if (verbose) {
    message(sprintf("Testing 1 to %d components", ncomp_max))
  }

  # Run CV
  result <- pls_cv_kfold(X, Y, ncomp_max, k, algorithm, scale, center,
                          seed, verbose)

  result$method <- method
  result$ncomp_max <- ncomp_max

  class(result) <- "pls_cv_result"
  return(result)
}


#' K-Fold Cross-Validation for PLS (Internal)
#'
#' @param X Predictor matrix
#' @param Y Response matrix
#' @param ncomp_max Maximum components
#' @param k Number of folds
#' @param algorithm PLS algorithm
#' @param scale Scale data
#' @param center Center data
#' @param seed Random seed
#' @param verbose Print progress
#'
#' @return List with CV results
#'
#' @keywords internal
pls_cv_kfold <- function(X, Y, ncomp_max, k, algorithm, scale, center,
                          seed, verbose) {
  n <- nrow(X)
  q <- ncol(Y)

  # Create folds
  if (!is.null(seed)) set.seed(seed)
  fold_ids <- sample(rep(1:k, length.out = n))

  # Storage for predictions
  Y_pred_all <- array(NA, dim = c(n, q, ncomp_max))

  for (fold in 1:k) {
    if (verbose && k <= 20) {
      message(sprintf("  Fold %d of %d", fold, k))
    } else if (verbose && fold %% 10 == 0) {
      message(sprintf("  Fold %d of %d", fold, k))
    }

    train_idx <- which(fold_ids != fold)
    test_idx <- which(fold_ids == fold)

    # Fit model on training data
    fit <- tryCatch({
      pls_regression(
        X = X[train_idx, , drop = FALSE],
        Y = Y[train_idx, , drop = FALSE],
        ncomp = ncomp_max,
        algorithm = algorithm,
        scale = scale,
        center = center,
        verbose = FALSE
      )
    }, error = function(e) {
      warning(sprintf("Fold %d failed: %s", fold, e$message))
      return(NULL)
    })

    if (is.null(fit)) next

    # Predict on test data for each ncomp
    for (a in 1:ncomp_max) {
      pred <- predict(fit, X[test_idx, , drop = FALSE], ncomp = a)
      if (is.vector(pred)) {
        Y_pred_all[test_idx, , a] <- matrix(pred, ncol = 1)
      } else {
        Y_pred_all[test_idx, , a] <- pred
      }
    }
  }

  # Calculate CV metrics
  Y_mean <- colMeans(Y)
  Y_centered <- sweep(Y, 2, Y_mean, "-")
  SST <- sum(Y_centered^2)

  PRESS <- numeric(ncomp_max)
  RMSEP <- numeric(ncomp_max)
  R2_cv <- numeric(ncomp_max)
  Q2 <- numeric(ncomp_max)

  for (a in 1:ncomp_max) {
    residuals <- Y - Y_pred_all[, , a, drop = FALSE]
    dim(residuals) <- c(n, q)

    PRESS[a] <- sum(residuals^2, na.rm = TRUE)
    RMSEP[a] <- sqrt(PRESS[a] / (n * q))
    R2_cv[a] <- 1 - PRESS[a] / SST
    Q2[a] <- 1 - PRESS[a] / SST
  }

  # Determine optimal ncomp (minimum RMSEP)
  optimal_ncomp <- which.min(RMSEP)

  list(
    PRESS = PRESS,
    RMSEP = RMSEP,
    R2_cv = R2_cv,
    Q2 = Q2,
    optimal_ncomp = optimal_ncomp,
    k = k,
    fold_ids = fold_ids,
    Y_pred = Y_pred_all
  )
}


#' Print Method for PLS CV Result
#'
#' @param x A pls_cv_result object
#' @param ... Additional arguments (ignored)
#'
#' @export
print.pls_cv_result <- function(x, ...) {
  method_name <- if (x$method == "loo") "Leave-One-Out" else sprintf("%d-Fold", x$k)

  cat("PLS Cross-Validation Result\n")
  cat("===========================\n\n")
  cat("Method:", method_name, "CV\n")
  cat("Components tested: 1 to", x$ncomp_max, "\n")
  cat("Optimal components:", x$optimal_ncomp, "\n\n")

  cat("CV Metrics:\n")
  cv_df <- data.frame(
    ncomp = 1:x$ncomp_max,
    RMSEP = round(x$RMSEP, 4),
    Q2 = round(x$Q2, 4)
  )

  # Mark optimal
  cv_df$Optimal <- ifelse(1:x$ncomp_max == x$optimal_ncomp, "*", "")

  print(cv_df, row.names = FALSE)

  cat("\n* = Optimal number of components\n")

  invisible(x)
}


#' Summary Method for PLS CV Result
#'
#' @param object A pls_cv_result object
#' @param ... Additional arguments (ignored)
#'
#' @export
summary.pls_cv_result <- function(object, ...) {
  print.pls_cv_result(object)

  cat("\nRecommendation:\n")
  cat(sprintf("  Use %d components (RMSEP = %.4f, Q2 = %.4f)\n",
              object$optimal_ncomp,
              object$RMSEP[object$optimal_ncomp],
              object$Q2[object$optimal_ncomp]))

  # Check for overfitting
  if (object$optimal_ncomp < object$ncomp_max) {
    q2_diff <- object$Q2[object$ncomp_max] - object$Q2[object$optimal_ncomp]
    if (q2_diff < 0) {
      cat(sprintf("\n  Note: Q2 decreases after %d components, suggesting overfitting.\n",
                  object$optimal_ncomp))
    }
  }

  invisible(object)
}


#' Plot Method for PLS CV Result
#'
#' Create plots of cross-validation results.
#'
#' @param x A pls_cv_result object
#' @param metric Which metric to plot: "RMSEP" (default), "Q2", or "both"
#' @param ... Additional arguments passed to ggplot2
#'
#' @return A ggplot2 object
#'
#' @export
plot.pls_cv_result <- function(x, metric = "RMSEP", ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plotting")
  }

  metric <- match.arg(metric, c("RMSEP", "Q2", "both"))

  ncomp_max <- x$ncomp_max

  if (metric == "both") {
    # Two-panel plot
    df <- data.frame(
      ncomp = rep(1:ncomp_max, 2),
      value = c(x$RMSEP, x$Q2),
      Metric = rep(c("RMSEP", "Q2"), each = ncomp_max)
    )

    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$ncomp, y = .data$value)) +
      ggplot2::geom_line(linewidth = 1) +
      ggplot2::geom_point(size = 3) +
      ggplot2::geom_vline(xintercept = x$optimal_ncomp, linetype = "dashed",
                          color = "red", linewidth = 0.8) +
      ggplot2::facet_wrap(~ .data$Metric, scales = "free_y") +
      ggplot2::scale_x_continuous(breaks = 1:ncomp_max) +
      ggplot2::labs(
        title = sprintf("PLS Cross-Validation (%s)",
                        if (x$method == "loo") "LOO" else sprintf("%d-fold", x$k)),
        x = "Number of Components",
        y = "Value"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        panel.grid.minor = ggplot2::element_blank(),
        strip.text = ggplot2::element_text(face = "bold")
      )

  } else {
    # Single metric plot
    df <- data.frame(
      ncomp = 1:ncomp_max,
      value = if (metric == "RMSEP") x$RMSEP else x$Q2
    )

    y_label <- if (metric == "RMSEP") {
      "Root Mean Squared Error of Prediction"
    } else {
      "Q-squared (Predictive R-squared)"
    }

    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$ncomp, y = .data$value)) +
      ggplot2::geom_line(linewidth = 1, color = "#2166AC") +
      ggplot2::geom_point(size = 3, color = "#2166AC") +
      ggplot2::geom_vline(xintercept = x$optimal_ncomp, linetype = "dashed",
                          color = "red", linewidth = 0.8) +
      ggplot2::scale_x_continuous(breaks = 1:ncomp_max) +
      ggplot2::labs(
        title = sprintf("PLS Cross-Validation (%s)",
                        if (x$method == "loo") "LOO" else sprintf("%d-fold", x$k)),
        x = "Number of Components",
        y = y_label
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(panel.grid.minor = ggplot2::element_blank()) +
      ggplot2::annotate("text", x = x$optimal_ncomp + 0.2,
                        y = df$value[x$optimal_ncomp],
                        label = sprintf("Optimal: %d", x$optimal_ncomp),
                        hjust = 0, color = "red", fontface = "bold")
  }

  p
}


#' One Standard Error Rule for Component Selection
#'
#' Select the simplest model within one standard error of the minimum.
#'
#' @param cv_result A pls_cv_result object
#'
#' @return Optimal number of components using 1SE rule
#'
#' @details
#' The one standard error rule selects the simplest model (fewest components)
#' whose error is within one standard error of the minimum error. This often
#' produces more parsimonious models.
#'
#' @export
pls_cv_1se <- function(cv_result) {
  if (!inherits(cv_result, "pls_cv_result")) {
    stop("cv_result must be a pls_cv_result object")
  }

  # Estimate SE from CV folds
  # This is a simplified version; for more accurate SE would need per-fold errors
  n <- length(cv_result$fold_ids)
  k <- cv_result$k

  # Approximate SE as RMSEP / sqrt(k)
  se <- cv_result$RMSEP / sqrt(k)

  # Find minimum and threshold
  min_idx <- which.min(cv_result$RMSEP)
  threshold <- cv_result$RMSEP[min_idx] + se[min_idx]

  # Find simplest model within threshold
  valid_idx <- which(cv_result$RMSEP <= threshold)
  optimal_1se <- min(valid_idx)

  return(optimal_1se)
}
