#' @title Missing Data Handling for PLS Analysis
#' @name imputation
#' @description Functions for checking, imputing, and handling missing data
#'   in PLS analysis with support for multiple imputation and Rubin's rules.
NULL


#' Check Missing Data in PLS Datasets
#'
#' Diagnose missing data patterns and generate warnings based on missingness
#' thresholds. This function should be called before PLS analysis to understand
#' the extent and pattern of missing data.
#'
#' @param datamat_lst List of data matrices, one per group.
#' @param behavdata Optional behavioral data matrix (for behavior PLS).
#' @param warn_threshold Proportion threshold for informational message (default 0.05 = 5%).
#' @param high_threshold Proportion threshold for warning (default 0.20 = 20%).
#' @param critical_threshold Proportion threshold for strong warning (default 0.50 = 50%).
#' @param verbose Print messages (default TRUE).
#'
#' @return A list with class "pls_missing_report" containing:
#' \describe{
#'   \item{has_missing}{Logical indicating if any missing data exists}
#'   \item{total_cells}{Total number of cells across all matrices}
#'   \item{total_missing}{Total number of missing cells}
#'   \item{overall_pct}{Overall percentage missing}
#'   \item{by_group}{List with per-group statistics}
#'   \item{by_variable}{Data frame with per-variable statistics}
#'   \item{by_subject}{Data frame with per-subject (row) statistics}
#'   \item{behavdata_stats}{Statistics for behavioral data (if provided)}
#'   \item{warnings}{Character vector of warning messages}
#'   \item{max_var_pct}{Maximum missing percentage for any variable}
#'   \item{max_subj_pct}{Maximum missing percentage for any subject}
#' }
#'
#' @details
#' Missing data thresholds are based on methodological literature:
#' \itemize{
#'   \item <5\%: Generally acceptable, minimal concern
#'   \item 5-20\%: Common range, imputation recommended
#'   \item 20-50\%: Requires careful consideration
#'   \item >50\%: High risk, consider variable removal
#' }
#'
#' @examples
#' \dontrun{
#' # Check missing data in brain matrices
#' data1 <- matrix(rnorm(300), 30, 10)
#' data1[sample(300, 15)] <- NA  # 5% missing
#'
#' report <- pls_check_missing(list(data1))
#' print(report)
#' }
#'
#' @export
pls_check_missing <- function(datamat_lst,
                              behavdata = NULL,
                              warn_threshold = 0.05,
                              high_threshold = 0.20,
                              critical_threshold = 0.50,
                              verbose = TRUE) {

  if (!is.list(datamat_lst)) {
    stop("datamat_lst must be a list of matrices")
  }

  num_groups <- length(datamat_lst)
  warnings_list <- character(0)

  # Initialize accumulators
  total_cells <- 0
  total_missing <- 0
  by_group <- list()

  # Track variable-level stats (assuming same columns across groups)
  num_vars <- ncol(datamat_lst[[1]])
  var_missing <- rep(0, num_vars)
  var_total <- rep(0, num_vars)

  # Track subject-level stats
  all_subj_pct <- numeric(0)
  all_subj_group <- integer(0)
  all_subj_row <- integer(0)

  # Process each group
  for (g in seq_len(num_groups)) {
    mat <- as.matrix(datamat_lst[[g]])

    if (ncol(mat) != num_vars) {
      stop(sprintf("Group %d has %d columns but group 1 has %d columns",
                   g, ncol(mat), num_vars))
    }

    n_cells <- length(mat)
    n_missing <- sum(is.na(mat))
    pct_missing <- n_missing / n_cells

    total_cells <- total_cells + n_cells
    total_missing <- total_missing + n_missing

    # Per-variable stats for this group
    var_missing_g <- colSums(is.na(mat))
    var_total_g <- rep(nrow(mat), ncol(mat))

    var_missing <- var_missing + var_missing_g
    var_total <- var_total + var_total_g

    # Per-subject stats
    subj_missing <- rowSums(is.na(mat))
    subj_pct <- subj_missing / ncol(mat)

    all_subj_pct <- c(all_subj_pct, subj_pct)
    all_subj_group <- c(all_subj_group, rep(g, nrow(mat)))
    all_subj_row <- c(all_subj_row, seq_len(nrow(mat)))

    by_group[[g]] <- list(
      n_rows = nrow(mat),
      n_cols = ncol(mat),
      n_cells = n_cells,
      n_missing = n_missing,
      pct_missing = pct_missing,
      vars_with_missing = sum(var_missing_g > 0),
      subjects_with_missing = sum(subj_missing > 0)
    )

    # Group-level warnings
    if (pct_missing > critical_threshold) {
      msg <- sprintf("Group %d: %.1f%% missing data (critical - >%.0f%%)",
                     g, pct_missing * 100, critical_threshold * 100)
      warnings_list <- c(warnings_list, msg)
      if (verbose) warning(msg, call. = FALSE)
    } else if (pct_missing > high_threshold) {
      msg <- sprintf("Group %d: %.1f%% missing data (high - >%.0f%%)",
                     g, pct_missing * 100, high_threshold * 100)
      warnings_list <- c(warnings_list, msg)
      if (verbose) warning(msg, call. = FALSE)
    } else if (pct_missing > warn_threshold) {
      msg <- sprintf("Group %d: %.1f%% missing data (moderate - >%.0f%%)",
                     g, pct_missing * 100, warn_threshold * 100)
      warnings_list <- c(warnings_list, msg)
      if (verbose) message(msg)
    }
  }

  # Overall stats
  overall_pct <- total_missing / total_cells
  has_missing <- total_missing > 0

  # Variable-level summary
  var_pct <- var_missing / var_total
  by_variable <- data.frame(
    variable = seq_len(num_vars),
    n_missing = var_missing,
    n_total = var_total,
    pct_missing = var_pct
  )

  # Flag problematic variables
  critical_vars <- which(var_pct > critical_threshold)
  high_vars <- which(var_pct > high_threshold & var_pct <= critical_threshold)

  if (length(critical_vars) > 0) {
    msg <- sprintf("%d variable(s) have >%.0f%% missing data: %s",
                   length(critical_vars), critical_threshold * 100,
                   paste(critical_vars, collapse = ", "))
    warnings_list <- c(warnings_list, msg)
    if (verbose) warning(msg, call. = FALSE)
  }

  if (length(high_vars) > 0) {
    msg <- sprintf("%d variable(s) have >%.0f%% missing data",
                   length(high_vars), high_threshold * 100)
    warnings_list <- c(warnings_list, msg)
    if (verbose) warning(msg, call. = FALSE)
  }

  # Subject-level summary
  by_subject <- data.frame(
    subject_row = all_subj_row,
    group = all_subj_group,
    pct_missing = all_subj_pct
  )

  # Flag problematic subjects
  critical_subj <- which(all_subj_pct > critical_threshold)
  if (length(critical_subj) > 0) {
    msg <- sprintf("%d subject(s) have >%.0f%% missing data",
                   length(critical_subj), critical_threshold * 100)
    warnings_list <- c(warnings_list, msg)
    if (verbose) warning(msg, call. = FALSE)
  }

  # Behavioral data stats
  behavdata_stats <- NULL
  if (!is.null(behavdata)) {
    behavdata <- as.matrix(behavdata)
    behav_missing <- sum(is.na(behavdata))
    behav_total <- length(behavdata)
    behav_pct <- behav_missing / behav_total

    behavdata_stats <- list(
      n_cells = behav_total,
      n_missing = behav_missing,
      pct_missing = behav_pct,
      vars_with_missing = sum(colSums(is.na(behavdata)) > 0)
    )

    if (behav_pct > critical_threshold) {
      msg <- sprintf("Behavioral data: %.1f%% missing (critical)", behav_pct * 100)
      warnings_list <- c(warnings_list, msg)
      if (verbose) warning(msg, call. = FALSE)
    } else if (behav_pct > high_threshold) {
      msg <- sprintf("Behavioral data: %.1f%% missing (high)", behav_pct * 100)
      warnings_list <- c(warnings_list, msg)
      if (verbose) warning(msg, call. = FALSE)
    } else if (behav_pct > warn_threshold && verbose) {
      message(sprintf("Behavioral data: %.1f%% missing", behav_pct * 100))
    }
  }

  result <- list(
    has_missing = has_missing,
    total_cells = total_cells,
    total_missing = total_missing,
    overall_pct = overall_pct,
    by_group = by_group,
    by_variable = by_variable,
    by_subject = by_subject,
    behavdata_stats = behavdata_stats,
    warnings = warnings_list,
    max_var_pct = max(var_pct),
    max_subj_pct = max(all_subj_pct),
    thresholds = list(
      warn = warn_threshold,
      high = high_threshold,
      critical = critical_threshold
    )
  )

  class(result) <- "pls_missing_report"
  return(result)
}


#' Print Missing Data Report
#'
#' @param x A pls_missing_report object
#' @param ... Additional arguments (ignored)
#'
#' @export
print.pls_missing_report <- function(x, ...) {
  cat("Missing Data Report for PLS Analysis\n")
  cat("=====================================\n\n")

  if (!x$has_missing) {
    cat("No missing data detected.\n")
    return(invisible(x))
  }

  cat(sprintf("Overall: %d missing of %d cells (%.2f%%)\n\n",
              x$total_missing, x$total_cells, x$overall_pct * 100))

  cat("By Group:\n")
  for (g in seq_along(x$by_group)) {
    grp <- x$by_group[[g]]
    cat(sprintf("  Group %d: %d missing (%.2f%%), %d vars affected, %d subjects affected\n",
                g, grp$n_missing, grp$pct_missing * 100,
                grp$vars_with_missing, grp$subjects_with_missing))
  }

  cat("\nVariable Summary:\n")
  vars_affected <- sum(x$by_variable$n_missing > 0)
  cat(sprintf("  %d of %d variables have missing data\n",
              vars_affected, nrow(x$by_variable)))
  cat(sprintf("  Max missing per variable: %.2f%%\n", x$max_var_pct * 100))

  cat("\nSubject Summary:\n")
  subj_affected <- sum(x$by_subject$pct_missing > 0)
  cat(sprintf("  %d of %d subjects have missing data\n",
              subj_affected, nrow(x$by_subject)))
  cat(sprintf("  Max missing per subject: %.2f%%\n", x$max_subj_pct * 100))

  if (!is.null(x$behavdata_stats)) {
    cat("\nBehavioral Data:\n")
    cat(sprintf("  %d missing of %d cells (%.2f%%)\n",
                x$behavdata_stats$n_missing,
                x$behavdata_stats$n_cells,
                x$behavdata_stats$pct_missing * 100))
  }

  if (length(x$warnings) > 0) {
    cat("\nWarnings:\n")
    for (w in x$warnings) {
      cat(sprintf("  - %s\n", w))
    }
  }

  invisible(x)
}


#' Impute Missing Data for PLS Analysis
#'
#' Perform single or multiple imputation on PLS data matrices using
#' mean imputation or random forest (via missRanger).
#'
#' @param datamat_lst List of data matrices, one per group.
#' @param behavdata Optional behavioral data matrix.
#' @param method Imputation method: "mean" or "rf" (random forest, default).
#' @param m Number of imputations. Use m=1 for single imputation,
#'   m>1 for multiple imputation (default 1).
#' @param maxiter Maximum iterations for random forest imputation (default 10).
#' @param pmm Use predictive mean matching (default TRUE). Only for method="rf".
#' @param num.trees Number of trees for random forest (default 100).
#' @param seed Random seed for reproducibility.
#' @param verbose Print progress messages (default TRUE).
#'
#' @return If m=1, returns a list with:
#' \describe{
#'   \item{datamat_lst}{List of imputed data matrices}
#'   \item{behavdata}{Imputed behavioral data (if provided)}
#'   \item{method}{Imputation method used}
#' }
#'
#' If m>1, returns a list with class "pls_mi_data" containing:
#' \describe{
#'   \item{m}{Number of imputations}
#'   \item{imputations}{List of m imputed datasets, each containing datamat_lst and behavdata}
#'   \item{method}{Imputation method used}
#' }
#'
#' @details
#' Random forest imputation (method="rf") uses the missRanger package, which
#' implements chained random forests with optional predictive mean matching.
#' This is generally preferred over mean imputation as it:
#' \itemize{
#'   \item Preserves relationships between variables
#'   \item Handles non-linear relationships
#'   \item Maintains appropriate variability
#' }
#'
#' For multiple imputation (m > 1), each imputation uses a different random seed
#' to generate variation across imputed datasets.
#'
#' @examples
#' \dontrun{
#' # Single RF imputation
#' data1 <- matrix(rnorm(300), 30, 10)
#' data1[sample(300, 15)] <- NA
#'
#' imputed <- pls_impute(list(data1), method = "rf", m = 1)
#'
#' # Multiple imputation with 5 datasets
#' mi_data <- pls_impute(list(data1), method = "rf", m = 5)
#' }
#'
#' @export
pls_impute <- function(datamat_lst,
                       behavdata = NULL,
                       method = "rf",
                       m = 1,
                       maxiter = 10,
                       pmm = TRUE,
                       num.trees = 100,
                       seed = NULL,
                       verbose = TRUE) {

  method <- tolower(method)
  if (!method %in% c("mean", "rf")) {
    stop("method must be 'mean' or 'rf'")
  }

  if (m < 1 || m != floor(m)) {
    stop("m must be a positive integer")
  }

  # Check for missRanger if using RF
  if (method == "rf") {
    if (!requireNamespace("missRanger", quietly = TRUE)) {
      stop("Package 'missRanger' is required for random forest imputation. ",
           "Install it with: install.packages('missRanger')")
    }
  }

  # Set seed if provided
  if (!is.null(seed)) {
    set.seed(seed)
  }

  # Check if any missing data exists
  has_missing_data <- any(sapply(datamat_lst, function(x) any(is.na(x))))
  has_missing_behav <- !is.null(behavdata) && any(is.na(behavdata))

  if (!has_missing_data && !has_missing_behav) {
    if (verbose) message("No missing data detected. Returning original data.")
    if (m == 1) {
      return(list(
        datamat_lst = datamat_lst,
        behavdata = behavdata,
        method = method,
        imputed = FALSE
      ))
    } else {
      imputations <- replicate(m, list(
        datamat_lst = datamat_lst,
        behavdata = behavdata
      ), simplify = FALSE)

      result <- list(
        m = m,
        imputations = imputations,
        method = method,
        imputed = FALSE
      )
      class(result) <- "pls_mi_data"
      return(result)
    }
  }

  if (verbose) {
    message(sprintf("Performing %s imputation (m=%d)...",
                    ifelse(method == "rf", "random forest", "mean"), m))
  }

  # Generate m imputations
  imputations <- vector("list", m)

  for (i in seq_len(m)) {
    if (verbose && m > 1) {
      message(sprintf("  Imputation %d of %d", i, m))
    }

    # Impute each group's data matrix
    imputed_datamat_lst <- vector("list", length(datamat_lst))

    for (g in seq_along(datamat_lst)) {
      mat <- as.matrix(datamat_lst[[g]])

      if (!any(is.na(mat))) {
        imputed_datamat_lst[[g]] <- mat
        next
      }

      if (method == "mean") {
        imputed_datamat_lst[[g]] <- impute_mean(mat)
      } else {
        # RF imputation - convert to data.frame for missRanger
        df <- as.data.frame(mat)

        # Use different seed for each imputation to get variation
        imp_seed <- if (!is.null(seed)) seed + i - 1 else NULL

        imputed_df <- missRanger::missRanger(
          df,
          maxiter = maxiter,
          pmm.k = if (pmm) 3 else 0,
          num.trees = num.trees,
          seed = imp_seed,
          verbose = 0
        )

        imputed_datamat_lst[[g]] <- as.matrix(imputed_df)
      }
    }

    # Impute behavioral data if provided
    imputed_behavdata <- NULL
    if (!is.null(behavdata)) {
      behavdata_mat <- as.matrix(behavdata)

      if (!any(is.na(behavdata_mat))) {
        imputed_behavdata <- behavdata_mat
      } else if (method == "mean") {
        imputed_behavdata <- impute_mean(behavdata_mat)
      } else {
        df <- as.data.frame(behavdata_mat)
        imp_seed <- if (!is.null(seed)) seed + i - 1 + 1000 else NULL

        imputed_df <- missRanger::missRanger(
          df,
          maxiter = maxiter,
          pmm.k = if (pmm) 3 else 0,
          num.trees = num.trees,
          seed = imp_seed,
          verbose = 0
        )

        imputed_behavdata <- as.matrix(imputed_df)
      }
    }

    imputations[[i]] <- list(
      datamat_lst = imputed_datamat_lst,
      behavdata = imputed_behavdata
    )
  }

  # Return appropriate structure

if (m == 1) {
    return(list(
      datamat_lst = imputations[[1]]$datamat_lst,
      behavdata = imputations[[1]]$behavdata,
      method = method,
      imputed = TRUE
    ))
  } else {
    result <- list(
      m = m,
      imputations = imputations,
      method = method,
      imputed = TRUE
    )
    class(result) <- "pls_mi_data"
    return(result)
  }
}


#' Mean Imputation (Internal)
#'
#' Replace NA values with column means.
#'
#' @param mat Matrix with missing values
#' @return Matrix with NA replaced by column means
#' @keywords internal
impute_mean <- function(mat) {
  col_means <- colMeans(mat, na.rm = TRUE)

  for (j in seq_len(ncol(mat))) {
    na_idx <- is.na(mat[, j])
    if (any(na_idx)) {
      mat[na_idx, j] <- col_means[j]
    }
  }

  # Handle columns that are all NA
  still_na <- which(is.na(mat), arr.ind = TRUE)
  if (nrow(still_na) > 0) {
    mat[still_na] <- 0
    warning("Some columns were entirely NA; replaced with 0")
  }

  return(mat)
}


#' Print Multiple Imputation Data
#'
#' @param x A pls_mi_data object
#' @param ... Additional arguments (ignored)
#'
#' @export
print.pls_mi_data <- function(x, ...) {
  cat("Multiple Imputation Data for PLS Analysis\n")
  cat("==========================================\n\n")

  cat(sprintf("Number of imputations: %d\n", x$m))
  cat(sprintf("Imputation method: %s\n", x$method))
  cat(sprintf("Data was imputed: %s\n", x$imputed))

  if (length(x$imputations) > 0) {
    first_imp <- x$imputations[[1]]
    cat(sprintf("Number of groups: %d\n", length(first_imp$datamat_lst)))

    for (g in seq_along(first_imp$datamat_lst)) {
      mat <- first_imp$datamat_lst[[g]]
      cat(sprintf("  Group %d: %d rows x %d columns\n", g, nrow(mat), ncol(mat)))
    }

    if (!is.null(first_imp$behavdata)) {
      cat(sprintf("Behavioral data: %d rows x %d columns\n",
                  nrow(first_imp$behavdata), ncol(first_imp$behavdata)))
    }
  }

  invisible(x)
}


#' Pool Estimates Using Rubin's Rules
#'
#' Combine parameter estimates from multiple imputations using Rubin's rules.
#'
#' @param estimates Matrix or array of estimates. If matrix, rows are parameters
#'   and columns are imputations. If array, third dimension is imputations.
#' @param variances Matrix or array of variance estimates (same structure as estimates).
#'   If NULL, only point estimates are pooled.
#'
#' @return A list containing:
#' \describe{
#'   \item{pooled}{Pooled point estimates}
#'   \item{within_var}{Within-imputation variance (if variances provided)}
#'   \item{between_var}{Between-imputation variance}
#'   \item{total_var}{Total variance (if variances provided)}
#'   \item{df}{Degrees of freedom (Barnard-Rubin)}
#'   \item{fmi}{Fraction of missing information}
#'   \item{m}{Number of imputations}
#' }
#'
#' @details
#' Rubin's rules:
#' \itemize{
#'   \item Point estimate: θ̄ = (1/m) Σ θᵢ
#'   \item Within-imputation variance: W̄ = (1/m) Σ Wᵢ
#'   \item Between-imputation variance: B = (1/(m-1)) Σ (θᵢ - θ̄)²
#'   \item Total variance: T = W̄ + (1 + 1/m) × B
#' }
#'
#' @export
pool_estimates <- function(estimates, variances = NULL) {

  # Handle matrix input (params x imputations)
  if (is.matrix(estimates)) {
    m <- ncol(estimates)

    # Pooled point estimate
    pooled <- rowMeans(estimates)

    # Between-imputation variance
    between_var <- apply(estimates, 1, var)

    if (!is.null(variances)) {
      if (!is.matrix(variances) || !identical(dim(variances), dim(estimates))) {
        stop("variances must have same dimensions as estimates")
      }

      # Within-imputation variance
      within_var <- rowMeans(variances)

      # Total variance with correction factor
      total_var <- within_var + (1 + 1/m) * between_var

      # Degrees of freedom (Barnard-Rubin)
      r <- (1 + 1/m) * between_var / within_var
      r[!is.finite(r)] <- 0
      df <- (m - 1) * (1 + 1/r)^2
      df[!is.finite(df)] <- Inf

      # Fraction of missing information
      fmi <- (between_var + between_var/m) / total_var
      fmi[!is.finite(fmi)] <- 0

      return(list(
        pooled = pooled,
        within_var = within_var,
        between_var = between_var,
        total_var = total_var,
        df = df,
        fmi = fmi,
        m = m
      ))
    } else {
      return(list(
        pooled = pooled,
        between_var = between_var,
        m = m
      ))
    }
  }

  # Handle array input (rows x cols x imputations)
  if (is.array(estimates) && length(dim(estimates)) == 3) {
    m <- dim(estimates)[3]

    # Pooled point estimate
    pooled <- apply(estimates, c(1, 2), mean)

    # Between-imputation variance
    between_var <- apply(estimates, c(1, 2), var)

    if (!is.null(variances)) {
      if (!is.array(variances) || !identical(dim(variances), dim(estimates))) {
        stop("variances must have same dimensions as estimates")
      }

      within_var <- apply(variances, c(1, 2), mean)
      total_var <- within_var + (1 + 1/m) * between_var

      r <- (1 + 1/m) * between_var / within_var
      r[!is.finite(r)] <- 0
      df <- (m - 1) * (1 + 1/r)^2
      df[!is.finite(df)] <- Inf

      fmi <- (between_var + between_var/m) / total_var
      fmi[!is.finite(fmi)] <- 0

      return(list(
        pooled = pooled,
        within_var = within_var,
        between_var = between_var,
        total_var = total_var,
        df = df,
        fmi = fmi,
        m = m
      ))
    } else {
      return(list(
        pooled = pooled,
        between_var = between_var,
        m = m
      ))
    }
  }

  stop("estimates must be a matrix or 3D array")
}


#' Pool P-Values Using Median Method
#'
#' Combine p-values from multiple imputations using the median p-value method.
#'
#' @param pvalues Matrix of p-values. Rows are tests, columns are imputations.
#'
#' @return Vector of pooled p-values (medians).
#'
#' @details
#' The median p-value method is a simple and robust approach for combining
#' p-values from multiply imputed data. It is particularly useful for
#' permutation test p-values where Rubin's rules cannot be directly applied.
#'
#' @references
#' Eekhout, I., et al. (2017). Methods for significance testing of categorical
#' covariates in logistic regression models after multiple imputation.
#' BMC Medical Research Methodology.
#'
#' @export
pool_pvalues <- function(pvalues) {
  if (!is.matrix(pvalues)) {
    pvalues <- matrix(pvalues, nrow = 1)
  }

  apply(pvalues, 1, median)
}


#' Pool PLS Results from Multiple Imputations
#'
#' Combine PLS analysis results from multiple imputations using Rubin's rules
#' for point estimates and median method for p-values.
#'
#' @param results_list List of pls_result objects from each imputation.
#' @param align_lvs Align latent variables across imputations using Procrustes
#'   rotation (default TRUE).
#' @param reference_idx Which imputation to use as reference for alignment (default 1).
#'
#' @return A list with class "pls_mi_result" containing pooled estimates.
#'
#' @details
#' Pooling strategy:
#' \itemize{
#'   \item \strong{u, v} (saliences): Procrustes alignment then average
#'   \item \strong{s} (singular values): Average across imputations
#'   \item \strong{scores}: Average across imputations
#'   \item \strong{p-values}: Median p-value method
#'   \item \strong{bootstrap ratios}: Rubin's rules with variance pooling
#' }
#'
#' @export
pool_pls_results <- function(results_list, align_lvs = TRUE, reference_idx = 1) {

  m <- length(results_list)
  if (m < 2) {
    stop("Need at least 2 results to pool")
  }

  # Validate all results have same structure
  ref <- results_list[[reference_idx]]
  method <- ref$method
  n_lvs <- length(ref$s)

  for (i in seq_len(m)) {
    if (results_list[[i]]$method != method) {
      stop("All results must use the same PLS method")
    }
    if (length(results_list[[i]]$s) != n_lvs) {
      stop("All results must have the same number of LVs")
    }
  }

  # Align LVs if requested
  if (align_lvs) {
    results_list <- align_pls_results(results_list, reference_idx)
  }

  # Pool singular values
  s_mat <- sapply(results_list, function(r) r$s)
  pooled_s <- rowMeans(s_mat)

  # Pool u (brain saliences)
  u_array <- array(NA, dim = c(nrow(ref$u), ncol(ref$u), m))
  for (i in seq_len(m)) {
    u_array[, , i] <- results_list[[i]]$u
  }
  u_pooled <- pool_estimates(u_array)

  # Pool v (design/behavior saliences)
  v_array <- array(NA, dim = c(nrow(ref$v), ncol(ref$v), m))
  for (i in seq_len(m)) {
    v_array[, , i] <- results_list[[i]]$v
  }
  v_pooled <- pool_estimates(v_array)

  # Pool scores
  usc_array <- array(NA, dim = c(nrow(ref$usc), ncol(ref$usc), m))
  for (i in seq_len(m)) {
    usc_array[, , i] <- results_list[[i]]$usc
  }
  usc_pooled <- pool_estimates(usc_array)

  # Initialize pooled result
  pooled <- list(
    method = method,
    u = u_pooled$pooled,
    s = pooled_s,
    v = v_pooled$pooled,
    usc = usc_pooled$pooled,
    u_between_var = u_pooled$between_var,
    v_between_var = v_pooled$between_var,
    m = m
  )

  # Pool vsc if present
  if (!is.null(ref$vsc)) {
    vsc_array <- array(NA, dim = c(nrow(ref$vsc), ncol(ref$vsc), m))
    for (i in seq_len(m)) {
      vsc_array[, , i] <- results_list[[i]]$vsc
    }
    vsc_pooled <- pool_estimates(vsc_array)
    pooled$vsc <- vsc_pooled$pooled
  }

  # Pool lvcorrs if present (behavior PLS)
  if (!is.null(ref$lvcorrs)) {
    lvcorrs_array <- array(NA, dim = c(nrow(ref$lvcorrs), ncol(ref$lvcorrs), m))
    for (i in seq_len(m)) {
      lvcorrs_array[, , i] <- results_list[[i]]$lvcorrs
    }
    lvcorrs_pooled <- pool_estimates(lvcorrs_array)
    pooled$lvcorrs <- lvcorrs_pooled$pooled
    pooled$lvcorrs_between_var <- lvcorrs_pooled$between_var
  }

  # Pool permutation results if present
  if (!is.null(ref$perm_result)) {
    pval_mat <- sapply(results_list, function(r) r$perm_result$sprob)
    if (is.vector(pval_mat)) pval_mat <- matrix(pval_mat, nrow = 1)

    pooled$perm_result <- list(
      num_perm = ref$perm_result$num_perm,
      sprob = pool_pvalues(pval_mat),
      sprob_all = pval_mat,  # Keep individual p-values for reference
      pooling_method = "median"
    )
  }

  # Pool bootstrap results if present
  if (!is.null(ref$boot_result)) {
    # Pool bootstrap ratios
    compare_u_array <- array(NA, dim = c(nrow(ref$boot_result$compare_u),
                                         ncol(ref$boot_result$compare_u), m))
    u_se_array <- array(NA, dim = c(nrow(ref$boot_result$u_se),
                                    ncol(ref$boot_result$u_se), m))

    for (i in seq_len(m)) {
      compare_u_array[, , i] <- results_list[[i]]$boot_result$compare_u
      u_se_array[, , i] <- results_list[[i]]$boot_result$u_se
    }

    # Use Rubin's rules for bootstrap ratios (treating SE^2 as variance)
    compare_u_pooled <- pool_estimates(compare_u_array, u_se_array^2)

    pooled$boot_result <- list(
      num_boot = ref$boot_result$num_boot,
      clim = ref$boot_result$clim,
      compare_u = compare_u_pooled$pooled,
      u_se = sqrt(compare_u_pooled$total_var),
      fmi = compare_u_pooled$fmi
    )

    # Pool correlation CIs if present (behavior PLS)
    if (!is.null(ref$boot_result$orig_corr)) {
      orig_corr_array <- array(NA, dim = c(nrow(ref$boot_result$orig_corr),
                                           ncol(ref$boot_result$orig_corr), m))
      for (i in seq_len(m)) {
        orig_corr_array[, , i] <- results_list[[i]]$boot_result$orig_corr
      }
      orig_corr_pooled <- pool_estimates(orig_corr_array)
      pooled$boot_result$orig_corr <- orig_corr_pooled$pooled

      # Pool confidence limits
      if (!is.null(ref$boot_result$ulcorr)) {
        ulcorr_array <- array(NA, dim = c(nrow(ref$boot_result$ulcorr),
                                          ncol(ref$boot_result$ulcorr), m))
        llcorr_array <- array(NA, dim = c(nrow(ref$boot_result$llcorr),
                                          ncol(ref$boot_result$llcorr), m))
        for (i in seq_len(m)) {
          ulcorr_array[, , i] <- results_list[[i]]$boot_result$ulcorr
          llcorr_array[, , i] <- results_list[[i]]$boot_result$llcorr
        }
        pooled$boot_result$ulcorr <- apply(ulcorr_array, c(1, 2), mean)
        pooled$boot_result$llcorr <- apply(llcorr_array, c(1, 2), mean)
      }
    }
  }

  # Store additional metadata
  pooled$num_subj_lst <- ref$num_subj_lst
  pooled$num_conditions <- ref$num_conditions
  pooled$other_input <- ref$other_input
  pooled$individual_results <- results_list

  class(pooled) <- c("pls_mi_result", "pls_result")
  return(pooled)
}


#' Align PLS Results Across Imputations
#'
#' Use Procrustes rotation to align LVs across imputations to a reference.
#'
#' @param results_list List of pls_result objects.
#' @param reference_idx Index of reference result (default 1).
#'
#' @return List of aligned pls_result objects.
#'
#' @keywords internal
align_pls_results <- function(results_list, reference_idx = 1) {
  m <- length(results_list)
  ref_v <- results_list[[reference_idx]]$v

  for (i in seq_len(m)) {
    if (i == reference_idx) next

    # Get rotation matrix using Procrustes
    rotatemat <- rri_bootprocrust(ref_v, results_list[[i]]$v)

    # Apply rotation to v
    s_i <- results_list[[i]]$s
    results_list[[i]]$v <- results_list[[i]]$v %*% diag(s_i, nrow = length(s_i)) %*% rotatemat

    # Recalculate s from rotated v
    results_list[[i]]$s <- sqrt(colSums(results_list[[i]]$v^2))

    # Normalize v
    for (lv in seq_len(ncol(results_list[[i]]$v))) {
      results_list[[i]]$v[, lv] <- results_list[[i]]$v[, lv] / results_list[[i]]$s[lv]
    }

    # Apply rotation to u
    results_list[[i]]$u <- results_list[[i]]$u %*% diag(s_i, nrow = length(s_i)) %*% rotatemat

    # Normalize u
    for (lv in seq_len(ncol(results_list[[i]]$u))) {
      norm_u <- sqrt(sum(results_list[[i]]$u[, lv]^2))
      if (norm_u > 0) {
        results_list[[i]]$u[, lv] <- results_list[[i]]$u[, lv] / norm_u
      }
    }

    # Rotate bootstrap results if present
    if (!is.null(results_list[[i]]$boot_result$compare_u)) {
      results_list[[i]]$boot_result$compare_u <-
        results_list[[i]]$boot_result$compare_u %*% rotatemat
    }
  }

  return(results_list)
}


#' PLS Analysis with Multiple Imputation
#'
#' Perform PLS analysis with multiple imputation, automatically handling
#' missing data and pooling results using Rubin's rules.
#'
#' @param datamat_lst List of data matrices, one per group. May contain NA values.
#' @param num_subj_lst Vector of number of subjects per group.
#' @param num_cond Number of conditions.
#' @param option PLS analysis options (see \code{\link{pls_analysis}}).
#' @param mi_option Multiple imputation options:
#' \describe{
#'   \item{method}{Imputation method: "rf" (default) or "mean"}
#'   \item{m}{Number of imputations (default 5)}
#'   \item{maxiter}{Max iterations for RF imputation (default 10)}
#'   \item{pmm}{Use predictive mean matching (default TRUE)}
#'   \item{seed}{Random seed for reproducibility}
#' }
#'
#' @return A pls_mi_result object containing pooled estimates.
#'
#' @details
#' This function:
#' \enumerate{
#'   \item Checks for missing data and generates warnings
#'   \item Creates m imputed datasets using random forest or mean imputation
#'   \item Runs PLS analysis on each imputed dataset
#'   \item Pools results using Rubin's rules for point estimates and
#'         median method for p-values
#' }
#'
#' @examples
#' \dontrun{
#' # Create data with missing values
#' data1 <- matrix(rnorm(300), 30, 10)
#' data1[sample(300, 15)] <- NA
#'
#' # Run PLS with multiple imputation
#' result <- pls_analysis_mi(
#'   datamat_lst = list(data1),
#'   num_subj_lst = c(10),
#'   num_cond = 3,
#'   option = list(method = 1, num_perm = 100, num_boot = 100),
#'   mi_option = list(m = 5, method = "rf")
#' )
#'
#' # View pooled p-values
#' print(result$perm_result$sprob)
#' }
#'
#' @export
pls_analysis_mi <- function(datamat_lst,
                            num_subj_lst,
                            num_cond,
                            option = NULL,
                            mi_option = NULL) {

  # Default MI options
  mi_defaults <- list(
    method = "rf",
    m = 5,
    maxiter = 10,
    pmm = TRUE,
    num.trees = 100,
    seed = NULL
  )

  if (!is.null(mi_option)) {
    for (nm in names(mi_option)) {
      mi_defaults[[nm]] <- mi_option[[nm]]
    }
  }
  mi_option <- mi_defaults

  # Extract behavioral data from option if present
  behavdata <- option$stacked_behavdata

  # Check missing data first
  missing_report <- pls_check_missing(
    datamat_lst,
    behavdata = behavdata,
    verbose = TRUE
  )

  if (!missing_report$has_missing) {
    message("No missing data detected. Running standard PLS analysis.")
    result <- pls_analysis(datamat_lst, num_subj_lst, num_cond, option)
    result$mi_info <- list(used_mi = FALSE, reason = "no_missing_data")
    return(result)
  }

  # Perform multiple imputation
  message(sprintf("\nPerforming multiple imputation (m=%d, method=%s)...",
                  mi_option$m, mi_option$method))

  mi_data <- pls_impute(
    datamat_lst,
    behavdata = behavdata,
    method = mi_option$method,
    m = mi_option$m,
    maxiter = mi_option$maxiter,
    pmm = mi_option$pmm,
    num.trees = mi_option$num.trees,
    seed = mi_option$seed,
    verbose = TRUE
  )

  # Run PLS on each imputed dataset
  message("\nRunning PLS analysis on each imputation...")
  results_list <- vector("list", mi_option$m)

  for (i in seq_len(mi_option$m)) {
    message(sprintf("  Analysis %d of %d", i, mi_option$m))

    # Update option with imputed behavioral data if present
    option_i <- option
    if (!is.null(behavdata)) {
      option_i$stacked_behavdata <- mi_data$imputations[[i]]$behavdata
    }

    # Suppress messages from individual analyses
    results_list[[i]] <- suppressMessages(
      pls_analysis(
        mi_data$imputations[[i]]$datamat_lst,
        num_subj_lst,
        num_cond,
        option_i
      )
    )
  }

  # Pool results
  message("\nPooling results using Rubin's rules...")
  pooled <- pool_pls_results(results_list, align_lvs = TRUE)

  # Add MI metadata
  pooled$mi_info <- list(
    used_mi = TRUE,
    m = mi_option$m,
    method = mi_option$method,
    missing_report = missing_report
  )

  message("Multiple imputation PLS analysis complete.")
  return(pooled)
}


#' Print Multiple Imputation PLS Result
#'
#' @param x A pls_mi_result object
#' @param ... Additional arguments (ignored)
#'
#' @export
print.pls_mi_result <- function(x, ...) {
  cat("PLS Analysis Result (Multiple Imputation)\n")
  cat("==========================================\n\n")

  method_names <- c(
    "Mean-Centering Task PLS",
    "Non-Rotated Task PLS",
    "Regular Behavior PLS",
    "Multiblock PLS",
    "Non-Rotated Behavior PLS",
    "Non-Rotated Multiblock PLS"
  )

  cat("Method:", method_names[x$method], "\n")
  cat("Number of LVs:", length(x$s), "\n")
  cat("Singular values:", round(x$s, 4), "\n\n")

  cat("Multiple Imputation:\n")
  cat("  Number of imputations:", x$m, "\n")
  if (!is.null(x$mi_info)) {
    cat("  Imputation method:", x$mi_info$method, "\n")
    if (!is.null(x$mi_info$missing_report)) {
      cat("  Original missing data:", sprintf("%.2f%%",
          x$mi_info$missing_report$overall_pct * 100), "\n")
    }
  }
  cat("\n")

  if (!is.null(x$perm_result)) {
    cat("Permutation Test Results (pooled using median method):\n")
    cat("  Number of permutations:", x$perm_result$num_perm, "\n")
    cat("  Pooled p-values:", round(x$perm_result$sprob, 4), "\n")

    if (!is.null(x$perm_result$sprob_all)) {
      cat("  P-value range across imputations:\n")
      for (lv in seq_along(x$perm_result$sprob)) {
        pvals <- x$perm_result$sprob_all[lv, ]
        cat(sprintf("    LV%d: [%.4f, %.4f]\n", lv, min(pvals), max(pvals)))
      }
    }
    cat("\n")
  }

  if (!is.null(x$boot_result)) {
    cat("Bootstrap Results (pooled using Rubin's rules):\n")
    cat("  Number of bootstraps:", x$boot_result$num_boot, "\n")
    cat("  Confidence level:", x$boot_result$clim, "%\n")

    if (!is.null(x$boot_result$fmi)) {
      avg_fmi <- mean(x$boot_result$fmi, na.rm = TRUE)
      cat("  Average fraction of missing information:", sprintf("%.2f%%", avg_fmi * 100), "\n")
    }
  }

  invisible(x)
}
