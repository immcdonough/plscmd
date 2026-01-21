#' Partial Least Squares Analysis
#'
#' Run PLS analysis on one or more data matrices with optional permutation
#' testing and bootstrap resampling for statistical inference.
#'
#' @param datamat_lst List of data matrices, one per group. Each matrix should
#'   have subjects in rows (organized by condition) and variables in columns.
#' @param num_subj_lst Vector of number of subjects per group. Can also be a
#'   list of vectors for split-subject-blocks (different subjects per condition).
#' @param num_cond Number of conditions.
#' @param option Optional list of analysis parameters:
#'   \describe{
#'     \item{method}{PLS method (1-6, default 1):
#'       \itemize{
#'         \item 1 = Mean-Centering Task PLS
#'         \item 2 = Non-Rotated Task PLS
#'         \item 3 = Regular Behavior PLS
#'         \item 4 = Multiblock PLS
#'         \item 5 = Non-Rotated Behavior PLS
#'         \item 6 = Non-Rotated Multiblock PLS
#'       }}
#'     \item{num_perm}{Number of permutations (default 0 = no permutation test)}
#'     \item{num_boot}{Number of bootstrap samples (default 0 = no bootstrap)}
#'     \item{num_split}{Number of split-half iterations (default 0)}
#'     \item{clim}{Confidence level 0-100 (default 95)}
#'     \item{bscan}{Subset of conditions for behavior block in multiblock PLS}
#'     \item{stacked_designdata}{Design contrast matrix for non-rotated PLS}
#'     \item{stacked_behavdata}{Stacked behavioral data for behavior PLS}
#'     \item{meancentering_type}{Mean-centering type 0-3 (default 0)}
#'     \item{cormode}{Correlation mode: 0=Pearson, 2=covariance, 4=cosine, 6=dot product}
#'     \item{boot_type}{"strat" (default) or "nonstrat" bootstrap}
#'     \item{is_struct}{Set to TRUE for structure PLS (default FALSE)}
#'     \item{verbose}{Print progress messages (default TRUE)}
#'     \item{robust_method}{Robust correlation method (only for cormode=0):
#'       "none" (default), "spearman", "winsorized", "biweight", "percentage_bend"}
#'     \item{robust_trim}{Trim proportion for winsorized correlation (default 0.1)}
#'     \item{robust_beta}{Bend constant for percentage bend correlation (default 0.2)}
#'     \item{impute_method}{Missing data imputation method: "none" (default, error on NA),
#'       "mean" (mean imputation), or "rf" (random forest via missRanger)}
#'     \item{check_missing}{Check for missing data and print report (default TRUE)}
#'     \item{missing_warn_threshold}{Threshold for missing data warning (default 0.05 = 5\%)}
#'   }
#'
#' @return A list with class "pls_result" containing:
#'   \describe{
#'     \item{method}{PLS method used}
#'     \item{u}{Brain latent variable (salience)}
#'     \item{s}{Singular values}
#'     \item{v}{Design/Behavior latent variable}
#'     \item{usc}{Brain scores}
#'     \item{vsc}{Design/Behavior scores}
#'     \item{lvcorrs}{LV correlations (behavior PLS only)}
#'     \item{datamatcorrs_lst}{Correlation maps (behavior PLS only)}
#'     \item{perm_result}{Permutation test results (if num_perm > 0)}
#'     \item{boot_result}{Bootstrap results (if num_boot > 0)}
#'     \item{perm_splithalf}{Split-half results (if num_split > 0)}
#'   }
#'
#' @details
#' PLS Methods:
#' \itemize{
#'   \item \strong{Method 1 (Mean-Centering Task PLS)}: Standard task PLS using
#'     SVD on mean-centered condition averages.
#'   \item \strong{Method 2 (Non-Rotated Task PLS)}: Uses design contrasts
#'     instead of SVD rotation.
#'   \item \strong{Method 3 (Regular Behavior PLS)}: Correlates brain activity
#'     with behavioral measures.
#'   \item \strong{Method 4 (Multiblock PLS)}: Combines task and behavior blocks.
#'   \item \strong{Method 5 (Non-Rotated Behavior PLS)}: Behavior PLS without rotation.
#'   \item \strong{Method 6 (Non-Rotated Multiblock PLS)}: Multiblock without rotation.
#' }
#'
#' @examples
#' \dontrun{
#' # Simple task PLS with 2 groups, 3 conditions, 10 subjects each
#' datamat1 <- matrix(rnorm(10 * 3 * 100), nrow = 30, ncol = 100)
#' datamat2 <- matrix(rnorm(10 * 3 * 100), nrow = 30, ncol = 100)
#'
#' result <- pls_analysis(
#'   datamat_lst = list(datamat1, datamat2),
#'   num_subj_lst = c(10, 10),
#'   num_cond = 3,
#'   option = list(
#'     method = 1,
#'     num_perm = 500,
#'     num_boot = 500
#'   )
#' )
#'
#' # View significant LVs
#' print(result$perm_result$sprob)
#' }
#'
#' @export
pls_analysis <- function(datamat_lst, num_subj_lst, num_cond, option = NULL) {

  # Check inputs
  if (missing(datamat_lst) || !is.list(datamat_lst)) {
    stop("datamat_lst should be a list of data matrices")
  }

  if (missing(num_subj_lst) || (!is.numeric(num_subj_lst) && !is.list(num_subj_lst))) {
    stop("num_subj_lst should be a numeric vector or list")
  }

  if (missing(num_cond) || !is.numeric(num_cond) || length(num_cond) != 1) {
    stop("num_cond should be a single number")
  }

  # Default options
  k <- num_cond
  method <- 1
  num_perm <- 0
  is_struct <- FALSE
  num_split <- 0
  num_boot <- 0
  clim <- 95
  bscan <- seq_len(num_cond)
  stacked_designdata <- NULL
  stacked_behavdata <- NULL
  meancentering_type <- 0
  cormode <- 0
  boot_type <- "strat"
  nonrotated_boot <- FALSE
  verbose <- TRUE
  robust_method <- "none"
  robust_trim <- 0.1
  robust_beta <- 0.2
  impute_method <- "none"
  check_missing <- TRUE
  missing_warn_threshold <- 0.05

  # Parse options
  if (!is.null(option)) {
    if (!is.list(option)) {
      stop("option argument should be a list")
    }

    if (!is.null(option$method)) {
      method <- option$method
      if (!method %in% 1:6) {
        stop("method should be 1, 2, 3, 4, 5, or 6")
      }
    }

    if (!is.null(option$num_perm)) {
      num_perm <- option$num_perm
      if (num_perm < 0 || num_perm != floor(num_perm)) {
        stop("num_perm should be a non-negative integer")
      }
    }

    if (!is.null(option$is_struct)) {
      is_struct <- option$is_struct
    }

    if (!is.null(option$num_split)) {
      num_split <- option$num_split
      if (num_split > 0) {
        nonrotated_boot <- TRUE
      }
    }

    if (!is.null(option$num_boot)) {
      num_boot <- option$num_boot
      if (num_boot < 0 || num_boot != floor(num_boot)) {
        stop("num_boot should be a non-negative integer")
      }
    }

    if (!is.null(option$clim)) {
      clim <- option$clim
      if (clim < 0 || clim > 100) {
        stop("clim should be between 0 and 100")
      }
    }

    if (!is.null(option$bscan)) {
      bscan <- option$bscan
      if (any(!bscan %in% seq_len(num_cond))) {
        stop("bscan should be a subset of 1:num_cond")
      }
    }

    if (!method %in% c(4, 6)) {
      bscan <- seq_len(num_cond)
    }

    if (!is.null(option$stacked_designdata)) {
      stacked_designdata <- as.matrix(option$stacked_designdata)
    }

    if (!is.null(option$stacked_behavdata)) {
      stacked_behavdata <- as.matrix(option$stacked_behavdata)
    }

    if (!is.null(option$meancentering_type)) {
      meancentering_type <- option$meancentering_type
      if (!meancentering_type %in% 0:3) {
        stop("meancentering_type should be 0, 1, 2, or 3")
      }
    }

    if (!is.null(option$cormode)) {
      cormode <- option$cormode
      if (!cormode %in% c(0, 2, 4, 6)) {
        stop("cormode should be 0, 2, 4, or 6")
      }
    }

    if (!is.null(option$boot_type)) {
      boot_type <- option$boot_type
      if (!boot_type %in% c("strat", "nonstrat")) {
        stop("boot_type should be 'strat' or 'nonstrat'")
      }
    }

    if (!is.null(option$verbose)) {
      verbose <- option$verbose
    }

    # Robust correlation options
    if (!is.null(option$robust_method)) {
      robust_method <- tolower(option$robust_method)
      if (!robust_method %in% c("none", "pearson", "spearman", "winsorized", "biweight", "percentage_bend")) {
        stop("robust_method should be 'none', 'pearson', 'spearman', 'winsorized', 'biweight', or 'percentage_bend'")
      }
    }

    if (!is.null(option$robust_trim)) {
      robust_trim <- option$robust_trim
      if (robust_trim <= 0 || robust_trim >= 0.5) {
        stop("robust_trim should be between 0 and 0.5 (exclusive)")
      }
    }

    if (!is.null(option$robust_beta)) {
      robust_beta <- option$robust_beta
      if (robust_beta <= 0 || robust_beta >= 0.5) {
        stop("robust_beta should be between 0 and 0.5 (exclusive)")
      }
    }

    # Missing data options
    if (!is.null(option$impute_method)) {
      impute_method <- tolower(option$impute_method)
      if (!impute_method %in% c("none", "mean", "rf")) {
        stop("impute_method should be 'none', 'mean', or 'rf'")
      }
    }

    if (!is.null(option$check_missing)) {
      check_missing <- option$check_missing
    }

    if (!is.null(option$missing_warn_threshold)) {
      missing_warn_threshold <- option$missing_warn_threshold
      if (missing_warn_threshold <= 0 || missing_warn_threshold >= 1) {
        stop("missing_warn_threshold should be between 0 and 1 (exclusive)")
      }
    }
  }

  # Initialize
  num_groups <- length(datamat_lst)
  total_rows <- 0

  for (g in seq_len(num_groups)) {
    total_rows <- total_rows + nrow(datamat_lst[[g]])
  }

  # Validate behavior data
  if (method %in% c(3, 4, 5, 6)) {
    if (is.null(stacked_behavdata)) {
      stop("stacked_behavdata is required for behavior PLS methods")
    }
    if (nrow(stacked_behavdata) != total_rows) {
      stop(sprintf("Wrong number of rows in behavior data, should be %d", total_rows))
    }
  }

  # Check for missing data and handle imputation
  has_missing_data <- any(sapply(datamat_lst, function(x) any(is.na(x))))
  has_missing_behav <- !is.null(stacked_behavdata) && any(is.na(stacked_behavdata))
  has_any_missing <- has_missing_data || has_missing_behav

  if (has_any_missing) {
    # Check and report missing data if requested
    if (check_missing && verbose) {
      missing_report <- pls_check_missing(
        datamat_lst,
        behavdata = stacked_behavdata,
        warn_threshold = missing_warn_threshold,
        verbose = TRUE
      )
    }

    if (impute_method == "none") {
      stop("Missing data (NA) detected in input matrices. ",
           "Either:\n",
           "  1. Remove or impute missing values before calling pls_analysis()\n",
           "  2. Set option$impute_method = 'mean' or 'rf' for automatic imputation\n",
           "  3. Use pls_analysis_mi() for multiple imputation with proper pooling\n",
           "Use pls_check_missing() to diagnose the extent of missing data.")
    }

    # Perform single imputation
    if (verbose) {
      message(sprintf("Imputing missing data using %s method...",
                      ifelse(impute_method == "rf", "random forest", "mean")))
    }

    imputed <- pls_impute(
      datamat_lst,
      behavdata = stacked_behavdata,
      method = impute_method,
      m = 1,
      verbose = FALSE
    )

    datamat_lst <- imputed$datamat_lst
    if (!is.null(stacked_behavdata)) {
      stacked_behavdata <- imputed$behavdata
    }

    if (verbose) {
      message("Imputation complete.")
    }
  }

  # Handle single condition case
  single_cond_lst <- NULL
  if (method %in% c(1, 2, 4, 6) && k == 1) {
    meancentering_type <- 1
    if (verbose) {
      message("Single condition detected, setting meancentering_type to 1")
    }
  }

  # Normalize design data for non-rotated PLS
  if (method %in% c(2, 5, 6)) {
    if (is.null(stacked_designdata)) {
      stop("stacked_designdata is required for non-rotated PLS methods")
    }
    stacked_designdata <- normalize(stacked_designdata)

    # Check rank
    if (qr(stacked_designdata)$rank != ncol(stacked_designdata)) {
      warning("Design contrast matrix is rank deficient")
    }
  }

  # Initialize result
  result <- list()
  result$method <- method
  result$is_struct <- is_struct
  result$robust_method <- robust_method

  # Stack data matrices
  if (verbose) {
    message("Stacking data matrices...")
  }
  stacked_datamat <- stacking_datamat(datamat_lst, single_cond_lst, verbose)

  # Calculate covariance/correlation
  if (verbose) {
    message("Calculating covariance/correlation data...")
  }

  datamat_reorder <- seq_len(nrow(stacked_datamat))
  behavdata_reorder <- if (method %in% c(3, 4, 5, 6)) seq_len(nrow(stacked_behavdata)) else NULL

  covcor_result <- rri_get_covcor(
    method, stacked_datamat, stacked_behavdata,
    num_groups, num_subj_lst, num_cond, bscan,
    meancentering_type, cormode, single_cond_lst,
    TRUE, num_boot, datamat_reorder, behavdata_reorder, NULL,
    robust_method = robust_method, trim = robust_trim, beta = robust_beta
  )

  datamatsvd <- covcor_result$datamatsvd
  datamatsvd_unnorm <- covcor_result$datamatsvd_unnorm
  datamatcorrs_lst <- covcor_result$datamatcorrs_lst
  stacked_smeanmat <- covcor_result$stacked_smeanmat

  # Save correlation maps for behavior PLS
  if (method %in% 3:6) {
    result$datamatcorrs_lst <- datamatcorrs_lst
  }

  # Calculate LVs
  if (verbose) {
    message("Calculating latent variables...")
  }

  if (method %in% c(2, 5, 6)) {
    # Non-rotated PLS
    crossblock <- t(stacked_designdata) %*% datamatsvd

    if (nonrotated_boot) {
      u <- normalize(t(crossblock))
    } else {
      u <- t(crossblock)
    }

    s <- sqrt(rowSums(crossblock^2))
    v <- stacked_designdata

    result$lvintercorrs <- t(normalize(u)) %*% normalize(u)

  } else {
    # SVD approach
    if (nrow(datamatsvd) <= ncol(datamatsvd)) {
      svd_result <- svd(t(datamatsvd))
      u <- svd_result$u
      s <- svd_result$d
      v <- svd_result$v
    } else {
      svd_result <- svd(datamatsvd)
      v <- svd_result$u
      s <- svd_result$d
      u <- svd_result$v
    }
  }

  org_s <- s
  org_v <- v
  original_u <- u %*% diag(s, nrow = length(s))
  original_v <- v %*% diag(s, nrow = length(s))

  kk <- length(bscan)

  # Adjust for multiblock PLS
  if (method %in% c(4, 6)) {
    result$bscan <- bscan

    total_s <- sum(datamatsvd_unnorm^2)
    per <- s^2 / sum(s^2)
    org_s <- sqrt(per * total_s)
    org_v <- v %*% diag(org_s, nrow = length(org_s))
  }

  # Save u, s, v
  result$u <- u
  result$s <- s
  result$v <- v

  # Calculate scores
  if (verbose) {
    message("Calculating scores...")
  }

  vsc <- NULL

  if (method %in% c(1, 2)) {
    if (method == 1) {
      usc <- stacked_datamat %*% u
      if (num_boot > 0) {
        usc2 <- stacked_smeanmat %*% u
      }
    } else {
      usc <- stacked_datamat %*% normalize(u)
      if (num_boot > 0) {
        usc2 <- stacked_smeanmat %*% u
      }
    }

    # Expand v to match subject structure
    for (g in seq_len(num_groups)) {
      if (!is.list(num_subj_lst)) {
        n <- num_subj_lst[g]
        v_g <- v[((g - 1) * k + 1):((g - 1) * k + k), , drop = FALSE]
        tmp <- v_g[rep(seq_len(k), each = n), , drop = FALSE]
      } else {
        n <- num_subj_lst[[g]]
        v_g <- v[((g - 1) * k + 1):((g - 1) * k + k), , drop = FALSE]
        tmp <- NULL
        for (k1 in seq_len(k)) {
          tmp <- rbind(tmp, v_g[rep(k1, n[k1]), , drop = FALSE])
        }
      }
      vsc <- rbind(vsc, tmp)
    }

  } else if (method %in% c(3, 5)) {
    # Behavior PLS scores
    if (!is.list(num_subj_lst)) {
      scores_result <- rri_get_behavscores(
        stacked_datamat, stacked_behavdata, u, v, k, num_subj_lst, cormode
      )
    } else {
      scores_result <- ssb_rri_get_behavscores(
        stacked_datamat, stacked_behavdata, u, v, k, num_subj_lst, cormode
      )
    }

    usc <- scores_result$scores
    vsc <- scores_result$fscores
    result$lvcorrs <- scores_result$lvcorrs

  } else if (method %in% c(4, 6)) {
    # Multiblock PLS scores
    if (num_boot > 0) {
      usc2 <- stacked_smeanmat %*% u
    }

    # Separate v into Task and Behavior parts
    t_behav <- ncol(stacked_behavdata)
    Tv <- NULL
    Bv <- NULL

    for (g in seq_len(num_groups)) {
      offset <- (g - 1) * k + (g - 1) * kk * t_behav
      Tv <- rbind(Tv, v[(offset + 1):(offset + k), , drop = FALSE])
      Bv <- rbind(Bv, v[(offset + k + 1):(offset + k + kk * t_behav), , drop = FALSE])
    }

    # Task scores
    Tusc <- stacked_datamat %*% normalize(u)
    Tvsc <- NULL

    row_idx <- NULL
    for (g in seq_len(num_groups)) {
      if (!is.list(num_subj_lst)) {
        n <- num_subj_lst[g]
        Tv_g <- Tv[((g - 1) * k + 1):((g - 1) * k + k), , drop = FALSE]
        tmp <- Tv_g[rep(seq_len(k), each = n), , drop = FALSE]
        Tvsc <- rbind(Tvsc, tmp)

        # Build row_idx for behavior
        offset <- if (g == 1) 0 else sum(num_subj_lst[1:(g-1)]) * k
        tmp_idx <- matrix(seq_len(n * k), nrow = n, ncol = k)
        row_idx <- c(row_idx, as.vector(tmp_idx[, bscan]) + offset)
      } else {
        n <- num_subj_lst[[g]]
        Tv_g <- Tv[((g - 1) * k + 1):((g - 1) * k + k), , drop = FALSE]
        tmp <- NULL
        for (k1 in seq_len(k)) {
          tmp <- rbind(tmp, Tv_g[rep(k1, n[k1]), , drop = FALSE])
        }
        Tvsc <- rbind(Tvsc, tmp)
      }
    }

    # Behavior scores
    if (!is.list(num_subj_lst)) {
      behav_result <- rri_get_behavscores(
        stacked_datamat[row_idx, , drop = FALSE],
        stacked_behavdata[row_idx, , drop = FALSE],
        u, Bv, kk, num_subj_lst, cormode
      )
    } else {
      num_subj_lst_4beh <- lapply(num_subj_lst, function(x) x[bscan])
      behav_result <- ssb_rri_get_behavscores(
        stacked_datamat[row_idx, , drop = FALSE],
        stacked_behavdata[row_idx, , drop = FALSE],
        u, Bv, kk, num_subj_lst_4beh, cormode
      )
    }

    Busc <- behav_result$scores
    Bvsc <- behav_result$fscores

    usc <- rbind(Tusc, Busc)
    vsc <- rbind(Tvsc, Bvsc)

    result$TBusc <- list(Tusc, Busc)
    result$TBvsc <- list(Tvsc, Bvsc)
    result$TBv <- list(Tv, Bv)
    result$lvcorrs <- behav_result$lvcorrs
  }

  result$usc <- usc
  result$vsc <- vsc
  result$num_subj_lst <- num_subj_lst
  result$num_conditions <- k

  if (method %in% c(3, 4, 5, 6)) {
    result$stacked_behavdata <- stacked_behavdata
  }
  if (method %in% c(2, 5, 6)) {
    result$stacked_designdata <- stacked_designdata
  }

  # ==================== PERMUTATION TEST ====================
  if (num_perm > 0 && num_split == 0) {
    if (verbose) {
      message(sprintf("Running permutation test (%d permutations)...", num_perm))
    }

    sp <- rep(0, length(s))

    # Generate permutation orders
    if (method %in% c(4, 6)) {
      if (!is.list(num_subj_lst)) {
        Treorder <- rri_perm_order(num_subj_lst, k, num_perm, is_struct)
        reorder <- matrix(0, nrow = total_rows, ncol = num_perm)
        for (p in seq_len(num_perm)) {
          reorder[, p] <- rri_randperm_notall(num_subj_lst, k, bscan)
        }
      } else {
        Treorder <- ssb_rri_perm_order(num_subj_lst, k, num_perm, is_struct)
        reorder <- Treorder
      }
    } else if (method %in% c(3, 5)) {
      reorder <- matrix(0, nrow = total_rows, ncol = num_perm)
      for (p in seq_len(num_perm)) {
        reorder[, p] <- sample(total_rows)
      }
    } else {
      if (!is.list(num_subj_lst)) {
        reorder <- rri_perm_order(num_subj_lst, k, num_perm, is_struct)
      } else {
        reorder <- ssb_rri_perm_order(num_subj_lst, k, num_perm, is_struct)
      }
    }

    # Run permutations
    for (p in seq_len(num_perm)) {
      if (verbose && p %% 100 == 0) {
        message(sprintf("  Permutation %d of %d", p, num_perm))
      }

      if (method %in% c(4, 6)) {
        datamat_reorder <- Treorder[, p]
        behavdata_reorder <- reorder[, p]
      } else if (method %in% c(3, 5)) {
        datamat_reorder <- seq_len(total_rows)
        behavdata_reorder <- reorder[, p]
      } else {
        datamat_reorder <- reorder[, p]
        behavdata_reorder <- NULL
      }

      perm_covcor <- rri_get_covcor(
        method, stacked_datamat, stacked_behavdata,
        num_groups, num_subj_lst, num_cond, bscan,
        meancentering_type, cormode, single_cond_lst,
        FALSE, 0, datamat_reorder, behavdata_reorder, NULL,
        robust_method = robust_method, trim = robust_trim, beta = robust_beta
      )

      datamatsvd_p <- perm_covcor$datamatsvd
      datamatsvd_unnorm_p <- perm_covcor$datamatsvd_unnorm

      if (method %in% c(2, 5, 6)) {
        crossblock_p <- t(normalize(stacked_designdata)) %*% datamatsvd_p
        sperm <- sqrt(rowSums(crossblock_p^2))
        sp <- sp + as.numeric(sperm >= s)
      } else {
        if (nrow(datamatsvd_p) <= ncol(datamatsvd_p)) {
          svd_p <- svd(t(datamatsvd_p))
          pu <- svd_p$u
          sperm <- svd_p$d
          pv <- svd_p$v
        } else {
          svd_p <- svd(datamatsvd_p)
          pv <- svd_p$u
          sperm <- svd_p$d
          pu <- svd_p$v
        }

        # Procrustes rotation
        rotatemat <- rri_bootprocrust(v, pv)
        pv <- pv %*% diag(sperm, nrow = length(sperm)) %*% rotatemat
        sperm <- sqrt(colSums(pv^2))

        if (method %in% c(4, 6)) {
          ptotal_s <- sum(datamatsvd_unnorm_p^2)
          per <- sperm^2 / sum(sperm^2)
          sperm <- sqrt(per * ptotal_s)
          sp <- sp + as.numeric(sperm >= org_s)
        } else {
          sp <- sp + as.numeric(sperm >= s)
        }
      }
    }

    result$perm_result <- list(
      num_perm = num_perm,
      sp = sp,
      sprob = sp / (num_perm + 1),
      permsamp = reorder,
      is_perm_splithalf = FALSE
    )

    if (method %in% c(4, 6)) {
      result$perm_result$Tpermsamp <- Treorder
    }
  }

  # ==================== BOOTSTRAP TEST ====================
  if (num_boot > 0) {
    if (verbose) {
      message(sprintf("Running bootstrap test (%d samples)...", num_boot))
    }

    incl_seq <- method == 1 || nonrotated_boot

    if (!is.list(num_subj_lst)) {
      boot_result <- rri_boot_order(num_subj_lst, k, num_boot, seq_len(k),
                                    incl_seq, boot_type)
      reorder <- boot_result$boot_order
      new_num_boot <- boot_result$new_num_boot

      if (method %in% c(4, 6)) {
        boot_result_4beh <- rri_boot_order(num_subj_lst, k, num_boot, bscan,
                                           incl_seq, boot_type)
        reorder_4beh <- boot_result_4beh$boot_order
        if (boot_result_4beh$new_num_boot < new_num_boot) {
          new_num_boot <- boot_result_4beh$new_num_boot
        }
      } else if (method %in% c(3, 5)) {
        reorder_4beh <- reorder
      }
    } else {
      # SSB version
      boot_result <- rri_boot_order(unlist(num_subj_lst), k, num_boot,
                                    seq_len(k), incl_seq, boot_type)
      reorder <- boot_result$boot_order
      new_num_boot <- boot_result$new_num_boot
      reorder_4beh <- reorder
    }

    num_boot <- new_num_boot

    # Initialize distributions
    if (method %in% c(3, 4, 5, 6)) {
      orig_corr <- result$lvcorrs
      distrib <- array(0, dim = c(nrow(orig_corr), ncol(orig_corr), num_boot + 1))
      distrib[, , 1] <- orig_corr
    }

    if (method %in% c(1, 2, 4, 6)) {
      # Calculate orig_usc
      orig_usc <- NULL
      first <- 1
      last <- 0

      for (g in seq_len(num_groups)) {
        if (!is.list(num_subj_lst)) {
          last <- last + k * num_subj_lst[g]
          if (method %in% c(2, 6)) {
            orig_usc <- rbind(orig_usc, rri_task_mean(usc[first:last, , drop = FALSE],
                                                      num_subj_lst[g]))
          } else {
            orig_usc <- rbind(orig_usc, rri_task_mean(usc2[first:last, , drop = FALSE],
                                                      num_subj_lst[g]))
          }
        } else {
          last <- last + sum(num_subj_lst[[g]])
          if (method %in% c(2, 6)) {
            orig_usc <- rbind(orig_usc, ssb_rri_task_mean(usc[first:last, , drop = FALSE],
                                                         num_subj_lst[[g]]))
          } else {
            orig_usc <- rbind(orig_usc, ssb_rri_task_mean(usc2[first:last, , drop = FALSE],
                                                         num_subj_lst[[g]]))
          }
        }
        first <- last + 1
      }

      if (method %in% c(4, 6)) {
        Tdistrib <- array(0, dim = c(nrow(orig_usc), ncol(orig_usc), num_boot + 1))
        Tdistrib[, , 1] <- orig_usc
      } else {
        distrib <- array(0, dim = c(nrow(orig_usc), ncol(orig_usc), num_boot + 1))
        distrib[, , 1] <- orig_usc
      }
    }

    # Initialize accumulators
    if (method == 1 || nonrotated_boot) {
      u_sum <- matrix(0, nrow = nrow(u), ncol = ncol(u))
    } else if (method %in% c(2, 5, 6)) {
      u_sum <- u
    } else {
      u_sum <- original_u
    }
    u_sq <- u_sum^2

    # Bootstrap loop
    for (p in seq_len(num_boot)) {
      if (verbose && p %% 100 == 0) {
        message(sprintf("  Bootstrap %d of %d", p, num_boot))
      }

      datamat_reorder <- reorder[, p]

      if (method %in% c(3, 4, 5, 6)) {
        datamat_reorder_4beh <- reorder_4beh[, p]
        behavdata_reorder <- reorder_4beh[, p]
      } else {
        datamat_reorder_4beh <- NULL
        behavdata_reorder <- NULL
      }

      boot_covcor <- rri_get_covcor(
        method, stacked_datamat, stacked_behavdata,
        num_groups, num_subj_lst, num_cond, bscan,
        meancentering_type, cormode, single_cond_lst,
        TRUE, num_boot, datamat_reorder, behavdata_reorder, datamat_reorder_4beh,
        robust_method = robust_method, trim = robust_trim, beta = robust_beta
      )

      datamatsvd_b <- boot_covcor$datamatsvd
      stacked_smeanmat_b <- boot_covcor$stacked_smeanmat

      if (method %in% c(2, 5, 6) && !nonrotated_boot) {
        crossblock_b <- t(normalize(stacked_designdata)) %*% datamatsvd_b
        u_sq <- u_sq + t(crossblock_b)^2
        u_sum <- u_sum + t(crossblock_b)

        if (method %in% c(5, 6)) {
          # Record distrib for behavior
        }

        if (method %in% c(2, 6)) {
          tmp_usc <- stacked_datamat %*% normalize(t(crossblock_b))
          tmp_orig_usc <- NULL
          first <- 1
          last <- 0

          for (g in seq_len(num_groups)) {
            if (!is.list(num_subj_lst)) {
              last <- last + k * num_subj_lst[g]
              tmp_orig_usc <- rbind(tmp_orig_usc,
                                    rri_task_mean(tmp_usc[first:last, , drop = FALSE],
                                                  num_subj_lst[g]))
            } else {
              last <- last + sum(num_subj_lst[[g]])
              tmp_orig_usc <- rbind(tmp_orig_usc,
                                    ssb_rri_task_mean(tmp_usc[first:last, , drop = FALSE],
                                                      num_subj_lst[[g]]))
            }
            first <- last + 1
          }

          if (method == 6) {
            Tdistrib[, , p + 1] <- tmp_orig_usc
          } else {
            distrib[, , p + 1] <- tmp_orig_usc
          }
        }

      } else {
        # SVD approach
        if (nrow(datamatsvd_b) <= ncol(datamatsvd_b)) {
          svd_b <- svd(t(datamatsvd_b))
          pu <- svd_b$u
          sboot <- svd_b$d
          pv <- svd_b$v
        } else {
          svd_b <- svd(datamatsvd_b)
          pv <- svd_b$u
          sboot <- svd_b$d
          pu <- svd_b$v
        }

        rotatemat <- rri_bootprocrust(v, pv)
        pu <- pu %*% diag(sboot, nrow = length(sboot)) %*% rotatemat
        pv <- pv %*% diag(sboot, nrow = length(sboot)) %*% rotatemat

        u_sum <- u_sum + pu
        u_sq <- u_sq + pu^2

        if (method == 1) {
          tmp_usc2 <- stacked_smeanmat_b %*% normalize(pu)
          tmp_orig_usc <- NULL
          first <- 1
          last <- 0

          for (g in seq_len(num_groups)) {
            if (!is.list(num_subj_lst)) {
              last <- last + k * num_subj_lst[g]
              tmp_orig_usc <- rbind(tmp_orig_usc,
                                    rri_task_mean(tmp_usc2[first:last, , drop = FALSE],
                                                  num_subj_lst[g]))
            } else {
              last <- last + sum(num_subj_lst[[g]])
              tmp_orig_usc <- rbind(tmp_orig_usc,
                                    ssb_rri_task_mean(tmp_usc2[first:last, , drop = FALSE],
                                                      num_subj_lst[[g]]))
            }
            first <- last + 1
          }

          distrib[, , p + 1] <- tmp_orig_usc
        }
      }
    }

    # Calculate standard errors
    if (method == 1 || nonrotated_boot) {
      u_sum2 <- u_sum^2 / num_boot
      u_se <- sqrt(abs(u_sq - u_sum2) / (num_boot - 1))
    } else {
      u_sum2 <- u_sum^2 / (num_boot + 1)
      u_se <- sqrt(abs(u_sq - u_sum2) / num_boot)
    }

    # Calculate confidence intervals
    ul <- clim
    ll <- 100 - clim
    climNi <- 0.5 * (1 - (clim * 0.01))

    if (method %in% c(1, 2)) {
      ci_result <- rri_distrib(distrib, ll, ul, num_boot, climNi, orig_usc)
      llusc <- ci_result$llcorr
      ulusc <- ci_result$ulcorr
      llusc_adj <- ci_result$llcorr_adj
      ulusc_adj <- ci_result$ulcorr_adj
      prop <- ci_result$prop
    } else if (method %in% c(3, 5)) {
      ci_result <- rri_distrib(distrib, ll, ul, num_boot, climNi, orig_corr)
      llcorr <- ci_result$llcorr
      ulcorr <- ci_result$ulcorr
      llcorr_adj <- ci_result$llcorr_adj
      ulcorr_adj <- ci_result$ulcorr_adj
      prop <- ci_result$prop
    }

    # Calculate bootstrap ratios
    test_zeros <- which(u_se <= 0)
    if (length(test_zeros) > 0) {
      u_se[test_zeros] <- 1
    }

    if (method %in% c(2, 5, 6)) {
      if (nonrotated_boot) {
        compare_u <- original_u / u_se
      } else {
        compare_u <- u / u_se
      }
    } else {
      compare_u <- original_u / u_se
    }

    if (length(test_zeros) > 0) {
      compare_u[test_zeros] <- 0
    }

    # Save bootstrap results
    result$boot_result <- list(
      num_boot = num_boot,
      clim = clim,
      boot_type = boot_type,
      nonrotated_boot = nonrotated_boot,
      bootsamp = reorder,
      compare_u = compare_u,
      u_se = u_se,
      zero_u_se = test_zeros
    )

    if (method %in% c(1, 2)) {
      result$boot_result$usc2 <- if (exists("usc2")) usc2 else NULL
      result$boot_result$orig_usc <- orig_usc
      result$boot_result$ulusc <- ulusc
      result$boot_result$llusc <- llusc
      result$boot_result$ulusc_adj <- ulusc_adj
      result$boot_result$llusc_adj <- llusc_adj
      result$boot_result$prop <- prop
      result$boot_result$distrib <- distrib
    }

    if (method %in% c(3, 4, 5, 6)) {
      result$boot_result$bootsamp_4beh <- reorder_4beh
      result$boot_result$orig_corr <- orig_corr
      result$boot_result$ulcorr <- if (exists("ulcorr")) ulcorr else NULL
      result$boot_result$llcorr <- if (exists("llcorr")) llcorr else NULL
      result$boot_result$ulcorr_adj <- if (exists("ulcorr_adj")) ulcorr_adj else NULL
      result$boot_result$llcorr_adj <- if (exists("llcorr_adj")) llcorr_adj else NULL
      result$boot_result$distrib <- distrib
      result$boot_result$prop <- if (exists("prop")) prop else NULL
    }

    if (method %in% c(4, 6)) {
      result$boot_result$Tdistrib <- if (exists("Tdistrib")) Tdistrib else NULL
    }
  }

  # Save other inputs
  result$other_input <- list(
    meancentering_type = meancentering_type,
    cormode = cormode
  )

  # Save imputation info if used
  if (exists("has_any_missing") && has_any_missing && impute_method != "none") {
    result$imputation_info <- list(
      was_imputed = TRUE,
      method = impute_method
    )
  }

  class(result) <- "pls_result"

  if (verbose) {
    message("PLS analysis complete.")
  }

  return(result)
}


#' Print PLS Result Summary
#'
#' @param x A pls_result object
#' @param ... Additional arguments (ignored)
#'
#' @export
print.pls_result <- function(x, ...) {
  cat("PLS Analysis Result\n")
  cat("==================\n\n")

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

  if (!is.null(x$perm_result)) {
    cat("Permutation Test Results:\n")
    cat("  Number of permutations:", x$perm_result$num_perm, "\n")
    cat("  p-values:", round(x$perm_result$sprob, 4), "\n\n")
  }

  if (!is.null(x$boot_result)) {
    cat("Bootstrap Results:\n")
    cat("  Number of bootstraps:", x$boot_result$num_boot, "\n")
    cat("  Confidence level:", x$boot_result$clim, "%\n")
  }

  invisible(x)
}
