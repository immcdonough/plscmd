#' Generate Bootstrap Sample Order
#'
#' Generate bootstrap resampling orders for all conditions and groups.
#' Sampling is done with replacement within each group.
#'
#' @param num_subj_lst Vector of number of subjects per group
#' @param num_cond Number of conditions
#' @param num_boot Number of bootstrap samples to generate
#' @param bscan Vector of conditions to include (default: all)
#' @param incl_seq Logical, include sequential (original) order as first sample
#' @param boot_type "strat" for stratified (default) or "nonstrat" for non-stratified
#'
#' @return List with:
#'   \itemize{
#'     \item boot_order: Matrix of bootstrap orders (rows x num_boot)
#'     \item new_num_boot: Actual number of bootstrap samples (may be reduced)
#'   }
#'
#' @details
#' For stratified bootstrap (default), subjects are resampled within each group
#' separately, maintaining the group structure. For non-stratified bootstrap,
#' all subjects are pooled and resampled together.
#'
#' Bootstrap requires at least 3 subjects per group for reliable estimates.
#'
#' @export
rri_boot_order <- function(num_subj_lst, num_cond, num_boot,
                           bscan = NULL, incl_seq = FALSE,
                           boot_type = "strat") {

  if (is.null(bscan)) {
    bscan <- seq_len(num_cond)
  }

  if (!boot_type %in% c("strat", "nonstrat")) {
    stop("boot_type must be 'strat' or 'nonstrat'")
  }

  total_subj <- sum(num_subj_lst)
  num_group <- length(num_subj_lst)
  num_cond0 <- num_cond
  num_cond_bscan <- length(bscan)
  total_rows0 <- num_cond0 * total_subj
  total_rows <- num_cond_bscan * total_subj

  new_num_boot <- num_boot

  if (boot_type == "nonstrat") {
    # Non-stratified bootstrap: treat all subjects as one group
    from_nonstrat <- list(
      total_rows0 = total_rows0,
      num_group = num_group,
      num_subj_lst = num_subj_lst,
      num_cond0 = num_cond0
    )

    # Create order as if all subjects are in one group
    result <- rri_boot_order(total_subj, num_cond, num_boot,
                             bscan, FALSE, "strat")
    tmp_order <- result$boot_order
    new_num_boot <- result$new_num_boot

    # Remap into original subj*cond*group ordering
    boot_order <- matrix(0, nrow = nrow(tmp_order), ncol = ncol(tmp_order))

    for (r in seq_len(nrow(tmp_order))) {
      for (p in seq_len(ncol(tmp_order))) {
        row <- tmp_order[r, p]
        ss <- ((row - 1) %% total_subj) + 1  # Subject index when all in one group

        # Find which group this subject belongs to
        g <- 1
        cumsum_subj <- cumsum(num_subj_lst)
        while (g <= num_group && ss > cumsum_subj[g]) {
          g <- g + 1
        }

        # Subject index within group
        if (g == 1) {
          s <- ss
        } else {
          s <- ss - cumsum_subj[g - 1]
        }

        # Find condition
        cond <- ceiling(row / total_subj)

        # Map to row of original datamat
        if (g == 1) {
          row1 <- num_subj_lst[g] * (cond - 1) + s
        } else {
          row1 <- sum(num_subj_lst[1:(g-1)]) * num_cond0 +
            num_subj_lst[g] * (cond - 1) + s
        }

        boot_order[r, p] <- row1
      }
    }

    return(list(boot_order = boot_order, new_num_boot = new_num_boot))
  }

  # Check bootstrap validity
  check_result <- rri_boot_check(num_subj_lst, num_cond, num_boot, incl_seq)
  min_subj_per_group <- check_result$min_subj_per_group
  is_boot_samples <- check_result$is_boot_samples
  boot_samples <- check_result$boot_samples
  new_num_boot <- check_result$new_num_boot

  # Generate bootstrap orders for each sample
  tmp_boot_order <- matrix(0, nrow = total_subj, ncol = new_num_boot)

  for (p in seq_len(new_num_boot)) {
    subj_order <- vector("list", num_group)
    not_done <- TRUE
    cnt <- 0

    while (not_done) {
      start_subj <- 1

      for (g in seq_len(num_group)) {
        num_subj <- num_subj_lst[g]
        all_samples_are_same <- TRUE

        while (all_samples_are_same) {
          if (is_boot_samples[g]) {
            # Get from pre-generated boot_samples
            new_subj_order <- boot_samples[[g]][p, ]
            all_samples_are_same <- FALSE
            not_done <- FALSE
          } else {
            not_done <- TRUE
            # Random resampling with replacement
            new_subj_order <- floor(runif(num_subj) * num_subj) + 1

            # Check for minimum diversity
            n_unique <- length(unique(new_subj_order))
            if (n_unique >= min_subj_per_group) {
              all_samples_are_same <- FALSE
            }
          }
        }

        subj_order[[g]] <- new_subj_order + start_subj - 1
        start_subj <- start_subj + num_subj
      }

      if (!all(is_boot_samples)) {
        # Check for duplicate orders
        not_done <- FALSE
        if (p > 1) {
          combined_order <- unlist(subj_order)
          for (i in seq_len(p - 1)) {
            if (all(tmp_boot_order[, i] == combined_order)) {
              not_done <- TRUE
              break
            }
          }
        }

        # Treat sequential order as duplicate
        if (!incl_seq && all(unlist(subj_order) == seq_len(total_subj))) {
          not_done <- TRUE
        }

        cnt <- cnt + 1
        if (cnt > 500) {
          not_done <- FALSE
          warning("Duplicated bootstrap orders are used!")
        }
      }
    }

    tmp_boot_order[, p] <- unlist(subj_order)
  }

  # Construct the resampling order matrix
  first <- 1
  last <- 0
  row_idx <- NULL

  for (g in seq_along(num_subj_lst)) {
    last <- last + num_cond * num_subj_lst[g]
    tmp <- matrix(seq(first, last), nrow = num_subj_lst[g], ncol = num_cond)
    row_idx <- cbind(row_idx, t(tmp))
    first <- last + 1
  }

  boot_order <- array(0, dim = c(num_cond, total_subj, new_num_boot))

  for (p in seq_len(new_num_boot)) {
    boot_order[, , p] <- row_idx[, tmp_boot_order[, p]]
  }

  # Reshape to final format
  b_order <- NULL
  for (g in seq_along(num_subj_lst)) {
    one_group <- NULL
    for (p in seq_len(new_num_boot)) {
      if (g == 1) {
        cols <- 1:num_subj_lst[g]
      } else {
        cols <- (sum(num_subj_lst[1:(g-1)]) + 1):sum(num_subj_lst[1:g])
      }
      tmp <- boot_order[, cols, p]
      tmp <- as.vector(t(tmp))
      one_group <- cbind(one_group, tmp)
    }
    b_order <- rbind(b_order, one_group)
  }

  boot_order <- b_order

  # Create template for bscan
  template <- rep(0, total_rows0)

  for (g in seq_len(num_group)) {
    n <- num_subj_lst[g]
    for (i in bscan) {
      if (g == 1) {
        offset <- 0
      } else {
        offset <- sum(num_cond0 * num_subj_lst[1:(g-1)])
      }
      start_idx <- offset + 1 + (n * (i - 1))
      end_idx <- offset + n * i
      template[start_idx:end_idx] <- 1
    }
  }

  template_idx <- which(template == 1)
  boot_order <- matrix(template_idx[boot_order], nrow = length(template_idx))

  # Add back non-selected conditions with original order
  boot_order0 <- matrix(rep(seq_len(total_rows0), new_num_boot),
                        nrow = total_rows0, ncol = new_num_boot)
  boot_order0[template_idx, ] <- boot_order
  boot_order <- boot_order0

  return(list(boot_order = boot_order, new_num_boot = new_num_boot))
}


#' Check Bootstrap Validity
#'
#' Validates bootstrap parameters and generates sample matrices for small groups.
#'
#' @param num_subj_lst Vector of number of subjects per group
#' @param num_cond Number of conditions
#' @param num_boot Requested number of bootstrap samples
#' @param incl_seq Include sequential order
#'
#' @return List with validation results and pre-generated samples
#'
#' @export
rri_boot_check <- function(num_subj_lst, num_cond, num_boot, incl_seq = FALSE) {
  total_subj <- sum(num_subj_lst)
  num_group <- length(num_subj_lst)

  new_num_boot <- num_boot

  # Minimum 3 subjects per group required
  if (min(num_subj_lst) < 3) {
    stop("Bootstrap analysis requires at least 3 subjects per group")
  }

  percentage <- 50
  min_subj_per_group <- ceiling(min(num_subj_lst) * percentage / 100)
  max_subj_per_group <- 8

  is_boot_samples <- rep(FALSE, num_group)
  boot_samples <- vector("list", num_group)

  # For small groups, pre-generate all possible bootstrap samples
  if (all(num_subj_lst <= max_subj_per_group)) {
    for (g in seq_len(num_group)) {
      num_subj <- num_subj_lst[g]
      boot_sample2 <- NULL

      if (incl_seq) {
        max_diff_subj <- num_subj
      } else {
        max_diff_subj <- num_subj - 1
      }

      for (diff_subj in min_subj_per_group:max_diff_subj) {
        boot_sample1 <- rri_boot_samples(num_subj, diff_subj)
        boot_sample2 <- rbind(boot_sample2, boot_sample1)
      }

      num_boot_samples <- nrow(boot_sample2)

      if (!is.null(boot_sample2) && num_boot_samples < new_num_boot) {
        warning(sprintf("%d subjects can only have %d different bootstrap samples. Setting num_boot to this maximum.",
                        num_subj, num_boot_samples))
        new_num_boot <- num_boot_samples
      }

      if (!is.null(boot_sample2)) {
        # Randomly shuffle
        boot_sample2 <- boot_sample2[sample(nrow(boot_sample2)), , drop = FALSE]
        boot_samples[[g]] <- boot_sample2
        is_boot_samples[g] <- TRUE
      }
    }
  }

  return(list(
    min_subj_per_group = min_subj_per_group,
    is_boot_samples = is_boot_samples,
    boot_samples = boot_samples,
    new_num_boot = new_num_boot
  ))
}


#' Generate All Bootstrap Samples with Given Diversity
#'
#' Generates all theoretically possible bootstrap samples with a minimum
#' number of unique subjects.
#'
#' @param n Number of subjects
#' @param diff_subjs Minimum number of unique subjects
#'
#' @return Matrix of bootstrap samples
#'
#' @export
rri_boot_samples <- function(n, diff_subjs) {
  B <- NULL
  bootstrap_num <- 0
  k <- n

  A <- rep(1, n)

  while (A[1] < n) {
    array_count <- rep(0, n)

    for (i in seq_len(n)) {
      array_count[A[i]] <- 1
    }

    sum_unique <- sum(array_count)

    if (sum_unique == diff_subjs) {
      B <- rbind(B, A)
      bootstrap_num <- bootstrap_num + 1
    }

    if (A[k] == n) {
      i <- k
      while (i > 0 && A[i] == n) {
        i <- i - 1
      }

      if (i > 0) {
        A[i] <- A[i] + 1
        for (j in (i + 1):n) {
          A[j] <- A[i]
        }
      } else {
        break
      }
    } else {
      A[k] <- A[k] + 1
    }
  }

  if (diff_subjs == 1) {
    B <- rbind(B, rep(n, n))
    bootstrap_num <- bootstrap_num + 1
  }

  return(B)
}


#' Calculate Bootstrap Confidence Intervals
#'
#' Calculates upper and lower confidence interval limits from bootstrap distribution.
#'
#' @param distrib 3D array of bootstrap distributions (rows x cols x samples)
#' @param ll Lower percentile (e.g., 5 for 95% CI)
#' @param ul Upper percentile (e.g., 95 for 95% CI)
#' @param num_boot Number of bootstrap samples
#' @param ClimNi Normalized confidence interval parameter
#' @param orig Original observed values
#'
#' @return List with:
#'   \itemize{
#'     \item llcorr: Lower CI limit
#'     \item ulcorr: Upper CI limit
#'     \item prop: Proportion of bootstrap samples <= original
#'     \item llcorr_adj: Adjusted lower CI limit (bias-corrected)
#'     \item ulcorr_adj: Adjusted upper CI limit (bias-corrected)
#'   }
#'
#' @export
rri_distrib <- function(distrib, ll, ul, num_boot, ClimNi, orig) {
  if (length(dim(distrib)) != 3) {
    stop("distrib must be a 3D array")
  }

  nrow_d <- dim(distrib)[1]
  ncol_d <- dim(distrib)[2]

  llcorr <- matrix(0, nrow = nrow_d, ncol = ncol_d)
  ulcorr <- matrix(0, nrow = nrow_d, ncol = ncol_d)
  prop <- matrix(0, nrow = nrow_d, ncol = ncol_d)
  llcorr_adj <- matrix(NA, nrow = nrow_d, ncol = ncol_d)
  ulcorr_adj <- matrix(NA, nrow = nrow_d, ncol = ncol_d)

  for (r in seq_len(nrow_d)) {
    for (c in seq_len(ncol_d)) {
      # Get bootstrap distribution (excluding first observation)
      boot_vals <- distrib[r, c, 2:(num_boot + 1)]

      ulcorr[r, c] <- percentile(boot_vals, ul)
      llcorr[r, c] <- percentile(boot_vals, ll)
      prop[r, c] <- sum(boot_vals <= orig[r, c]) / num_boot

      if (prop[r, c] == 1 || prop[r, c] == 0) {
        llcorr_adj[r, c] <- NA
        ulcorr_adj[r, c] <- NA
      } else {
        # Adjusted confidence intervals for skewed distributions
        ni <- cumulative_gaussian_inv(prop[r, c])

        # Recalculate confidence interval limits
        uli <- (2 * ni) + cumulative_gaussian_inv(1 - ClimNi)
        lli <- (2 * ni) + cumulative_gaussian_inv(ClimNi)

        ncdf_lli <- cumulative_gaussian(lli) * 100
        ncdf_uli <- cumulative_gaussian(uli) * 100

        llcorr_adj[r, c] <- percentile(boot_vals, ncdf_lli)
        ulcorr_adj[r, c] <- percentile(boot_vals, ncdf_uli)
      }
    }
  }

  return(list(
    llcorr = llcorr,
    ulcorr = ulcorr,
    prop = prop,
    llcorr_adj = llcorr_adj,
    ulcorr_adj = ulcorr_adj
  ))
}
