#' Generate Permutation Order
#'
#' Generates permutation orders for significance testing. The permutation
#' is done in two steps: first permuting conditions within each subject,
#' then permuting subjects across groups.
#'
#' @param num_subj_lst Vector of number of subjects per group
#' @param num_cond Number of conditions
#' @param num_perm Number of permutations to generate
#' @param not_in_cond Logical, if TRUE do not permute conditions within subjects
#'   (for structure PLS)
#'
#' @return Matrix of permutation orders (total_rows x num_perm)
#'
#' @details
#' The permutation procedure:
#' 1. For each subject, randomly permute the condition labels
#' 2. Randomly permute subjects across groups
#' 3. Ensure the permuted average is different from the original
#' 4. Ensure no duplicate permutations
#'
#' @export
rri_perm_order <- function(num_subj_lst, num_cond, num_perm, not_in_cond = FALSE) {
  num_subj_grp <- sum(num_subj_lst)
  total_rows <- num_subj_grp * num_cond
  num_group <- length(num_subj_lst)

  perm_order <- matrix(0, nrow = total_rows, ncol = num_perm)

  # Set random seed
  set.seed(NULL)

  for (p in seq_len(num_perm)) {
    cnt <- 0
    duplicated <- TRUE

    while (duplicated) {
      cnt <- cnt + 1

      # Create task_group matrix
      first <- 1
      last <- 0
      task_group <- NULL

      for (g in seq_along(num_subj_lst)) {
        last <- last + num_cond * num_subj_lst[g]
        tmp <- matrix(seq(first, last), nrow = num_subj_lst[g], ncol = num_cond)
        task_group <- cbind(task_group, t(tmp))
        first <- last + 1
      }

      origin_task_group <- task_group

      # Permute conditions within each subject (unless not_in_cond)
      if (!not_in_cond) {
        for (i in seq_len(num_subj_grp)) {
          task_perm <- sample(num_cond)
          task_group[, i] <- task_group[task_perm, i]
        }
      }

      # Permute subjects across groups
      group_perm <- sample(num_subj_grp)
      task_group <- task_group[, group_perm, drop = FALSE]

      # Check that average is different from original
      duplicated <- FALSE

      for (c in seq_len(num_cond)) {
        accum <- 0
        for (g in seq_along(num_subj_lst)) {
          cols <- (accum + 1):(accum + num_subj_lst[g])
          if (all(sort(task_group[c, cols]) == origin_task_group[c, cols])) {
            duplicated <- TRUE
          }
          accum <- accum + num_subj_lst[g]
        }
      }

      # Convert to order vector
      new_perm_order <- NULL
      for (g in seq_along(num_subj_lst)) {
        if (g == 1) {
          cols <- 1:num_subj_lst[g]
        } else {
          cols <- (sum(num_subj_lst[1:(g-1)]) + 1):sum(num_subj_lst[1:g])
        }
        tmp <- task_group[, cols, drop = FALSE]
        tmp <- as.vector(t(tmp))
        new_perm_order <- c(new_perm_order, tmp)
      }

      # Check for duplicate permutations
      if (p > 1) {
        for (i in seq_len(p - 1)) {
          if (all(perm_order[, i] == new_perm_order)) {
            duplicated <- TRUE
            break
          }
        }
      }

      # Check for sequential order
      if (all(new_perm_order == seq_len(total_rows))) {
        duplicated <- TRUE
      }

      if (cnt > 500) {
        duplicated <- FALSE
        warning("Duplicated permutation orders are used!")
      }
    }

    perm_order[, p] <- new_perm_order
  }

  return(perm_order)
}


#' Generate Permutation Order for Split-Subject-Blocks
#'
#' Version of rri_perm_order for variable sample sizes per condition.
#'
#' @param num_subj_lst List of vectors, each containing number of subjects
#'   per condition for a group
#' @param num_cond Number of conditions
#' @param num_perm Number of permutations
#' @param not_in_cond Do not permute conditions
#'
#' @return Matrix of permutation orders
#'
#' @keywords internal
ssb_rri_perm_order <- function(num_subj_lst, num_cond, num_perm, not_in_cond = FALSE) {
  # For split-subject-blocks, num_subj_lst is a list of vectors
  if (!is.list(num_subj_lst)) {
    return(rri_perm_order(num_subj_lst, num_cond, num_perm, not_in_cond))
  }

  total_rows <- sum(unlist(num_subj_lst))
  num_group <- length(num_subj_lst)

  perm_order <- matrix(0, nrow = total_rows, ncol = num_perm)

  for (p in seq_len(num_perm)) {
    cnt <- 0
    duplicated <- TRUE

    while (duplicated) {
      cnt <- cnt + 1

      new_perm_order <- NULL
      offset <- 0

      for (g in seq_len(num_group)) {
        n <- num_subj_lst[[g]]
        group_size <- sum(n)

        # Create index for this group
        group_idx <- seq(offset + 1, offset + group_size)

        if (!not_in_cond) {
          # Permute within conditions
          for (k in seq_len(num_cond)) {
            if (k == 1) {
              start_idx <- 1
            } else {
              start_idx <- sum(n[1:(k-1)]) + 1
            }
            end_idx <- sum(n[1:k])

            cond_idx <- group_idx[start_idx:end_idx]
            cond_idx <- sample(cond_idx)
            new_perm_order <- c(new_perm_order, cond_idx)
          }
        } else {
          new_perm_order <- c(new_perm_order, sample(group_idx))
        }

        offset <- offset + group_size
      }

      # Check for sequential order
      duplicated <- all(new_perm_order == seq_len(total_rows))

      # Check for duplicates
      if (!duplicated && p > 1) {
        for (i in seq_len(p - 1)) {
          if (all(perm_order[, i] == new_perm_order)) {
            duplicated <- TRUE
            break
          }
        }
      }

      if (cnt > 500) {
        duplicated <- FALSE
        warning("Duplicated permutation orders are used!")
      }
    }

    perm_order[, p] <- new_perm_order
  }

  return(perm_order)
}


#' Generate Permutation Order Including Original
#'
#' Same as rri_perm_order but includes the original (unpermuted) order
#' as the first column.
#'
#' @param num_subj_lst Vector of number of subjects per group
#' @param num_cond Number of conditions
#' @param num_perm Number of permutations
#' @param not_in_cond Do not permute conditions
#'
#' @return Matrix of permutation orders with original as first column
#'
#' @keywords internal
missnk_rri_perm_order <- function(num_subj_lst, num_cond, num_perm, not_in_cond = FALSE) {
  num_subj_grp <- sum(num_subj_lst)
  total_rows <- num_subj_grp * num_cond

  # First column is original order
  perm_order <- matrix(0, nrow = total_rows, ncol = num_perm)
  perm_order[, 1] <- seq_len(total_rows)

  if (num_perm > 1) {
    remaining <- rri_perm_order(num_subj_lst, num_cond, num_perm - 1, not_in_cond)
    perm_order[, 2:num_perm] <- remaining
  }

  return(perm_order)
}


#' Random Permutation Excluding Certain Conditions
#'
#' Generates a random permutation that only includes specified conditions.
#'
#' @param num_subj_lst Vector of number of subjects per group
#' @param num_cond Number of conditions
#' @param bscan Vector of conditions to include
#'
#' @return Vector of permuted indices
#'
#' @keywords internal
rri_randperm_notall <- function(num_subj_lst, num_cond, bscan) {
  total_subj <- sum(num_subj_lst)
  total_rows <- total_subj * num_cond

  # Create mask for selected conditions
  mask <- rep(0, total_rows)

  for (g in seq_along(num_subj_lst)) {
    n <- num_subj_lst[g]
    if (g == 1) {
      offset <- 0
    } else {
      offset <- sum(num_subj_lst[1:(g-1)]) * num_cond
    }

    for (cond in bscan) {
      start_idx <- offset + n * (cond - 1) + 1
      end_idx <- offset + n * cond
      mask[start_idx:end_idx] <- 1
    }
  }

  selected_idx <- which(mask == 1)

  # Create result with original order
  result <- seq_len(total_rows)

  # Permute only selected indices
  result[selected_idx] <- sample(selected_idx)

  return(result)
}


#' Split-Half Permutation Test
#'
#' Performs split-half reliability testing with permutation.
#'
#' @param sop Count of singular values exceeding original
#' @param ucorr_distrib Distribution of u correlations
#' @param vcorr_distrib Distribution of v correlations
#' @param num_perm Number of outer permutations
#' @param num_split Number of split-half iterations
#' @param num_lvs Number of latent variables
#' @param num_groups Number of groups
#' @param num_cond Number of conditions
#' @param num_subj_lst Number of subjects per group
#' @param num_subj_lst1 First half subject counts
#' @param num_subj_lst2 Second half subject counts
#' @param rows1 Row indices for first half
#' @param rows2 Row indices for second half
#' @param meancentering_type Mean-centering type
#' @param cormode Correlation mode
#' @param single_cond_lst Single condition list (if applicable)
#' @param method PLS method
#' @param s Singular values
#' @param org_s Original singular values
#' @param org_v Original v matrix
#' @param bscan Conditions for behavior block
#' @param stacked_datamat Stacked data matrix
#' @param opt Options list
#'
#' @return List with updated sop and correlation distributions
#'
#' @export
splithalf_perm <- function(sop, ucorr_distrib, vcorr_distrib,
                           num_perm, num_split, num_lvs,
                           num_groups, num_cond, num_subj_lst,
                           num_subj_lst1, num_subj_lst2,
                           rows1, rows2,
                           meancentering_type, cormode, single_cond_lst,
                           method, s, org_s, org_v, bscan,
                           stacked_datamat, opt) {

  # Extract options
  stacked_behavdata <- opt$stacked_behavdata
  stacked_designdata <- opt$stacked_designdata

  if (method %in% c(4, 6)) {
    Treorder <- if (!is.null(opt$Treorder)) opt$Treorder else opt$reorder
    Breorder <- if (!is.null(opt$Breorder)) opt$Breorder else opt$reorder
  } else {
    reorder <- opt$reorder
  }

  # Outer permutation loop
  for (op in seq_len(num_perm)) {
    message(sprintf("Outer permutation: %d of %d", op, num_perm))

    # Set up reorder indices for this permutation
    if (method %in% c(4, 6)) {
      datamat_reorder <- Treorder[, op]
      behavdata_reorder <- seq_len(nrow(stacked_behavdata))
      datamat_reorder_4beh <- Breorder[, op]
    } else if (method %in% c(3, 5)) {
      datamat_reorder <- reorder[, op]
      behavdata_reorder <- seq_len(nrow(stacked_behavdata))
      datamat_reorder_4beh <- NULL
    } else {
      datamat_reorder <- reorder[, op]
      behavdata_reorder <- NULL
      datamat_reorder_4beh <- NULL
    }

    # Get permuted covariance matrix
    covcor_result <- rri_get_covcor(
      method, stacked_datamat, stacked_behavdata,
      num_groups, num_subj_lst, num_cond, bscan,
      meancentering_type, cormode, single_cond_lst,
      FALSE, 0, datamat_reorder, behavdata_reorder, datamat_reorder_4beh
    )

    datamatsvd_op <- covcor_result$datamatsvd
    datamatsvd_unnorm_op <- covcor_result$datamatsvd_unnorm

    # Compute LVs for this permutation
    if (method %in% c(2, 5, 6)) {
      # Non-rotated PLS
      v_op <- stacked_designdata
      crossblock <- t(normalize(stacked_designdata)) %*% datamatsvd_op
      u_op <- t(crossblock)
      sperm <- sqrt(rowSums(crossblock^2))
    } else {
      # SVD
      if (nrow(datamatsvd_op) <= ncol(datamatsvd_op)) {
        svd_result <- misssvd(t(datamatsvd_op), 0)
        u_op <- svd_result$u
        sperm <- diag(svd_result$d)
        v_op <- svd_result$v
      } else {
        svd_result <- misssvd(datamatsvd_op, 0)
        v_op <- svd_result$u
        sperm <- diag(svd_result$d)
        u_op <- svd_result$v
      }
    }

    # Update sop count (skip first = original)
    if (op > 1) {
      if (method %in% c(4, 6)) {
        ptotal_s <- sum(datamatsvd_unnorm_op^2)
        per <- sperm^2 / sum(sperm^2)
        sperm <- sqrt(per * ptotal_s)
      }

      sop <- sop + as.numeric(sperm >= org_s)
    }

    # Inner split-half loop
    for (p in seq_len(num_split)) {
      # Create split-half by randomly permuting subjects within each group
      in_reorder <- NULL

      for (g in seq_len(num_groups)) {
        gperm <- sample(num_subj_lst[g])
        offset <- if (g == 1) 0 else sum(num_subj_lst[1:(g-1)]) * num_cond

        for (cond in seq_len(num_cond)) {
          in_reorder <- c(in_reorder, offset + num_subj_lst[g] * (cond - 1) + gperm)
        }
      }

      # Calculate datamatsvd for each half
      if (method %in% c(4, 6)) {
        datamat_reorder1 <- datamat_reorder[in_reorder[rows1]]
        datamat_reorder2 <- datamat_reorder[in_reorder[rows2]]
        behavdata_reorder1 <- behavdata_reorder[in_reorder[rows1]]
        behavdata_reorder2 <- behavdata_reorder[in_reorder[rows2]]
        datamat_reorder_4beh1 <- datamat_reorder_4beh[in_reorder[rows1]]
        datamat_reorder_4beh2 <- datamat_reorder_4beh[in_reorder[rows2]]
      } else if (method %in% c(3, 5)) {
        datamat_reorder1 <- datamat_reorder[in_reorder[rows1]]
        datamat_reorder2 <- datamat_reorder[in_reorder[rows2]]
        behavdata_reorder1 <- behavdata_reorder[in_reorder[rows1]]
        behavdata_reorder2 <- behavdata_reorder[in_reorder[rows2]]
        datamat_reorder_4beh1 <- NULL
        datamat_reorder_4beh2 <- NULL
      } else {
        datamat_reorder1 <- datamat_reorder[in_reorder[rows1]]
        datamat_reorder2 <- datamat_reorder[in_reorder[rows2]]
        behavdata_reorder1 <- NULL
        behavdata_reorder2 <- NULL
        datamat_reorder_4beh1 <- NULL
        datamat_reorder_4beh2 <- NULL
      }

      # Get covariance for each half
      covcor1 <- rri_get_covcor(
        method, stacked_datamat, stacked_behavdata,
        num_groups, num_subj_lst1, num_cond, bscan,
        meancentering_type, cormode, single_cond_lst,
        FALSE, 0, datamat_reorder1, behavdata_reorder1, datamat_reorder_4beh1
      )

      covcor2 <- rri_get_covcor(
        method, stacked_datamat, stacked_behavdata,
        num_groups, num_subj_lst2, num_cond, bscan,
        meancentering_type, cormode, single_cond_lst,
        FALSE, 0, datamat_reorder2, behavdata_reorder2, datamat_reorder_4beh2
      )

      datamatsvd_p1 <- covcor1$datamatsvd
      datamatsvd_p2 <- covcor2$datamatsvd

      # Project design onto both halves
      u_p1 <- t(datamatsvd_p1) %*% v_op  # brain factor score
      u_p2 <- t(datamatsvd_p2) %*% v_op

      # Project brain pattern onto both halves
      v_p1 <- datamatsvd_p1 %*% u_op  # effect factor score
      v_p2 <- datamatsvd_p2 %*% u_op

      # Correlate and update mean correlations
      for (lv in seq_len(num_lvs)) {
        ucorr_distrib[op, lv] <- ucorr_distrib[op, lv] +
          rri_xcor(u_p1[, lv, drop = FALSE], u_p2[, lv, drop = FALSE], cormode)
        vcorr_distrib[op, lv] <- vcorr_distrib[op, lv] +
          rri_xcor(v_p1[, lv, drop = FALSE], v_p2[, lv, drop = FALSE], cormode)
      }
    }
  }

  return(list(
    sop = sop,
    ucorr_distrib = ucorr_distrib,
    vcorr_distrib = vcorr_distrib
  ))
}
