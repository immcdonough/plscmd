#' Calculate Covariance/Correlation Data for PLS
#'
#' Prepares covariance or correlation data for SVD or crossblock calculation.
#' This is a core function that implements the various mean-centering strategies.
#'
#' @param method PLS method (1-6)
#' @param stacked_datamat Stacked data matrix from all groups
#' @param stacked_behavdata Stacked behavioral data (for behavior PLS)
#' @param num_groups Number of groups
#' @param num_subj_lst Vector of number of subjects per group
#' @param num_cond Number of conditions
#' @param bscan Conditions for behavior block in multiblock PLS
#' @param meancentering_type Mean-centering type:
#'   \itemize{
#'     \item 0 = Remove group condition means (boost condition differences)
#'     \item 1 = Remove grand condition means (boost group differences)
#'     \item 2 = Remove grand mean
#'     \item 3 = Remove all main effects (interaction only)
#'   }
#' @param cormode Correlation mode (0, 2, 4, or 6)
#' @param single_cond_lst Single condition list (special case)
#' @param is_observation Whether this is the observed (not resampled) data
#' @param num_boot Number of bootstrap samples
#' @param datamat_reorder Reorder indices for data matrix
#' @param behavdata_reorder Reorder indices for behavior data
#' @param datamat_reorder_4beh Reorder indices for multiblock behavior
#' @param robust_method Robust correlation method: "none", "spearman",
#'   "winsorized", "biweight", or "percentage_bend"
#' @param trim Trim proportion for winsorized correlation (default: 0.1)
#' @param beta Bend constant for percentage bend correlation (default: 0.2)
#'
#' @return List containing:
#'   \itemize{
#'     \item datamatsvd: Stacked covariance/correlation data for SVD
#'     \item datamatsvd_unnorm: Unnormalized version (for multiblock)
#'     \item datamatcorrs_lst: Correlation list (behavior PLS)
#'     \item stacked_smeanmat: Mean-centered data for bootstrap CI
#'   }
#'
#' @export
rri_get_covcor <- function(method, stacked_datamat, stacked_behavdata,
                           num_groups, num_subj_lst, num_cond, bscan,
                           meancentering_type, cormode, single_cond_lst = NULL,
                           is_observation = TRUE, num_boot = 0,
                           datamat_reorder = NULL, behavdata_reorder = NULL,
                           datamat_reorder_4beh = NULL,
                           robust_method = "none", trim = 0.1, beta = 0.2) {

  if (!is.matrix(stacked_datamat)) {
    stacked_datamat <- as.matrix(stacked_datamat)
  }

  if (is.null(datamat_reorder)) {
    datamat_reorder <- seq_len(nrow(stacked_datamat))
  }

  if (is.null(datamat_reorder_4beh)) {
    datamat_reorder_4beh <- seq_len(nrow(stacked_datamat))
  }

  # Initialize outputs
  datamatsvd <- NULL
  datamatsvd_unnorm <- NULL
  datamatcorrs_lst <- list()
  stacked_smeanmat <- NULL

  k <- num_cond

  # Handle different mean-centering types
  if (meancentering_type == 1) {
    # Calculate grand condition means across all groups
    grand_mean <- matrix(0, nrow = num_cond, ncol = ncol(stacked_datamat))

    for (g in seq_len(num_groups)) {
      if (!is.list(num_subj_lst)) {
        n <- num_subj_lst[g]
        span <- if (g == 1) 0 else sum(num_subj_lst[1:(g-1)]) * k
        datamat <- stacked_datamat[datamat_reorder, , drop = FALSE]
        datamat <- datamat[(1 + span):(n * k + span), , drop = FALSE]
        grand_mean <- grand_mean + rri_task_mean(datamat, n)
      } else {
        n <- num_subj_lst[[g]]
        span <- if (g == 1) 0 else sum(unlist(num_subj_lst[1:(g-1)]))
        datamat <- stacked_datamat[datamat_reorder, , drop = FALSE]
        datamat <- datamat[(1 + span):(sum(n) + span), , drop = FALSE]
        grand_mean <- grand_mean + ssb_rri_task_mean(datamat, n)
      }
    }

    grand_mean <- grand_mean / num_groups

  } else if (meancentering_type == 2) {
    # Calculate grand mean over all subjects and conditions
    datamat <- stacked_datamat[datamat_reorder, , drop = FALSE]
    grand_mean <- colMeans(datamat)

  } else if (meancentering_type == 3) {
    # Calculate condition and group means for interaction
    cond_group_mean <- array(0, dim = c(num_cond, num_groups, ncol(stacked_datamat)))

    for (g in seq_len(num_groups)) {
      if (!is.list(num_subj_lst)) {
        n <- num_subj_lst[g]
        span <- if (g == 1) 0 else sum(num_subj_lst[1:(g-1)]) * k
        datamat <- stacked_datamat[datamat_reorder, , drop = FALSE]
        datamat <- datamat[(1 + span):(n * k + span), , drop = FALSE]
        cond_group_mean[, g, ] <- rri_task_mean(datamat, n)
      } else {
        n <- num_subj_lst[[g]]
        span <- if (g == 1) 0 else sum(unlist(num_subj_lst[1:(g-1)]))
        datamat <- stacked_datamat[datamat_reorder, , drop = FALSE]
        datamat <- datamat[(1 + span):(sum(n) + span), , drop = FALSE]
        cond_group_mean[, g, ] <- ssb_rri_task_mean(datamat, n)
      }
    }

    cond_mean <- apply(cond_group_mean, c(2, 3), mean)  # Mean over conditions
    group_mean <- apply(cond_group_mean, c(1, 3), mean)  # Mean over groups
    grand_mean <- colMeans(matrix(cond_group_mean, ncol = ncol(stacked_datamat)))
  }

  # Loop across groups
  for (g in seq_len(num_groups)) {
    if (!is.list(num_subj_lst)) {
      n <- num_subj_lst[g]
      span <- if (g == 1) 0 else sum(num_subj_lst[1:(g-1)]) * k
      datamat <- stacked_datamat[datamat_reorder, , drop = FALSE]

      if (is.null(single_cond_lst)) {
        datamat <- datamat[(1 + span):(n * k + span), , drop = FALSE]
      }

      if (method %in% c(4, 6)) {
        datamat_4beh <- stacked_datamat[datamat_reorder_4beh, , drop = FALSE]
        datamat_4beh <- datamat_4beh[(1 + span):(n * k + span), , drop = FALSE]
      }

      if (method %in% c(3, 4, 5, 6)) {
        behavdata <- stacked_behavdata[behavdata_reorder, , drop = FALSE]
        behavdata <- behavdata[(1 + span):(n * k + span), , drop = FALSE]
      }

      # Calculate correlation or covariance based on method
      if (method == 1) {
        # Mean-Centering Task PLS
        if (!is.null(single_cond_lst) && g == 1) {
          datamat <- single_cond_lst[[1]]
          datamatcorrs <- ssb_rri_task_mean(datamat, num_subj_lst) -
            matrix(1, nrow = num_groups, ncol = 1) %*% matrix(colMeans(datamat), nrow = 1)

          if (is_observation && num_boot > 0) {
            smeanmat <- datamat - matrix(1, nrow = nrow(datamat), ncol = 1) %*%
              matrix(colMeans(datamat), nrow = 1)
          }
        } else if (is.null(single_cond_lst) && meancentering_type == 0) {
          datamatcorrs <- rri_task_mean(datamat, n) -
            matrix(1, nrow = k, ncol = 1) %*% matrix(colMeans(datamat), nrow = 1)

          if (is_observation && num_boot > 0) {
            smeanmat <- datamat - matrix(1, nrow = nrow(datamat), ncol = 1) %*%
              matrix(colMeans(datamat), nrow = 1)
          }
        } else if (is.null(single_cond_lst) && meancentering_type == 1) {
          datamatcorrs <- rri_task_mean(datamat, n) - grand_mean

          if (is_observation && num_boot > 0) {
            smeanmat <- NULL
            for (cond_idx in seq_len(num_cond)) {
              start_row <- (cond_idx - 1) * num_subj_lst[g] + 1
              end_row <- cond_idx * num_subj_lst[g]
              tmp <- datamat[start_row:end_row, , drop = FALSE] -
                matrix(1, nrow = num_subj_lst[g], ncol = 1) %*%
                matrix(grand_mean[cond_idx, ], nrow = 1)
              smeanmat <- rbind(smeanmat, tmp)
            }
          }
        } else if (is.null(single_cond_lst) && meancentering_type == 2) {
          datamatcorrs <- rri_task_mean(datamat, n) -
            matrix(1, nrow = k, ncol = 1) %*% matrix(grand_mean, nrow = 1)

          if (is_observation && num_boot > 0) {
            smeanmat <- datamat - matrix(1, nrow = nrow(datamat), ncol = 1) %*%
              matrix(grand_mean, nrow = 1)
          }
        } else if (is.null(single_cond_lst) && meancentering_type == 3) {
          datamatcorrs <- rri_task_mean(datamat, n) - group_mean -
            matrix(1, nrow = k, ncol = 1) %*% matrix(cond_mean[g, ], nrow = 1) +
            matrix(1, nrow = k, ncol = 1) %*% matrix(grand_mean, nrow = 1)

          if (is_observation && num_boot > 0) {
            smeanmat <- NULL
            for (cond_idx in seq_len(num_cond)) {
              start_row <- (cond_idx - 1) * num_subj_lst[g] + 1
              end_row <- cond_idx * num_subj_lst[g]
              tmp <- datamat[start_row:end_row, , drop = FALSE] -
                matrix(1, nrow = num_subj_lst[g], ncol = 1) %*%
                matrix(group_mean[cond_idx, ], nrow = 1) -
                matrix(1, nrow = num_subj_lst[g], ncol = 1) %*%
                matrix(cond_mean[g, ], nrow = 1) +
                matrix(1, nrow = num_subj_lst[g], ncol = 1) %*%
                matrix(grand_mean, nrow = 1)
              smeanmat <- rbind(smeanmat, tmp)
            }
          }
        }

        TBdatamatcorrs <- NULL

      } else if (method %in% c(2, 6)) {
        # Non-Rotated Task PLS
        datamatcorrs <- rri_task_mean(datamat, n)

        if (is_observation && num_boot > 0) {
          if (meancentering_type == 0) {
            smeanmat <- datamat - matrix(1, nrow = nrow(datamat), ncol = 1) %*%
              matrix(colMeans(datamat), nrow = 1)
          } else if (meancentering_type == 1) {
            smeanmat <- NULL
            for (cond_idx in seq_len(num_cond)) {
              start_row <- (cond_idx - 1) * num_subj_lst[g] + 1
              end_row <- cond_idx * num_subj_lst[g]
              tmp <- datamat[start_row:end_row, , drop = FALSE] -
                matrix(1, nrow = num_subj_lst[g], ncol = 1) %*%
                matrix(grand_mean[cond_idx, ], nrow = 1)
              smeanmat <- rbind(smeanmat, tmp)
            }
          } else if (meancentering_type == 2) {
            smeanmat <- datamat - matrix(1, nrow = nrow(datamat), ncol = 1) %*%
              matrix(grand_mean, nrow = 1)
          } else if (meancentering_type == 3) {
            smeanmat <- NULL
            for (cond_idx in seq_len(num_cond)) {
              start_row <- (cond_idx - 1) * num_subj_lst[g] + 1
              end_row <- cond_idx * num_subj_lst[g]
              tmp <- datamat[start_row:end_row, , drop = FALSE] -
                matrix(1, nrow = num_subj_lst[g], ncol = 1) %*%
                matrix(group_mean[cond_idx, ], nrow = 1) -
                matrix(1, nrow = num_subj_lst[g], ncol = 1) %*%
                matrix(cond_mean[g, ], nrow = 1) +
                matrix(1, nrow = num_subj_lst[g], ncol = 1) %*%
                matrix(grand_mean, nrow = 1)
              smeanmat <- rbind(smeanmat, tmp)
            }
          }
        }

        if (method == 6) {
          Tdatamatcorrs <- datamatcorrs
          Bdatamatcorrs <- rri_corr_maps_notall(behavdata, datamat_4beh, n, bscan, cormode,
                                                 robust_method = robust_method, trim = trim, beta = beta)
          datamatcorrs_lst <- c(datamatcorrs_lst, list(Bdatamatcorrs))

          TBdatamatcorrs <- rbind(Tdatamatcorrs, Bdatamatcorrs)
          datamatcorrs <- rbind(normalize(Tdatamatcorrs, 2), normalize(Bdatamatcorrs, 2))
        } else {
          TBdatamatcorrs <- NULL
        }

      } else if (method %in% c(3, 5)) {
        # Behavior PLS
        datamatcorrs <- rri_corr_maps(behavdata, datamat, n, k, cormode,
                                      robust_method = robust_method, trim = trim, beta = beta)
        datamatcorrs_lst <- c(datamatcorrs_lst, list(datamatcorrs))
        TBdatamatcorrs <- NULL

      } else if (method == 4) {
        # Multiblock PLS
        if (meancentering_type == 0) {
          Tdatamatcorrs <- rri_task_mean(datamat, n) -
            matrix(1, nrow = k, ncol = 1) %*% matrix(colMeans(datamat), nrow = 1)

          if (is_observation && num_boot > 0) {
            smeanmat <- datamat - matrix(1, nrow = nrow(datamat), ncol = 1) %*%
              matrix(colMeans(datamat), nrow = 1)
          }
        } else if (meancentering_type == 1) {
          Tdatamatcorrs <- rri_task_mean(datamat, n) - grand_mean

          if (is_observation && num_boot > 0) {
            smeanmat <- NULL
            for (cond_idx in seq_len(num_cond)) {
              start_row <- (cond_idx - 1) * num_subj_lst[g] + 1
              end_row <- cond_idx * num_subj_lst[g]
              tmp <- datamat[start_row:end_row, , drop = FALSE] -
                matrix(1, nrow = num_subj_lst[g], ncol = 1) %*%
                matrix(grand_mean[cond_idx, ], nrow = 1)
              smeanmat <- rbind(smeanmat, tmp)
            }
          }
        } else if (meancentering_type == 2) {
          Tdatamatcorrs <- rri_task_mean(datamat, n) -
            matrix(1, nrow = k, ncol = 1) %*% matrix(grand_mean, nrow = 1)

          if (is_observation && num_boot > 0) {
            smeanmat <- datamat - matrix(1, nrow = nrow(datamat), ncol = 1) %*%
              matrix(grand_mean, nrow = 1)
          }
        } else if (meancentering_type == 3) {
          Tdatamatcorrs <- rri_task_mean(datamat, n) - group_mean -
            matrix(1, nrow = k, ncol = 1) %*% matrix(cond_mean[g, ], nrow = 1) +
            matrix(1, nrow = k, ncol = 1) %*% matrix(grand_mean, nrow = 1)

          if (is_observation && num_boot > 0) {
            smeanmat <- NULL
            for (cond_idx in seq_len(num_cond)) {
              start_row <- (cond_idx - 1) * num_subj_lst[g] + 1
              end_row <- cond_idx * num_subj_lst[g]
              tmp <- datamat[start_row:end_row, , drop = FALSE] -
                matrix(1, nrow = num_subj_lst[g], ncol = 1) %*%
                matrix(group_mean[cond_idx, ], nrow = 1) -
                matrix(1, nrow = num_subj_lst[g], ncol = 1) %*%
                matrix(cond_mean[g, ], nrow = 1) +
                matrix(1, nrow = num_subj_lst[g], ncol = 1) %*%
                matrix(grand_mean, nrow = 1)
              smeanmat <- rbind(smeanmat, tmp)
            }
          }
        }

        Bdatamatcorrs <- rri_corr_maps_notall(behavdata, datamat_4beh, n, bscan, cormode,
                                               robust_method = robust_method, trim = trim, beta = beta)
        datamatcorrs_lst <- c(datamatcorrs_lst, list(Bdatamatcorrs))

        TBdatamatcorrs <- rbind(Tdatamatcorrs, Bdatamatcorrs)
        datamatcorrs <- rbind(normalize(Tdatamatcorrs, 2), normalize(Bdatamatcorrs, 2))
      }

      # Accumulate results
      if (is.null(single_cond_lst) || g == 1) {
        datamatsvd_unnorm <- rbind(datamatsvd_unnorm, TBdatamatcorrs)
        datamatsvd <- rbind(datamatsvd, datamatcorrs)

        if (is_observation && num_boot > 0 && method %in% c(1, 2, 4, 6)) {
          if (exists("smeanmat")) {
            stacked_smeanmat <- rbind(stacked_smeanmat, smeanmat)
          }
        }
      }

    } else {
      # Handle list-based num_subj_lst (split-subject-blocks)
      # Similar logic but with variable sample sizes per condition
      # Implementation follows the same pattern as above
      n <- num_subj_lst[[g]]
      span <- if (g == 1) 0 else sum(unlist(num_subj_lst[1:(g-1)]))
      datamat <- stacked_datamat[datamat_reorder, , drop = FALSE]
      datamat <- datamat[(1 + span):(sum(n) + span), , drop = FALSE]

      # Simplified handling for SSB
      if (method == 1) {
        if (meancentering_type == 0) {
          datamatcorrs <- ssb_rri_task_mean(datamat, n) -
            matrix(1, nrow = k, ncol = 1) %*% matrix(colMeans(datamat), nrow = 1)
        } else {
          datamatcorrs <- ssb_rri_task_mean(datamat, n)
        }
        TBdatamatcorrs <- NULL
      } else {
        datamatcorrs <- ssb_rri_task_mean(datamat, n)
        TBdatamatcorrs <- NULL
      }

      datamatsvd <- rbind(datamatsvd, datamatcorrs)
    }
  }

  return(list(
    datamatsvd = datamatsvd,
    datamatsvd_unnorm = datamatsvd_unnorm,
    datamatcorrs_lst = datamatcorrs_lst,
    stacked_smeanmat = stacked_smeanmat
  ))
}


#' Correlation Maps Excluding Some Conditions
#'
#' Creates correlation maps for selected conditions only.
#'
#' @param behav Behavioral data
#' @param datamat Data matrix
#' @param n Number of subjects
#' @param bscan Conditions to include
#' @param cormode Correlation mode
#' @param robust_method Robust correlation method (see rri_xcor)
#' @param trim Trim proportion for winsorized correlation
#' @param beta Bend constant for percentage bend correlation
#'
#' @return Correlation maps for selected conditions
#'
#' @keywords internal
rri_corr_maps_notall <- function(behav, datamat, n, bscan, cormode = 0,
                                  robust_method = "none", trim = 0.1, beta = 0.2) {
  if (!is.matrix(behav)) behav <- as.matrix(behav)
  if (!is.matrix(datamat)) datamat <- as.matrix(datamat)

  maps <- NULL

  for (i in bscan) {
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
