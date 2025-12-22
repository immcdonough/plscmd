#' Calculate Behavioral PLS Scores
#'
#' Computes brain scores, behavior scores, and LV correlations from PLS analysis.
#'
#' @param stacked_datamat Stacked data matrix
#' @param stacked_behavdata Stacked behavioral data
#' @param brainlv Brain latent variable (salience)
#' @param behavlv Behavior latent variable
#' @param k Number of conditions
#' @param num_subj_lst Vector of number of subjects per group
#' @param cormode Correlation mode (see rri_xcor)
#'
#' @return List containing:
#'   \itemize{
#'     \item scores: Brain scores
#'     \item fscores: Behavior scores
#'     \item lvcorrs: LV correlations (orig_corr)
#'   }
#'
#' @export
rri_get_behavscores <- function(stacked_datamat, stacked_behavdata,
                                brainlv, behavlv, k, num_subj_lst,
                                cormode = 0) {

  if (!is.matrix(stacked_datamat)) {
    stacked_datamat <- as.matrix(stacked_datamat)
  }
  if (!is.matrix(stacked_behavdata)) {
    stacked_behavdata <- as.matrix(stacked_behavdata)
  }
  if (!is.matrix(brainlv)) {
    brainlv <- as.matrix(brainlv)
  }
  if (!is.matrix(behavlv)) {
    behavlv <- as.matrix(behavlv)
  }

  # Brain scores
  scores <- stacked_datamat %*% brainlv

  fscores <- NULL
  lvcorrs <- NULL

  num_groups <- length(num_subj_lst)

  for (g in seq_len(num_groups)) {
    n <- num_subj_lst[g]
    t <- ncol(stacked_behavdata)  # Number of behavior measures
    tmp <- NULL

    span <- if (g == 1) 0 else sum(num_subj_lst[1:(g-1)]) * k

    for (i in seq_len(k)) {
      # Calculate behavior scores for this condition
      behav_idx <- (1 + n * (i - 1) + span):(n * i + span)
      behavlv_idx <- (1 + t * (i - 1) + (g - 1) * t * k):(t * i + (g - 1) * t * k)

      tmp_k <- stacked_behavdata[behav_idx, , drop = FALSE] %*%
        behavlv[behavlv_idx, , drop = FALSE]
      tmp <- rbind(tmp, tmp_k)
    }

    fscores <- rbind(fscores, tmp)

    # Calculate LV correlations
    behav_range <- (1 + span):(n * k + span)
    tmp <- rri_corr_maps(
      stacked_behavdata[behav_range, , drop = FALSE],
      scores[behav_range, , drop = FALSE],
      n, k, cormode
    )
    lvcorrs <- rbind(lvcorrs, tmp)
  }

  return(list(
    scores = scores,
    fscores = fscores,
    lvcorrs = lvcorrs
  ))
}


#' Calculate Behavioral PLS Scores for Split-Subject-Blocks
#'
#' Version of rri_get_behavscores for variable sample sizes per condition.
#'
#' @param stacked_datamat Stacked data matrix
#' @param stacked_behavdata Stacked behavioral data
#' @param brainlv Brain latent variable
#' @param behavlv Behavior latent variable
#' @param k Number of conditions
#' @param num_subj_lst List of vectors with number of subjects per condition
#' @param cormode Correlation mode
#'
#' @return List with scores, fscores, and lvcorrs
#'
#' @export
ssb_rri_get_behavscores <- function(stacked_datamat, stacked_behavdata,
                                    brainlv, behavlv, k, num_subj_lst,
                                    cormode = 0) {

  if (!is.list(num_subj_lst)) {
    return(rri_get_behavscores(stacked_datamat, stacked_behavdata,
                               brainlv, behavlv, k, num_subj_lst, cormode))
  }

  if (!is.matrix(stacked_datamat)) {
    stacked_datamat <- as.matrix(stacked_datamat)
  }
  if (!is.matrix(stacked_behavdata)) {
    stacked_behavdata <- as.matrix(stacked_behavdata)
  }
  if (!is.matrix(brainlv)) {
    brainlv <- as.matrix(brainlv)
  }
  if (!is.matrix(behavlv)) {
    behavlv <- as.matrix(behavlv)
  }

  # Brain scores
  scores <- stacked_datamat %*% brainlv

  fscores <- NULL
  lvcorrs <- NULL

  num_groups <- length(num_subj_lst)

  for (g in seq_len(num_groups)) {
    n <- num_subj_lst[[g]]
    t <- ncol(stacked_behavdata)
    tmp <- NULL

    span <- if (g == 1) 0 else sum(unlist(num_subj_lst[1:(g-1)]))

    step <- 0
    for (i in seq_len(k)) {
      behav_idx <- (1 + step + span):(n[i] + step + span)
      behavlv_idx <- (1 + t * (i - 1) + (g - 1) * t * k):(t * i + (g - 1) * t * k)

      tmp_k <- stacked_behavdata[behav_idx, , drop = FALSE] %*%
        behavlv[behavlv_idx, , drop = FALSE]
      tmp <- rbind(tmp, tmp_k)
      step <- step + n[i]
    }

    fscores <- rbind(fscores, tmp)

    # Calculate LV correlations using SSB version
    behav_range <- (1 + span):(sum(n) + span)
    tmp <- ssb_rri_corr_maps(
      stacked_behavdata[behav_range, , drop = FALSE],
      scores[behav_range, , drop = FALSE],
      n, k, cormode
    )
    lvcorrs <- rbind(lvcorrs, tmp)
  }

  return(list(
    scores = scores,
    fscores = fscores,
    lvcorrs = lvcorrs
  ))
}


#' Correlation Maps for Split-Subject-Blocks
#'
#' Creates correlation maps with variable sample sizes per condition.
#'
#' @param behav Behavioral data
#' @param datamat Data matrix
#' @param n Vector of subjects per condition
#' @param k Number of conditions
#' @param cormode Correlation mode
#'
#' @return Stacked correlation maps
#'
#' @keywords internal
ssb_rri_corr_maps <- function(behav, datamat, n, k, cormode = 0) {
  if (!is.matrix(behav)) behav <- as.matrix(behav)
  if (!is.matrix(datamat)) datamat <- as.matrix(datamat)

  maps <- NULL
  step <- 0

  for (i in seq_len(k)) {
    start_row <- 1 + step
    end_row <- step + n[i]

    behav_subset <- behav[start_row:end_row, , drop = FALSE]
    data_subset <- datamat[start_row:end_row, , drop = FALSE]

    temp <- rri_xcor(behav_subset, data_subset, cormode)
    maps <- rbind(maps, temp)

    step <- step + n[i]
  }

  return(maps)
}


#' Correlation Maps for Split-Subject-Blocks (Selected Conditions)
#'
#' Creates correlation maps for selected conditions with variable sample sizes.
#'
#' @param behav Behavioral data
#' @param datamat Data matrix
#' @param n Vector of subjects per condition
#' @param bscan Conditions to include
#' @param cormode Correlation mode
#'
#' @return Correlation maps
#'
#' @keywords internal
ssb_rri_corr_maps_notall <- function(behav, datamat, n, bscan, cormode = 0) {
  if (!is.matrix(behav)) behav <- as.matrix(behav)
  if (!is.matrix(datamat)) datamat <- as.matrix(datamat)

  maps <- NULL
  step <- 0

  for (i in seq_along(n)) {
    if (i %in% bscan) {
      start_row <- 1 + step
      end_row <- step + n[i]

      behav_subset <- behav[start_row:end_row, , drop = FALSE]
      data_subset <- datamat[start_row:end_row, , drop = FALSE]

      temp <- rri_xcor(behav_subset, data_subset, cormode)
      maps <- rbind(maps, temp)
    }

    step <- step + n[i]
  }

  return(maps)
}
