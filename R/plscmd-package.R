#' plscmd: Partial Least Squares Analysis with Permutation and Bootstrap Testing
#'
#' Comprehensive Partial Least Squares (PLS) analysis toolkit for neuroimaging
#' data. This package provides six different PLS methods with support for
#' permutation testing, bootstrap resampling, and split-half validation.
#'
#' @section PLS Methods:
#' The package supports six PLS analysis methods:
#' \describe{
#'   \item{Method 1: Mean-Centering Task PLS}{Standard task PLS that decomposes
#'     mean-centered condition averages using SVD.}
#'   \item{Method 2: Non-Rotated Task PLS}{Uses user-specified design contrasts
#'     instead of data-driven rotation.}
#'   \item{Method 3: Regular Behavior PLS}{Identifies brain patterns that
#'     correlate with behavioral measures.}
#'   \item{Method 4: Multiblock PLS}{Combines task and behavior blocks in a
#'     single analysis.}
#'   \item{Method 5: Non-Rotated Behavior PLS}{Behavior PLS with user-specified
#'     contrasts.}
#'   \item{Method 6: Non-Rotated Multiblock PLS}{Multiblock PLS with user-specified
#'     contrasts.}
#' }
#'
#' @section Statistical Inference:
#' \describe{
#'   \item{Permutation Testing}{Assesses statistical significance of latent
#'     variables by comparing observed singular values to a null distribution
#'     generated through permutation.}
#'   \item{Bootstrap Resampling}{Provides confidence intervals for brain
#'     saliences and LV correlations through resampling with replacement.}
#'   \item{Split-Half Validation}{Assesses reliability of latent variable
#'     structure by splitting subjects and correlating results.}
#' }
#'
#' @section Main Functions:
#' \describe{
#'   \item{\code{\link{pls_analysis}}}{Main function for running PLS analysis}
#'   \item{\code{\link{normalize}}}{Normalize vectors to unit length}
#'   \item{\code{\link{percentile}}}{Calculate percentiles}
#'   \item{\code{\link{rri_boot_order}}}{Generate bootstrap resampling orders}
#'   \item{\code{\link{rri_perm_order}}}{Generate permutation orders}
#' }
#'
#' @section Data Organization:
#' Data should be organized with subjects in rows (stacked by condition within
#' each group) and variables (e.g., voxels, electrodes) in columns:
#' \preformatted{
#'   Group 1:
#'     Condition 1: Subject 1, Subject 2, ..., Subject n
#'     Condition 2: Subject 1, Subject 2, ..., Subject n
#'     ...
#'   Group 2:
#'     Condition 1: Subject 1, Subject 2, ..., Subject m
#'     ...
#' }
#'
#' @section References:
#' McIntosh, A.R., & Lobaugh, N.J. (2004). Partial least squares analysis of
#' neuroimaging data: applications and advances. NeuroImage, 23, S250-S263.
#'
#' @docType package
#' @name plscmd-package
#' @aliases plscmd
NULL
