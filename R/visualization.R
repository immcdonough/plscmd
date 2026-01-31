#' PLS Visualization Functions
#'
#' Functions for creating publication-quality visualizations of PLS results.
#' Requires ggplot2 and patchwork packages.
#'
#' @name visualization
#' @keywords internal
NULL

#' Create PLS Summary Table
#'
#' Creates a summary table of latent variable statistics including singular values,
#' variance explained, and permutation p-values.
#'
#' @param result A PLS result object from \code{\link{pls_analysis}}
#' @param digits Integer; number of decimal places for rounding (default: 4)
#'
#' @return A data.frame with columns: LV, Singular_Value, Variance_Explained_Pct,
#'   Cumulative_Pct, P_Value, and Significant
#'
#' @export
#'
#' @examples
#' \dontrun{
#' result <- pls_analysis(...)
#' summary_table <- pls_summary_table(result)
#' print(summary_table)
#' }
pls_summary_table <- function(result, digits = 4) {
  s <- result$s
  prop_explained <- (s^2) / sum(s^2) * 100
  pvals <- result$perm_result$sprob

  summary_df <- data.frame(
    LV = seq_along(s),
    Singular_Value = round(s, digits),
    Variance_Explained_Pct = round(prop_explained, 2),
    Cumulative_Pct = round(cumsum(prop_explained), 2),
    P_Value = round(pvals, digits),
    Significant = ifelse(pvals < 0.05, "*", "")
  )

  class(summary_df) <- c("pls_summary", "data.frame")
  return(summary_df)
}

#' Default PLS Plot Theme
#'
#' Creates a clean ggplot2 theme for PLS visualizations with customizable
#' font sizes.
#'
#' @param base_size Base font size (default: 12)
#' @param title_size Title font size (default: 14)
#' @param axis_title_size Axis title font size (default: 12)
#' @param axis_text_size Axis text font size (default: 10)
#'
#' @return A ggplot2 theme object
#'
#' @export
#'
#' @examples
#' \dontrun{
#' library(ggplot2)
#' p + theme_pls(base_size = 14, title_size = 16)
#' }
theme_pls <- function(base_size = 12, title_size = 14,
                      axis_title_size = 12, axis_text_size = 10) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plotting. Please install it.")
  }

  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = title_size, hjust = 0.5),
      plot.subtitle = ggplot2::element_text(size = base_size, hjust = 0.5),
      axis.title = ggplot2::element_text(face = "bold", size = axis_title_size),
      axis.text = ggplot2::element_text(size = axis_text_size),
      panel.grid.minor = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 1),
      legend.position = "bottom",
      legend.title = ggplot2::element_text(face = "bold"),
      strip.text = ggplot2::element_text(face = "bold", size = axis_title_size)
    )
}

#' Default PLS Color Palette
#'
#' Returns the default color palette used in PLS visualizations.
#' Can be customized by passing to plotting functions.
#'
#' @return A named list of colors
#'
#' @export
#'
#' @examples
#' colors <- pls_colors()
#' colors$positive <- "red"  # Customize
pls_colors <- function() {
  list(
    positive = "#d73027",      # Red for positive values
    negative = "#4575b4",      # Blue for negative values
    significant = "gray40",    # Dark gray for significant
    nonsignificant = "gray80", # Light gray for non-significant
    points = "#4575b4",        # Blue for scatter points
    line = "black",            # Black for regression lines
    ci_fill = "gray80",        # Gray for confidence bands
    bar_fill = "gray50",       # Medium gray for bars
    bar_outline = "black"      # Black for bar outlines
  )
}

#' Plot Brain-Behavior Correlations
#'
#' Creates a bar plot showing correlations between brain scores and behavioral
#' measures with bootstrap confidence intervals.
#'
#' @param result A PLS result object from \code{\link{pls_analysis}}
#' @param lv Integer; which latent variable to plot (default: 1)
#' @param behav_names Character vector of behavior variable names (optional)
#' @param n_boot Integer; number of bootstrap samples for CIs (default: 1000)
#' @param colors Named list of colors (see \code{\link{pls_colors}})
#' @param font_sizes Named list with base_size, title_size, axis_title_size, axis_text_size
#'
#' @return A ggplot2 object
#'
#' @export
#'
#' @examples
#' \dontrun{
#' result <- pls_analysis(...)
#' p <- plot_behavior_correlations(result, lv = 1,
#'        behav_names = c("RT", "Accuracy"))
#' print(p)
#' }
plot_behavior_correlations <- function(result, lv = 1, behav_names = NULL,
                                        n_boot = 1000, colors = NULL,
                                        font_sizes = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plotting. Please install it.")
  }

  # Set defaults
  if (is.null(colors)) colors <- pls_colors()
  if (is.null(font_sizes)) {
    font_sizes <- list(base_size = 12, title_size = 14,
                       axis_title_size = 12, axis_text_size = 10)
  }

  # Get brain scores and behavior data
  usc <- result$usc[, lv]
  behav <- result$stacked_behavdata
  n_behav <- ncol(behav)

  # Set default behavior names
  if (is.null(behav_names)) {
    behav_names <- paste0("Behavior ", 1:n_behav)
  }

  # Calculate correlations
  correlations <- sapply(1:n_behav, function(b) {
    stats::cor(usc, behav[, b], use = "complete.obs")
  })

  # Bootstrap for CIs
  boot_corrs <- matrix(NA, n_boot, n_behav)
  n <- length(usc)

  for (i in 1:n_boot) {
    idx <- sample(1:n, n, replace = TRUE)
    for (b in 1:n_behav) {
      boot_corrs[i, b] <- stats::cor(usc[idx], behav[idx, b], use = "complete.obs")
    }
  }

  ci_low <- apply(boot_corrs, 2, stats::quantile, probs = 0.025, na.rm = TRUE)
  ci_high <- apply(boot_corrs, 2, stats::quantile, probs = 0.975, na.rm = TRUE)

  # Determine significance (CI doesn't cross zero)
  significant <- sign(ci_low) == sign(ci_high)

  # Create data frame for plotting
  plot_df <- data.frame(
    Behavior = factor(behav_names, levels = behav_names),
    Correlation = correlations,
    CI_low = ci_low,
    CI_high = ci_high,
    Significant = significant
  )

  # Create plot
  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = .data$Behavior, y = .data$Correlation,
                                              fill = .data$Significant)) +
    ggplot2::geom_bar(stat = "identity", color = colors$bar_outline, width = 0.7) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = .data$CI_low, ymax = .data$CI_high),
                           width = 0.2, linewidth = 0.8) +
    ggplot2::geom_hline(yintercept = 0, linetype = "solid", linewidth = 0.5) +
    ggplot2::scale_fill_manual(
      values = c("TRUE" = colors$significant, "FALSE" = colors$nonsignificant),
      labels = c("TRUE" = "p < .05", "FALSE" = "n.s."),
      name = ""
    ) +
    ggplot2::labs(
      title = sprintf("Brain-Behavior Correlations (LV%d)", lv),
      x = "", y = "Correlation (r)"
    ) +
    theme_pls(font_sizes$base_size, font_sizes$title_size,
              font_sizes$axis_title_size, font_sizes$axis_text_size) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      axis.title.y = ggplot2::element_text(vjust = 1),
      plot.margin = ggplot2::margin(t = 5, r = 5, b = 5, l = 5)
    )

  return(p)
}

#' Plot Brain vs Behavior Scores Scatter
#'
#' Creates a scatter plot of brain scores versus behavior scores with
#' a regression line and correlation statistics.
#'
#' @param result A PLS result object from \code{\link{pls_analysis}}
#' @param lv Integer; which latent variable to plot (default: 1)
#' @param colors Named list of colors (see \code{\link{pls_colors}})
#' @param font_sizes Named list with base_size, title_size, axis_title_size, axis_text_size
#' @param point_size Numeric; size of scatter points (default: 3)
#' @param point_alpha Numeric; transparency of points 0-1 (default: 0.7)
#' @param p_digits Integer; number of decimal places for p-value display.
#'   If NULL (default), uses scientific notation for very small values.
#'   Set to 3 for fixed 3 decimal places.
#'
#' @return A ggplot2 object
#'
#' @export
#'
#' @examples
#' \dontrun{
#' result <- pls_analysis(...)
#' p <- plot_brain_behavior_scatter(result, lv = 1)
#' print(p)
#' }
plot_brain_behavior_scatter <- function(result, lv = 1, colors = NULL,
                                         font_sizes = NULL, point_size = 3,
                                         point_alpha = 0.7, p_digits = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plotting. Please install it.")
  }

  # Set defaults
  if (is.null(colors)) colors <- pls_colors()
  if (is.null(font_sizes)) {
    font_sizes <- list(base_size = 12, title_size = 14,
                       axis_title_size = 12, axis_text_size = 10)
  }

  usc <- result$usc[, lv]  # Brain scores
  vsc <- result$vsc[, lv]  # Behavior scores

  # Remove NAs
  valid <- !is.na(usc) & !is.na(vsc)
  usc <- usc[valid]
  vsc <- vsc[valid]

  # Calculate correlation
  r_val <- stats::cor(vsc, usc)
  p_val <- stats::cor.test(vsc, usc)$p.value

  # Format p-value
  if (is.null(p_digits)) {
    p_label <- sprintf("r = %.2f, p = %.3g", r_val, p_val)
  } else {
    p_format <- paste0("%.", p_digits, "f")
    p_label <- sprintf(paste0("r = %.2f, p = ", p_format), r_val, p_val)
  }

  # Create data frame
  plot_df <- data.frame(
    Behavior_Score = vsc,
    Brain_Score = usc
  )

  # Create plot
  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = .data$Behavior_Score,
                                              y = .data$Brain_Score)) +
    ggplot2::geom_point(size = point_size, alpha = point_alpha,
                        color = colors$points, shape = 16) +
    ggplot2::geom_smooth(method = "lm", se = TRUE, color = colors$line,
                         fill = colors$ci_fill, linewidth = 1) +
    ggplot2::annotate("text",
                      x = min(vsc) + 0.05 * diff(range(vsc)),
                      y = max(usc) - 0.05 * diff(range(usc)),
                      label = p_label,
                      hjust = 0, fontface = "bold",
                      size = font_sizes$axis_text_size / 2.5) +
    ggplot2::labs(
      title = sprintf("LV%d: Brain vs Behavior Scores", lv),
      x = "Behavior Score", y = "Brain Score"
    ) +
    theme_pls(font_sizes$base_size, font_sizes$title_size,
              font_sizes$axis_title_size, font_sizes$axis_text_size) +
    ggplot2::theme(
      axis.title.x = ggplot2::element_text(vjust = 1),
      plot.margin = ggplot2::margin(t = 5, r = 5, b = 5, l = 5)
    )

  return(p)
}

#' Plot Bootstrap Ratios for Brain Regions
#'
#' Creates a bar plot of bootstrap ratios for brain regions, showing only
#' regions where the confidence interval lower bound exceeds the significance
#' threshold (i.e., reliably significant).
#'
#' @param result A PLS result object from \code{\link{pls_analysis}}
#' @param lv Integer; which latent variable to plot (default: 1)
#' @param region_names Character vector of region names (optional)
#' @param threshold Numeric; bootstrap ratio threshold for significance (default: 1.96)
#' @param top_n Integer; maximum number of regions to show if none significant (default: 30)
#' @param colors Named list of colors (see \code{\link{pls_colors}})
#' @param font_sizes Named list with base_size, title_size, axis_title_size, axis_text_size
#' @param horizontal Logical; if TRUE, bars are horizontal (default: TRUE)
#' @param show_all Logical; if TRUE, show all regions with non-significant ones in
#'   light gray. If FALSE (default), only show significant regions.
#'
#' @return A ggplot2 object
#'
#' @details
#' A region is considered reliably significant only if the lower bound of its
#' bootstrap ratio confidence interval (|BSR| - SE) exceeds the threshold.
#' This is more conservative than simply checking if |BSR| > threshold.
#'
#' When \code{show_all = TRUE}, all regions are displayed with:
#' \itemize{
#'   \item Significant positive BSR: red (colors$positive)
#'   \item Significant negative BSR: blue (colors$negative)
#'   \item Non-significant: light gray (colors$nonsignificant)
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' result <- pls_analysis(...)
#' p <- plot_bootstrap_ratios(result, lv = 1,
#'        region_names = paste0("ROI_", 1:100))
#' print(p)
#' }
plot_bootstrap_ratios <- function(result, lv = 1, region_names = NULL,
                                   threshold = 1.96, top_n = 30,
                                   colors = NULL, font_sizes = NULL,
                                   horizontal = TRUE, show_all = FALSE) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plotting. Please install it.")
  }


  # Set defaults
  if (is.null(colors)) colors <- pls_colors()
  if (is.null(font_sizes)) {
    font_sizes <- list(base_size = 12, title_size = 14,
                       axis_title_size = 12, axis_text_size = 10)
  }

  bsr <- result$boot_result$compare_u[, lv]
  se <- result$boot_result$u_se[, lv]

  n_regions <- length(bsr)
  if (is.null(region_names)) {
    region_names <- paste0("Region ", 1:n_regions)
  }

  # Calculate lower bound of CI for significance determination
  # A region is significant only if the LOWER BOUND of CI exceeds threshold
  lower_bound <- abs(bsr) - se
  significant <- lower_bound > threshold

  # Create data frame
  plot_df <- data.frame(
    Region = region_names,
    BSR = bsr,
    SE = se,
    Lower_Bound = lower_bound,
    Significant = significant
  )

  if (show_all) {
    # Show all regions with non-significant in gray
    # Create fill category: "pos_sig", "neg_sig", or "nonsig"
    plot_df$FillCategory <- ifelse(!plot_df$Significant, "nonsig",
                                    ifelse(plot_df$BSR > 0, "pos_sig", "neg_sig"))

    # Sort by BSR value
    plot_df <- plot_df[order(plot_df$BSR, decreasing = TRUE), ]
    plot_df$Region <- factor(plot_df$Region, levels = plot_df$Region)

    # Create plot with three-color scheme
    p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = .data$Region, y = .data$BSR,
                                                fill = .data$FillCategory)) +
      ggplot2::geom_bar(stat = "identity", color = colors$bar_outline, width = 0.7) +
      ggplot2::geom_errorbar(ggplot2::aes(ymin = .data$BSR - .data$SE,
                                           ymax = .data$BSR + .data$SE),
                             width = 0.2, linewidth = 0.6) +
      ggplot2::geom_hline(yintercept = c(-threshold, threshold),
                          linetype = "dashed", color = "gray40", linewidth = 0.8) +
      ggplot2::geom_hline(yintercept = c(-3.3, 3.3),
                          linetype = "dotted", color = "gray40", linewidth = 0.8) +
      ggplot2::scale_fill_manual(
        values = c("pos_sig" = colors$positive,
                   "neg_sig" = colors$negative,
                   "nonsig" = colors$nonsignificant),
        labels = c("pos_sig" = "Positive (sig.)",
                   "neg_sig" = "Negative (sig.)",
                   "nonsig" = "Non-significant"),
        name = ""
      ) +
      ggplot2::labs(
        title = sprintf("Bootstrap Ratios (LV%d)", lv),
        x = "", y = "Bootstrap Ratio"
      ) +
      theme_pls(font_sizes$base_size, font_sizes$title_size,
                font_sizes$axis_title_size, font_sizes$axis_text_size)

  } else {
    # Original behavior: filter to significant regions only
    if (sum(plot_df$Significant) > 0) {
      plot_df <- plot_df[plot_df$Significant, ]
    } else {
      # If none significant, show top N by lower bound
      plot_df <- plot_df[order(plot_df$Lower_Bound, decreasing = TRUE)[1:min(top_n, n_regions)], ]
      message(sprintf("No regions with lower bound CI > %.2f. Showing top %d by |BSR| - SE.",
                      threshold, min(top_n, n_regions)))
    }

    # Sort by BSR value
    plot_df <- plot_df[order(plot_df$BSR, decreasing = TRUE), ]
    plot_df$Region <- factor(plot_df$Region, levels = plot_df$Region)

    # Create plot
    p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = .data$Region, y = .data$BSR,
                                                fill = .data$BSR > 0)) +
      ggplot2::geom_bar(stat = "identity", color = colors$bar_outline, width = 0.7) +
      ggplot2::geom_errorbar(ggplot2::aes(ymin = .data$BSR - .data$SE,
                                           ymax = .data$BSR + .data$SE),
                             width = 0.2, linewidth = 0.6) +
      ggplot2::geom_hline(yintercept = c(-threshold, threshold),
                          linetype = "dashed", color = "gray40", linewidth = 0.8) +
      ggplot2::geom_hline(yintercept = c(-3.3, 3.3),
                          linetype = "dotted", color = "gray40", linewidth = 0.8) +
      ggplot2::scale_fill_manual(
        values = c("TRUE" = colors$positive, "FALSE" = colors$negative),
        guide = "none"
      ) +
      ggplot2::labs(
        title = sprintf("Bootstrap Ratios (LV%d)", lv),
        x = "", y = "Bootstrap Ratio"
      ) +
      theme_pls(font_sizes$base_size, font_sizes$title_size,
                font_sizes$axis_title_size, font_sizes$axis_text_size)
  }

  if (horizontal) {
    p <- p + ggplot2::coord_flip()
  } else {
    p <- p + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  }

  return(p)
}

#' Plot LV Correlations with Bootstrap CIs
#'
#' Creates a bar plot of behavior correlations for a latent variable using
#' bootstrap confidence intervals from the PLS analysis.
#'
#' @param result A PLS result object from \code{\link{pls_analysis}}
#' @param lv Integer; which latent variable to plot (default: 1)
#' @param behav_names Character vector of behavior variable names (optional)
#' @param colors Named list of colors (see \code{\link{pls_colors}})
#' @param font_sizes Named list with base_size, title_size, axis_title_size, axis_text_size
#' @param horizontal Logical; if TRUE, bars are horizontal (default: TRUE)
#'
#' @return A ggplot2 object
#'
#' @export
#'
#' @examples
#' \dontrun{
#' result <- pls_analysis(...)
#' p <- plot_lv_correlations(result, lv = 1,
#'        behav_names = c("Age", "IQ", "Memory"))
#' print(p)
#' }
plot_lv_correlations <- function(result, lv = 1, behav_names = NULL,
                                  colors = NULL, font_sizes = NULL,
                                  horizontal = TRUE) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plotting. Please install it.")
  }

  # Set defaults
  if (is.null(colors)) colors <- pls_colors()
  if (is.null(font_sizes)) {
    font_sizes <- list(base_size = 12, title_size = 14,
                       axis_title_size = 12, axis_text_size = 10)
  }

  # Original correlations
  orig_corr <- result$boot_result$orig_corr[, lv]
  ul_corr <- result$boot_result$ulcorr[, lv]
  ll_corr <- result$boot_result$llcorr[, lv]

  n_behav <- length(orig_corr)
  if (is.null(behav_names)) {
    behav_names <- paste0("Behavior ", 1:n_behav)
  }

  # Determine significance (CI doesn't cross zero)
  significant <- sign(ul_corr) == sign(ll_corr)

  # Create data frame
  plot_df <- data.frame(
    Behavior = factor(behav_names, levels = behav_names),
    Correlation = orig_corr,
    CI_low = ll_corr,
    CI_high = ul_corr,
    Significant = significant
  )

  # Sort by correlation
  plot_df <- plot_df[order(plot_df$Correlation, decreasing = TRUE), ]
  plot_df$Behavior <- factor(plot_df$Behavior, levels = plot_df$Behavior)

  # Create plot
  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = .data$Behavior, y = .data$Correlation,
                                              fill = .data$Significant)) +
    ggplot2::geom_bar(stat = "identity", color = colors$bar_outline, width = 0.7) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = .data$CI_low, ymax = .data$CI_high),
                           width = 0.2, linewidth = 0.6) +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.5) +
    ggplot2::scale_fill_manual(
      values = c("TRUE" = colors$significant, "FALSE" = "white"),
      labels = c("TRUE" = "p < .05", "FALSE" = "n.s."),
      name = ""
    ) +
    ggplot2::labs(
      title = sprintf("LV%d Behavior Correlations", lv),
      x = "", y = "Correlation (r)"
    ) +
    theme_pls(font_sizes$base_size, font_sizes$title_size,
              font_sizes$axis_title_size, font_sizes$axis_text_size)

  if (horizontal) {
    p <- p + ggplot2::coord_flip()
  } else {
    p <- p + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  }

  return(p)
}

#' Create Combined LV Summary Figure
#'
#' Creates a multi-panel figure combining behavior correlations, brain-behavior
#' scatter plot, and bootstrap ratios for a single latent variable.
#'
#' @param result A PLS result object from \code{\link{pls_analysis}}
#' @param lv Integer; which latent variable to plot (default: 1)
#' @param behav_names Character vector of behavior variable names (optional)
#' @param region_names Character vector of region names (optional)
#' @param colors Named list of colors (see \code{\link{pls_colors}})
#' @param font_sizes Named list with base_size, title_size, axis_title_size, axis_text_size
#' @param save_path Character; file path to save the figure (optional)
#' @param width Numeric; figure width in inches (default: 12)
#' @param height Numeric; figure height in inches (default: 10)
#' @param dpi Numeric; resolution for saved figure (default: 300)
#' @param show_all_regions Logical; if TRUE, show all brain regions in bootstrap
#'   ratio plot with non-significant ones in light gray (default: FALSE)
#' @param p_digits Integer; number of decimal places for p-value display.
#'   If NULL (default), uses 4 decimal places for title and scientific notation
#'   for scatter plot. Set to 3 for fixed 3 decimal places everywhere.
#' @param tag_levels Character; panel tag style. Use "A" for letters (A, B, C),
#'   "1" for numbers (1, 2, 3), or NULL for no tags (default).
#'
#' @return A patchwork object combining the plots, or NULL if LV not significant
#'
#' @export
#'
#' @examples
#' \dontrun{
#' result <- pls_analysis(...)
#' fig <- plot_lv_summary(result, lv = 1,
#'          behav_names = c("RT", "Accuracy"),
#'          save_path = "LV1_summary.png",
#'          tag_levels = "A")  # Add A, B, C tags
#' }
plot_lv_summary <- function(result, lv = 1, behav_names = NULL,
                             region_names = NULL, colors = NULL,
                             font_sizes = NULL, save_path = NULL,
                             width = 12, height = 10, dpi = 300,
                             show_all_regions = FALSE, p_digits = NULL,
                             tag_levels = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plotting. Please install it.")
  }
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop("Package 'patchwork' is required for combined plots. Please install it.")
  }

  # Check if LV is significant
  if (result$perm_result$sprob[lv] >= 0.05) {
    message(sprintf("LV%d is not significant (p = %.3f)", lv, result$perm_result$sprob[lv]))
    return(NULL)
  }

  # Create individual plots
  p1 <- plot_behavior_correlations(result, lv, behav_names, colors = colors,
                                    font_sizes = font_sizes)
  p2 <- plot_brain_behavior_scatter(result, lv, colors = colors,
                                     font_sizes = font_sizes, p_digits = p_digits)
  p3 <- plot_bootstrap_ratios(result, lv, region_names, colors = colors,
                               font_sizes = font_sizes, show_all = show_all_regions)

  # Format p-value for title
  p_val <- result$perm_result$sprob[lv]
  if (is.null(p_digits)) {
    title_text <- sprintf("Latent Variable %d Summary (p = %.4f)", lv, p_val)
  } else {
    p_format <- paste0("%.", p_digits, "f")
    title_text <- sprintf(paste0("Latent Variable %d Summary (p = ", p_format, ")"), lv, p_val)
  }

  # Combine with patchwork
  combined <- (p1 | p2) / p3 +
    patchwork::plot_annotation(
      title = title_text,
      tag_levels = tag_levels,
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold", size = 16, hjust = 0.5),
        plot.tag = ggplot2::element_text(face = "bold", size = 14)
      )
    )

  # Save if path provided
  if (!is.null(save_path)) {
    ggplot2::ggsave(save_path, combined, width = width, height = height, dpi = dpi)
    message(sprintf("Figure saved to: %s", save_path))
  }

  return(combined)
}

#' Plot All Significant LVs
#'
#' Creates and optionally saves summary figures for all statistically
#' significant latent variables.
#'
#' @param result A PLS result object from \code{\link{pls_analysis}}
#' @param behav_names Character vector of behavior variable names (optional)
#' @param region_names Character vector of region names (optional)
#' @param colors Named list of colors (see \code{\link{pls_colors}})
#' @param font_sizes Named list with base_size, title_size, axis_title_size, axis_text_size
#' @param output_dir Character; directory to save figures (default: current directory)
#' @param prefix Character; prefix for saved file names (default: "LV")
#' @param width Numeric; figure width in inches (default: 12)
#' @param height Numeric; figure height in inches (default: 10)
#' @param dpi Numeric; resolution for saved figures (default: 300)
#' @param show_all_regions Logical; if TRUE, show all brain regions in bootstrap
#'   ratio plots with non-significant ones in light gray (default: FALSE)
#' @param p_digits Integer; number of decimal places for p-value display.
#'   If NULL (default), uses 4 decimal places for title and scientific notation
#'   for scatter plot. Set to 3 for fixed 3 decimal places everywhere.
#' @param tag_levels Character; panel tag style. Use "A" for letters (A, B, C),
#'   "1" for numbers (1, 2, 3), or NULL for no tags (default).
#'
#' @return A list of patchwork objects, one per significant LV
#'
#' @export
#'
#' @examples
#' \dontrun{
#' result <- pls_analysis(...)
#' figures <- plot_all_significant_lvs(result,
#'              behav_names = c("Age", "IQ"),
#'              output_dir = "figures",
#'              tag_levels = "A")  # Add A, B, C tags
#' }
plot_all_significant_lvs <- function(result, behav_names = NULL,
                                      region_names = NULL, colors = NULL,
                                      font_sizes = NULL, output_dir = ".",
                                      prefix = "LV", width = 12,
                                      height = 10, dpi = 300,
                                      show_all_regions = FALSE,
                                      p_digits = NULL,
                                      tag_levels = NULL) {
  sig_lvs <- which(result$perm_result$sprob < 0.05)

  if (length(sig_lvs) == 0) {
    message("No significant LVs found (p < 0.05)")
    return(NULL)
  }

  message(sprintf("Found %d significant LV(s): %s",
                  length(sig_lvs), paste(sig_lvs, collapse = ", ")))

  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  figures <- list()
  for (lv in sig_lvs) {
    save_path <- file.path(output_dir, sprintf("%s%d_summary.png", prefix, lv))
    figures[[lv]] <- plot_lv_summary(
      result, lv, behav_names, region_names,
      colors = colors, font_sizes = font_sizes,
      save_path = save_path, width = width, height = height, dpi = dpi,
      show_all_regions = show_all_regions, p_digits = p_digits,
      tag_levels = tag_levels
    )
  }

  return(figures)
}


# ===========================================================================
# PLS Regression Visualization Functions
# ===========================================================================

#' Plot VIP Scores for PLS Regression
#'
#' Creates a bar plot of Variable Importance in Projection (VIP) scores
#' with an optional threshold line.
#'
#' @param result A pls_regression_result object
#' @param comp Component number to use for VIP (default: final component)
#' @param threshold VIP threshold line (default: 1.0)
#' @param var_names Optional character vector of variable names
#' @param top_n Number of top variables to show (default: 30)
#' @param colors Named list of colors (see \code{\link{pls_colors}})
#' @param font_sizes Named list with base_size, title_size, axis_title_size, axis_text_size
#' @param show_ci Show bootstrap confidence intervals if available (default: TRUE)
#'
#' @return A ggplot2 object
#'
#' @export
#'
#' @examples
#' \dontrun{
#' fit <- pls_regression(X, Y, ncomp = 5)
#' p <- plot_vip(fit, top_n = 20)
#' print(p)
#' }
plot_vip <- function(result, comp = NULL, threshold = 1.0,
                      var_names = NULL, top_n = 30,
                      colors = NULL, font_sizes = NULL,
                      show_ci = TRUE) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plotting")
  }

  if (!inherits(result, "pls_regression_result")) {
    stop("result must be a pls_regression_result object")
  }

  if (is.null(colors)) colors <- pls_colors()
  if (is.null(font_sizes)) {
    font_sizes <- list(base_size = 12, title_size = 14,
                       axis_title_size = 12, axis_text_size = 10)
  }

  if (is.null(comp)) comp <- result$ncomp
  p <- result$p

  vip <- result$vip[, comp]

  if (is.null(var_names)) {
    var_names <- paste0("X", seq_len(p))
  }

  # Create data frame
  df <- data.frame(
    Variable = var_names,
    VIP = vip,
    Important = vip > threshold
  )

  # Add CIs if available
  has_ci <- !is.null(result$boot_result) && show_ci
  if (has_ci) {
    df$CI_low <- result$boot_result$vip_ci_low
    df$CI_high <- result$boot_result$vip_ci_high
  }

  # Sort and filter
  df <- df[order(df$VIP, decreasing = TRUE), ]
  if (!is.null(top_n) && top_n < nrow(df)) {
    df <- df[1:top_n, ]
  }
  df$Variable <- factor(df$Variable, levels = rev(df$Variable))

  # Create plot
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$Variable, y = .data$VIP,
                                         fill = .data$Important)) +
    ggplot2::geom_bar(stat = "identity", color = colors$bar_outline, width = 0.7) +
    ggplot2::geom_hline(yintercept = threshold, linetype = "dashed",
                        color = colors$positive, linewidth = 0.8) +
    ggplot2::scale_fill_manual(
      values = c("TRUE" = colors$positive, "FALSE" = colors$nonsignificant),
      labels = c("TRUE" = paste0("VIP > ", threshold), "FALSE" = paste0("VIP <= ", threshold)),
      name = "Importance"
    ) +
    ggplot2::coord_flip() +
    ggplot2::labs(
      title = sprintf("Variable Importance in Projection (Component %d)", comp),
      x = "",
      y = "VIP Score"
    ) +
    theme_pls(base_size = font_sizes$base_size, title_size = font_sizes$title_size,
              axis_title_size = font_sizes$axis_title_size,
              axis_text_size = font_sizes$axis_text_size)

  # Add CIs if available
  if (has_ci) {
    p <- p + ggplot2::geom_errorbar(
      ggplot2::aes(ymin = .data$CI_low, ymax = .data$CI_high),
      width = 0.2, linewidth = 0.5
    )
  }

  p
}


#' Plot PLS Regression Loadings
#'
#' Creates a bar plot of X or Y loadings for a specified component.
#'
#' @param result A pls_regression_result object
#' @param comp Component number (default: 1)
#' @param type Which loadings: "X", "Y", or "both" (default: "X")
#' @param var_names Optional character vector of variable names
#' @param top_n Number of top variables to show (default: 30 for X, all for Y)
#' @param colors Named list of colors
#' @param font_sizes Named list of font sizes
#'
#' @return A ggplot2 object
#'
#' @export
plot_loadings <- function(result, comp = 1, type = "X",
                           var_names = NULL, top_n = NULL,
                           colors = NULL, font_sizes = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plotting")
  }

  if (!inherits(result, "pls_regression_result")) {
    stop("result must be a pls_regression_result object")
  }

  if (is.null(colors)) colors <- pls_colors()
  if (is.null(font_sizes)) {
    font_sizes <- list(base_size = 12, title_size = 14,
                       axis_title_size = 12, axis_text_size = 10)
  }

  type <- match.arg(type, c("X", "Y", "both"))

  if (type == "both") {
    if (!requireNamespace("patchwork", quietly = TRUE)) {
      stop("Package 'patchwork' is required for combined plots")
    }
    p1 <- plot_loadings(result, comp, "X", var_names, top_n, colors, font_sizes)
    p2 <- plot_loadings(result, comp, "Y", NULL, NULL, colors, font_sizes)
    return(p1 / p2)
  }

  if (type == "X") {
    loadings <- result$X_loadings[, comp]
    title_prefix <- "X Loadings"
    if (is.null(var_names)) var_names <- paste0("X", seq_len(result$p))
    if (is.null(top_n)) top_n <- 30
  } else {
    loadings <- result$Y_loadings[, comp]
    title_prefix <- "Y Loadings"
    if (is.null(var_names)) var_names <- paste0("Y", seq_len(result$q))
    if (is.null(top_n)) top_n <- result$q
  }

  # Create data frame
  df <- data.frame(
    Variable = var_names,
    Loading = loadings,
    Direction = ifelse(loadings >= 0, "Positive", "Negative")
  )

  # Sort by absolute value and filter
  df <- df[order(abs(df$Loading), decreasing = TRUE), ]
  if (top_n < nrow(df)) {
    df <- df[1:top_n, ]
  }
  df$Variable <- factor(df$Variable, levels = rev(df$Variable))

  # Create plot
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$Variable, y = .data$Loading,
                                         fill = .data$Direction)) +
    ggplot2::geom_bar(stat = "identity", color = colors$bar_outline, width = 0.7) +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.5) +
    ggplot2::scale_fill_manual(
      values = c("Positive" = colors$positive, "Negative" = colors$negative),
      name = "Direction"
    ) +
    ggplot2::coord_flip() +
    ggplot2::labs(
      title = sprintf("%s (Component %d)", title_prefix, comp),
      x = "",
      y = "Loading"
    ) +
    theme_pls(base_size = font_sizes$base_size, title_size = font_sizes$title_size,
              axis_title_size = font_sizes$axis_title_size,
              axis_text_size = font_sizes$axis_text_size)

  p
}


#' Plot Predicted vs Observed Values
#'
#' Creates a scatter plot of predicted vs observed Y values with regression line.
#'
#' @param result A pls_regression_result object
#' @param response_idx Which response variable to plot (default: 1)
#' @param colors Named list of colors
#' @param font_sizes Named list of font sizes
#'
#' @return A ggplot2 object
#'
#' @export
plot_predicted_vs_observed <- function(result, response_idx = 1,
                                        colors = NULL, font_sizes = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plotting")
  }

  if (!inherits(result, "pls_regression_result")) {
    stop("result must be a pls_regression_result object")
  }

  if (is.null(colors)) colors <- pls_colors()
  if (is.null(font_sizes)) {
    font_sizes <- list(base_size = 12, title_size = 14,
                       axis_title_size = 12, axis_text_size = 10)
  }

  # Get observed and predicted
  if (result$q == 1) {
    observed <- as.vector(result$Y_fitted - result$Y_residuals *
                           ifelse(is.null(result$Y_scale), 1, result$Y_scale))
    predicted <- as.vector(result$Y_fitted)
  } else {
    observed <- result$Y_fitted[, response_idx] -
                result$Y_residuals[, response_idx] *
                ifelse(is.null(result$Y_scale), 1, result$Y_scale[response_idx])
    predicted <- result$Y_fitted[, response_idx]
  }

  df <- data.frame(Observed = observed, Predicted = predicted)

  # Calculate R2
  r2 <- result$R2Y_cum[result$ncomp]

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$Observed, y = .data$Predicted)) +
    ggplot2::geom_point(color = colors$points, alpha = 0.6, size = 2) +
    ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed",
                         color = colors$line, linewidth = 0.8) +
    ggplot2::geom_smooth(method = "lm", se = TRUE, color = colors$positive,
                         fill = colors$ci_fill, alpha = 0.3) +
    ggplot2::labs(
      title = "Predicted vs Observed",
      subtitle = sprintf("R-squared = %.3f", r2),
      x = "Observed",
      y = "Predicted"
    ) +
    theme_pls(base_size = font_sizes$base_size, title_size = font_sizes$title_size,
              axis_title_size = font_sizes$axis_title_size,
              axis_text_size = font_sizes$axis_text_size)

  p
}


#' Plot Variance Explained
#'
#' Creates a bar plot showing variance explained per component.
#'
#' @param result A pls_regression_result object
#' @param type "Y" (default), "X", or "both"
#' @param colors Named list of colors
#' @param font_sizes Named list of font sizes
#'
#' @return A ggplot2 object
#'
#' @export
plot_variance_explained <- function(result, type = "Y",
                                     colors = NULL, font_sizes = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plotting")
  }

  if (!inherits(result, "pls_regression_result")) {
    stop("result must be a pls_regression_result object")
  }

  if (is.null(colors)) colors <- pls_colors()
  if (is.null(font_sizes)) {
    font_sizes <- list(base_size = 12, title_size = 14,
                       axis_title_size = 12, axis_text_size = 10)
  }

  type <- match.arg(type, c("Y", "X", "both"))

  ncomp <- result$ncomp

  if (type == "both") {
    df <- data.frame(
      Component = rep(1:ncomp, 2),
      Variance = c(result$R2X * 100, result$R2Y * 100),
      Type = rep(c("X", "Y"), each = ncomp)
    )

    p <- ggplot2::ggplot(df, ggplot2::aes(x = factor(.data$Component),
                                           y = .data$Variance,
                                           fill = .data$Type)) +
      ggplot2::geom_bar(stat = "identity", position = "dodge", width = 0.7) +
      ggplot2::scale_fill_manual(values = c("X" = colors$negative, "Y" = colors$positive)) +
      ggplot2::labs(
        title = "Variance Explained per Component",
        x = "Component",
        y = "Variance Explained (%)"
      )
  } else {
    var_exp <- if (type == "Y") result$R2Y * 100 else result$R2X * 100
    cum_var <- cumsum(var_exp)

    df <- data.frame(
      Component = factor(1:ncomp),
      Variance = var_exp,
      Cumulative = cum_var
    )

    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$Component, y = .data$Variance)) +
      ggplot2::geom_bar(stat = "identity", fill = colors$bar_fill,
                        color = colors$bar_outline, width = 0.7) +
      ggplot2::geom_line(ggplot2::aes(y = .data$Cumulative, group = 1),
                         color = colors$positive, linewidth = 1) +
      ggplot2::geom_point(ggplot2::aes(y = .data$Cumulative),
                          color = colors$positive, size = 3) +
      ggplot2::labs(
        title = sprintf("%s Variance Explained", type),
        subtitle = sprintf("Total: %.1f%%", cum_var[ncomp]),
        x = "Component",
        y = "Variance Explained (%)"
      )
  }

  p + theme_pls(base_size = font_sizes$base_size, title_size = font_sizes$title_size,
                axis_title_size = font_sizes$axis_title_size,
                axis_text_size = font_sizes$axis_text_size)
}


#' Plot PLS Regression Summary
#'
#' Creates a combined multi-panel figure summarizing PLS regression results.
#'
#' @param result A pls_regression_result object
#' @param var_names Optional predictor variable names
#' @param colors Named list of colors
#' @param font_sizes Named list of font sizes
#' @param save_path Optional path to save figure
#' @param width Figure width in inches (default: 12)
#' @param height Figure height in inches (default: 10)
#' @param dpi Resolution for saved figure (default: 300)
#'
#' @return A patchwork combined plot
#'
#' @export
#'
#' @examples
#' \dontrun{
#' fit <- pls_regression(X, Y, ncomp = 5)
#' p <- plot_regression_summary(fit)
#' print(p)
#' }
plot_regression_summary <- function(result, var_names = NULL,
                                     colors = NULL, font_sizes = NULL,
                                     save_path = NULL, width = 12,
                                     height = 10, dpi = 300) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plotting")
  }
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop("Package 'patchwork' is required for combined plots")
  }

  if (!inherits(result, "pls_regression_result")) {
    stop("result must be a pls_regression_result object")
  }

  if (is.null(colors)) colors <- pls_colors()
  if (is.null(font_sizes)) {
    font_sizes <- list(base_size = 10, title_size = 12,
                       axis_title_size = 10, axis_text_size = 8)
  }

  # Create panels
  p1 <- plot_variance_explained(result, type = "both", colors = colors, font_sizes = font_sizes)
  p2 <- plot_predicted_vs_observed(result, colors = colors, font_sizes = font_sizes)
  p3 <- plot_vip(result, var_names = var_names, top_n = 20,
                  colors = colors, font_sizes = font_sizes)
  p4 <- plot_loadings(result, comp = 1, type = "X", var_names = var_names,
                       top_n = 20, colors = colors, font_sizes = font_sizes)

  # Combine with patchwork
  combined <- (p1 | p2) / (p3 | p4) +
    patchwork::plot_annotation(
      title = sprintf("PLS Regression Summary (%s, %d components)",
                      toupper(result$algorithm), result$ncomp),
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold", size = 16, hjust = 0.5)
      )
    )

  if (!is.null(save_path)) {
    ggplot2::ggsave(save_path, combined, width = width, height = height, dpi = dpi)
    message(sprintf("Saved to %s", save_path))
  }

  combined
}
