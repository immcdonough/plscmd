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
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

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
                                         point_alpha = 0.7) {
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
                      label = sprintf("r = %.2f, p = %.3g", r_val, p_val),
                      hjust = 0, fontface = "bold",
                      size = font_sizes$axis_text_size / 2.5) +
    ggplot2::labs(
      title = sprintf("LV%d: Brain vs Behavior Scores", lv),
      x = "Behavior Score", y = "Brain Score"
    ) +
    theme_pls(font_sizes$base_size, font_sizes$title_size,
              font_sizes$axis_title_size, font_sizes$axis_text_size)

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
        subtitle = "Dashed: p < .05, Dotted: p < .001 (CI lower bound must exceed threshold)",
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
        subtitle = "Dashed: p < .05, Dotted: p < .001 (CI lower bound must exceed threshold)",
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
#'          save_path = "LV1_summary.png")
#' }
plot_lv_summary <- function(result, lv = 1, behav_names = NULL,
                             region_names = NULL, colors = NULL,
                             font_sizes = NULL, save_path = NULL,
                             width = 12, height = 10, dpi = 300,
                             show_all_regions = FALSE) {
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
                                     font_sizes = font_sizes)
  p3 <- plot_bootstrap_ratios(result, lv, region_names, colors = colors,
                               font_sizes = font_sizes, show_all = show_all_regions)

  # Combine with patchwork
  combined <- (p1 | p2) / p3 +
    patchwork::plot_annotation(
      title = sprintf("Latent Variable %d Summary (p = %.4f)",
                      lv, result$perm_result$sprob[lv]),
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold", size = 16, hjust = 0.5)
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
#'              output_dir = "figures")
#' }
plot_all_significant_lvs <- function(result, behav_names = NULL,
                                      region_names = NULL, colors = NULL,
                                      font_sizes = NULL, output_dir = ".",
                                      prefix = "LV", width = 12,
                                      height = 10, dpi = 300,
                                      show_all_regions = FALSE) {
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
      show_all_regions = show_all_regions
    )
  }

  return(figures)
}
