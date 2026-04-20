library(ggplot2)
library(ggsignif)

OKABE_ITO <- c(
  "#E69F00", "#56B4E9", "#009E73", "#F0E442",
  "#0072B2", "#D55E00", "#CC79A7", "#000000"
)

format_pvalue <- function(p) {
  ifelse(p < 0.0001, "p < 0.0001", sprintf("p = %.4f", p))
}

make_barplot <- function(full_df, stats_df, gene, error_type = "SD",
                         y_min = NULL, y_max = NULL,
                         plot_type = "scatter",
                         show_sig = TRUE,
                         sig_comparisons = NULL) {
  sub <- full_df[full_df$Gene == gene, ]
  sub$Group <- factor(sub$Group, levels = unique(sub$Group))
  n_groups <- nlevels(sub$Group)

  err_fun <- if (error_type == "SEM") {
    function(x) {
      x <- x[!is.na(x)]
      m <- mean(x); s <- sd(x) / sqrt(length(x))
      data.frame(y = m, ymin = m - s, ymax = m + s)
    }
  } else {
    function(x) {
      x <- x[!is.na(x)]
      m <- mean(x); s <- sd(x)
      data.frame(y = m, ymin = m - s, ymax = m + s)
    }
  }

  p <- ggplot(sub, aes(x = Group, y = log2_fold_change, fill = Group)) +
    geom_hline(yintercept = 0, linetype = "dotted", linewidth = 0.8, color = "gray40")

  if (plot_type == "column") {
    # Prism "column + dots" style: bar for mean, error bars, dots overlaid
    p <- p +
      stat_summary(fun = mean, geom = "col", color = "black",
                   linewidth = 0.6, width = 0.65, alpha = 0.85) +
      stat_summary(fun.data = err_fun, geom = "errorbar",
                   width = 0.25, linewidth = 0.8, color = "black") +
      geom_jitter(shape = 21, color = "black", stroke = 0.5,
                  size = 3, width = 0.12, alpha = 0.95)
  } else {
    # "Scatter + error bars" style
    p <- p +
      stat_summary(fun.data = err_fun, geom = "errorbar",
                   width = 0.25, linewidth = 0.8, color = "black") +
      stat_summary(fun = mean, geom = "errorbar",
                   fun.min = mean, fun.max = mean,
                   width = 0.4, linewidth = 1.2, color = "black") +
      geom_jitter(shape = 21, color = "black", stroke = 0.5,
                  size = 3, width = 0.12, alpha = 0.95)
  }

  p <- p +
    scale_fill_manual(values = rep_len(OKABE_ITO, n_groups)) +
    labs(
      title = gene,
      y     = expression(Log[2]~"(Fold Change)"),
      x     = NULL
    ) +
    theme_classic(base_size = 14) +
    theme(
      legend.position   = "none",
      axis.line         = element_line(linewidth = 1, color = "black"),
      axis.ticks        = element_line(linewidth = 0.8, color = "black"),
      axis.ticks.length = unit(0.2, "cm"),
      axis.title        = element_text(face = "bold", color = "black"),
      axis.text         = element_text(face = "bold", color = "black"),
      plot.title        = element_text(face = "bold", size = 18, hjust = 0.5),
      aspect.ratio      = 1
    )

  if (!is.null(y_min) || !is.null(y_max))
    p <- p + coord_cartesian(ylim = c(y_min, y_max))

  # Significance bars
  if (isTRUE(show_sig)) {
    sig <- stats_df[stats_df$Gene == gene & stats_df$Significant == "Yes", ]
    # Skip two-way ANOVA rows (no pair to bracket)
    sig <- sig[grepl(" vs ", sig$Comparison, fixed = TRUE), ]
    # Filter to user-selected comparisons if provided
    if (!is.null(sig_comparisons))
      sig <- sig[sig$Comparison %in% sig_comparisons, ]
    if (nrow(sig) > 0) {
      comparisons_list <- lapply(sig$Comparison, function(comp)
        trimws(strsplit(comp, " vs ", fixed = TRUE)[[1]]))
      p <- p + geom_signif(
        comparisons      = comparisons_list,
        annotations      = format_pvalue(sig$p_value),
        map_signif_level = FALSE,
        step_increase    = 0.1,
        textsize         = 4.5,
        fontface         = "bold",
        vjust            = -0.4,
        tip_length       = 0.025,
        size             = 0.8,
        color            = "black"
      )
    }
  }

  p
}
