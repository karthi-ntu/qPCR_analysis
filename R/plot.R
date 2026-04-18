library(ggplot2)
library(ggsignif)

OKABE_ITO <- c(
  "#E69F00", "#56B4E9", "#009E73", "#F0E442",
  "#0072B2", "#D55E00", "#CC79A7", "#000000"
)

make_barplot <- function(full_df, stats_df, gene, error_type = "SD",
                         y_min = NULL, y_max = NULL) {
  sub <- full_df[full_df$Gene == gene, ]
  sub$Group <- factor(sub$Group, levels = unique(sub$Group))
  n_groups <- nlevels(sub$Group)

  p <- ggplot(sub, aes(x = Group, y = log2_fold_change, color = Group)) +
    geom_hline(yintercept = 0, linetype = "dotted", linewidth = 0.8, color = "gray40") +
    geom_jitter(width = 0.12, size = 2.5, alpha = 0.9) +
    stat_summary(fun = mean, geom = "errorbar",
                 fun.min = mean, fun.max = mean,
                 width = 0.4, linewidth = 1.2, color = "black") +
    scale_color_manual(values = rep_len(OKABE_ITO, n_groups)) +
    labs(
      title = gene,
      y     = expression(Log[2]~"(Fold Change)"),
      x     = NULL
    ) +
    theme_classic(base_size = 14) +
    theme(legend.position = "none")

  if (!is.null(y_min) || !is.null(y_max))
    p <- p + coord_cartesian(ylim = c(y_min, y_max))

  sig <- stats_df[stats_df$Gene == gene & stats_df$Significant == "Yes", ]
  if (nrow(sig) > 0) {
    comparisons_list <- lapply(sig$Comparison, function(comp)
      as.list(trimws(strsplit(comp, " vs ", fixed = TRUE)[[1]])))
    p <- p + geom_signif(
      comparisons      = comparisons_list,
      annotations      = paste0("p = ", sig$p_value),
      map_signif_level = FALSE,
      textsize         = 3.5,
      vjust            = -0.2,
      color            = "black"
    )
  }

  p
}
