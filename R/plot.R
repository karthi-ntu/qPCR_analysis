library(ggplot2)
library(ggsignif)

OKABE_ITO <- c(
  "#E69F00", "#56B4E9", "#009E73", "#F0E442",
  "#0072B2", "#D55E00", "#CC79A7", "#000000"
)

make_barplot <- function(summary_df, stats_df, gene, error_type = "SD") {
  sub <- summary_df[summary_df$Gene == gene, ]
  sub$error <- if (error_type == "SD") sub$sd_fc else sub$sem_fc

  sub$Group <- factor(sub$Group, levels = unique(sub$Group))

  p <- ggplot(sub, aes(x = Group, y = mean_fc, fill = Group)) +
    geom_bar(stat = "identity", color = "black", width = 0.6) +
    geom_errorbar(
      aes(ymin = mean_fc - error, ymax = mean_fc + error),
      width = 0.2, linewidth = 0.7
    ) +
    scale_fill_manual(values = rep_len(OKABE_ITO, nrow(sub))) +
    labs(
      title = gene,
      y     = expression("Fold Change (2"^{-Delta*Delta*Ct}*")"),
      x     = NULL
    ) +
    theme_classic(base_size = 14) +
    theme(legend.position = "none")

  sig <- stats_df[stats_df$Gene == gene & stats_df$Significant == "Yes", ]

  if (nrow(sig) > 0) {
    comparisons_list <- lapply(sig$Comparison, function(comp) {
      parts <- trimws(strsplit(comp, " vs ", fixed = TRUE)[[1]])
      as.list(parts)
    })
    annotations <- paste0("p = ", sig$p_value)

    p <- p + geom_signif(
      comparisons      = comparisons_list,
      annotations      = annotations,
      map_signif_level = FALSE,
      textsize         = 3.5,
      vjust            = -0.2
    )
  }

  p
}
