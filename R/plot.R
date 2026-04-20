library(ggplot2)
library(ggsignif)

OKABE_ITO <- c(
  "#E69F00", "#56B4E9", "#009E73", "#F0E442",
  "#0072B2", "#D55E00", "#CC79A7", "#000000"
)

format_pvalue <- function(p) {
  ifelse(p < 0.0001, "p < 0.0001", sprintf("p = %.4f", p))
}

# Error-bar summary function factory
make_err_fun <- function(error_type) {
  if (error_type == "SEM") {
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
}

prism_theme <- function(rotate_x = FALSE) {
  theme_classic(base_size = 14) +
    theme(
      legend.position   = "none",
      axis.line         = element_line(linewidth = 1, color = "black"),
      axis.ticks        = element_line(linewidth = 0.8, color = "black"),
      axis.ticks.length = unit(0.2, "cm"),
      axis.title        = element_text(face = "bold", color = "black"),
      axis.text.y       = element_text(face = "bold", color = "black"),
      axis.text.x       = element_text(
        face  = "bold", color = "black",
        angle = if (rotate_x) 45 else 0,
        hjust = if (rotate_x) 1  else 0.5
      ),
      plot.title        = element_text(face = "bold.italic", size = 18, hjust = 0.5),
      strip.text        = element_text(face = "bold.italic", size = 14),
      strip.background  = element_blank()
    )
}

make_barplot <- function(full_df, stats_df, gene,
                         error_type    = "SD",
                         y_min         = NULL, y_max = NULL,
                         plot_type     = "scatter",
                         show_sig      = TRUE,
                         sig_comparisons = NULL,
                         control_group = "Control",
                         rotate_x      = FALSE,
                         aspect_ratio  = 1) {
  sub <- full_df[full_df$Gene == gene, ]
  sub$Group <- factor(sub$Group, levels = unique(sub$Group))
  n_groups <- nlevels(sub$Group)
  err_fun  <- make_err_fun(error_type)

  p <- ggplot(sub, aes(x = Group, y = log2_fold_change, fill = Group)) +
    geom_hline(yintercept = 0, linetype = "dotted", linewidth = 0.8, color = "gray40")

  if (plot_type == "column") {
    p <- p +
      stat_summary(fun = mean, geom = "col", color = "black",
                   linewidth = 0.6, width = 0.65, alpha = 0.85) +
      stat_summary(fun.data = err_fun, geom = "errorbar",
                   width = 0.25, linewidth = 0.8, color = "black") +
      geom_jitter(shape = 21, color = "black", stroke = 0.5,
                  size = 3, width = 0.12, alpha = 0.95)
  } else {
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
      y     = bquote(Log[2] ~ "(FC to " * .(control_group) * ")"),
      x     = NULL
    ) +
    prism_theme(rotate_x = rotate_x) +
    theme(aspect.ratio = aspect_ratio)

  if (!is.null(y_min) || !is.null(y_max))
    p <- p + coord_cartesian(ylim = c(if (is.null(y_min)) NA else y_min,
                                       if (is.null(y_max)) NA else y_max))

  if (isTRUE(show_sig)) {
    sig <- stats_df[stats_df$Gene == gene & stats_df$Significant == "Yes", ]
    sig <- sig[grepl(" vs ", sig$Comparison, fixed = TRUE), ]
    if (!is.null(sig_comparisons))
      sig <- sig[sig$Comparison %in% sig_comparisons, ]
    if (nrow(sig) > 0) {
      comparisons_list <- lapply(sig$Comparison, function(comp)
        trimws(strsplit(comp, " vs ", fixed = TRUE)[[1]]))
      p <- p + geom_signif(
        comparisons      = comparisons_list,
        annotations      = format_pvalue(sig$p_value),
        map_signif_level = FALSE,
        step_increase    = 0.13,
        textsize         = 3.5,
        fontface         = "bold",
        vjust            = -0.6,
        tip_length       = 0.02,
        size             = 0.7,
        color            = "black"
      )
    }
  }

  p
}

# Combined multi-gene figure using facet_wrap (g-tibo style)
make_combined_plot <- function(full_df, stats_df,
                               error_type    = "SD",
                               y_min         = NULL, y_max = NULL,
                               plot_type     = "scatter",
                               show_sig      = TRUE,
                               sig_comparisons_all = NULL,
                               control_group = "Control",
                               rotate_x      = FALSE,
                               aspect_ratio  = 1,
                               ncol          = 3) {
  full_df$Gene  <- factor(full_df$Gene,  levels = unique(full_df$Gene))
  full_df$Group <- factor(full_df$Group, levels = unique(full_df$Group))
  n_groups <- nlevels(full_df$Group)
  err_fun  <- make_err_fun(error_type)

  p <- ggplot(full_df, aes(x = Group, y = log2_fold_change, fill = Group)) +
    geom_hline(yintercept = 0, linetype = "dotted", linewidth = 0.8, color = "gray40")

  if (plot_type == "column") {
    p <- p +
      stat_summary(fun = mean, geom = "col", color = "black",
                   linewidth = 0.6, width = 0.65, alpha = 0.85) +
      stat_summary(fun.data = err_fun, geom = "errorbar",
                   width = 0.25, linewidth = 0.8, color = "black") +
      geom_jitter(shape = 21, color = "black", stroke = 0.5,
                  size = 2.5, width = 0.12, alpha = 0.95)
  } else {
    p <- p +
      stat_summary(fun.data = err_fun, geom = "errorbar",
                   width = 0.25, linewidth = 0.8, color = "black") +
      stat_summary(fun = mean, geom = "errorbar",
                   fun.min = mean, fun.max = mean,
                   width = 0.4, linewidth = 1.2, color = "black") +
      geom_jitter(shape = 21, color = "black", stroke = 0.5,
                  size = 2.5, width = 0.12, alpha = 0.95)
  }

  p <- p +
    scale_fill_manual(values = rep_len(OKABE_ITO, n_groups)) +
    facet_wrap(~ Gene, scales = "free_y", ncol = ncol) +
    labs(
      y = bquote(Log[2] ~ "(FC to " * .(control_group) * ")"),
      x = NULL
    ) +
    prism_theme(rotate_x = rotate_x) +
    theme(aspect.ratio = aspect_ratio)

  if (!is.null(y_min) || !is.null(y_max))
    p <- p + coord_cartesian(ylim = c(if (is.null(y_min)) NA else y_min,
                                       if (is.null(y_max)) NA else y_max))

  # Per-facet significance bars using manual positioning
  if (isTRUE(show_sig) && !is.null(stats_df) && nrow(stats_df) > 0) {
    sig <- stats_df[stats_df$Significant == "Yes" &
                      grepl(" vs ", stats_df$Comparison, fixed = TRUE), ]
    if (!is.null(sig_comparisons_all)) {
      keep <- paste0(sig$Gene, ": ", sig$Comparison) %in% sig_comparisons_all
      sig  <- sig[keep, ]
    }
    if (nrow(sig) > 0) {
      # Build per-gene bracket data with computed y positions
      rows <- do.call(rbind, lapply(split(sig, sig$Gene), function(gsub) {
        gene_data <- full_df[full_df$Gene == as.character(gsub$Gene[1]), ]
        y_hi  <- max(gene_data$log2_fold_change, na.rm = TRUE)
        y_lo  <- min(gene_data$log2_fold_change, na.rm = TRUE)
        rng   <- max(y_hi - y_lo, 1)
        step  <- 0.16 * rng
        pairs <- strsplit(gsub$Comparison, " vs ", fixed = TRUE)
        data.frame(
          Gene         = gsub$Gene,
          start        = sapply(pairs, `[`, 1),
          end          = sapply(pairs, `[`, 2),
          y_position   = y_hi + step * seq_len(nrow(gsub)),
          annotations  = format_pvalue(gsub$p_value),
          stringsAsFactors = FALSE
        )
      }))
      p <- p + geom_signif(
        data    = rows,
        aes(xmin = start, xmax = end, annotations = annotations,
            y_position = y_position),
        manual           = TRUE,
        textsize         = 3.5,
        fontface         = "bold",
        vjust            = -0.4,
        tip_length       = 0.02,
        size             = 0.7,
        color            = "black",
        inherit.aes      = FALSE
      )
    }
  }

  p
}

# --- Paired / repeated-measures before-after plot --------------------------
make_paired_plot <- function(full_df, stats_df,
                             error_type    = "SD",
                             y_min         = NULL, y_max = NULL,
                             show_sig      = TRUE,
                             sig_comparisons_all = NULL,
                             control_group = "Control",
                             rotate_x      = FALSE,
                             aspect_ratio  = 1,
                             ncol          = 3) {
  full_df$Gene   <- factor(full_df$Gene,   levels = unique(full_df$Gene))
  full_df$Group  <- factor(full_df$Group,  levels = unique(full_df$Group))
  full_df$Sample <- factor(full_df$Sample, levels = unique(full_df$Sample))

  n_samples <- nlevels(full_df$Sample)
  palette   <- rep_len(OKABE_ITO, n_samples)

  p <- ggplot(full_df, aes(x = Group, y = log2_fold_change,
                           group = Sample, color = Sample)) +
    geom_hline(yintercept = 0, linetype = "dotted",
               linewidth = 0.8, color = "gray40") +
    geom_line(linewidth = 0.6, alpha = 0.7) +
    geom_point(size = 3.5, alpha = 0.95) +
    scale_color_manual(values = palette) +
    facet_wrap(~ Gene, scales = "free_y", ncol = ncol) +
    labs(
      y     = bquote(Log[2] ~ "(FC to " * .(control_group) * ")"),
      x     = NULL,
      color = "Sample"
    ) +
    prism_theme(rotate_x = rotate_x) +
    theme(
      aspect.ratio    = aspect_ratio,
      legend.position = "right",
      legend.title    = element_text(face = "bold")
    )

  if (!is.null(y_min) || !is.null(y_max))
    p <- p + coord_cartesian(ylim = c(if (is.null(y_min)) NA else y_min,
                                       if (is.null(y_max)) NA else y_max))

  if (isTRUE(show_sig) && !is.null(stats_df) && nrow(stats_df) > 0) {
    sig <- stats_df[stats_df$Significant == "Yes" &
                      grepl(" vs ", stats_df$Comparison, fixed = TRUE), ]
    if (!is.null(sig_comparisons_all)) {
      keep <- paste0(sig$Gene, ": ", sig$Comparison) %in% sig_comparisons_all
      sig  <- sig[keep, ]
    }
    if (nrow(sig) > 0) {
      rows <- do.call(rbind, lapply(split(sig, sig$Gene), function(gsub) {
        gene_data <- full_df[full_df$Gene == as.character(gsub$Gene[1]), ]
        y_hi  <- max(gene_data$log2_fold_change, na.rm = TRUE)
        y_lo  <- min(gene_data$log2_fold_change, na.rm = TRUE)
        rng   <- max(y_hi - y_lo, 1)
        step  <- 0.16 * rng
        pairs <- strsplit(gsub$Comparison, " vs ", fixed = TRUE)
        data.frame(
          Gene        = gsub$Gene,
          start       = sapply(pairs, `[`, 1),
          end         = sapply(pairs, `[`, 2),
          y_position  = y_hi + step * seq_len(nrow(gsub)),
          annotations = format_pvalue(gsub$p_value),
          stringsAsFactors = FALSE
        )
      }))
      p <- p + geom_signif(
        data    = rows,
        aes(xmin = start, xmax = end, annotations = annotations,
            y_position = y_position),
        manual      = TRUE,
        textsize    = 3.5,
        fontface    = "bold",
        vjust       = -0.4,
        tip_length  = 0.02,
        size        = 0.7,
        color       = "black",
        inherit.aes = FALSE
      )
    }
  }

  p
}
