library(ggplot2)
library(ggsignif)

OKABE_ITO <- c(
  "#E69F00", "#56B4E9", "#009E73", "#F0E442",
  "#0072B2", "#D55E00", "#CC79A7", "#000000"
)

# Default text sizes used across all plots. Overridable via the 'text_sizes'
# argument to prism_theme().
DEFAULT_TEXT_SIZES <- list(
  title      = 18,
  axis_title = 14,
  axis_text  = 12,
  legend     = 12,
  facet      = 14,
  sig_bar    = 3.5
)

format_pvalue <- function(p) {
  ifelse(p < 0.0001, "p < 0.0001", sprintf("p = %.4f", p))
}

# GraphPad-Prism-style significance stars:
#   ns  p >= 0.05
#   *    p < 0.05
#   **   p < 0.01
#   ***  p < 0.001
#   **** p < 0.0001
format_pvalue_stars <- function(p) {
  out <- rep("ns", length(p))
  out[!is.na(p) & p < 0.05]   <- "*"
  out[!is.na(p) & p < 0.01]   <- "**"
  out[!is.na(p) & p < 0.001]  <- "***"
  out[!is.na(p) & p < 0.0001] <- "****"
  out[is.na(p)] <- "NA"
  out
}

# Pick the annotation formatter based on user preference.
pick_sig_formatter <- function(sig_format = "exact") {
  if (identical(sig_format, "stars")) format_pvalue_stars else format_pvalue
}

# Error-bar summary function factory. Supports SD, SEM, and CI95 (95 % CI
# via the t-distribution — matches Prism's default).
make_err_fun <- function(error_type) {
  if (error_type == "SEM") {
    function(x) {
      x <- x[!is.na(x)]
      m <- mean(x); s <- sd(x) / sqrt(length(x))
      data.frame(y = m, ymin = m - s, ymax = m + s)
    }
  } else if (error_type == "CI95") {
    function(x) {
      x <- x[!is.na(x)]
      n <- length(x)
      if (n < 2) return(data.frame(y = mean(x), ymin = mean(x), ymax = mean(x)))
      m <- mean(x); se <- sd(x) / sqrt(n); t_crit <- qt(0.975, df = n - 1)
      data.frame(y = m, ymin = m - t_crit * se, ymax = m + t_crit * se)
    }
  } else {
    function(x) {
      x <- x[!is.na(x)]
      m <- mean(x); s <- sd(x)
      data.frame(y = m, ymin = m - s, ymax = m + s)
    }
  }
}

# Builds nested-axis labels for the grouped (2-factor) layout: the inner
# factor (subgroup, under each dot) and the outer factor (group, centered
# under each dodge cluster). Returns a list of ggplot layers to add.
# Uses y = -Inf with vjust/clip=off to place labels below the axis.
#
# Note: for this to render, the caller MUST set
#   theme(plot.margin = margin(..., bottom = >= 60 pt))
# and use coord_cartesian(clip = "off"). The helper also returns the
# minimum bottom margin needed so callers can query it via attr().
nested_axis_labels <- function(x_levels, col_levels, dodge_w,
                               text_sizes = DEFAULT_TEXT_SIZES,
                               font_family = "") {
  ts <- modifyList(DEFAULT_TEXT_SIZES, as.list(text_sizes))
  n_col <- length(col_levels)
  inner_size <- ts$axis_text * 0.82           # pt — smaller to avoid overlap
  outer_size <- ts$axis_text * 1.15           # pt
  # vjust is in units of text height; position labels just below the axis.
  inner_vjust <- 1.8
  outer_vjust <- 4.0
  # Inner: one label per (outer, inner) cell at the dodged x position.
  inner_df <- do.call(rbind, lapply(seq_along(x_levels), function(i) {
    data.frame(
      x     = i + (seq_along(col_levels) - (n_col + 1) / 2) * (dodge_w / n_col),
      label = col_levels,
      stringsAsFactors = FALSE
    )
  }))
  outer_df <- data.frame(x = seq_along(x_levels),
                         label = x_levels,
                         stringsAsFactors = FALSE)
  fam <- if (nzchar(font_family)) font_family else ""
  layers <- list(
    geom_text(data = inner_df,
              aes(x = x, label = label), y = -Inf,
              vjust = inner_vjust, size = inner_size / .pt,
              fontface = "bold", family = fam,
              inherit.aes = FALSE),
    geom_text(data = outer_df,
              aes(x = x, label = label), y = -Inf,
              vjust = outer_vjust, size = outer_size / .pt,
              fontface = "bold", family = fam,
              inherit.aes = FALSE)
  )
  # Tell the caller the bottom margin required to avoid clipping.
  # outer_vjust text-heights below the axis, scaled by text size + safety.
  needed_bottom_pt <- ceiling(outer_vjust * outer_size + 8)
  attr(layers, "bottom_margin_pt") <- needed_bottom_pt
  layers
}

prism_theme <- function(rotate_x    = FALSE,
                        text_sizes  = DEFAULT_TEXT_SIZES,
                        font_family = "") {
  ts <- modifyList(DEFAULT_TEXT_SIZES, as.list(text_sizes))
  base_size <- as.numeric(ts$axis_text)
  base <- theme_classic(base_size = base_size)
  if (nzchar(font_family)) base <- base %+replace% theme(text = element_text(family = font_family))
  base +
    theme(
      legend.position   = "none",
      legend.text       = element_text(size = ts$legend, face = "bold"),
      legend.title      = element_text(size = ts$legend, face = "bold"),
      axis.line         = element_line(linewidth = 1, color = "black"),
      axis.ticks        = element_line(linewidth = 0.8, color = "black"),
      axis.ticks.length = unit(0.2, "cm"),
      axis.title        = element_text(face = "bold", color = "black",
                                       size = ts$axis_title),
      axis.text.y       = element_text(face = "bold", color = "black",
                                       size = ts$axis_text),
      axis.text.x       = element_text(
        face  = "bold", color = "black",
        size  = ts$axis_text,
        angle = if (rotate_x) 45 else 0,
        hjust = if (rotate_x) 1  else 0.5
      ),
      plot.title        = element_text(face = "bold.italic", size = ts$title, hjust = 0.5),
      strip.text        = element_text(face = "bold.italic", size = ts$facet),
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
                         aspect_ratio  = 1,
                         has_subgroup  = FALSE,
                         fill_override = NULL,
                         text_sizes    = DEFAULT_TEXT_SIZES,
                         font_family   = "",
                         sig_format    = "exact",
                         x_axis_var    = "Group") {
  sub <- full_df[full_df$Gene == gene, ]
  sub$Group <- factor(sub$Group, levels = unique(sub$Group))
  err_fun  <- make_err_fun(error_type)

  use_sub <- isTRUE(has_subgroup) && "Subgroup" %in% names(sub) &&
             length(unique(sub$Subgroup)) >= 2
  if (use_sub) sub$Subgroup <- factor(sub$Subgroup, levels = unique(sub$Subgroup))

  dodge_w <- 0.85
  if (use_sub) {
    # Grouped (Prism) layout with nested axis labels.
    if (identical(x_axis_var, "Subgroup")) {
      x_fac <- sub$Subgroup; col_fac <- sub$Group
      x_name <- "Subgroup";  col_name <- "Group"
    } else {
      x_fac <- sub$Group;    col_fac <- sub$Subgroup
      x_name <- "Group";     col_name <- "Subgroup"
    }
    sub$XVar <- x_fac; sub$ColVar <- col_fac
    pd <- position_dodge(width = dodge_w)
    pj <- position_jitterdodge(jitter.width = 0.12, dodge.width = dodge_w, seed = 1)
    p <- ggplot(sub, aes(x = XVar, y = log2_fold_change, fill = ColVar)) +
      geom_hline(yintercept = 0, linetype = "dotted",
                 linewidth = 0.8, color = "gray40")
    if (plot_type == "column") {
      p <- p +
        stat_summary(fun = mean, geom = "col", color = "black",
                     linewidth = 0.6, width = dodge_w * 0.85, alpha = 0.85,
                     position = pd) +
        stat_summary(fun.data = err_fun, geom = "errorbar",
                     width = 0.25, linewidth = 0.8, color = "black",
                     position = pd) +
        geom_jitter(shape = 21, color = "black", stroke = 0.5,
                    size = 3, alpha = 0.95, position = pj)
    } else {
      p <- p +
        stat_summary(fun.data = err_fun, geom = "errorbar",
                     width = 0.25, linewidth = 0.8, color = "black",
                     position = pd) +
        stat_summary(fun = mean, geom = "errorbar",
                     fun.min = mean, fun.max = mean,
                     width = 0.4, linewidth = 1.2, color = "black",
                     position = pd) +
        geom_jitter(shape = 21, color = "black", stroke = 0.5,
                    size = 3, alpha = 0.95, position = pj)
    }
    fill_levels <- levels(col_fac)
    fill_values <- resolve_fill(fill_override, fill_levels)
  } else {
    # Single-factor (no subgroup)
    p <- ggplot(sub, aes(x = Group, y = log2_fold_change, fill = Group)) +
      geom_hline(yintercept = 0, linetype = "dotted",
                 linewidth = 0.8, color = "gray40")
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
    fill_levels <- levels(sub$Group)
    fill_values <- resolve_fill(fill_override, fill_levels)
  }

  fill_label <- if (use_sub) col_name else NULL

  p <- p +
    scale_fill_manual(values = fill_values) +
    labs(
      title = gene,
      y     = bquote(Log[2] ~ "(FC to " * .(control_group) * ")"),
      x     = NULL,
      fill  = fill_label
    ) +
    prism_theme(rotate_x    = rotate_x,
                text_sizes  = text_sizes,
                font_family = font_family) +
    theme(aspect.ratio = aspect_ratio)

  if (use_sub) {
    # Add nested Group + Subgroup axis labels below the axis line, hide the
    # default axis text, and make room with an expanded bottom margin sized
    # to actually fit both rows (clipping bug fix).
    nal <- nested_axis_labels(levels(x_fac), levels(col_fac), dodge_w,
                              text_sizes = text_sizes, font_family = font_family)
    p <- p + nal + theme(
      legend.position = "top",
      legend.title    = element_text(face = "bold"),
      axis.text.x     = element_blank(),
      axis.ticks.x    = element_blank(),
      plot.margin     = margin(5.5, 5.5,
                               attr(nal, "bottom_margin_pt"), 5.5)
    )
  }

  # Combine ylim + clip=off into one coord_cartesian (required by ggplot2).
  p <- p + coord_cartesian(
    ylim = c(if (is.null(y_min)) NA else y_min,
             if (is.null(y_max)) NA else y_max),
    clip = "off"
  )

  if (isTRUE(show_sig)) {
    sig <- stats_df[stats_df$Gene == gene & stats_df$Significant == "Yes", ]
    sig <- sig[grepl(" vs ", sig$Comparison, fixed = TRUE), ]
    if (!is.null(sig_comparisons))
      sig <- sig[sig$Comparison %in% sig_comparisons, ]
    if (nrow(sig) > 0) {
      sig_fmt  <- pick_sig_formatter(sig_format)
      sig_size <- as.numeric(modifyList(DEFAULT_TEXT_SIZES,
                                        as.list(text_sizes))$sig_bar)

      if (use_sub) {
        # Compute numeric x positions for the dodged layout.
        x_levels   <- levels(x_fac)
        col_levels <- levels(col_fac)
        n_col <- length(col_levels)
        cell_pos <- function(x_lab, col_lab) {
          xi <- match(x_lab, x_levels)
          ci <- match(col_lab, col_levels)
          xi + (ci - (n_col + 1) / 2) * (dodge_w / n_col)
        }
        parse_cell <- function(s) {
          bits <- trimws(strsplit(s, "|", fixed = TRUE)[[1]])
          list(Group = bits[1], Subgroup = bits[2])
        }
        pairs <- strsplit(sig$Comparison, " vs ", fixed = TRUE)
        valid <- vapply(pairs, function(p) length(p) == 2 &&
                          grepl("|", p[1], fixed = TRUE) &&
                          grepl("|", p[2], fixed = TRUE), logical(1))
        sig   <- sig[valid, ]; pairs <- pairs[valid]
        if (nrow(sig) > 0) {
          y_hi <- max(sub$log2_fold_change, na.rm = TRUE)
          y_lo <- min(sub$log2_fold_change, na.rm = TRUE)
          step <- 0.13 * max(y_hi - y_lo, 1)
          xmin <- numeric(nrow(sig)); xmax <- numeric(nrow(sig))
          for (i in seq_len(nrow(sig))) {
            c1 <- parse_cell(pairs[[i]][1]); c2 <- parse_cell(pairs[[i]][2])
            x1 <- cell_pos(c1[[x_name]], c1[[col_name]])
            x2 <- cell_pos(c2[[x_name]], c2[[col_name]])
            xmin[i] <- min(x1, x2); xmax[i] <- max(x1, x2)
          }
          rows <- data.frame(
            xmin       = xmin,
            xmax       = xmax,
            y_position = y_hi + step * seq_len(nrow(sig)),
            annotations = sig_fmt(sig$p_value),
            stringsAsFactors = FALSE
          )
          p <- p + geom_signif(
            data    = rows,
            aes(xmin = xmin, xmax = xmax,
                annotations = annotations, y_position = y_position),
            manual      = TRUE,
            textsize    = sig_size,
            fontface    = "bold",
            family      = if (nzchar(font_family)) font_family else "",
            vjust       = -0.4,
            tip_length  = 0.02,
            size        = 0.7,
            color       = "black",
            inherit.aes = FALSE
          )
        }
      } else {
        comparisons_list <- lapply(sig$Comparison, function(comp)
          trimws(strsplit(comp, " vs ", fixed = TRUE)[[1]]))
        p <- p + geom_signif(
          comparisons      = comparisons_list,
          annotations      = sig_fmt(sig$p_value),
          map_signif_level = FALSE,
          step_increase    = 0.13,
          textsize         = sig_size,
          fontface         = "bold",
          family           = if (nzchar(font_family)) font_family else "",
          vjust            = -0.6,
          tip_length       = 0.02,
          size             = 0.7,
          color            = "black"
        )
      }
    }
  }

  p
}

# Helper: merge user's per-level color override with the default Okabe-Ito
# palette. `override` is a named character vector (names are factor levels,
# values are hex codes). Missing levels fall back to the palette.
resolve_fill <- function(override, levels) {
  default <- setNames(rep_len(OKABE_ITO, length(levels)), levels)
  if (is.null(override) || length(override) == 0) return(default)
  # Accept only hex codes like #RGB or #RRGGBB; drop invalid entries.
  ok <- grepl("^#[0-9A-Fa-f]{3}([0-9A-Fa-f]{3})?$", override)
  override <- override[ok]
  if (length(override) == 0) return(default)
  default[names(override)] <- override
  default
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
                               ncol          = 3,
                               has_subgroup  = FALSE,
                               fill_override = NULL,
                               text_sizes    = DEFAULT_TEXT_SIZES,
                               font_family   = "",
                               sig_format    = "exact",
                               x_axis_var    = "Group") {
  full_df$Gene  <- factor(full_df$Gene,  levels = unique(full_df$Gene))
  full_df$Group <- factor(full_df$Group, levels = unique(full_df$Group))
  err_fun  <- make_err_fun(error_type)

  use_sub <- isTRUE(has_subgroup) && "Subgroup" %in% names(full_df) &&
             length(unique(full_df$Subgroup)) >= 2
  if (use_sub) full_df$Subgroup <- factor(full_df$Subgroup,
                                          levels = unique(full_df$Subgroup))

  dodge_w <- 0.85
  if (use_sub) {
    if (identical(x_axis_var, "Subgroup")) {
      full_df$XVar   <- full_df$Subgroup; full_df$ColVar <- full_df$Group
      x_levels   <- levels(full_df$Subgroup); col_levels <- levels(full_df$Group)
      x_name <- "Subgroup"; col_name <- "Group"
    } else {
      full_df$XVar   <- full_df$Group;    full_df$ColVar <- full_df$Subgroup
      x_levels   <- levels(full_df$Group);    col_levels <- levels(full_df$Subgroup)
      x_name <- "Group";    col_name <- "Subgroup"
    }
    n_fills <- length(col_levels)
    p <- ggplot(full_df, aes(x = XVar, y = log2_fold_change, fill = ColVar))
  } else {
    n_fills <- nlevels(full_df$Group)
    p <- ggplot(full_df, aes(x = Group, y = log2_fold_change, fill = Group))
  }
  p <- p +
    geom_hline(yintercept = 0, linetype = "dotted", linewidth = 0.8, color = "gray40")

  if (use_sub) {
    pd <- position_dodge(width = dodge_w)
    pj <- position_jitterdodge(jitter.width = 0.12, dodge.width = dodge_w, seed = 1)
    if (plot_type == "column") {
      p <- p +
        stat_summary(fun = mean, geom = "col", color = "black",
                     linewidth = 0.6, width = dodge_w * 0.85, alpha = 0.85,
                     position = pd) +
        stat_summary(fun.data = err_fun, geom = "errorbar",
                     width = 0.25, linewidth = 0.8, color = "black",
                     position = pd) +
        geom_jitter(shape = 21, color = "black", stroke = 0.5,
                    size = 2.5, alpha = 0.95, position = pj)
    } else {
      p <- p +
        stat_summary(fun.data = err_fun, geom = "errorbar",
                     width = 0.25, linewidth = 0.8, color = "black",
                     position = pd) +
        stat_summary(fun = mean, geom = "errorbar",
                     fun.min = mean, fun.max = mean,
                     width = 0.4, linewidth = 1.2, color = "black",
                     position = pd) +
        geom_jitter(shape = 21, color = "black", stroke = 0.5,
                    size = 2.5, alpha = 0.95, position = pj)
    }
  } else if (plot_type == "column") {
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

  fill_levels <- if (use_sub) col_levels
                 else levels(full_df$Group)
  fill_values <- resolve_fill(fill_override, fill_levels)

  fill_label <- if (use_sub) col_name else NULL

  p <- p +
    scale_fill_manual(values = fill_values) +
    facet_wrap(~ Gene, scales = "free_y", ncol = ncol) +
    labs(
      y    = bquote(Log[2] ~ "(FC to " * .(control_group) * ")"),
      x    = NULL,
      fill = fill_label
    ) +
    prism_theme(rotate_x    = rotate_x,
                text_sizes  = text_sizes,
                font_family = font_family) +
    theme(aspect.ratio = aspect_ratio)

  if (use_sub) {
    nal <- nested_axis_labels(x_levels, col_levels, dodge_w,
                              text_sizes = text_sizes, font_family = font_family)
    p <- p + nal + theme(
      legend.position = "top",
      legend.title    = element_text(face = "bold"),
      axis.text.x     = element_blank(),
      axis.ticks.x    = element_blank(),
      plot.margin     = margin(5.5, 5.5,
                               attr(nal, "bottom_margin_pt"), 5.5)
    )
  }

  # Combine ylim + clip="off" in one coord_cartesian call.
  p <- p + coord_cartesian(
    ylim = c(if (is.null(y_min)) NA else y_min,
             if (is.null(y_max)) NA else y_max),
    clip = "off"
  )

  # Per-facet significance bars using manual positioning
  if (isTRUE(show_sig) && !is.null(stats_df) && nrow(stats_df) > 0) {
    sig <- stats_df[stats_df$Significant == "Yes" &
                      grepl(" vs ", stats_df$Comparison, fixed = TRUE), ]
    if (!is.null(sig_comparisons_all)) {
      keep <- paste0(sig$Gene, ": ", sig$Comparison) %in% sig_comparisons_all
      sig  <- sig[keep, ]
    }
    if (nrow(sig) > 0) {
      sig_fmt <- pick_sig_formatter(sig_format)
      if (use_sub) {
        n_col <- length(col_levels)
        cell_pos <- function(x_lab, col_lab) {
          xi <- match(x_lab, x_levels); ci <- match(col_lab, col_levels)
          xi + (ci - (n_col + 1) / 2) * (dodge_w / n_col)
        }
        parse_cell <- function(s) {
          bits <- trimws(strsplit(s, "|", fixed = TRUE)[[1]])
          list(Group = bits[1], Subgroup = bits[2])
        }
        rows <- do.call(rbind, lapply(split(sig, sig$Gene), function(gsub) {
          gene_data <- full_df[full_df$Gene == as.character(gsub$Gene[1]), ]
          y_hi <- max(gene_data$log2_fold_change, na.rm = TRUE)
          y_lo <- min(gene_data$log2_fold_change, na.rm = TRUE)
          step <- 0.13 * max(y_hi - y_lo, 1)
          pairs <- strsplit(gsub$Comparison, " vs ", fixed = TRUE)
          keep  <- vapply(pairs, function(p) length(p) == 2 &&
                            grepl("|", p[1], fixed = TRUE) &&
                            grepl("|", p[2], fixed = TRUE), logical(1))
          if (!any(keep)) return(NULL)
          pairs <- pairs[keep]; gs2 <- gsub[keep, , drop = FALSE]
          xmin <- numeric(length(pairs)); xmax <- numeric(length(pairs))
          for (i in seq_along(pairs)) {
            c1 <- parse_cell(pairs[[i]][1]); c2 <- parse_cell(pairs[[i]][2])
            x1 <- cell_pos(c1[[x_name]], c1[[col_name]])
            x2 <- cell_pos(c2[[x_name]], c2[[col_name]])
            xmin[i] <- min(x1, x2); xmax[i] <- max(x1, x2)
          }
          data.frame(
            Gene        = gs2$Gene,
            xmin        = xmin,
            xmax        = xmax,
            y_position  = y_hi + step * seq_len(nrow(gs2)),
            annotations = sig_fmt(gs2$p_value),
            stringsAsFactors = FALSE
          )
        }))
        if (!is.null(rows) && nrow(rows) > 0) {
          p <- p + geom_signif(
            data    = rows,
            aes(xmin = xmin, xmax = xmax,
                annotations = annotations, y_position = y_position),
            manual      = TRUE,
            textsize    = as.numeric(modifyList(DEFAULT_TEXT_SIZES,
                                                as.list(text_sizes))$sig_bar),
            fontface    = "bold",
            family      = if (nzchar(font_family)) font_family else "",
            vjust       = -0.4,
            tip_length  = 0.02,
            size        = 0.7,
            color       = "black",
            inherit.aes = FALSE
          )
        }
      } else {
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
            annotations  = sig_fmt(gsub$p_value),
            stringsAsFactors = FALSE
          )
        }))
        p <- p + geom_signif(
          data    = rows,
          aes(xmin = start, xmax = end, annotations = annotations,
              y_position = y_position),
          manual           = TRUE,
          textsize         = as.numeric(modifyList(DEFAULT_TEXT_SIZES,
                                                   as.list(text_sizes))$sig_bar),
          fontface         = "bold",
          family           = if (nzchar(font_family)) font_family else "",
          vjust            = -0.4,
          tip_length       = 0.02,
          size             = 0.7,
          color            = "black",
          inherit.aes      = FALSE
        )
      }
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
                             ncol          = 3,
                             text_sizes    = DEFAULT_TEXT_SIZES,
                             font_family   = "",
                             sig_format    = "exact") {
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
    prism_theme(rotate_x    = rotate_x,
                text_sizes  = text_sizes,
                font_family = font_family) +
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
      sig_fmt <- pick_sig_formatter(sig_format)
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
          annotations = sig_fmt(gsub$p_value),
          stringsAsFactors = FALSE
        )
      }))
      p <- p + geom_signif(
        data    = rows,
        aes(xmin = start, xmax = end, annotations = annotations,
            y_position = y_position),
        manual      = TRUE,
        textsize    = as.numeric(modifyList(DEFAULT_TEXT_SIZES,
                                            as.list(text_sizes))$sig_bar),
        fontface    = "bold",
        family      = if (nzchar(font_family)) font_family else "",
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
