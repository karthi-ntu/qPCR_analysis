library(ggplot2)
library(ggsignif)

# Okabe-Ito colourblind-safe palette (grey instead of black as the 8th colour
# so that black-outlined points stay visible).
OKABE_ITO <- c(
  "#E69F00", "#56B4E9", "#009E73", "#F0E442",
  "#0072B2", "#D55E00", "#CC79A7", "#999999"
)

DEFAULT_FONT <- "Arial"

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
  ifelse(is.na(p), "p = NA",
         ifelse(p < 0.0001, "p < 0.0001", sprintf("p = %.4f", p)))
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

# Error-bar summary function factory. Supports SD, SEM, and CI95 (95% CI via
# the t-distribution, which matches Prism's default).
make_err_fun <- function(error_type) {
  if (error_type == "SEM") {
    function(x) {
      x <- x[!is.na(x)]
      m <- mean(x); s <- if (length(x) > 1) sd(x) / sqrt(length(x)) else 0
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
      m <- mean(x); s <- if (length(x) > 1) sd(x) else 0
      data.frame(y = m, ymin = m - s, ymax = m + s)
    }
  }
}

# Font family for the plot; "" means the device default.
resolve_font <- function(font_family) {
  if (is.null(font_family) || is.na(font_family)) return("")
  font_family
}

# Builds nested-axis labels for the grouped (2-factor) layout: the inner
# factor (subgroup, under each dot) and the outer factor (group, centered
# under each dodge cluster). Returns a list of ggplot layers to add.
# Uses y = -Inf with vjust/clip=off to place labels below the axis.
#
# The caller must set theme(plot.margin = margin(..., bottom = >= the value
# in attr(layers, "bottom_margin_pt"))) and use coord_cartesian(clip = "off").
nested_axis_labels <- function(x_levels, col_levels, dodge_w,
                               text_sizes = DEFAULT_TEXT_SIZES,
                               font_family = "") {
  ts <- modifyList(DEFAULT_TEXT_SIZES, as.list(text_sizes))
  n_col <- length(col_levels)
  inner_size <- ts$axis_text * 0.78
  outer_size <- ts$axis_text * 1.15
  inner_vjust <- 1.8
  outer_vjust <- 4.0
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
  fam <- resolve_font(font_family)
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
  attr(layers, "bottom_margin_pt") <- ceiling(outer_vjust * outer_size + 8)
  layers
}

prism_theme <- function(rotate_x    = FALSE,
                        text_sizes  = DEFAULT_TEXT_SIZES,
                        font_family = "") {
  ts <- modifyList(DEFAULT_TEXT_SIZES, as.list(text_sizes))
  base_size <- as.numeric(ts$axis_text)
  # base_family sets the family on the root text element without wiping its
  # other properties (size, colour, ...).
  base <- theme_classic(base_size = base_size, base_family = resolve_font(font_family))
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

# Helper: merge user's per-level colour override with the default Okabe-Ito
# palette. `override` is a named character vector (names are factor levels,
# values are hex codes). Missing levels fall back to the palette.
resolve_fill <- function(override, levels) {
  default <- setNames(rep_len(OKABE_ITO, length(levels)), levels)
  if (is.null(override) || length(override) == 0) return(default)
  ok <- grepl("^#[0-9A-Fa-f]{3}([0-9A-Fa-f]{3})?$", override) & names(override) %in% levels
  override <- override[ok]
  if (length(override) == 0) return(default)
  default[names(override)] <- override
  default
}

# y-axis variable, reference line and label for the chosen scale.
y_scale_spec <- function(y_scale, control_group) {
  if (identical(y_scale, "linear")) {
    list(var = "fold_change", ref = 1,
         label = bquote("FC to" ~ .(control_group)))
  } else {
    list(var = "log2_fold_change", ref = 0,
         label = bquote(Log[2] * "FC to" ~ .(control_group)))
  }
}

# Select which comparisons to draw as brackets.
#   sig_comparisons = NULL          -> all significant comparisons
#   sig_comparisons = character(0)  -> none
#   otherwise                       -> exactly those listed ("Gene: A vs B" or "A vs B")
select_bracket_rows <- function(stats_df, sig_comparisons, gene = NULL) {
  if (is.null(stats_df) || nrow(stats_df) == 0) return(NULL)
  sig <- stats_df[grepl(" vs ", stats_df$Comparison, fixed = TRUE), , drop = FALSE]
  if (!is.null(gene)) sig <- sig[sig$Gene == gene, , drop = FALSE]
  if (is.null(sig_comparisons)) {
    sig <- sig[sig$Significant == "Yes", , drop = FALSE]
  } else {
    keys_full <- paste0(sig$Gene, ": ", sig$Comparison)
    sig <- sig[keys_full %in% sig_comparisons | sig$Comparison %in% sig_comparisons, , drop = FALSE]
  }
  if (nrow(sig) == 0) NULL else sig
}

# Top of the tallest error bar (or point) per gene, used as the base height
# for significance brackets so they never overlap the bars.
bracket_base <- function(gene_data, yvar, group_key, err_fun) {
  vals <- gene_data[[yvar]]
  top  <- max(vals, na.rm = TRUE)
  bottom <- min(vals, na.rm = TRUE)
  err_top <- tapply(vals, group_key, function(v) if (length(v[!is.na(v)])) err_fun(v)$ymax else NA)
  top <- max(top, err_top, na.rm = TRUE)
  list(hi = top, lo = bottom)
}

# --- Core builder ----------------------------------------------------------
# Draws one or several genes. `facet = TRUE` gives the combined grid
# (facet_wrap by Gene); `facet = FALSE` expects a single gene and uses its
# name as the plot title.
build_expression_plot <- function(full_df, stats_df,
                                  facet         = TRUE,
                                  error_type    = "SD",
                                  y_min         = NULL, y_max = NULL,
                                  plot_type     = "scatter",
                                  show_sig      = TRUE,
                                  sig_comparisons = NULL,
                                  control_group = "Control",
                                  rotate_x      = FALSE,
                                  aspect_ratio  = 1,
                                  ncol          = 3,
                                  has_subgroup  = FALSE,
                                  fill_override = NULL,
                                  text_sizes    = DEFAULT_TEXT_SIZES,
                                  font_family   = DEFAULT_FONT,
                                  sig_format    = "exact",
                                  x_axis_var    = "Group",
                                  y_scale       = "log2",
                                  point_size    = 3) {
  full_df$Gene  <- factor(full_df$Gene,  levels = unique(full_df$Gene))
  full_df$Group <- factor(full_df$Group, levels = unique(full_df$Group))
  err_fun <- make_err_fun(error_type)
  ys      <- y_scale_spec(y_scale, control_group)
  yvar    <- ys$var
  fam     <- resolve_font(font_family)
  ts      <- modifyList(DEFAULT_TEXT_SIZES, as.list(text_sizes))

  use_sub <- isTRUE(has_subgroup) && "Subgroup" %in% names(full_df) &&
             length(unique(full_df$Subgroup)) >= 2
  if (use_sub) full_df$Subgroup <- factor(full_df$Subgroup, levels = unique(full_df$Subgroup))

  dodge_w <- 0.85
  if (use_sub) {
    if (identical(x_axis_var, "Subgroup")) {
      full_df$XVar <- full_df$Subgroup; full_df$ColVar <- full_df$Group
      x_name <- "Subgroup"; col_name <- "Group"
    } else {
      full_df$XVar <- full_df$Group;    full_df$ColVar <- full_df$Subgroup
      x_name <- "Group";    col_name <- "Subgroup"
    }
    x_levels   <- levels(full_df$XVar)
    col_levels <- levels(full_df$ColVar)
    p <- ggplot(full_df, aes(x = XVar, y = .data[[yvar]], fill = ColVar))
    pd <- position_dodge(width = dodge_w)
    pj <- position_jitterdodge(jitter.width = 0.12, dodge.width = dodge_w, seed = 1)
    bar_w <- dodge_w * 0.85
  } else {
    p <- ggplot(full_df, aes(x = Group, y = .data[[yvar]], fill = Group))
    pd <- position_identity()
    pj <- position_jitter(width = 0.12, height = 0, seed = 1)
    bar_w <- 0.65
    col_levels <- levels(full_df$Group)
  }

  p <- p + geom_hline(yintercept = ys$ref, linetype = "dotted", linewidth = 0.8, color = "gray40")

  if (plot_type == "column") {
    p <- p +
      stat_summary(fun = mean, geom = "col", color = "black",
                   linewidth = 0.6, width = bar_w, alpha = 0.85, position = pd) +
      stat_summary(fun.data = err_fun, geom = "errorbar",
                   width = 0.25, linewidth = 0.8, color = "black", position = pd) +
      geom_point(shape = 21, color = "black", stroke = 0.5,
                 size = point_size, alpha = 0.95, position = pj)
  } else {
    p <- p +
      stat_summary(fun.data = err_fun, geom = "errorbar",
                   width = 0.25, linewidth = 0.8, color = "black", position = pd) +
      stat_summary(fun = mean, geom = "errorbar", fun.min = mean, fun.max = mean,
                   width = 0.4, linewidth = 1.2, color = "black", position = pd) +
      geom_point(shape = 21, color = "black", stroke = 0.5,
                 size = point_size, alpha = 0.95, position = pj)
  }

  p <- p +
    scale_fill_manual(values = resolve_fill(fill_override, col_levels)) +
    labs(y = ys$label, x = NULL, fill = if (use_sub) col_name else NULL) +
    prism_theme(rotate_x = rotate_x, text_sizes = text_sizes, font_family = fam) +
    theme(aspect.ratio = aspect_ratio)

  if (facet) {
    p <- p + facet_wrap(~ Gene, scales = "free_y", ncol = ncol)
  } else {
    p <- p + labs(title = as.character(full_df$Gene[1]))
  }

  if (use_sub) {
    nal <- nested_axis_labels(x_levels, col_levels, dodge_w,
                              text_sizes = text_sizes, font_family = fam)
    p <- p + nal + theme(
      legend.position = "top",
      legend.title    = element_text(face = "bold"),
      axis.text.x     = element_blank(),
      axis.ticks.x    = element_blank(),
      plot.margin     = margin(5.5, 5.5, attr(nal, "bottom_margin_pt"), 5.5)
    )
  }

  p <- p + coord_cartesian(
    ylim = c(if (is.null(y_min)) NA else y_min, if (is.null(y_max)) NA else y_max),
    clip = "off"
  )

  # --- Significance brackets ----------------------------------------------
  if (isTRUE(show_sig)) {
    sig <- select_bracket_rows(stats_df, sig_comparisons)
    sig <- if (is.null(sig)) NULL else sig[sig$Gene %in% levels(full_df$Gene), , drop = FALSE]
    if (!is.null(sig) && nrow(sig) > 0) {
      sig_fmt <- pick_sig_formatter(sig_format)
      rows <- do.call(rbind, lapply(split(sig, sig$Gene, drop = TRUE), function(gs) {
        gene_data <- full_df[full_df$Gene == as.character(gs$Gene[1]), ]
        if (nrow(gene_data) == 0) return(NULL)
        key  <- if (use_sub) paste(gene_data$XVar, gene_data$ColVar) else gene_data$Group
        base <- bracket_base(gene_data, yvar, key, err_fun)
        # Leave room for the label text above each bracket line: exact p
        # strings are wider and taller than stars.
        step_frac <- if (identical(sig_format, "stars")) 0.13 else 0.2
        step <- step_frac * max(base$hi - base$lo, if (yvar == "fold_change") 0.5 else 1)
        pairs <- strsplit(gs$Comparison, " vs ", fixed = TRUE)
        pairs <- lapply(pairs, trimws)
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
          keep <- vapply(pairs, function(pp) length(pp) == 2 &&
                           grepl("|", pp[1], fixed = TRUE) && grepl("|", pp[2], fixed = TRUE),
                         logical(1))
          if (!any(keep)) return(NULL)
          pairs <- pairs[keep]; gs <- gs[keep, , drop = FALSE]
          xs <- t(vapply(pairs, function(pp) {
            c1 <- parse_cell(pp[1]); c2 <- parse_cell(pp[2])
            x1 <- cell_pos(c1[[x_name]], c1[[col_name]])
            x2 <- cell_pos(c2[[x_name]], c2[[col_name]])
            c(min(x1, x2), max(x1, x2))
          }, numeric(2)))
          ok <- !is.na(xs[, 1]) & !is.na(xs[, 2])
          if (!any(ok)) return(NULL)
          data.frame(Gene = gs$Gene[ok], xmin = xs[ok, 1], xmax = xs[ok, 2],
                     y_position = base$hi + step * seq_len(sum(ok)),
                     annotations = sig_fmt(gs$p_value[ok]), stringsAsFactors = FALSE)
        } else {
          keep <- vapply(pairs, function(pp) length(pp) == 2 &&
                           all(pp %in% levels(full_df$Group)), logical(1))
          if (!any(keep)) return(NULL)
          pairs <- pairs[keep]; gs <- gs[keep, , drop = FALSE]
          xs <- t(vapply(pairs, function(pp) {
            pos <- match(pp, levels(full_df$Group)); c(min(pos), max(pos))
          }, numeric(2)))
          data.frame(Gene = gs$Gene, xmin = xs[, 1], xmax = xs[, 2],
                     y_position = base$hi + step * seq_len(nrow(gs)),
                     annotations = sig_fmt(gs$p_value), stringsAsFactors = FALSE)
        }
      }))
      if (!is.null(rows) && nrow(rows) > 0) {
        rows$Gene <- factor(rows$Gene, levels = levels(full_df$Gene))
        rows$grp  <- seq_len(nrow(rows))
        p <- p + suppressWarnings(geom_signif(
          data = rows,
          aes(xmin = xmin, xmax = xmax, annotations = annotations,
              y_position = y_position, group = grp),
          manual      = TRUE,
          textsize    = as.numeric(ts$sig_bar),
          fontface    = "bold",
          family      = fam,
          vjust       = -0.4,
          tip_length  = 0.02,
          size        = 0.7,
          color       = "black",
          inherit.aes = FALSE
        ))
      }
    }
  }
  p
}

# Single-gene plot (title = gene name).
make_barplot <- function(full_df, stats_df, gene, ..., sig_comparisons = NULL) {
  sub <- full_df[full_df$Gene == gene, , drop = FALSE]
  build_expression_plot(sub, stats_df, facet = FALSE,
                        sig_comparisons = sig_comparisons, ...)
}

# Combined multi-gene figure using facet_wrap.
make_combined_plot <- function(full_df, stats_df, ..., sig_comparisons_all = NULL) {
  build_expression_plot(full_df, stats_df, facet = TRUE,
                        sig_comparisons = sig_comparisons_all, ...)
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
                             font_family   = DEFAULT_FONT,
                             sig_format    = "exact",
                             y_scale       = "log2") {
  full_df$Gene   <- factor(full_df$Gene,   levels = unique(full_df$Gene))
  full_df$Group  <- factor(full_df$Group,  levels = unique(full_df$Group))
  full_df$Sample <- factor(full_df$Sample, levels = unique(full_df$Sample))
  err_fun <- make_err_fun(error_type)
  ys   <- y_scale_spec(y_scale, control_group)
  yvar <- ys$var
  fam  <- resolve_font(font_family)
  ts   <- modifyList(DEFAULT_TEXT_SIZES, as.list(text_sizes))

  n_samples <- nlevels(full_df$Sample)
  palette   <- if (n_samples <= length(OKABE_ITO)) OKABE_ITO[seq_len(n_samples)]
               else grDevices::hcl.colors(n_samples, "Dark 3")

  p <- ggplot(full_df, aes(x = Group, y = .data[[yvar]], group = Sample, color = Sample)) +
    geom_hline(yintercept = ys$ref, linetype = "dotted", linewidth = 0.8, color = "gray40") +
    geom_line(linewidth = 0.6, alpha = 0.7) +
    geom_point(size = 3.5, alpha = 0.95) +
    stat_summary(aes(group = 1), fun.data = err_fun, geom = "errorbar",
                 width = 0.18, linewidth = 0.8, color = "black") +
    stat_summary(aes(group = 1), fun = mean, geom = "errorbar", fun.min = mean, fun.max = mean,
                 width = 0.3, linewidth = 1.2, color = "black") +
    scale_color_manual(values = palette) +
    facet_wrap(~ Gene, scales = "free_y", ncol = ncol) +
    labs(y = ys$label, x = NULL, color = "Sample") +
    prism_theme(rotate_x = rotate_x, text_sizes = text_sizes, font_family = fam) +
    theme(aspect.ratio = aspect_ratio, legend.position = "right",
          legend.title = element_text(face = "bold")) +
    coord_cartesian(ylim = c(if (is.null(y_min)) NA else y_min,
                             if (is.null(y_max)) NA else y_max),
                    clip = "off")

  if (isTRUE(show_sig)) {
    sig <- select_bracket_rows(stats_df, sig_comparisons_all)
    if (!is.null(sig)) {
      sig_fmt <- pick_sig_formatter(sig_format)
      rows <- do.call(rbind, lapply(split(sig, sig$Gene, drop = TRUE), function(gs) {
        gene_data <- full_df[full_df$Gene == as.character(gs$Gene[1]), ]
        if (nrow(gene_data) == 0) return(NULL)
        base  <- bracket_base(gene_data, yvar, gene_data$Group, err_fun)
        step  <- (if (identical(sig_format, "stars")) 0.13 else 0.2) * max(base$hi - base$lo, 1)
        pairs <- lapply(strsplit(gs$Comparison, " vs ", fixed = TRUE), trimws)
        keep  <- vapply(pairs, function(pp) length(pp) == 2 &&
                          all(pp %in% levels(full_df$Group)), logical(1))
        if (!any(keep)) return(NULL)
        pairs <- pairs[keep]; gs <- gs[keep, , drop = FALSE]
        xs <- t(vapply(pairs, function(pp) {
          pos <- match(pp, levels(full_df$Group)); c(min(pos), max(pos))
        }, numeric(2)))
        data.frame(Gene = gs$Gene, xmin = xs[, 1], xmax = xs[, 2],
                   y_position = base$hi + step * seq_len(nrow(gs)),
                   annotations = sig_fmt(gs$p_value), stringsAsFactors = FALSE)
      }))
      if (!is.null(rows) && nrow(rows) > 0) {
        rows$Gene <- factor(rows$Gene, levels = levels(full_df$Gene))
        rows$grp  <- seq_len(nrow(rows))
        p <- p + suppressWarnings(geom_signif(
          data = rows,
          aes(xmin = xmin, xmax = xmax, annotations = annotations,
              y_position = y_position, group = grp),
          manual = TRUE, textsize = as.numeric(ts$sig_bar), fontface = "bold",
          family = fam, vjust = -0.4, tip_length = 0.02, size = 0.7,
          color = "black", inherit.aes = FALSE
        ))
      }
    }
  }
  p
}

# --- Saving ----------------------------------------------------------------
# Save a ggplot to file in the requested format. PDF uses cairo_pdf when
# available so that system fonts such as Arial embed correctly; otherwise the
# base pdf device with a matching standard family is used.
save_plot_file <- function(file, plot, format, width_mm, height_mm, dpi = 300,
                           font_family = DEFAULT_FONT) {
  format <- tolower(format)
  args <- list(filename = file, plot = plot, width = width_mm, height = height_mm,
               units = "mm", bg = "white", limitsize = FALSE)
  if (format == "png") {
    do.call(ggsave, c(args, list(device = "png", dpi = dpi)))
  } else if (format == "tiff") {
    do.call(ggsave, c(args, list(device = "tiff", dpi = dpi, compression = "lzw")))
  } else if (format == "svg") {
    dev <- if (requireNamespace("svglite", quietly = TRUE)) svglite::svglite
           else grDevices::svg
    do.call(ggsave, c(args, list(device = dev)))
  } else if (format == "pdf") {
    has_cairo <- isTRUE(tryCatch(capabilities("cairo"), error = function(e) FALSE))
    if (has_cairo) {
      ok <- tryCatch({ do.call(ggsave, c(args, list(device = grDevices::cairo_pdf))); TRUE },
                     error = function(e) FALSE)
      if (ok) return(invisible(file))
    }
    fam <- resolve_font(font_family)
    pdf_family <- if (fam %in% c("", "sans", "Arial", "Helvetica")) "Helvetica"
                  else if (fam %in% c("serif", "Times New Roman", "Times")) "Times"
                  else if (fam %in% c("mono", "Courier New", "Courier")) "Courier"
                  else "Helvetica"
    # Remap the plot's font to a family the base pdf device knows.
    plot <- plot + theme(text = element_text(family = pdf_family))
    args$plot <- plot
    do.call(ggsave, c(args, list(device = "pdf", family = pdf_family)))
  } else {
    stop("Unknown format: ", format)
  }
  invisible(file)
}
