# qPCR analysis pure functions (no Shiny dependencies).
#
# All statistics are computed on delta Ct (Ct_target - Ct_reference), which is
# the log2-scale quantity. Effect sizes and confidence intervals are reported
# on the log2 fold-change scale (log2FC = -ddCt), so a positive difference for
# "A vs B" means A is expressed more highly than B.

ALPHA <- 0.05

# Column layout shared by every stats table the app produces.
STATS_COLS <- c("Gene", "Comparison", "Test", "Statistic", "df",
                "Diff_log2FC", "CI95_low", "CI95_high", "p_value", "Significant")

compute_delta_ct <- function(df) {
  df$delta_ct <- df$Ct_target - df$Ct_reference
  df
}

# --- Technical replicate averaging ----------------------------------------
# Collapses rows that share the same (Sample, Group, [Subgroup], Gene) to a
# single row by averaging Ct_target and Ct_reference independently. Original
# row order is preserved. attr(out, "tech_counts") gives the number of rows
# that went into each averaged row.
average_tech_reps <- function(df) {
  if (nrow(df) == 0) return(df)
  has_sub <- "Subgroup" %in% names(df)
  keys <- if (has_sub)
    paste(df$Sample, df$Group, df$Subgroup, df$Gene, sep = "\u001f")
  else
    paste(df$Sample, df$Group, df$Gene, sep = "\u001f")
  split_df <- split(df, factor(keys, levels = unique(keys)))
  rep_counts <- vapply(split_df, nrow, integer(1))
  out <- do.call(rbind, lapply(split_df, function(rows) {
    rows[1, "Ct_target"]    <- mean(rows$Ct_target,    na.rm = TRUE)
    rows[1, "Ct_reference"] <- mean(rows$Ct_reference, na.rm = TRUE)
    rows[1, , drop = FALSE]
  }))
  rownames(out) <- NULL
  attr(out, "tech_counts") <- rep_counts
  out
}

# Returns a string describing the replicate structure for methods text / UI.
replicate_summary <- function(df, rep_mode = "biological") {
  if (nrow(df) == 0) return("")
  has_sub <- "Subgroup" %in% names(df)
  keys <- if (has_sub)
    paste(df$Sample, df$Group, df$Subgroup, sep = "\u001f")
  else
    paste(df$Sample, df$Group, sep = "\u001f")
  per_group_keys <- if (has_sub) paste(df$Group, df$Subgroup, sep = " | ")
                    else df$Group
  bio_per_group <- tapply(keys, per_group_keys, function(x) length(unique(x)))
  bio_range <- range(bio_per_group, na.rm = TRUE)
  fmt_range <- function(r) if (r[1] == r[2]) as.character(r[1]) else paste0(r[1], " to ", r[2])
  bio_txt <- paste0("n = ", fmt_range(bio_range), " biological replicates per group")
  if (rep_mode == "technical") {
    tc <- attr(df, "tech_counts")
    if (!is.null(tc) && length(tc) > 0) {
      return(paste0(bio_txt, ", each the mean of ", fmt_range(range(tc)),
                    " technical replicates"))
    }
  }
  bio_txt
}

# --- Fold change -----------------------------------------------------------
# control_subgroup (optional): restrict the baseline to one Group x Subgroup
# cell. Errors with a clear message if the baseline has no data for a gene.
compute_fold_change <- function(df, control_group, control_subgroup = NULL) {
  if (!control_group %in% df$Group) {
    stop(sprintf("Control group '%s' not found in the Group column.", control_group))
  }
  use_sub <- !is.null(control_subgroup) &&
             nzchar(control_subgroup) &&
             "Subgroup" %in% names(df)
  if (use_sub && !any(df$Group == control_group & df$Subgroup == control_subgroup)) {
    stop(sprintf("No rows with Group '%s' and Subgroup '%s' to use as the baseline.",
                 control_group, control_subgroup))
  }
  genes <- unique(df$Gene)
  result_list <- lapply(genes, function(gene) {
    sub <- df[df$Gene == gene, ]
    idx <- if (use_sub) sub$Group == control_group & sub$Subgroup == control_subgroup
           else sub$Group == control_group
    ctrl_vals <- sub$delta_ct[idx]
    ctrl_vals <- ctrl_vals[!is.na(ctrl_vals)]
    if (length(ctrl_vals) == 0) {
      stop(sprintf("Gene '%s' has no usable Ct values in the control baseline.", gene))
    }
    control_mean <- mean(ctrl_vals)
    sub$delta_delta_ct    <- sub$delta_ct - control_mean
    sub$fold_change       <- 2^(-sub$delta_delta_ct)
    sub$log2_fold_change  <- -sub$delta_delta_ct
    sub
  })
  out <- do.call(rbind, result_list)
  rownames(out) <- NULL
  out
}

# --- Internal helpers ------------------------------------------------------
fmt_num <- function(x, digits = 3) {
  ifelse(is.na(x), NA_character_, as.character(signif(x, digits)))
}

stat_rows <- function(gene, comparison, test, statistic = NA_character_,
                      df = NA_character_, diff = NA_real_, ci_low = NA_real_,
                      ci_high = NA_real_, p) {
  n <- length(comparison)
  data.frame(
    Gene        = rep(gene, n),
    Comparison  = comparison,
    Test        = rep_len(test, n),
    Statistic   = rep_len(as.character(statistic), n),
    df          = rep_len(as.character(df), n),
    Diff_log2FC = rep_len(signif(diff, 4), n),
    CI95_low    = rep_len(signif(ci_low, 4), n),
    CI95_high   = rep_len(signif(ci_high, 4), n),
    p_value     = signif(p, 4),
    Significant = ifelse(!is.na(p) & p < ALPHA, "Yes", "No"),
    stringsAsFactors = FALSE
  )
}

empty_stats <- function() {
  stat_rows(character(0), character(0), character(0), p = numeric(0))
}

# Split TukeyHSD / pairwise rownames of the form "<levelA>-<levelB>" using
# the known level set, so level names that contain "-" are handled.
split_pair_name <- function(rn, levels) {
  for (L in levels) {
    suffix <- paste0("-", L)
    if (endsWith(rn, suffix)) {
      pre <- substr(rn, 1, nchar(rn) - nchar(suffix))
      if (pre %in% levels) return(c(pre, L))
    }
  }
  c(NA_character_, NA_character_)
}

# Cell label used for Group x Subgroup combinations. Must match the label the
# plot code parses (split on " | ").
cell_label <- function(group, subgroup) paste(group, subgroup, sep = " | ")

cell_factor <- function(sub) {
  lab <- cell_label(sub$Group, sub$Subgroup)
  factor(lab, levels = unique(lab))
}

# --- Main dispatcher -------------------------------------------------------
# Returns a data.frame with STATS_COLS. attr(, "notes") carries messages about
# settings that could not be applied (shown in the UI and used to keep the
# methods text honest).
run_stats <- function(df, paired = FALSE, test_type = "parametric",
                      has_subgroup = FALSE, control_group = NULL,
                      control_subgroup = NULL, var_equal = FALSE) {
  genes <- unique(df$Gene)
  notes <- character(0)
  vs_control    <- test_type %in% c("parametric_vs_control", "nonparametric_vs_control")
  nonparametric <- test_type %in% c("nonparametric", "nonparametric_vs_control")
  use_sub <- isTRUE(has_subgroup) && "Subgroup" %in% names(df) &&
             length(unique(df$Subgroup)) >= 2

  if (isTRUE(paired) && use_sub) {
    notes <- c(notes, "Paired analysis is not available for two-factor designs; an unpaired analysis was run.")
    paired <- FALSE
  }
  if (isTRUE(paired) && vs_control) {
    notes <- c(notes, "Paired analysis does not apply to the vs-control tests; an unpaired analysis was run.")
    paired <- FALSE
  }
  if (isTRUE(paired) && !use_sub && length(unique(df$Group)) > 2) {
    notes <- c(notes, "Paired analysis is only available for exactly 2 groups here; an unpaired analysis was run. Use the Repeated Measures tab for 3+ paired conditions.")
    paired <- FALSE
  }
  if (use_sub && nonparametric && !vs_control) {
    notes <- c(notes, "There is no nonparametric two-way ANOVA. Cells (Group x Subgroup) were compared with a Kruskal-Wallis test and pairwise Mann-Whitney tests (Holm).")
  }
  if (use_sub && vs_control && (is.null(control_subgroup) || !nzchar(control_subgroup))) {
    notes <- c(notes, "No control subgroup selected: every cell was compared against the pooled control group.")
  }

  results <- lapply(genes, function(gene) {
    sub <- df[df$Gene == gene, ]
    sub <- sub[!is.na(sub$delta_ct), ]
    if (use_sub) {
      if (vs_control) {
        if (is.null(control_group) || !control_group %in% sub$Group) return(NULL)
        return(run_vs_control_cells(sub, gene, control_group, control_subgroup, nonparametric))
      }
      if (nonparametric) return(run_cells_nonparametric(sub, gene))
      return(run_two_way(sub, gene))
    }
    groups   <- unique(sub$Group)
    n_groups <- length(groups)
    if (n_groups < 2) return(NULL)
    if (vs_control) {
      if (is.null(control_group) || !control_group %in% groups) return(NULL)
      run_vs_control(sub, gene, control_group, nonparametric)
    } else if (n_groups == 2) {
      run_two_group(sub, gene, groups, paired, nonparametric, var_equal)
    } else {
      run_multi_group(sub, gene, nonparametric)
    }
  })

  out <- do.call(rbind, Filter(Negate(is.null), results))
  if (is.null(out)) out <- empty_stats()
  rownames(out) <- NULL
  attr(out, "notes") <- notes
  out
}

# --- vs control (Dunnett-type) --------------------------------------------
# Parametric: pooled-variance t-tests (the one-way ANOVA error term, as in
# Dunnett's procedure) with Holm-Bonferroni correction across the k-1
# comparisons against the control. Nonparametric: Mann-Whitney U with Holm.
vs_control_core <- function(values, groups, control, nonparametric, gene, test_suffix = "") {
  groups <- as.character(groups)
  levels <- unique(groups)
  others <- setdiff(levels, control)
  if (length(others) == 0) return(NULL)
  ctrl_vals <- values[groups == control]
  if (nonparametric) {
    res <- lapply(others, function(g) {
      x <- values[groups == g]
      w <- tryCatch(suppressWarnings(wilcox.test(x, ctrl_vals)), error = function(e) NULL)
      if (is.null(w)) return(list(stat = NA_character_, p = NA_real_, diff = NA_real_))
      list(stat = paste0("W = ", w$statistic), p = w$p.value,
           diff = -(median(x) - median(ctrl_vals)))
    })
    p_adj <- p.adjust(vapply(res, `[[`, numeric(1), "p"), method = "holm")
    return(stat_rows(gene, paste(others, "vs", control),
                     paste0("Mann-Whitney U vs control (Holm)", test_suffix),
                     statistic = vapply(res, `[[`, character(1), "stat"),
                     diff = vapply(res, `[[`, numeric(1), "diff"),
                     p = p_adj))
  }
  # Pooled SD across all groups (ANOVA MSE).
  ns  <- tapply(values, groups, length)
  vars <- tapply(values, groups, var)
  df_res <- sum(ns) - length(ns)
  if (df_res < 1 || any(ns < 2)) return(NULL)
  s_pooled <- sqrt(sum((ns - 1) * vars, na.rm = TRUE) / df_res)
  means <- tapply(values, groups, mean)
  res <- lapply(others, function(g) {
    se <- s_pooled * sqrt(1 / ns[[g]] + 1 / ns[[control]])
    d  <- means[[g]] - means[[control]]              # delta Ct scale
    t  <- d / se
    p  <- 2 * pt(-abs(t), df_res)
    half <- qt(0.975, df_res) * se
    list(t = t, p = p, diff = -d, lo = -(d + half), hi = -(d - half))
  })
  p_adj <- p.adjust(vapply(res, `[[`, numeric(1), "p"), method = "holm")
  stat_rows(gene, paste(others, "vs", control),
            paste0("Pooled-variance t-test vs control (Holm)", test_suffix),
            statistic = paste0("t = ", fmt_num(vapply(res, `[[`, numeric(1), "t"))),
            df = df_res,
            diff = vapply(res, `[[`, numeric(1), "diff"),
            ci_low = vapply(res, `[[`, numeric(1), "lo"),
            ci_high = vapply(res, `[[`, numeric(1), "hi"),
            p = p_adj)
}

run_vs_control <- function(sub, gene, control_group, nonparametric = FALSE) {
  vs_control_core(sub$delta_ct, sub$Group, control_group, nonparametric, gene)
}

# Two-factor version: every Group | Subgroup cell is compared against the
# control cell (control_group | control_subgroup). Without a control subgroup
# the control group's cells are pooled into one baseline.
run_vs_control_cells <- function(sub, gene, control_group, control_subgroup, nonparametric) {
  has_ctrl_sub <- !is.null(control_subgroup) && nzchar(control_subgroup) &&
                  any(sub$Group == control_group & sub$Subgroup == control_subgroup)
  cells <- as.character(cell_factor(sub))
  if (has_ctrl_sub) {
    control <- cell_label(control_group, control_subgroup)
  } else {
    control <- control_group
    cells[sub$Group == control_group] <- control_group
  }
  vs_control_core(sub$delta_ct, cells, control, nonparametric, gene)
}

# --- Two groups ------------------------------------------------------------
run_two_group <- function(sub, gene, groups, paired, nonparametric, var_equal = FALSE) {
  g1_vals <- sub$delta_ct[sub$Group == groups[1]]
  g2_vals <- sub$delta_ct[sub$Group == groups[2]]
  use_paired <- isTRUE(paired)
  if (use_paired) {
    sub1 <- sub[sub$Group == groups[1], ]
    sub2 <- sub[sub$Group == groups[2], ]
    common <- intersect(sub1$Sample, sub2$Sample)
    dup1 <- any(duplicated(sub1$Sample)); dup2 <- any(duplicated(sub2$Sample))
    if (length(common) >= 2 && !dup1 && !dup2) {
      g1_vals <- sub1$delta_ct[match(common, sub1$Sample)]
      g2_vals <- sub2$delta_ct[match(common, sub2$Sample)]
    } else {
      use_paired <- FALSE
    }
  }
  comparison <- paste(groups[1], "vs", groups[2])
  if (nonparametric) {
    test_name <- if (use_paired) "Wilcoxon signed-rank test" else "Mann-Whitney U test"
    w <- tryCatch(suppressWarnings(wilcox.test(g1_vals, g2_vals, paired = use_paired)),
                  error = function(e) NULL)
    if (is.null(w)) return(stat_rows(gene, comparison, test_name, p = NA_real_))
    diff <- if (use_paired) -median(g1_vals - g2_vals) else -(median(g1_vals) - median(g2_vals))
    return(stat_rows(gene, comparison, test_name,
                     statistic = paste0(names(w$statistic), " = ", w$statistic),
                     diff = diff, p = w$p.value))
  }
  test_name <- if (use_paired) "Paired t-test"
               else if (isTRUE(var_equal)) "Student's t-test (equal variances)"
               else "Welch's t-test"
  tt <- tryCatch(t.test(g1_vals, g2_vals, paired = use_paired, var.equal = isTRUE(var_equal)),
                 error = function(e) NULL)
  if (is.null(tt)) return(stat_rows(gene, comparison, test_name, p = NA_real_))
  est <- if (use_paired) unname(tt$estimate) else unname(tt$estimate[1] - tt$estimate[2])
  stat_rows(gene, comparison, test_name,
            statistic = paste0("t = ", fmt_num(unname(tt$statistic))),
            df = fmt_num(unname(tt$parameter), 4),
            diff = -est, ci_low = -tt$conf.int[2], ci_high = -tt$conf.int[1],
            p = tt$p.value)
}

# --- Three or more groups --------------------------------------------------
run_multi_group <- function(sub, gene, nonparametric) {
  sub$Group <- factor(sub$Group, levels = unique(sub$Group))
  if (nonparametric) {
    kw <- tryCatch(kruskal.test(delta_ct ~ Group, data = sub), error = function(e) NULL)
    omnibus <- if (is.null(kw)) NULL else
      stat_rows(gene, "Kruskal-Wallis (omnibus)", "Kruskal-Wallis test",
                statistic = paste0("H = ", fmt_num(unname(kw$statistic))),
                df = unname(kw$parameter), p = kw$p.value)
    pw <- tryCatch(suppressWarnings(
      pairwise.wilcox.test(sub$delta_ct, sub$Group, p.adjust.method = "holm")$p.value),
      error = function(e) NULL)
    if (is.null(pw)) return(omnibus)
    comps <- c(); pvals <- c(); diffs <- c()
    meds <- tapply(sub$delta_ct, sub$Group, median)
    for (r in rownames(pw)) for (cc in colnames(pw)) {
      v <- pw[r, cc]
      if (!is.na(v)) {
        comps <- c(comps, paste(r, "vs", cc)); pvals <- c(pvals, v)
        diffs <- c(diffs, -(meds[[r]] - meds[[cc]]))
      }
    }
    return(rbind(omnibus, stat_rows(gene, comps, "Pairwise Mann-Whitney U (Holm)",
                                    diff = diffs, p = pvals)))
  }
  aov_fit <- tryCatch(aov(delta_ct ~ Group, data = sub), error = function(e) NULL)
  if (is.null(aov_fit)) return(NULL)
  s <- summary(aov_fit)[[1]]
  omnibus <- stat_rows(gene, "One-way ANOVA (omnibus)", "One-way ANOVA",
                       statistic = paste0("F = ", fmt_num(s[1, "F value"])),
                       df = paste(s[1, "Df"], s[2, "Df"], sep = ", "),
                       p = s[1, "Pr(>F)"])
  tukey <- tryCatch(TukeyHSD(aov_fit)$Group, error = function(e) NULL)
  if (is.null(tukey) || nrow(tukey) == 0) return(omnibus)
  pairs <- t(vapply(rownames(tukey), split_pair_name, character(2), levels = levels(sub$Group)))
  keep  <- !is.na(pairs[, 1])
  # TukeyHSD row "B-A" is B minus A on the delta Ct scale.
  rbind(omnibus,
        stat_rows(gene, paste(pairs[keep, 1], "vs", pairs[keep, 2]), "Tukey HSD",
                  diff = -tukey[keep, "diff"],
                  ci_low = -tukey[keep, "upr"], ci_high = -tukey[keep, "lwr"],
                  p = tukey[keep, "p adj"]))
}

# --- Two-way ANOVA (Type II sums of squares) ------------------------------
# Type II SS do not depend on the order of the factors, which matters for
# unbalanced designs (Prism and car::Anova default). Falls back gracefully if
# a term cannot be estimated.
run_two_way <- function(sub, gene) {
  sub$Group    <- factor(sub$Group,    levels = unique(sub$Group))
  sub$Subgroup <- factor(sub$Subgroup, levels = unique(sub$Subgroup))
  rss  <- function(f) sum(residuals(f)^2)
  fits <- tryCatch(list(
    g   = lm(delta_ct ~ Group,            data = sub),
    s   = lm(delta_ct ~ Subgroup,         data = sub),
    gs  = lm(delta_ct ~ Group + Subgroup, data = sub),
    full = lm(delta_ct ~ Group * Subgroup, data = sub)
  ), error = function(e) NULL)
  if (is.null(fits)) return(NULL)
  df_res <- df.residual(fits$full)
  main_df <- NULL
  if (df_res > 0) {
    mse <- rss(fits$full) / df_res
    term <- function(reduced, fuller) {
      d_df <- fuller$rank - reduced$rank
      if (d_df < 1) return(list(F = NA_real_, df = NA_integer_, p = NA_real_))
      ss <- rss(reduced) - rss(fuller)
      F  <- (ss / d_df) / mse
      list(F = F, df = d_df, p = pf(F, d_df, df_res, lower.tail = FALSE))
    }
    t_g <- term(fits$s,  fits$gs)
    t_s <- term(fits$g,  fits$gs)
    t_i <- term(fits$gs, fits$full)
    main_df <- stat_rows(
      gene,
      c("Group (main effect)", "Subgroup (main effect)", "Group x Subgroup (interaction)"),
      "Two-way ANOVA (Type II SS)",
      statistic = paste0("F = ", fmt_num(c(t_g$F, t_s$F, t_i$F))),
      df = paste(c(t_g$df, t_s$df, t_i$df), df_res, sep = ", "),
      p = c(t_g$p, t_s$p, t_i$p))
  }

  # Pairwise Tukey HSD on the Group x Subgroup cells (same error term).
  sub$Cell <- cell_factor(sub)
  fit_c <- tryCatch(aov(delta_ct ~ Cell, data = sub), error = function(e) NULL)
  pair_df <- NULL
  if (!is.null(fit_c)) {
    tk <- tryCatch(TukeyHSD(fit_c)$Cell, error = function(e) NULL)
    if (!is.null(tk) && nrow(tk) > 0) {
      pairs <- t(vapply(rownames(tk), split_pair_name, character(2), levels = levels(sub$Cell)))
      keep  <- !is.na(pairs[, 1])
      if (any(keep)) {
        pair_df <- stat_rows(gene, paste(pairs[keep, 1], "vs", pairs[keep, 2]),
                             "Tukey HSD (Group x Subgroup cells)",
                             diff = -tk[keep, "diff"],
                             ci_low = -tk[keep, "upr"], ci_high = -tk[keep, "lwr"],
                             p = tk[keep, "p adj"])
      }
    }
  }
  rbind(main_df, pair_df)
}

# Nonparametric fallback for two-factor designs: Kruskal-Wallis across cells
# plus pairwise Mann-Whitney U tests with Holm correction.
run_cells_nonparametric <- function(sub, gene) {
  sub$Cell <- cell_factor(sub)
  kw <- tryCatch(kruskal.test(delta_ct ~ Cell, data = sub), error = function(e) NULL)
  omnibus <- if (is.null(kw)) NULL else
    stat_rows(gene, "Kruskal-Wallis across cells (omnibus)", "Kruskal-Wallis test",
              statistic = paste0("H = ", fmt_num(unname(kw$statistic))),
              df = unname(kw$parameter), p = kw$p.value)
  pw <- tryCatch(suppressWarnings(
    pairwise.wilcox.test(sub$delta_ct, sub$Cell, p.adjust.method = "holm")$p.value),
    error = function(e) NULL)
  if (is.null(pw)) return(omnibus)
  meds <- tapply(sub$delta_ct, sub$Cell, median)
  comps <- c(); pvals <- c(); diffs <- c()
  for (r in rownames(pw)) for (cc in colnames(pw)) {
    v <- pw[r, cc]
    if (!is.na(v)) {
      comps <- c(comps, paste(r, "vs", cc)); pvals <- c(pvals, v)
      diffs <- c(diffs, -(meds[[r]] - meds[[cc]]))
    }
  }
  rbind(omnibus, stat_rows(gene, comps, "Pairwise Mann-Whitney U between cells (Holm)",
                           diff = diffs, p = pvals))
}

# --- Per-group summary (for the results export) ---------------------------
summarize_groups <- function(df) {
  has_sub <- "Subgroup" %in% names(df)
  key <- if (has_sub) paste(df$Gene, df$Group, df$Subgroup, sep = "\u001f")
         else paste(df$Gene, df$Group, sep = "\u001f")
  parts <- split(df, factor(key, levels = unique(key)))
  rows <- lapply(parts, function(d) {
    l2 <- d$log2_fold_change[!is.na(d$log2_fold_change)]
    fc <- d$fold_change[!is.na(d$fold_change)]
    n  <- length(l2)
    out <- data.frame(Gene = d$Gene[1], Group = d$Group[1], stringsAsFactors = FALSE)
    if (has_sub) out$Subgroup <- d$Subgroup[1]
    out$n               <- n
    out$mean_log2FC     <- if (n) mean(l2) else NA_real_
    out$sd_log2FC       <- if (n > 1) sd(l2) else NA_real_
    out$sem_log2FC      <- if (n > 1) sd(l2) / sqrt(n) else NA_real_
    out$mean_FC         <- if (n) mean(fc) else NA_real_
    out$sd_FC           <- if (n > 1) sd(fc) else NA_real_
    out$sem_FC          <- if (n > 1) sd(fc) / sqrt(n) else NA_real_
    out$geo_mean_FC     <- if (n) 2^mean(l2) else NA_real_
    out
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

# --- Paired / repeated-measures analysis -----------------------------------
is_repeated_measures <- function(df) {
  if (nrow(df) == 0) return(FALSE)
  counts <- tapply(df$Group, df$Sample, function(x) length(unique(x)))
  any(counts >= 2, na.rm = TRUE)
}

# Per-gene paired t-tests between every pair of groups, Holm-corrected within
# each gene when there is more than one comparison.
run_paired_stats <- function(df) {
  genes <- unique(df$Gene)
  results <- lapply(genes, function(gene) {
    sub    <- df[df$Gene == gene & !is.na(df$delta_ct), ]
    groups <- unique(sub$Group)
    if (length(groups) < 2) return(NULL)
    pair_idx <- combn(length(groups), 2, simplify = FALSE)
    rows <- lapply(pair_idx, function(ij) {
      g1 <- groups[ij[1]]; g2 <- groups[ij[2]]
      d1 <- sub[sub$Group == g1, ]
      d2 <- sub[sub$Group == g2, ]
      common <- intersect(d1$Sample, d2$Sample)
      if (length(common) < 2) return(NULL)
      v1 <- d1$delta_ct[match(common, d1$Sample)]
      v2 <- d2$delta_ct[match(common, d2$Sample)]
      test_name <- paste0("Paired t-test (n = ", length(common), ")")
      tt <- tryCatch(t.test(v1, v2, paired = TRUE), error = function(e) NULL)
      if (is.null(tt)) return(stat_rows(gene, paste(g1, "vs", g2), test_name, p = NA_real_))
      stat_rows(gene, paste(g1, "vs", g2), test_name,
                statistic = paste0("t = ", fmt_num(unname(tt$statistic))),
                df = unname(tt$parameter),
                diff = -unname(tt$estimate), ci_low = -tt$conf.int[2], ci_high = -tt$conf.int[1],
                p = tt$p.value)
    })
    out <- do.call(rbind, Filter(Negate(is.null), rows))
    if (!is.null(out) && nrow(out) > 1) {
      out$p_value     <- signif(p.adjust(out$p_value, method = "holm"), 4)
      out$Significant <- ifelse(!is.na(out$p_value) & out$p_value < ALPHA, "Yes", "No")
      out$Test        <- paste0(out$Test, ", Holm-corrected")
    }
    out
  })
  out <- do.call(rbind, Filter(Negate(is.null), results))
  if (is.null(out)) out <- empty_stats()
  rownames(out) <- NULL
  out
}

# --- Methods text generator ------------------------------------------------
# Describes what was actually run, based on the Test column of stats_df.
generate_methods_text <- function(df, stats_df, test_type = NULL,
                                  paired = FALSE, control_group, has_subgroup = FALSE,
                                  rep_mode = "biological", sig_format = "exact",
                                  var_equal = FALSE) {
  n_txt <- replicate_summary(df, rep_mode = rep_mode)
  genes <- unique(df$Gene)
  tests_used <- unique(stats_df$Test)
  has <- function(pattern) any(grepl(pattern, tests_used, fixed = TRUE))

  baseline_phrase <- if (grepl(" | ", control_group, fixed = TRUE))
    paste0("the '", sub(" | ", " / ", control_group, fixed = TRUE), "' group (Group / Subgroup)")
  else
    paste0("the '", control_group, "' group")

  test_sentences <- c(
    if (has("Welch")) "Two groups were compared with Welch's unpaired t-test on delta Ct values.",
    if (has("Student")) "Two groups were compared with Student's unpaired t-test (equal variances assumed) on delta Ct values.",
    if (has("Paired t-test") && !has("Holm-corrected")) "Matched samples were compared with a paired t-test on delta Ct values.",
    if (has("Holm-corrected")) "Matched samples were compared with paired t-tests on delta Ct values; p values were adjusted for multiple comparisons with the Holm-Bonferroni method.",
    if (has("Mann-Whitney U test")) "Two groups were compared with the Mann-Whitney U test on delta Ct values.",
    if (has("Wilcoxon signed-rank")) "Matched samples were compared with the Wilcoxon signed-rank test on delta Ct values.",
    if (has("One-way ANOVA")) "Groups were compared by one-way ANOVA on delta Ct values followed by Tukey's HSD test for all pairwise comparisons.",
    if (has("Kruskal-Wallis") && !has("between cells")) "Groups were compared with the Kruskal-Wallis test on delta Ct values followed by pairwise Mann-Whitney U tests with Holm-Bonferroni correction.",
    if (has("Two-way ANOVA")) "The two-factor design was analysed by two-way ANOVA (Type II sums of squares) on delta Ct values with Group, Subgroup and their interaction as terms, followed by Tukey's HSD test between all Group x Subgroup cells.",
    if (has("between cells")) "Group x Subgroup cells were compared with the Kruskal-Wallis test on delta Ct values followed by pairwise Mann-Whitney U tests with Holm-Bonferroni correction.",
    if (has("Pooled-variance t-test vs control")) "Each group was compared with the control using pooled-variance t-tests (one-way ANOVA error term) with Holm-Bonferroni correction for the comparisons against the control (a Dunnett-type procedure).",
    if (has("Mann-Whitney U vs control")) "Each group was compared with the control using Mann-Whitney U tests with Holm-Bonferroni correction for the comparisons against the control."
  )
  if (length(test_sentences) == 0) test_sentences <- "No statistical test could be run on these data."

  p_sentence <- if (identical(sig_format, "stars"))
    "Significance is shown as * p < 0.05, ** p < 0.01, *** p < 0.001, **** p < 0.0001."
  else
    "Exact p values are reported; values below 0.0001 are shown as p < 0.0001."

  parts <- c(
    paste0("Relative gene expression was calculated with the delta-delta Ct method ",
           "(Livak and Schmittgen, 2001): delta Ct = Ct(target) - Ct(reference), and fold change = 2^(-delta-delta Ct) ",
           "relative to the mean delta Ct of ", baseline_phrase, "."),
    paste0("Expression of ", length(genes), " gene", if (length(genes) == 1) "" else "s", " (",
           paste(genes, collapse = ", "), ") was analysed with ", n_txt, "."),
    test_sentences,
    paste0("Statistical significance was defined as p < ", ALPHA, ". ", p_sentence),
    "Analyses were performed in R."
  )
  paste(parts, collapse = " ")
}
