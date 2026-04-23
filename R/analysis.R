# qPCR analysis pure functions

compute_delta_ct <- function(df) {
  df$delta_ct <- df$Ct_target - df$Ct_reference
  df
}

# --- Technical replicate averaging ----------------------------------------
# Collapses rows that share the same (Sample, Group, [Subgroup], Gene) to a
# single row by averaging Ct_target and Ct_reference independently. Returns
# a data frame with one row per biological replicate per gene.
# Also returns attr(out, "tech_counts") — vector of rep-counts per sample.
average_tech_reps <- function(df) {
  if (nrow(df) == 0) return(df)
  has_sub <- "Subgroup" %in% names(df)
  keys <- if (has_sub)
    paste(df$Sample, df$Group, df$Subgroup, df$Gene, sep = "\u001f")
  else
    paste(df$Sample, df$Group, df$Gene, sep = "\u001f")
  split_df <- split(df, keys)
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
# e.g. "n = 3 biological replicates" or
#      "n = 3 biological replicates (averaged from 3 technical replicates per well)"
replicate_summary <- function(df, rep_mode = "biological") {
  if (nrow(df) == 0) return("")
  has_sub <- "Subgroup" %in% names(df)
  keys <- if (has_sub)
    paste(df$Sample, df$Group, df$Subgroup, sep = "\u001f")
  else
    paste(df$Sample, df$Group, sep = "\u001f")
  n_bio <- length(unique(keys))
  per_group_keys <- if (has_sub) paste(df$Group, df$Subgroup, sep = " | ")
                    else df$Group
  bio_per_group <- tapply(keys, per_group_keys,
                          function(x) length(unique(x)))
  bio_range <- range(bio_per_group, na.rm = TRUE)
  bio_txt <- if (bio_range[1] == bio_range[2])
    paste0("n = ", bio_range[1]) else paste0("n = ", bio_range[1], "-", bio_range[2])
  if (rep_mode == "technical") {
    tc <- attr(df, "tech_counts")
    if (!is.null(tc) && length(tc) > 0) {
      tc_range <- range(tc)
      tc_txt <- if (tc_range[1] == tc_range[2])
        paste0(tc_range[1], " technical replicates")
      else
        paste0(tc_range[1], "-", tc_range[2], " technical replicates")
      return(paste0(bio_txt, " biological replicates (averaged from ",
                    tc_txt, " per sample)"))
    }
  }
  paste0(bio_txt, " biological replicates")
}

compute_fold_change <- function(df, control_group, control_subgroup = NULL) {
  if (!control_group %in% df$Group) {
    stop(sprintf("control_group '%s' not found in df$Group", control_group))
  }
  use_sub <- !is.null(control_subgroup) &&
             nzchar(control_subgroup) &&
             "Subgroup" %in% names(df) &&
             control_subgroup %in% df$Subgroup
  genes <- unique(df$Gene)
  result_list <- lapply(genes, function(gene) {
    sub <- df[df$Gene == gene, ]
    if (use_sub) {
      idx <- sub$Group == control_group & sub$Subgroup == control_subgroup
    } else {
      idx <- sub$Group == control_group
    }
    control_mean <- mean(sub$delta_ct[idx], na.rm = TRUE)
    sub$delta_delta_ct    <- sub$delta_ct - control_mean
    sub$fold_change       <- 2^(-sub$delta_delta_ct)
    sub$log2_fold_change  <- log2(sub$fold_change)
    sub
  })
  do.call(rbind, result_list)
}

# --- Main stats (parametric / nonparametric / vs-control) ------------------
run_stats <- function(df, paired = FALSE, test_type = "parametric",
                      has_subgroup = FALSE, control_group = NULL) {
  genes <- unique(df$Gene)

  vs_control <- test_type %in% c("parametric_vs_control", "nonparametric_vs_control")

  results <- lapply(genes, function(gene) {
    sub <- df[df$Gene == gene, ]

    # Two-way ANOVA branch (only when NOT in vs-control mode; the vs-control
    # flow ignores Subgroup and runs direct comparisons against the control).
    if (!vs_control && has_subgroup && "Subgroup" %in% names(sub) &&
        length(unique(sub$Subgroup)) >= 2) {
      return(run_two_way(sub, gene))
    }

    groups   <- unique(sub$Group)
    n_groups <- length(groups)
    if (n_groups < 2) return(NULL)

    if (vs_control) {
      if (is.null(control_group) || !control_group %in% groups) return(NULL)
      run_vs_control(sub, gene, control_group, test_type)
    } else if (n_groups == 2) {
      run_two_group(sub, gene, groups, paired, test_type)
    } else {
      run_multi_group(sub, gene, test_type)
    }
  })

  do.call(rbind, Filter(Negate(is.null), results))
}

# Dunnett-style: each non-control group vs control, Holm-Bonferroni adjusted.
# Uses pooled-SD pairwise t-test (closer to true Dunnett's method) in the
# parametric branch, and Mann-Whitney U + Holm in the nonparametric branch.
run_vs_control <- function(sub, gene, control_group, test_type) {
  groups <- unique(sub$Group)
  others <- setdiff(groups, control_group)
  if (length(others) == 0) return(NULL)

  if (test_type == "nonparametric_vs_control") {
    raw <- vapply(others, function(g) {
      tryCatch(
        suppressWarnings(wilcox.test(
          sub$delta_ct[sub$Group == g],
          sub$delta_ct[sub$Group == control_group],
          exact = FALSE)$p.value),
        error = function(e) NA_real_)
    }, numeric(1))
    adj <- p.adjust(raw, method = "holm")
    test_name <- "Mann-Whitney vs control (Holm)"
  } else {
    # Pooled-SD pairwise t-test is the standard Dunnett-approximating approach
    # when multcomp is unavailable (which is the case in webR/Shinylive).
    pw <- tryCatch(
      suppressWarnings(
        pairwise.t.test(sub$delta_ct, sub$Group,
                        p.adjust.method = "holm", pool.sd = TRUE)$p.value),
      error = function(e) NULL
    )
    if (is.null(pw)) return(NULL)
    adj <- vapply(others, function(g) {
      v1 <- if (g %in% rownames(pw) && control_group %in% colnames(pw))
              pw[g, control_group] else NA_real_
      v2 <- if (control_group %in% rownames(pw) && g %in% colnames(pw))
              pw[control_group, g] else NA_real_
      if (!is.na(v1)) v1 else v2
    }, numeric(1))
    test_name <- "Dunnett-style t-test vs control (Holm, pooled SD)"
  }

  data.frame(
    Gene        = rep(gene, length(others)),
    Comparison  = paste(others, "vs", control_group),
    Test        = rep(test_name, length(others)),
    p_value     = signif(adj, 4),
    Significant = ifelse(!is.na(adj) & adj < 0.05, "Yes", "No"),
    stringsAsFactors = FALSE
  )
}

run_two_group <- function(sub, gene, groups, paired, test_type) {
  g1_vals <- sub$delta_ct[sub$Group == groups[1]]
  g2_vals <- sub$delta_ct[sub$Group == groups[2]]

  use_paired <- paired
  if (paired) {
    sub1 <- sub[sub$Group == groups[1], ]
    sub2 <- sub[sub$Group == groups[2], ]
    common_samples <- intersect(sub1$Sample, sub2$Sample)
    if (length(common_samples) >= 2) {
      sub1 <- sub1[sub1$Sample %in% common_samples, ]
      sub2 <- sub2[sub2$Sample %in% common_samples, ]
      sub1 <- sub1[order(sub1$Sample), ]
      sub2 <- sub2[order(sub2$Sample), ]
      g1_vals <- sub1$delta_ct
      g2_vals <- sub2$delta_ct
    } else {
      use_paired <- FALSE
    }
  }

  if (test_type == "nonparametric") {
    test_name <- if (use_paired) "Wilcoxon signed-rank" else "Mann-Whitney U"
    p <- tryCatch(
      wilcox.test(g1_vals, g2_vals, paired = use_paired, exact = FALSE)$p.value,
      error = function(e) NA_real_
    )
  } else {
    test_name <- if (use_paired) "Paired t-test" else "Welch's t-test"
    p <- tryCatch(
      t.test(g1_vals, g2_vals, paired = use_paired)$p.value,
      error = function(e) NA_real_
    )
  }

  data.frame(
    Gene        = gene,
    Comparison  = paste(groups[1], "vs", groups[2]),
    Test        = test_name,
    p_value     = signif(p, 4),
    Significant = ifelse(!is.na(p) & p < 0.05, "Yes", "No"),
    stringsAsFactors = FALSE
  )
}

run_multi_group <- function(sub, gene, test_type) {
  if (test_type == "nonparametric") {
    pw <- tryCatch(
      pairwise.wilcox.test(sub$delta_ct, sub$Group,
                           p.adjust.method = "bonferroni", exact = FALSE)$p.value,
      error = function(e) NULL
    )
    if (is.null(pw)) return(NULL)
    # pw is a lower triangular matrix; extract pairs
    comps <- c(); pvals <- c()
    for (r in rownames(pw)) {
      for (cc in colnames(pw)) {
        v <- pw[r, cc]
        if (!is.na(v)) {
          comps <- c(comps, paste(r, "vs", cc))
          pvals <- c(pvals, v)
        }
      }
    }
    data.frame(
      Gene        = gene,
      Comparison  = comps,
      Test        = "Kruskal-Wallis + Dunn (Bonferroni)",
      p_value     = signif(pvals, 4),
      Significant = ifelse(pvals < 0.05, "Yes", "No"),
      stringsAsFactors = FALSE
    )
  } else {
    aov_fit <- aov(delta_ct ~ Group, data = sub)
    tukey   <- TukeyHSD(aov_fit)$Group
    comps   <- rownames(tukey)
    data.frame(
      Gene        = gene,
      Comparison  = sub("^(.+)-([^-]+)$", "\\1 vs \\2", comps),
      Test        = "One-way ANOVA + Tukey HSD",
      p_value     = signif(tukey[, "p adj"], 4),
      Significant = ifelse(tukey[, "p adj"] < 0.05, "Yes", "No"),
      stringsAsFactors = FALSE
    )
  }
}

run_two_way <- function(sub, gene) {
  fit <- tryCatch(aov(delta_ct ~ Group * Subgroup, data = sub),
                  error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  s <- summary(fit)[[1]]
  # Extract main effect and interaction p-values
  terms <- trimws(rownames(s))
  pick  <- function(name) {
    i <- which(terms == name)
    if (length(i) == 0) return(NA_real_)
    s[i, "Pr(>F)"]
  }
  p_group     <- pick("Group")
  p_subgroup  <- pick("Subgroup")
  p_inter     <- pick("Group:Subgroup")
  main_df <- data.frame(
    Gene        = rep(gene, 3),
    Comparison  = c("Group (main effect)", "Subgroup (main effect)",
                    "Group x Subgroup (interaction)"),
    Test        = rep("Two-way ANOVA", 3),
    p_value     = signif(c(p_group, p_subgroup, p_inter), 4),
    Significant = ifelse(!is.na(c(p_group, p_subgroup, p_inter)) &
                           c(p_group, p_subgroup, p_inter) < 0.05, "Yes", "No"),
    stringsAsFactors = FALSE
  )

  # --- Pairwise Tukey HSD on the 4 (or N) interaction cells ---------------
  # Build a combined factor whose levels match the plot x-axis exactly.
  sub$Combined <- interaction(sub$Group, sub$Subgroup,
                              sep = " | ", drop = TRUE, lex.order = FALSE)
  fit_c <- tryCatch(aov(delta_ct ~ Combined, data = sub),
                    error = function(e) NULL)
  pair_df <- NULL
  if (!is.null(fit_c)) {
    tk <- tryCatch(TukeyHSD(fit_c)$Combined, error = function(e) NULL)
    if (!is.null(tk) && nrow(tk) > 0) {
      # Parse rownames safely: each rowname is "<levelA>-<levelB>" where both
      # levels are known. Search the level set to find the split point (levels
      # may themselves contain "-").
      levs <- levels(sub$Combined)
      parse_pair <- function(rn) {
        for (L in levs) {
          suffix <- paste0("-", L)
          if (endsWith(rn, suffix)) {
            pre <- substr(rn, 1, nchar(rn) - nchar(suffix))
            if (pre %in% levs) return(c(pre, L))
          }
        }
        c(NA_character_, NA_character_)
      }
      pairs  <- t(vapply(rownames(tk), parse_pair, character(2)))
      pvals  <- tk[, "p adj"]
      keep   <- !is.na(pairs[, 1])
      if (any(keep)) {
        pair_df <- data.frame(
          Gene        = rep(gene, sum(keep)),
          Comparison  = paste(pairs[keep, 1], "vs", pairs[keep, 2]),
          Test        = rep("Two-way ANOVA + Tukey HSD (pairwise)", sum(keep)),
          p_value     = signif(pvals[keep], 4),
          Significant = ifelse(!is.na(pvals[keep]) & pvals[keep] < 0.05,
                               "Yes", "No"),
          stringsAsFactors = FALSE
        )
      }
    }
  }

  if (!is.null(pair_df)) rbind(main_df, pair_df) else main_df
}

summarize_groups <- function(df) {
  genes  <- unique(df$Gene)
  groups <- unique(df$Group)

  rows <- lapply(genes, function(gene) {
    lapply(groups, function(grp) {
      vals <- df$fold_change[df$Gene == gene & df$Group == grp]
      if (length(vals) == 0) return(NULL)
      data.frame(
        Gene    = gene,
        Group   = grp,
        mean_fc = mean(vals, na.rm = TRUE),
        sd_fc   = sd(vals, na.rm = TRUE),
        sem_fc  = sd(vals, na.rm = TRUE) / sqrt(length(vals)),
        stringsAsFactors = FALSE
      )
    })
  })

  do.call(rbind, Filter(Negate(is.null), unlist(rows, recursive = FALSE)))
}

# --- Paired / repeated-measures analysis -----------------------------------
# Detect whether data is repeated-measures: same Sample in 2+ Groups
is_repeated_measures <- function(df) {
  if (nrow(df) == 0) return(FALSE)
  counts <- tapply(df$Group, df$Sample, function(x) length(unique(x)))
  any(counts >= 2, na.rm = TRUE)
}

# Per-gene paired tests between each pair of groups
run_paired_stats <- function(df) {
  genes <- unique(df$Gene)
  results <- lapply(genes, function(gene) {
    sub    <- df[df$Gene == gene, ]
    groups <- unique(sub$Group)
    if (length(groups) < 2) return(NULL)
    # All unordered pairs of groups
    pair_idx <- combn(length(groups), 2, simplify = FALSE)
    rows <- lapply(pair_idx, function(ij) {
      g1 <- groups[ij[1]]; g2 <- groups[ij[2]]
      d1 <- sub[sub$Group == g1, ]
      d2 <- sub[sub$Group == g2, ]
      common <- intersect(d1$Sample, d2$Sample)
      if (length(common) < 2) return(NULL)
      v1 <- d1$delta_ct[match(common, d1$Sample)]
      v2 <- d2$delta_ct[match(common, d2$Sample)]
      p  <- tryCatch(t.test(v1, v2, paired = TRUE)$p.value,
                     error = function(e) NA_real_)
      data.frame(
        Gene        = gene,
        Comparison  = paste(g1, "vs", g2),
        Test        = paste0("Paired t-test (n=", length(common), ")"),
        p_value     = signif(p, 4),
        Significant = ifelse(!is.na(p) & p < 0.05, "Yes", "No"),
        stringsAsFactors = FALSE
      )
    })
    do.call(rbind, Filter(Negate(is.null), rows))
  })
  do.call(rbind, Filter(Negate(is.null), results))
}

# --- Methods text generator ------------------------------------------------
generate_methods_text <- function(df, stats_df, test_type,
                                  paired, control_group, has_subgroup = FALSE,
                                  rep_mode = "biological") {
  n_txt <- replicate_summary(df, rep_mode = rep_mode)

  genes <- unique(df$Gene)
  tests_used <- unique(stats_df$Test)
  test_desc <- paste(tests_used, collapse = " / ")

  baseline_phrase <- if (grepl(" | ", control_group, fixed = TRUE))
    paste0("'", control_group, "' cell (Group | Subgroup combination)")
  else
    paste0("'", control_group, "' group")
  parts <- c(
    paste0("Relative gene expression was calculated using the \u0394\u0394Ct method ",
           "(Livak & Schmittgen, 2001), with 2^(-\u0394\u0394Ct) as fold change ",
           "relative to the ", baseline_phrase, "."),
    paste0("Data for ", length(genes), " gene(s) (",
           paste(genes, collapse = ", "),
           ") were analyzed (", n_txt, " per group)."),
    paste0("Statistical testing: ", test_desc, "."),
    if (paired) "Paired design was applied where sample pairing was available." else NULL,
    if (has_subgroup) "A two-factor design (Group x Subgroup) was analyzed by two-way ANOVA." else NULL,
    if (any(grepl("vs control", tests_used, ignore.case = TRUE)))
      "Multiple comparisons against the control group were adjusted using the Holm-Bonferroni method (Dunnett-style pooled-variance t-test in the parametric branch)." else NULL,
    "Significance threshold: \u03b1 = 0.05. P-values < 0.0001 are displayed as 'p < 0.0001'."
  )
  paste(Filter(Negate(is.null), parts), collapse = " ")
}
