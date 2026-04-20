# qPCR analysis pure functions

compute_delta_ct <- function(df) {
  df$delta_ct <- df$Ct_target - df$Ct_reference
  df
}

compute_fold_change <- function(df, control_group) {
  if (!control_group %in% df$Group) {
    stop(sprintf("control_group '%s' not found in df$Group", control_group))
  }
  genes <- unique(df$Gene)
  result_list <- lapply(genes, function(gene) {
    sub <- df[df$Gene == gene, ]
    control_mean <- mean(sub$delta_ct[sub$Group == control_group], na.rm = TRUE)
    sub$delta_delta_ct    <- sub$delta_ct - control_mean
    sub$fold_change       <- 2^(-sub$delta_delta_ct)
    sub$log2_fold_change  <- log2(sub$fold_change)
    sub
  })
  do.call(rbind, result_list)
}

# --- Normality test (Shapiro-Wilk per group per gene) ----------------------
run_normality <- function(df) {
  genes <- unique(df$Gene)
  rows <- lapply(genes, function(gene) {
    sub    <- df[df$Gene == gene, ]
    groups <- unique(sub$Group)
    do.call(rbind, lapply(groups, function(grp) {
      vals <- sub$delta_ct[sub$Group == grp]
      vals <- vals[!is.na(vals)]
      if (length(vals) < 3) {
        data.frame(Gene = gene, Group = grp, n = length(vals),
                   shapiro_p = NA_real_, Normal = NA_character_,
                   stringsAsFactors = FALSE)
      } else {
        p <- tryCatch(shapiro.test(vals)$p.value, error = function(e) NA_real_)
        data.frame(Gene = gene, Group = grp, n = length(vals),
                   shapiro_p = round(p, 4),
                   Normal    = ifelse(is.na(p), "NA", ifelse(p > 0.05, "Yes", "No")),
                   stringsAsFactors = FALSE)
      }
    }))
  })
  do.call(rbind, rows)
}

# Is all data normal per gene? (TRUE if all Shapiro p > 0.05 or too few samples)
is_all_normal <- function(norm_df, gene) {
  sub <- norm_df[norm_df$Gene == gene, ]
  ps  <- sub$shapiro_p[!is.na(sub$shapiro_p)]
  if (length(ps) == 0) return(TRUE)
  all(ps > 0.05)
}

# --- Main stats (parametric / nonparametric / auto) ------------------------
run_stats <- function(df, paired = FALSE, test_type = "parametric",
                      norm_df = NULL, has_subgroup = FALSE) {
  genes <- unique(df$Gene)

  results <- lapply(genes, function(gene) {
    sub <- df[df$Gene == gene, ]

    # Two-way ANOVA branch
    if (has_subgroup && "Subgroup" %in% names(sub) &&
        length(unique(sub$Subgroup)) >= 2) {
      return(run_two_way(sub, gene))
    }

    # Decide parametric vs nonparametric for this gene
    use_type <- test_type
    if (test_type == "auto" && !is.null(norm_df))
      use_type <- if (is_all_normal(norm_df, gene)) "parametric" else "nonparametric"

    groups   <- unique(sub$Group)
    n_groups <- length(groups)
    if (n_groups < 2) return(NULL)

    if (n_groups == 2) {
      run_two_group(sub, gene, groups, paired, use_type)
    } else {
      run_multi_group(sub, gene, use_type)
    }
  })

  do.call(rbind, Filter(Negate(is.null), results))
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
  terms <- rownames(s)
  terms <- trimws(terms)
  pick  <- function(name) {
    i <- which(terms == name)
    if (length(i) == 0) return(NA_real_)
    s[i, "Pr(>F)"]
  }
  p_group     <- pick("Group")
  p_subgroup  <- pick("Subgroup")
  p_inter     <- pick("Group:Subgroup")
  data.frame(
    Gene        = rep(gene, 3),
    Comparison  = c("Group (main effect)", "Subgroup (main effect)",
                    "Group x Subgroup (interaction)"),
    Test        = rep("Two-way ANOVA", 3),
    p_value     = signif(c(p_group, p_subgroup, p_inter), 4),
    Significant = ifelse(c(p_group, p_subgroup, p_inter) < 0.05, "Yes", "No"),
    stringsAsFactors = FALSE
  )
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

# --- Methods text generator ------------------------------------------------
generate_methods_text <- function(df, stats_df, norm_df, test_type,
                                  paired, control_group, has_subgroup = FALSE) {
  n_per_grp <- table(df$Sample[!duplicated(df[, c("Sample", "Group")])],
                     df$Group[!duplicated(df[, c("Sample", "Group")])])
  n_range <- range(colSums(n_per_grp > 0))
  n_txt <- if (n_range[1] == n_range[2]) paste0("n = ", n_range[1])
           else paste0("n = ", n_range[1], "-", n_range[2])

  genes <- unique(df$Gene)

  # Determine which test was actually used
  tests_used <- unique(stats_df$Test)
  test_desc <- paste(tests_used, collapse = " / ")

  parts <- c(
    paste0("Relative gene expression was calculated using the \u0394\u0394Ct method ",
           "(Livak & Schmittgen, 2001), with 2^(-\u0394\u0394Ct) as fold change ",
           "relative to the '", control_group, "' group."),
    paste0("Data for ", length(genes), " gene(s) (",
           paste(genes, collapse = ", "),
           ") were analyzed (", n_txt, " biological replicates per group)."),
    if (!is.null(norm_df) && any(!is.na(norm_df$shapiro_p)))
      paste0("Normality was assessed with the Shapiro-Wilk test."),
    paste0("Statistical testing: ", test_desc, "."),
    if (paired) "Paired design was applied where sample pairing was available." else NULL,
    if (has_subgroup) "A two-factor design (Group x Subgroup) was analyzed by two-way ANOVA." else NULL,
    "Significance threshold: \u03b1 = 0.05. P-values < 0.0001 are displayed as 'p < 0.0001'."
  )
  paste(Filter(Negate(is.null), parts), collapse = " ")
}
