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
    sub$delta_delta_ct <- sub$delta_ct - control_mean
    sub$fold_change <- 2^(-sub$delta_delta_ct)
    sub
  })
  do.call(rbind, result_list)
}

run_stats <- function(df, paired = FALSE) {
  genes <- unique(df$Gene)

  results <- lapply(genes, function(gene) {
    sub <- df[df$Gene == gene, ]
    groups <- unique(sub$Group)
    n_groups <- length(groups)

    if (n_groups < 2) return(NULL)

    if (n_groups == 2) {
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

      test <- t.test(g1_vals, g2_vals, paired = use_paired)
      data.frame(
        Gene        = gene,
        Comparison  = paste(groups[1], "vs", groups[2]),
        p_value     = round(test$p.value, 4),
        Significant = ifelse(test$p.value < 0.05, "Yes", "No"),
        stringsAsFactors = FALSE
      )
    } else {
      aov_fit <- aov(delta_ct ~ Group, data = sub)
      tukey   <- TukeyHSD(aov_fit)$Group
      comps   <- rownames(tukey)
      data.frame(
        Gene        = gene,
        Comparison  = sub("^(.+)-([^-]+)$", "\\1 vs \\2", comps),
        p_value     = round(tukey[, "p adj"], 4),
        Significant = ifelse(tukey[, "p adj"] < 0.05, "Yes", "No"),
        stringsAsFactors = FALSE
      )
    }
  })

  do.call(rbind, Filter(Negate(is.null), results))
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
