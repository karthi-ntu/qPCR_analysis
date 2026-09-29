library(testthat)
source("../../R/analysis.R")

test_that("compute_delta_ct subtracts reference from target", {
  df <- data.frame(
    Sample = c("S1", "S2"),
    Group = c("Control", "Treatment"),
    Gene = c("ACTB", "ACTB"),
    Ct_target = c(20.0, 22.0),
    Ct_reference = c(18.0, 18.0)
  )
  result <- compute_delta_ct(df)
  expect_equal(result$delta_ct, c(2.0, 4.0))
})

test_that("compute_fold_change sets control group to 1.0", {
  df <- data.frame(
    Sample = c("S1", "S2", "S3"),
    Group = c("Control", "Control", "Treatment"),
    Gene = c("ACTB", "ACTB", "ACTB"),
    delta_ct = c(2.0, 2.0, 4.0)
  )
  result <- compute_fold_change(df, control_group = "Control")
  expect_equal(result$fold_change[result$Group == "Control"], c(1.0, 1.0))
  expect_equal(result$fold_change[result$Group == "Treatment"], 2^(-2.0))
  expect_equal(result$log2_fold_change[result$Group == "Treatment"], -2.0)
})

test_that("compute_fold_change applies per-gene control mean independently", {
  df <- data.frame(
    Sample   = c("S1", "S2", "S3", "S4"),
    Group    = c("Control", "Treatment", "Control", "Treatment"),
    Gene     = c("ACTB", "ACTB", "MYC", "MYC"),
    delta_ct = c(2.0, 4.0, 5.0, 6.0)
  )
  result <- compute_fold_change(df, control_group = "Control")
  expect_equal(result$fold_change[result$Gene == "ACTB" & result$Group == "Treatment"], 2^(-2.0))
  expect_equal(result$fold_change[result$Gene == "MYC" & result$Group == "Treatment"], 2^(-1.0))
})

test_that("compute_fold_change uses the control subgroup cell as baseline", {
  df <- data.frame(
    Sample = paste0("S", 1:4), Group = c("Y", "Y", "A", "A"),
    Subgroup = c("Ctl", "Trt", "Ctl", "Trt"), Gene = "G",
    delta_ct = c(2, 1, 3, 0)
  )
  r <- compute_fold_change(df, "Y", "Ctl")
  expect_equal(r$log2_fold_change, c(0, 1, -1, 2))
  expect_error(compute_fold_change(df, "Y", "Nope"), "Subgroup")
  expect_error(compute_fold_change(df, "Missing"), "Control group")
})

test_that("average_tech_reps averages both Ct columns and keeps order", {
  df <- data.frame(
    Sample = c("b", "b", "a", "a"), Group = "C", Gene = "G",
    Ct_target = c(20, 22, 30, 32), Ct_reference = c(10, 12, 15, 17),
    stringsAsFactors = FALSE
  )
  out <- average_tech_reps(df)
  expect_equal(out$Sample, c("b", "a"))
  expect_equal(out$Ct_target, c(21, 31))
  expect_equal(out$Ct_reference, c(11, 16))
  expect_equal(unname(attr(out, "tech_counts")), c(2L, 2L))
})

two_group_df <- function() {
  data.frame(
    Sample = paste0("S", 1:6),
    Group = rep(c("Control", "Treatment"), each = 3),
    Gene = "ACTB",
    delta_ct = c(2.0, 2.1, 1.9, 4.0, 4.1, 3.9)
  )
}

test_that("run_stats returns Welch t-test with statistic, df and CI for 2 groups", {
  result <- run_stats(two_group_df(), paired = FALSE)
  expect_equal(nrow(result), 1)
  expect_equal(result$Test, "Welch's t-test")
  expect_true(result$p_value < 0.05)
  expect_equal(result$Significant, "Yes")
  expect_match(result$Statistic, "^t = ")
  # Control delta Ct is 2 lower than Treatment, so Control is 4x higher:
  # "Control vs Treatment" difference = +2 on the log2 scale.
  expect_equal(result$Diff_log2FC, 2, tolerance = 1e-6)
  expect_true(result$CI95_low < 2 && result$CI95_high > 2)
  expect_true(all(STATS_COLS %in% names(result)))
})

test_that("run_stats can use Student's t-test", {
  result <- run_stats(two_group_df(), var_equal = TRUE)
  expect_match(result$Test, "Student")
  expect_equal(result$df, "4")
})

test_that("run_stats paired t-test for 2 groups matched by Sample", {
  df <- data.frame(
    Sample = c("S1", "S2", "S3", "S1", "S2", "S3"),
    Group = c(rep("Control", 3), rep("Treatment", 3)),
    Gene = "ACTB",
    delta_ct = c(2.0, 2.3, 1.8, 4.1, 4.0, 3.9)
  )
  result <- run_stats(df, paired = TRUE)
  expect_equal(result$Test, "Paired t-test")
  expect_true(result$p_value < 0.05)
  expect_equal(result$df, "2")
})

test_that("nonparametric two-group test is Mann-Whitney / signed-rank", {
  r <- run_stats(two_group_df(), test_type = "nonparametric")
  expect_equal(r$Test, "Mann-Whitney U test")
  expect_equal(r$p_value, 0.1, tolerance = 1e-6)   # exact p for 3 vs 3, no ties
})

test_that("run_stats returns omnibus ANOVA + Tukey rows and parses hyphenated names", {
  df <- data.frame(
    Sample = paste0("S", 1:9),
    Group = rep(c("si-Ctrl", "si-A", "si-B"), each = 3),
    Gene = "ACTB",
    delta_ct = c(0, .1, -.1, 1, 1.1, .9, 2, 2.1, 1.9)
  )
  result <- run_stats(df)
  expect_equal(result$Test[1], "One-way ANOVA")
  expect_setequal(result$Comparison[-1],
                  c("si-A vs si-Ctrl", "si-B vs si-Ctrl", "si-B vs si-A"))
  expect_equal(result$Diff_log2FC[result$Comparison == "si-B vs si-Ctrl"], -2, tolerance = 1e-6)
})

test_that("nonparametric multi-group runs Kruskal-Wallis and pairwise Mann-Whitney", {
  df <- data.frame(Sample = paste0("S", 1:9), Group = rep(c("A", "B", "C"), each = 3),
                   Gene = "G", delta_ct = c(0, .1, -.1, 1, 1.1, .9, 2, 2.1, 1.9))
  r <- run_stats(df, test_type = "nonparametric")
  expect_equal(r$Test[1], "Kruskal-Wallis test")
  expect_true(all(grepl("Mann-Whitney", r$Test[-1])))
  expect_equal(nrow(r), 4)
})

test_that("vs-control uses Holm across the k-1 control comparisons only", {
  set.seed(1)
  df <- data.frame(Gene = "G", Group = rep(c("C", "A", "B", "D"), each = 4),
                   delta_ct = c(rnorm(4, 0, .3), rnorm(4, 1, .3), rnorm(4, .4, .3), rnorm(4, .2, .3)))
  r <- run_stats(df, test_type = "parametric_vs_control", control_group = "C")
  expect_equal(r$Comparison, c("A vs C", "B vs C", "D vs C"))
  pw  <- pairwise.t.test(df$delta_ct, df$Group, p.adjust.method = "none", pool.sd = TRUE)$p.value
  raw <- c(pw["C", "A"], pw["C", "B"], pw["D", "C"])
  expect_equal(r$p_value, signif(p.adjust(raw, "holm"), 4))
  expect_equal(r$df, rep("12", 3))
  rn <- run_stats(df, test_type = "nonparametric_vs_control", control_group = "C")
  expect_true(all(grepl("Mann-Whitney", rn$Test)))
})

two_factor_df <- function(seed = 2) {
  set.seed(seed)
  d <- expand.grid(rep = 1:3, Subgroup = c("Control", "Treated"), Group = c("Young", "Aged"),
                   stringsAsFactors = FALSE)
  d$Sample <- paste0("s", seq_len(nrow(d))); d$Gene <- "G1"
  d$delta_ct <- 5 - ifelse(d$Subgroup == "Treated", 1, 0) -
    ifelse(d$Subgroup == "Treated" & d$Group == "Aged", 1.5, 0) + rnorm(12, 0, .2)
  d
}

test_that("two-way ANOVA (Type II) matches aov on balanced data and adds Tukey cells", {
  d <- two_factor_df()
  r <- run_stats(d, has_subgroup = TRUE)
  expect_equal(r$Comparison[1:3],
               c("Group (main effect)", "Subgroup (main effect)", "Group x Subgroup (interaction)"))
  p_aov <- summary(aov(delta_ct ~ Group * Subgroup, d))[[1]][1:3, "Pr(>F)"]
  expect_equal(r$p_value[1:3], signif(p_aov, 4))
  expect_equal(sum(grepl("Tukey", r$Test)), 6)
  expect_true(all(grepl(" | ", r$Comparison[-(1:3)], fixed = TRUE)))
})

test_that("two-factor vs-control compares cells against the control cell", {
  d <- two_factor_df()
  r <- run_stats(d, has_subgroup = TRUE, test_type = "parametric_vs_control",
                 control_group = "Young", control_subgroup = "Control")
  expect_equal(nrow(r), 3)
  expect_true(all(endsWith(r$Comparison, "vs Young | Control")))
})

test_that("settings that do not apply produce notes instead of silent changes", {
  d <- two_factor_df()
  r <- run_stats(d, has_subgroup = TRUE, test_type = "nonparametric", paired = TRUE)
  notes <- attr(r, "notes")
  expect_true(any(grepl("Paired", notes)))
  expect_true(any(grepl("nonparametric two-way", notes)))
  expect_true(all(grepl("Kruskal|Mann-Whitney", r$Test)))
  d3 <- data.frame(Sample = paste0("S", 1:9), Group = rep(c("A", "B", "C"), each = 3),
                   Gene = "G", delta_ct = c(0, .1, -.1, 1, 1.1, .9, 2, 2.1, 1.9))
  r3 <- run_stats(d3, paired = TRUE)
  expect_true(any(grepl("Repeated Measures", attr(r3, "notes"))))
})

test_that("run_stats returns an empty table with a single group", {
  df <- data.frame(Sample = c("a", "b"), Group = "A", Gene = "G", delta_ct = c(1, 2))
  r <- run_stats(df)
  expect_equal(nrow(r), 0)
  expect_true(all(STATS_COLS %in% names(r)))
})

test_that("run_paired_stats applies Holm within gene", {
  set.seed(3)
  dp <- data.frame(Sample = rep(paste0("m", 1:5), 3), Group = rep(c("Base", "T1", "T2"), each = 5),
                   Gene = "G", delta_ct = c(rnorm(5, 5, .3), rnorm(5, 4, .3), rnorm(5, 3, .3)))
  r <- run_paired_stats(dp)
  expect_equal(nrow(r), 3)
  expect_true(all(grepl("Holm", r$Test)))
  raw <- sapply(list(c("Base", "T1"), c("Base", "T2"), c("T1", "T2")), function(g)
    t.test(dp$delta_ct[dp$Group == g[1]], dp$delta_ct[dp$Group == g[2]], paired = TRUE)$p.value)
  expect_equal(r$p_value, signif(p.adjust(raw, "holm"), 4))
  expect_true(is_repeated_measures(dp))
})

test_that("summarize_groups reports n and means per group (and subgroup)", {
  d <- compute_fold_change(two_factor_df(), "Young", "Control")
  s <- summarize_groups(d)
  expect_equal(nrow(s), 4)
  expect_equal(s$n, rep(3, 4))
  expect_equal(s$mean_log2FC[1], 0, tolerance = 1e-8)
  expect_true(all(c("Subgroup", "mean_FC", "sem_log2FC", "geo_mean_FC") %in% names(s)))
})

test_that("methods text describes the tests that were actually run", {
  d <- two_factor_df()
  fc <- compute_fold_change(d, "Young", "Control")
  r  <- run_stats(fc, has_subgroup = TRUE, paired = TRUE)
  txt <- generate_methods_text(fc, r, paired = TRUE, control_group = "Young | Control",
                               has_subgroup = TRUE, sig_format = "stars")
  expect_match(txt, "two-way ANOVA")
  expect_false(grepl("[Pp]aired", txt))
  expect_match(txt, "\\*\\*\\*\\*")
  expect_false(grepl("—|–", txt))
  r2 <- run_stats(two_group_df(), test_type = "parametric_vs_control", control_group = "Control")
  txt2 <- generate_methods_text(two_group_df(), r2, control_group = "Control")
  expect_match(txt2, "Dunnett-type")
  expect_match(txt2, "Holm")
})
