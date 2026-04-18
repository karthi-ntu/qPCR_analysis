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
})

test_that("compute_fold_change applies per-gene control mean independently", {
  df <- data.frame(
    Sample   = c("S1", "S2", "S3", "S4"),
    Group    = c("Control", "Treatment", "Control", "Treatment"),
    Gene     = c("ACTB", "ACTB", "MYC", "MYC"),
    delta_ct = c(2.0, 4.0, 5.0, 6.0)
  )
  result <- compute_fold_change(df, control_group = "Control")
  # ACTB: control mean = 2.0, Treatment ΔΔCt = 2, FC = 2^(-2) = 0.25
  expect_equal(result$fold_change[result$Gene == "ACTB" & result$Group == "Control"], 1.0)
  expect_equal(result$fold_change[result$Gene == "ACTB" & result$Group == "Treatment"], 2^(-2.0))
  # MYC: control mean = 5.0, Treatment ΔΔCt = 1, FC = 2^(-1) = 0.5
  expect_equal(result$fold_change[result$Gene == "MYC" & result$Group == "Control"], 1.0)
  expect_equal(result$fold_change[result$Gene == "MYC" & result$Group == "Treatment"], 2^(-1.0))
})

test_that("compute_delta_ct preserves all input columns", {
  df <- data.frame(
    Sample = "S1", Group = "Control", Gene = "ACTB",
    Ct_target = 20.0, Ct_reference = 18.0
  )
  result <- compute_delta_ct(df)
  expect_true(all(c("Sample", "Group", "Gene", "Ct_target", "Ct_reference", "delta_ct") %in% names(result)))
})

test_that("run_stats returns unpaired t-test for 2 groups", {
  set.seed(42)
  df <- data.frame(
    Sample = paste0("S", 1:6),
    Group = rep(c("Control", "Treatment"), each = 3),
    Gene = "ACTB",
    delta_ct = c(2.0, 2.1, 1.9, 4.0, 4.1, 3.9)
  )
  result <- run_stats(df, paired = FALSE)
  expect_equal(nrow(result), 1)
  expect_true("p_value" %in% names(result))
  expect_true(result$p_value < 0.05)
  expect_equal(result$Significant, "Yes")
})

test_that("run_stats returns paired t-test result for 2 groups paired=TRUE", {
  df <- data.frame(
    Sample = c("S1", "S2", "S3", "S1", "S2", "S3"),
    Group = c(rep("Control", 3), rep("Treatment", 3)),
    Gene = "ACTB",
    delta_ct = c(2.0, 2.1, 1.9, 4.0, 4.1, 3.9)
  )
  result <- run_stats(df, paired = TRUE)
  expect_equal(nrow(result), 1)
  expect_true(result$p_value < 0.05)
})

test_that("run_stats returns ANOVA + Tukey for 3+ groups", {
  set.seed(42)
  df <- data.frame(
    Sample = paste0("S", 1:9),
    Group = rep(c("Control", "TreatA", "TreatB"), each = 3),
    Gene = "ACTB",
    delta_ct = c(2.0, 2.1, 1.9, 4.0, 4.1, 3.9, 6.0, 6.1, 5.9)
  )
  result <- run_stats(df, paired = FALSE)
  expect_true(nrow(result) >= 2)
  expect_true(all(c("Gene", "Comparison", "p_value", "Significant") %in% names(result)))
})
