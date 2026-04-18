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
