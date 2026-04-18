library(testthat)
library(ggplot2)
source("../../R/analysis.R")
source("../../R/plot.R")

test_that("make_barplot returns a ggplot object", {
  summary_df <- data.frame(
    Gene    = "ACTB",
    Group   = c("Control", "Treatment"),
    mean_fc = c(1.0, 2.0),
    sd_fc   = c(0.1, 0.2),
    sem_fc  = c(0.05, 0.1),
    stringsAsFactors = FALSE
  )
  stats_df <- data.frame(
    Gene        = "ACTB",
    Comparison  = "Control vs Treatment",
    p_value     = 0.02,
    Significant = "Yes",
    stringsAsFactors = FALSE
  )
  p <- make_barplot(summary_df, stats_df, gene = "ACTB", error_type = "SD")
  expect_s3_class(p, "gg")
})

test_that("make_barplot works with no significant comparisons", {
  summary_df <- data.frame(
    Gene    = "ACTB",
    Group   = c("Control", "Treatment"),
    mean_fc = c(1.0, 1.1),
    sd_fc   = c(0.1, 0.1),
    sem_fc  = c(0.05, 0.05),
    stringsAsFactors = FALSE
  )
  stats_df <- data.frame(
    Gene        = "ACTB",
    Comparison  = "Control vs Treatment",
    p_value     = 0.8,
    Significant = "No",
    stringsAsFactors = FALSE
  )
  p <- make_barplot(summary_df, stats_df, gene = "ACTB", error_type = "SEM")
  expect_s3_class(p, "gg")
})
