library(testthat)
source("../../R/data.R")

test_that("parse_ct handles decimal commas and undetermined tokens", {
  expect_equal(parse_ct(c("25.3", "25,3", " 24 ")), c(25.3, 25.3, 24))
  expect_true(all(is.na(parse_ct(c("Undetermined", "N/A", "", "NA", "abc")))))
})

test_that("to_wide_df keeps technical replicate rows separate", {
  long <- data.frame(Sample = rep("s1", 6), Group = "C",
                     Gene = rep(c("A", "B"), each = 3),
                     Ct_target = 1:6, Ct_reference = c(10, 11, 12, 10, 11, 12),
                     stringsAsFactors = FALSE)
  wide <- to_wide_df(long)
  expect_equal(nrow(wide), 3)
  expect_equal(wide$Ct_reference, c(10, 11, 12))
  expect_equal(wide$A, c(1, 2, 3))
  expect_equal(wide$B, c(4, 5, 6))
})

test_that("to_wide_df / to_long_df round trip preserves order and Subgroup", {
  wide <- data.frame(Sample = c("a", "b", "c"), Group = c("Y", "Y", "A"),
                     Subgroup = c("Ctl", "Trt", "Ctl"), Ct_reference = c(13, 13.5, 14),
                     G1 = c(20, 21, 22), G2 = c(25, NA, 27), stringsAsFactors = FALSE)
  long <- to_long_df(wide)
  expect_equal(nrow(long), 6)
  expect_equal(long$Gene[1:2], c("G1", "G2"))
  expect_true(is.na(long$Ct_target[long$Sample == "b" & long$Gene == "G2"]))
  back <- to_wide_df(long)
  expect_equal(back$Sample, wide$Sample)
  expect_equal(back$Subgroup, wide$Subgroup)
  expect_equal(back$G1, wide$G1)
  expect_equal(back$G2, wide$G2)
})

test_that("parse_pasted_text reads tab separated Excel data", {
  txt <- "Sample\tGroup\tCt_ref\tG1\tG2\nm1\tC\t13,5\t20.1\tUndetermined\nm2\tC\t13.6\t20.2\t25\nm3\tT\t13.4\t18\t23\n"
  res <- parse_pasted_text(txt)
  expect_true(res$ok)
  expect_equal(names(res$wide), c("Sample", "Group", "Ct_reference", "G1", "G2"))
  expect_equal(res$wide$Ct_reference, c(13.5, 13.6, 13.4))
  expect_true(is.na(res$wide$G2[1]))
  expect_true(length(res$warnings) == 1)
  expect_false(res$has_subgroup)
})

test_that("parse_pasted_text recognises Subgroup aliases and comma separation", {
  txt <- "Sample,Group,Treatment,Reference gene,G1\nm1,Y,Ctl,13.5,20\nm2,Y,Trt,13.5,19\n"
  res <- parse_pasted_text(txt)
  expect_true(res$ok)
  expect_true(res$has_subgroup)
  expect_equal(names(res$wide), c("Sample", "Group", "Subgroup", "Ct_reference", "G1"))
  expect_equal(res$wide$Subgroup, c("Ctl", "Trt"))
})

test_that("parse_pasted_text rejects data without a header row", {
  res <- parse_pasted_text("m1\tA\t13\t20\nm2\tA\t13\t21\n")
  expect_false(res$ok)
  expect_match(res$error, "header")
})

test_that("parse_pasted_text rejects too few columns and duplicate genes", {
  expect_false(parse_pasted_text("Sample\tGroup\tCt\nm1\tA\t13\n")$ok)
  expect_false(parse_pasted_text("Sample\tGroup\tCt\tG\tG\nm1\tA\t13\t1\t2\n")$ok)
})

test_that("missing_ct_cells lists the sample and gene", {
  long <- data.frame(Sample = c("a", "b"), Group = "C", Gene = "G1",
                     Ct_target = c(NA, 20), Ct_reference = c(13, NA), stringsAsFactors = FALSE)
  m <- missing_ct_cells(long)
  expect_equal(length(m), 2)
  expect_true(any(grepl("reference", m)))
})
