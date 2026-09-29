library(testthat)
library(ggplot2)
source("../../R/analysis.R")
source("../../R/plot.R")

# Rendering (ggplot_build) is what catches theme/aesthetic errors; building
# the ggplot object alone does not.
renders <- function(p) {
  expect_s3_class(p, "gg")
  expect_error(suppressWarnings(ggplot_build(p)), NA)
  invisible(p)
}

one_factor <- function() {
  d <- data.frame(Sample = paste0("s", 1:9), Group = rep(c("si-Ctrl", "si-A", "si-B"), each = 3),
                  Gene = "G", delta_ct = c(0, .1, -.1, 1, 1.1, .9, 2, 2.1, 1.9),
                  stringsAsFactors = FALSE)
  fc <- compute_fold_change(d, "si-Ctrl")
  list(fc = fc, st = run_stats(fc))
}

two_factor <- function() {
  set.seed(2)
  d <- expand.grid(rep = 1:3, Subgroup = c("Control", "Treated"), Group = c("Young", "Aged"),
                   stringsAsFactors = FALSE)
  d$Sample <- paste0("s", seq_len(nrow(d)))
  d <- rbind(transform(d, Gene = "G1"), transform(d, Gene = "G2"))
  d$delta_ct <- 5 - ifelse(d$Subgroup == "Treated", 1, 0) -
    ifelse(d$Subgroup == "Treated" & d$Group == "Aged", 1.5, 0) + rnorm(24, 0, .2)
  fc <- compute_fold_change(d, "Young", "Control")
  list(fc = fc, st = run_stats(fc, has_subgroup = TRUE))
}

test_that("single-gene plot renders for column and scatter, with and without stats", {
  x <- one_factor()
  renders(make_barplot(x$fc, x$st, gene = "G", plot_type = "column", error_type = "CI95"))
  renders(make_barplot(x$fc, x$st, gene = "G", plot_type = "scatter", error_type = "SEM",
                       sig_format = "stars", font_family = "Arial"))
  renders(make_barplot(x$fc, NULL, gene = "G"))
  renders(make_barplot(x$fc, x$st[0, ], gene = "G"))
})

test_that("brackets are drawn for every requested comparison, including ns", {
  x <- one_factor()
  x$st$Significant[x$st$Comparison == "si-B vs si-A"] <- "No"
  count_brackets <- function(p) {
    b <- ggplot_build(p)
    layer <- which(vapply(p$layers, function(l) inherits(l$geom, "GeomSignif"), logical(1)))
    if (length(layer) == 0) return(0)
    length(unique(b$data[[layer]]$group))
  }
  expect_equal(count_brackets(make_barplot(x$fc, x$st, "G")), 2)                      # significant only
  expect_equal(count_brackets(make_barplot(x$fc, x$st, "G", sig_comparisons = character(0))), 0)
  expect_equal(count_brackets(make_barplot(x$fc, x$st, "G",
                                           sig_comparisons = x$st$Comparison[-1])), 3)  # incl. ns
})

test_that("combined plot renders in both layouts, both scales and with all comparisons", {
  x <- two_factor()
  renders(make_combined_plot(x$fc, x$st, has_subgroup = TRUE, control_group = "Young | Control"))
  renders(make_combined_plot(x$fc, x$st, has_subgroup = TRUE, x_axis_var = "Subgroup",
                             plot_type = "column"))
  renders(make_combined_plot(x$fc, x$st, has_subgroup = TRUE, y_scale = "linear",
                             sig_comparisons_all = paste0(x$st$Gene, ": ", x$st$Comparison)))
  renders(make_combined_plot(x$fc, x$st, has_subgroup = TRUE, sig_comparisons_all = character(0)))
  renders(make_barplot(x$fc, x$st, "G1", has_subgroup = TRUE, font_family = "Arial"))
})

test_that("colour overrides apply only to valid hex codes", {
  f <- resolve_fill(c(A = "#ff0000", B = "red"), c("A", "B", "C"))
  expect_equal(unname(f["A"]), "#ff0000")
  expect_equal(unname(f["B"]), OKABE_ITO[2])
  expect_equal(length(f), 3)
})

test_that("paired plot renders with more than 8 samples and brackets", {
  set.seed(4)
  dp <- data.frame(Sample = rep(paste0("m", 1:10), 3), Group = rep(c("Base", "T1", "T2"), each = 10),
                   Gene = "G", delta_ct = c(rnorm(10, 5, .3), rnorm(10, 4, .3), rnorm(10, 3, .3)))
  fc <- compute_fold_change(dp, "Base")
  renders(make_paired_plot(fc, run_paired_stats(fc), control_group = "Base", error_type = "SEM"))
})

test_that("p-value formatters follow Prism thresholds", {
  expect_equal(format_pvalue_stars(c(0.2, 0.04, 0.009, 0.0009, 0.00009, NA)),
               c("ns", "*", "**", "***", "****", "NA"))
  expect_equal(format_pvalue(c(0.5, 0.00001)), c("p = 0.5000", "p < 0.0001"))
})

test_that("plots save to png and pdf at the requested size", {
  x <- one_factor()
  p <- make_barplot(x$fc, x$st, "G", font_family = "Arial")
  tmp <- tempfile(fileext = ".png")
  suppressWarnings(save_plot_file(tmp, p, "png", 85, 85, dpi = 150))
  expect_true(file.exists(tmp) && file.size(tmp) > 1000)
  tmp2 <- tempfile(fileext = ".pdf")
  suppressWarnings(save_plot_file(tmp2, p, "pdf", 85, 85))
  expect_true(file.exists(tmp2) && file.size(tmp2) > 1000)
})
