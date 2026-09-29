# qPCR Analysis (delta-delta Ct)

A browser app for relative gene expression from qPCR Ct values: paste from Excel, get fold changes, statistics, publication-style figures and a methods paragraph. Built with R Shiny and deployed as a static Shinylive site, so it runs in the browser with no installation and no data leaving your computer.

Live app: https://karthi-ntu.github.io/qPCR_analysis/

## What it does

- Paste raw Ct values straight from Excel (tab, comma or semicolon separated; decimal commas and "Undetermined" handled).
- delta Ct, delta-delta Ct and fold change per sample against a chosen control group (or a control Group x Subgroup cell in two-factor designs).
- Technical replicate averaging (rows sharing a Sample name) before delta Ct.
- Statistics on delta Ct with test statistic, degrees of freedom, difference and 95% CI on the log2 fold-change scale:
  - 2 groups: Welch's or Student's t-test, paired t-test, Mann-Whitney U, Wilcoxon signed-rank.
  - 3+ groups: one-way ANOVA + Tukey HSD, or Kruskal-Wallis + pairwise Mann-Whitney (Holm).
  - vs control only: pooled-variance t-tests or Mann-Whitney with Holm correction (Dunnett-type).
  - Two-factor designs: two-way ANOVA (Type II SS) + Tukey HSD between cells.
  - Repeated measures tab: paired t-tests between conditions, Holm-corrected.
- Prism-style plots: scatter or column with individual points, SD / SEM / 95% CI error bars, significance brackets with exact p or stars, grouped two-factor layout with nested axis labels, log2 or linear fold-change axis, Okabe-Ito colourblind-safe palette, Arial by default.
- Exports: PNG, TIFF (chosen DPI), PDF, SVG at exact millimetre sizes; per-sample results, per-group summary and statistics as CSV; methods paragraph as text.

Settings that cannot be applied to the current data (for example "Paired" with three groups) are reported above the plot rather than silently ignored, and the methods text always describes what was actually run.

## Input format

First row is the header. Columns: `Sample`, `Group`, optional `Subgroup` (also accepted: Factor2, Treatment, Condition), reference gene Ct (`Ct_ref`, `Reference gene`, ...), then one column per target gene.

```
Sample   Group    Ct_ref   Gene1    Gene2
mouse1   Control  13.50    20.51    25.40
mouse2   Control  13.56    20.62    25.23
mouse3   Treated  17.60    20.54    25.24
```

See the Help tab in the app for two-factor and repeated-measures layouts.

## Running locally

Requires R (4.1 or later) with `shiny`, `bslib`, `ggplot2`, `ggsignif`, `DT` and, for tests, `testthat`.

```r
install.packages(c("shiny", "bslib", "ggplot2", "ggsignif", "DT", "testthat"))
shiny::runApp(".")
```

On Windows, `run_local.bat` does the same.

## Tests

```r
testthat::test_dir("tests/testthat")
```

## Deploying the static site

The `docs/` folder is the Shinylive export served by GitHub Pages. After changing `app.R` or anything in `R/`, re-export and commit:

```r
install.packages("shinylive", repos = "https://posit-dev.r-universe.dev")
shinylive::export(appdir = ".", destdir = "docs")
```

`update_app.bat` runs the tests, exports, commits and pushes in one go.

## Project layout

```
app.R                 Shiny UI and server
R/data.R              paste parsing, wide/long conversion
R/analysis.R          delta-delta Ct, statistics, methods text
R/plot.R              ggplot builders, theme, file export
tests/testthat/       unit tests (data, analysis, plot rendering)
docs/                 Shinylive export (generated, do not edit by hand)
```

## Statistical notes

- All tests run on delta Ct (log2 scale), which is closer to normal than linear fold change.
- The "vs control" option is a Dunnett-type procedure using the pooled ANOVA variance and Holm-Bonferroni correction over the k-1 comparisons against the control; it is slightly more conservative than Prism's Dunnett test.
- The nonparametric multi-group option runs pairwise Mann-Whitney tests with Holm correction after Kruskal-Wallis, not Dunn's test.
- Two-way ANOVA uses Type II sums of squares, so main-effect p values do not depend on factor order in unbalanced designs.

## References

- Livak KJ, Schmittgen TD (2001). Analysis of relative gene expression data using real-time quantitative PCR and the 2^(-delta delta Ct) method. Methods 25(4):402-408.
- Bustin SA et al. (2009). The MIQE guidelines. Clin Chem 55(4):611-622.
