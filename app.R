library(shiny)
library(ggplot2)
library(DT)
library(ggsignif)
library(bslib)

# Shinylive workaround: strip `download` attribute from downloadButton, otherwise
# Chromium browsers save the file as .htm/.xml in WebAssembly.
# https://posit-dev.github.io/r-shinylive/ (Chromium issue 468227 workaround)
downloadButton <- function(...) {
  tag <- shiny::downloadButton(...)
  tag$attribs$download <- NULL
  tag
}

source("R/analysis.R")
source("R/plot.R")

`%||%` <- function(a, b) if (is.null(a)) b else a

# --- Data helpers ----------------------------------------------------------
empty_wide <- function(with_subgroup = FALSE) {
  if (with_subgroup) {
    data.frame(Sample = character(0), Group = character(0),
               Subgroup = character(0), Ct_reference = numeric(0),
               stringsAsFactors = FALSE)
  } else {
    data.frame(Sample = character(0), Group = character(0),
               Ct_reference = numeric(0), stringsAsFactors = FALSE)
  }
}

empty_long <- function() {
  data.frame(Sample = character(0), Group = character(0),
             Gene = character(0), Ct_target = numeric(0),
             Ct_reference = numeric(0), Subgroup = character(0),
             stringsAsFactors = FALSE)
}

meta_cols <- function(wide) {
  intersect(c("Sample", "Group", "Subgroup", "Ct_reference"), names(wide))
}

to_long_df <- function(wide) {
  mcols     <- meta_cols(wide)
  gene_cols <- setdiff(names(wide), mcols)
  if (length(gene_cols) == 0 || nrow(wide) == 0) return(empty_long())
  has_sub <- "Subgroup" %in% names(wide)
  rows <- lapply(seq_len(nrow(wide)), function(i) {
    lapply(gene_cols, function(g) {
      r <- data.frame(
        Sample       = wide$Sample[i],
        Group        = wide$Group[i],
        Gene         = g,
        Ct_target    = suppressWarnings(as.numeric(wide[[g]][i])),
        Ct_reference = suppressWarnings(as.numeric(wide$Ct_reference[i])),
        stringsAsFactors = FALSE
      )
      if (has_sub) r$Subgroup <- wide$Subgroup[i]
      r
    })
  })
  do.call(rbind, Filter(Negate(is.null), unlist(rows, recursive = FALSE)))
}

to_wide_df <- function(long) {
  genes <- unique(long$Gene[nchar(trimws(long$Gene)) > 0])
  if (length(genes) == 0 || nrow(long) == 0) return(empty_wide())
  has_sub <- "Subgroup" %in% names(long)
  key_cols <- if (has_sub) c("Sample", "Group", "Subgroup") else c("Sample", "Group")
  unique_idx <- !duplicated(do.call(paste, c(long[key_cols], sep = "\t")))
  wide <- long[unique_idx, c(key_cols, "Ct_reference")]
  for (g in genes) {
    gene_sub <- long[long$Gene == g, c(key_cols, "Ct_target")]
    names(gene_sub)[length(key_cols) + 1] <- g
    wide <- merge(wide, gene_sub, by = key_cols, all.x = TRUE, sort = FALSE)
  }
  wide
}

# --- Example data generators (for "Load example" buttons) ----------------
example_data_1factor <- function() {
  data.frame(
    Sample       = as.character(1:8),
    Group        = rep(c("Control", "Treated"), each = 4),
    Ct_reference = c(13.5, 13.6, 13.4, 13.5, 13.5, 13.4, 13.6, 13.5),
    HERPUD1      = c(25.1, 25.3, 25.2, 25.0, 22.1, 22.3, 22.0, 22.4),
    ATF4         = c(24.5, 24.7, 24.6, 24.4, 21.8, 22.0, 21.9, 22.1),
    stringsAsFactors = FALSE
  )
}

example_data_2factor <- function() {
  data.frame(
    Sample       = as.character(1:12),
    Group        = c(rep("Young", 6), rep("Aged", 6)),
    Subgroup     = c(rep("Control", 3), rep("Thaps", 3),
                     rep("Control", 3), rep("Thaps", 3)),
    Ct_reference = c(13.50, 13.58, 13.42, 13.55, 13.62, 13.44,
                     13.50, 13.46, 13.52, 13.48, 13.44, 13.36),
    HERPUD1      = c(25.1, 25.3, 25.2, 22.1, 22.3, 22.0,
                     24.8, 24.9, 24.7, 21.5, 21.6, 21.4),
    stringsAsFactors = FALSE
  )
}

# --- App-specific CSS -----------------------------------------------------
# Visible paste dropzone + copy-button + data-summary-card styling. Scoped
# via class names so the defaults don't leak into other Bootstrap elements.
app_css <- "
.paste-dropzone {
  border: 2px dashed #00A8A8;
  border-radius: 8px;
  padding: 14px;
  background: rgba(0, 168, 168, 0.05);
  text-align: center;
  margin-bottom: 10px;
  transition: background 0.2s;
}
.paste-dropzone:hover { background: rgba(0, 168, 168, 0.12); }
.paste-dropzone .pz-title { font-weight: 600; color: #007373; margin-bottom: 4px; }
.paste-dropzone .pz-sub   { font-size: 12px; color: #666; }
.data-summary {
  padding: 8px 12px;
  border-radius: 6px;
  margin: 8px 0 12px;
  font-size: 13px;
  background: rgba(0, 168, 168, 0.08);
  border-left: 3px solid #00A8A8;
}
.data-summary.warn {
  background: rgba(217, 119, 6, 0.08);
  border-left-color: #d97706;
}
.data-summary.empty {
  background: rgba(120, 120, 120, 0.08);
  border-left-color: #aaa;
  color: #777;
}
.data-badge {
  font-size: 12px; padding: 2px 8px; border-radius: 10px;
  background: rgba(0, 168, 168, 0.15); color: #007373;
  margin-right: 10px;
}
.data-badge.empty { background: rgba(120,120,120,0.12); color: #777; }
.data-badge.warn  { background: rgba(217, 119, 6, 0.15); color: #b45309; }
.methods-copy-wrap { position: relative; }
.methods-copy-btn {
  position: absolute; top: 6px; right: 6px; z-index: 5;
  font-size: 11px; padding: 2px 8px;
}
.plot-card .card-header {
  font-size: 13px; color: #555; padding: 6px 12px;
  background: transparent; border-bottom: 1px solid #eee;
}
"

# Small helper: generate the navbar data-status badge from reactive state.
data_status_badge_html <- function(wide, has_sub) {
  if (nrow(wide) == 0)
    return(span(class = "data-badge empty", icon("circle"), " no data"))
  n_samples <- nrow(wide)
  gene_cols <- setdiff(names(wide), c("Sample", "Group", "Subgroup", "Ct_reference"))
  n_genes   <- length(gene_cols)
  n_groups  <- length(unique(wide$Group[nzchar(wide$Group)]))
  txt <- paste(n_samples, "samples \u00b7",
               n_genes, if (n_genes == 1) "gene" else "genes", "\u00b7",
               n_groups, "groups")
  if (has_sub) txt <- paste(txt, "\u00b7 2-factor")
  span(class = "data-badge", icon("check-circle"), " ", txt)
}

# Sidebar tooltip helper: defaults placement to "top" so tooltips stay
# inside the sidebar column and never overlap the main panel content.
sidebar_tip <- function(trigger, msg, placement = "top", ...) {
  bslib::tooltip(trigger, msg, placement = placement, ...)
}

# --- Sidebar (shared controls) ---------------------------------------------
# Reorganized into 6 numbered accordion sections. Section 1 (Data) is open
# by default; others are collapsed until needed. All input IDs preserved
# unchanged so server logic does not need to be touched.
shared_sidebar <- function() {
  tagList(
    # Global paste listener — moved up so it runs on any page.
    tags$script(HTML("
      document.addEventListener('paste', function(e) {
        var tag = document.activeElement ? document.activeElement.tagName.toLowerCase() : '';
        if (tag === 'input' || tag === 'textarea') return;
        var text = e.clipboardData.getData('text/plain');
        if (text && text.trim().length > 0) {
          Shiny.setInputValue('pasted_clipboard', text, {priority: 'event'});
          e.preventDefault();
        }
      });
    ")),

    # Guided vs Expert mode: guided opens Step 1 only (auto-collapse others);
    # expert opens all 6 steps so power users see everything at once.
    sidebar_tip(
      checkboxInput("guided_mode",
                    tagList(icon("graduation-cap"),
                            " Guided mode (walk through sections step by step)"),
                    value = TRUE),
      "Off = all sections open at once (expert view). On = start with Data only, expand others as you go."
    ),

    bslib::accordion(
      id = "sidebar_accordion",
      open = c("data"),  # only Step 1 open by default
      multiple = TRUE,

      # --- Step 1: Data -------------------------------------------------------
      bslib::accordion_panel(
        value = "data",
        title = tagList(icon("database"), " 1. Data"),
        div(class = "paste-dropzone",
            div(class = "pz-title", icon("clipboard"), " Paste your Excel data here"),
            div(class = "pz-sub",
                "Press ", tags$kbd("Ctrl"), "+", tags$kbd("V"),
                " anywhere on the page. Header: ",
                code("Sample | Group | [Subgroup] | Ct_ref | Gene1 | ..."))
        ),
        fluidRow(
          column(6, actionButton("load_ex_1factor",
                                 "Load 1-factor example",
                                 icon = icon("flask"),
                                 style = "width:100%;")),
          column(6, actionButton("load_ex_2factor",
                                 "Load 2-factor example",
                                 icon = icon("flask-vial"),
                                 style = "width:100%;"))
        ),
        uiOutput("data_summary_ui"),
        DTOutput("table"),
        br(),
        fluidRow(
          column(6, actionButton("add_row",   "Add row",
                                 icon = icon("plus"), style = "width:100%;")),
          column(6, actionButton("clear_tbl", "Clear table",
                                 icon = icon("trash"), style = "width:100%;"))
        ),
        uiOutput("error_msg")
      ),

      # --- Step 2: Experiment -------------------------------------------------
      bslib::accordion_panel(
        value = "experiment",
        title = tagList(icon("vial"), " 2. Experiment"),
        sidebar_tip(
          uiOutput("ctrl_group_ui"),
          "The group used as the reference for log\u2082 fold change. Choose from your pasted data."
        ),
        sidebar_tip(
          uiOutput("ctrl_subgroup_ui"),
          "Optional. Pin the baseline to a specific Group\u00d7Subgroup cell (e.g. Young+Control)."
        ),
        sidebar_tip(
          radioButtons("rep_mode", "Replicate type",
                       choices = c("Biological (one row per sample)" = "biological",
                                   "Technical (average rows sharing Sample ID)" = "technical"),
                       selected = "biological"),
          "Biological: each row is a distinct sample. Technical: rows sharing Sample are averaged first (MIQE-compliant)."
        ),
        uiOutput("layout_question_ui")
      ),

      # --- Step 3: Statistics -------------------------------------------------
      bslib::accordion_panel(
        value = "stats",
        title = tagList(icon("square-root-variable"), " 3. Statistics"),
        sidebar_tip(
          radioButtons("test_type", "Statistical test",
                       choices = c("Parametric (t-test / ANOVA)"              = "parametric",
                                   "Non-parametric (Wilcoxon / Kruskal)"      = "nonparametric",
                                   "Dunnett-style: vs control only (Holm)"    = "parametric_vs_control",
                                   "Mann-Whitney vs control only (Holm)"      = "nonparametric_vs_control"),
                       selected = "parametric"),
          "Dunnett/Holm options compare every non-control group to the control only \u2014 Prism's default for dose-response."
        ),
        sidebar_tip(
          radioButtons("sig_format", "P-value display",
                       choices = c("Exact (p = 0.0123)" = "exact",
                                   "Stars (*, **, ***, ****)" = "stars"),
                       selected = "exact"),
          "Some journals require star notation instead of exact values."
        ),
        checkboxInput("paired", "Paired samples (2 groups only)", value = FALSE)
      ),

      # --- Step 4: Plot -------------------------------------------------------
      bslib::accordion_panel(
        value = "plot",
        title = tagList(icon("chart-column"), " 4. Plot"),
        radioButtons("plot_type", "Plot type",
                     choices = c("Scatter + error bars" = "scatter",
                                 "Column + dots" = "column"),
                     inline = TRUE),
        radioButtons("plot_layout", "Plot layout",
                     choices = c("Individual" = "individual",
                                 "Combined grid" = "combined"),
                     selected = "individual"),
        conditionalPanel(
          "input.plot_layout == 'combined'",
          numericInput("facet_ncol", "Columns in grid",
                       value = 3, min = 1, max = 6, step = 1)
        ),
        sidebar_tip(
          radioButtons("error_type", "Error bars",
                       choices = c("SD" = "SD", "SEM" = "SEM", "95% CI" = "CI95"),
                       inline = TRUE),
          "95 % CI uses the t-distribution (Prism default)."
        ),
        checkboxInput("show_sig", "Show significance bars", value = TRUE),
        checkboxInput("rotate_x", "Rotate x-axis labels 45\u00b0", value = FALSE),
        sliderInput("aspect_ratio", "Plot aspect ratio (height/width)",
                    min = 0.5, max = 2.5, value = 1.0, step = 0.1),
        fluidRow(
          column(6, numericInput("y_min", "Y-axis min", value = NA, step = 0.5)),
          column(6, numericInput("y_max", "Y-axis max", value = NA, step = 0.5))
        ),
        uiOutput("hide_comps_ui")
      ),

      # --- Step 5: Styling ----------------------------------------------------
      bslib::accordion_panel(
        value = "styling",
        title = tagList(icon("palette"), " 5. Styling"),
        sidebar_tip(
          selectInput("font_family", "Font family",
                      choices = c("Default (sans)" = "",
                                  "Serif"          = "serif",
                                  "Mono"           = "mono",
                                  "Arial"          = "Arial",
                                  "Helvetica"      = "Helvetica",
                                  "Times New Roman"= "Times New Roman",
                                  "Courier New"    = "Courier New"),
                      selected = ""),
          "Generic families always work; named fonts depend on OS availability."
        ),
        sliderInput("ts_title",      "Plot title size",     6, 28, 18, 1),
        sliderInput("ts_axis_title", "Axis title size",     6, 24, 14, 1),
        sliderInput("ts_axis_text",  "Axis tick text size", 6, 22, 12, 1),
        sliderInput("ts_legend",     "Legend text size",    6, 20, 12, 1),
        sliderInput("ts_facet",      "Facet (gene) title size", 6, 24, 14, 1),
        sliderInput("ts_sig_bar",    "Sig-bar text size",   2, 10, 3.5, 0.5),
        uiOutput("color_picker_ui"),
        actionButton("reset_colors", "Reset colors to Okabe-Ito",
                     icon = icon("rotate-left"),
                     style = "width:100%; margin-top:6px;")
      ),

      # --- Step 6: Export -----------------------------------------------------
      bslib::accordion_panel(
        value = "export",
        title = tagList(icon("download"), " 6. Export size"),
        p(style = "font-size:11px; color:#666; margin-bottom:4px;",
          "Exact dimensions for manuscript figures. Common widths: ",
          strong("85 mm"), " (single column), ", strong("174 mm"),
          " (double column), ", strong("89/183 mm"), " (Nature)."),
        fluidRow(
          column(6, numericInput("export_w_mm", "Width (mm)",
                                 value = 174, min = 30, max = 400, step = 1)),
          column(6, numericInput("export_h_mm", "Height (mm)",
                                 value = 120, min = 30, max = 400, step = 1))
        ),
        radioButtons("export_dpi", "DPI (PNG/TIFF)",
                     choices = c("300" = 300, "600" = 600, "1200" = 1200),
                     selected = 300, inline = TRUE),
        p(style = "font-size:11px; color:#666;",
          em("SVG and PDF are vector \u2014 DPI ignored."))
      )
    )
  )
}

# --- Help tab content ------------------------------------------------------
help_content <- function() {
  tagList(
    h2("qPCR Analysis \u2014 Help & Documentation"),
    h3("Which tab should I use?"),
    tags$ul(
      tags$li(strong("Analysis"), " \u2014 standard qPCR. Each sample belongs to one group. ",
              "Different animals per condition."),
      tags$li(strong("Repeated Measures / Animal"), " \u2014 same animal measured across ",
              "multiple conditions (paired design). Samples repeat in 2+ groups."),
      tags$li(strong("Help"), " \u2014 this page.")
    ),
    h3("Expected input format \u2014 Analysis tab"),
    p("Paste from Excel (tab-separated) with this header row:"),
    tags$pre(
      "Sample   Group    Reference gene   Gene1    Gene2\n",
      "mouse1   Control  13.50            20.51    25.40\n",
      "mouse2   Control  13.56            20.62    25.23\n",
      "mouse3   Control  13.64            20.53    25.28\n",
      "mouse4   Treated  17.60            20.54    25.24\n",
      "mouse5   Treated  17.58            20.51    26.25\n",
      "mouse6   Treated  18.38            21.02    26.39"
    ),
    p(em("Optional \u2014 two-factor design (e.g. Young/Aged \u00d7 Carrier/Thapsigargin):"),
      " add a ", code("Subgroup"), " column as the ", strong("3rd column"),
      " (between ", code("Group"), " and ", code("Reference gene"),
      "). The header must be exactly ", code("Subgroup"),
      " (or ", code("Factor2"), "). Example:"),
    tags$pre(
      "Sample   Group   Subgroup  Reference gene   Gene1   Gene2\n",
      "cell1    Young   Carrier   13.50            20.51   25.40\n",
      "cell2    Young   Carrier   13.56            20.62   25.23\n",
      "cell3    Young   Carrier   13.64            20.53   25.28\n",
      "cell4    Young   Thaps     17.60            17.20   22.10\n",
      "cell5    Young   Thaps     17.58            17.15   22.30\n",
      "cell6    Young   Thaps     18.38            17.40   22.40\n",
      "cell7    Aged    Carrier   13.55            20.35   25.15\n",
      "cell8    Aged    Carrier   13.60            20.40   25.20\n",
      "cell9    Aged    Carrier   13.65            20.30   25.10\n",
      "cell10   Aged    Thaps     17.62            16.10   21.25\n",
      "cell11   Aged    Thaps     17.55            16.25   21.40\n",
      "cell12   Aged    Thaps     18.40            16.35   21.55"
    ),
    p("With a ", code("Subgroup"), " column present, the app automatically runs ",
      strong("two-way ANOVA"), " (main effects of Group, Subgroup, and the Group\u00d7Subgroup interaction) ",
      "plus ", strong("Tukey HSD"), " pairwise tests between all 4 interaction cells. ",
      "The plot x-axis shows the 4 combinations (", code("Young | Carrier"),
      ", ", code("Young | Thaps"), ", ...) and dots are colored by Subgroup."),
    p(strong("Control subgroup (optional):"), " when a Subgroup column is detected, a ",
      code("Control subgroup"), " selector appears below the Control group. ",
      "Leave it at ", em("'(use Group mean - mix all subgroups)'"),
      " to compute log\u2082 fold change against the mean of all samples in the control group. ",
      "Pick a specific subgroup (e.g. ", code("Carrier"), ") to use ",
      em("only the Young+Carrier cells"), " as the fold-change baseline \u2014 the biologically ",
      "typical choice for Young/Aged \u00d7 vehicle/treatment designs."),
    h3("Expected input format \u2014 Repeated Measures tab"),
    p("Same header, but ", strong("each Sample name must appear in 2 or more Groups"), ":"),
    tags$pre(
      "Sample   Group      Reference gene   Gene1    Gene2\n",
      "mouse1   Baseline   13.50            20.51    25.40\n",
      "mouse1   Treated    13.56            17.70    23.19\n",
      "mouse2   Baseline   13.64            20.53    25.28\n",
      "mouse2   Treated    14.22            17.67    23.12\n",
      "mouse3   Baseline   13.60            20.55    25.30\n",
      "mouse3   Treated    14.12            17.67    23.10"
    ),
    p("Each animal gets a unique color. Lines connect the same animal across conditions (before-after style)."),
    h3("Statistical methods"),
    tags$ul(
      tags$li(strong("\u0394\u0394Ct method"), " (Livak & Schmittgen, 2001): ",
              code("\u0394Ct = Ct_target \u2212 Ct_reference"), "; ",
              code("\u0394\u0394Ct = \u0394Ct \u2212 mean(\u0394Ct_control)"), "; ",
              code("Fold change = 2^(\u2212\u0394\u0394Ct)"), "."),
      tags$li(strong("Welch's t-test"), " \u2014 2 groups, assumes unequal variances (default)."),
      tags$li(strong("Paired t-test"), " \u2014 2 groups, matched samples."),
      tags$li(strong("Mann-Whitney U"), " \u2014 2 groups, non-parametric."),
      tags$li(strong("Wilcoxon signed-rank"), " \u2014 2 groups, paired non-parametric."),
      tags$li(strong("One-way ANOVA + Tukey HSD"), " \u2014 3+ groups, all pairwise comparisons."),
      tags$li(strong("Kruskal-Wallis + Dunn (Bonferroni)"), " \u2014 3+ groups, non-parametric."),
      tags$li(strong("Two-way ANOVA"), " \u2014 when Subgroup column present; tests Group, Subgroup, and interaction."),
      tags$li(strong("Dunnett-style vs-control (Holm-Bonferroni)"), " \u2014 compares every non-control group to the control only (like Prism's default for dose-response). Uses pooled-variance t-test with Holm adjustment; statistically similar to a true Dunnett's test but built-in (no external package)."),
      tags$li(strong("Mann-Whitney vs control (Holm)"), " \u2014 non-parametric version of the above.")
    ),
    h3("Replicate type (biological vs technical)"),
    p("At the top of the sidebar, pick how each row in your paste should be interpreted:"),
    tags$ul(
      tags$li(strong("Biological replicates (default)"), " \u2014 each row is a distinct biological sample (animal, well, culture). ",
              "The app uses all rows as independent observations for stats."),
      tags$li(strong("Technical replicates"), " \u2014 rows sharing the same ",
              code("Sample"), " + ", code("Group"), " (+ ", code("Subgroup"),
              " if present) are treated as technical repeats of one biological sample. ",
              "The app averages their ", code("Ct_target"), " and ", code("Ct_reference"),
              " before computing \u0394Ct, giving one observation per biological sample (MIQE-compliant).")
    ),
    p(em("Tip:"), " most users pre-average tech reps in Excel before pasting, so ",
      strong("Biological"), " is usually the right choice. Switch to ",
      strong("Technical"), " only when your paste contains 2-3 rows per sample."),

    h3("Two-factor layout (Prism grouped style)"),
    p("When a ", code("Subgroup"), " column is detected, the plot uses the ",
      strong("Prism grouped layout"), ": one factor drives the x-axis and the other ",
      "appears as ", em("dodged colored clusters"), " within each x-axis cluster. ",
      "The axis shows ", strong("nested labels"),
      " \u2014 Subgroup names directly under each dot, Group names centered under each cluster. ",
      "A blue panel appears at the top of the sidebar asking which factor drives the x-axis ",
      "(like Prism's new-graph wizard). Example for Young/Aged \u00d7 Carrier/Thapsigargin ",
      "with Group on x-axis:"),
    tags$pre(
      "     [•]  [•]      [•]  [•]            <- dots\n",
      "  Carrier Thaps  Carrier Thaps         <- Subgroup labels\n",
      "     Young          Aged               <- Group labels"
    ),
    p("Significance brackets span the correct dodged positions (within-cluster, ",
      "across-cluster same-subgroup, or diagonal comparisons)."),

    h3("Manuscript-ready export"),
    p("Expand the ", strong("Export size (for downloads)"), " panel in the sidebar to set exact figure dimensions in ", strong("millimeters"), " and choose ", strong("DPI"), " (for PNG/TIFF):"),
    tags$ul(
      tags$li(code("85 mm"), " \u2014 most journals' single column"),
      tags$li(code("174 mm"), " \u2014 most journals' double column"),
      tags$li(code("89 mm"), " / ", code("183 mm"), " \u2014 Nature family"),
      tags$li(code("300 dpi"), " minimum for most journals; ", code("600 dpi"), " for Cell/Nature; ",
              code("1200 dpi"), " for line art in some cases."),
      tags$li(em("SVG / PDF are vector formats \u2014 DPI is ignored; size scales without loss."))
    ),

    h3("Plot styling"),
    tags$ul(
      tags$li(strong("Font family"), " \u2014 pick Arial, Helvetica, Times New Roman, or generic sans/serif/mono."),
      tags$li(strong("Text size sliders"), " \u2014 independent control over plot title, axis titles, tick labels, legend, facet (gene) titles, and significance-bar text."),
      tags$li(strong("Per-group hex colors"), " \u2014 one ", code("#RRGGBB"), " textInput per detected Group (or Subgroup in two-factor designs). Invalid hex falls back to the default Okabe-Ito palette."),
      tags$li(strong("P-value display"), " \u2014 switch between exact p-values (", code("p = 0.0042"), ") and stars (", code("ns / * / ** / *** / ****"), "); affects both the plot and the stats table."),
      tags$li(strong("Error bars"), " \u2014 SD, SEM, or 95 % CI (t-distribution)."),
      tags$li(strong("Export formats"), " \u2014 PNG (300 dpi), PDF (vector), SVG (vector), and TIFF (LZW-compressed, 300 dpi). Use SVG or PDF for journals that require vector figures.")
    ),
    h3("Tips"),
    tags$ul(
      tags$li("Set the ", code("Control group name"), " to match exactly what's in your Group column."),
      tags$li("Use the ", strong("aspect ratio slider"), " to change plot shape: 1.0 = square (Prism-like), <1 = wider, >1 = taller."),
      tags$li("Toggle individual significance bars via the ", strong("Significance bars to show"), " checkboxes. Use Select/Deselect all for bulk changes."),
      tags$li("For publication, use ", strong("Combined grid"), " layout to produce one multi-panel figure."),
      tags$li("PNG export is 600 dpi; PDF is vector (infinite resolution).")
    ),
    h3("References"),
    tags$ul(
      tags$li("Livak KJ, Schmittgen TD (2001). Analysis of relative gene expression data using real-time quantitative PCR and the 2^(\u2212\u0394\u0394Ct) method. ", em("Methods"), " 25(4):402-408."),
      tags$li("Bustin SA et al. (2009). The MIQE guidelines: minimum information for publication of quantitative real-time PCR experiments. ", em("Clin Chem"), " 55(4):611-622.")
    ),
    hr(),
    p(style = "color:#888;",
      "Built with R Shiny + Shinylive. Source: ",
      a("github.com/karthi-ntu/qPCR_analysis", href = "https://github.com/karthi-ntu/qPCR_analysis", target = "_blank"))
  )
}

# ---- UI -------------------------------------------------------------------
# bslib theme: Bootstrap 5, lab-teal primary color, Inter font from Google.
# Supports light/dark mode via `input_dark_mode()`.
app_theme <- bslib::bs_theme(
  version     = 5,
  primary     = "#00A8A8",
  "navbar-bg" = "#00A8A8",
  base_font   = bslib::font_google("Inter", local = FALSE),
  heading_font = bslib::font_google("Inter", local = FALSE)
)

# Panel used by the Analysis tab — main content area, reused in the layout.
analysis_main_panel <- function() {
  tagList(
    # Plot card (title + content slot)
    bslib::card(
      class = "plot-card",
      bslib::card_header(
        div(style = "display:flex; justify-content:space-between; align-items:center;",
            div(icon("chart-column"), strong(" Plot")),
            uiOutput("plot_metadata_chip", inline = TRUE))
      ),
      bslib::card_body(
        min_height = "420px",
        uiOutput("plots_ui")
      )
    ),

    # Statistics card
    bslib::card(
      bslib::card_header(icon("square-root-variable"), strong(" Statistics")),
      bslib::card_body(DTOutput("stats_table"))
    ),

    # Methods card with copy-to-clipboard button
    bslib::card(
      bslib::card_header(icon("file-lines"), strong(" Methods (copy-paste for papers)")),
      bslib::card_body(
        div(class = "methods-copy-wrap",
            actionButton("copy_methods", "Copy",
                         icon = icon("copy"),
                         class = "btn-sm btn-outline-primary methods-copy-btn"),
            verbatimTextOutput("methods_text")
        )
      )
    ),

    # Downloads card
    bslib::card(
      bslib::card_header(icon("download"), strong(" Downloads")),
      bslib::card_body(
        fluidRow(
          column(2, downloadButton("dl_png",     "PNG",     icon = icon("image"))),
          column(2, downloadButton("dl_pdf",     "PDF",     icon = icon("file-pdf"))),
          column(2, downloadButton("dl_svg",     "SVG",     icon = icon("bezier-curve"))),
          column(2, downloadButton("dl_tiff",    "TIFF",    icon = icon("image"))),
          column(2, downloadButton("dl_csv",     "Stats",   icon = icon("table"))),
          column(2, downloadButton("dl_methods", "Methods", icon = icon("file-lines")))
        )
      )
    )
  )
}

# Panel used by the Repeated Measures tab.
rm_main_panel <- function() {
  tagList(
    bslib::card(
      class = "plot-card",
      bslib::card_header(icon("link"), strong(" Paired plot (animal-matched)")),
      bslib::card_body(
        min_height = "420px",
        uiOutput("rm_status"),
        uiOutput("rm_plot_ui")
      )
    ),
    bslib::card(
      bslib::card_header(icon("square-root-variable"), strong(" Paired statistics")),
      bslib::card_body(DTOutput("rm_stats_table"))
    ),
    bslib::card(
      bslib::card_header(icon("download"), strong(" Downloads")),
      bslib::card_body(
        fluidRow(
          column(2, downloadButton("rm_dl_png",  "PNG",   icon = icon("image"))),
          column(2, downloadButton("rm_dl_pdf",  "PDF",   icon = icon("file-pdf"))),
          column(2, downloadButton("rm_dl_svg",  "SVG",   icon = icon("bezier-curve"))),
          column(2, downloadButton("rm_dl_tiff", "TIFF",  icon = icon("image"))),
          column(3, downloadButton("rm_dl_csv",  "Stats (CSV)",
                                   icon = icon("table")))
        )
      )
    )
  )
}

rm_sidebar <- function() {
  tagList(
    p(style = "font-size:13px;",
      icon("info-circle"), " This tab uses the same data pasted in the ",
      strong("Analysis"), " tab, but treats ",
      strong("each Sample name as one animal"),
      ". Each Sample must appear in 2+ Groups."),
    hr(),
    h6(icon("link"), " RM plot style"),
    p(style = "font-size:12px; color:#888;",
      "Dots are colored by animal; thin lines connect the same animal across groups.")
  )
}

ui <- bslib::page_navbar(
  title = tagList(icon("dna"), " qPCR Analysis"),
  window_title = "qPCR Analysis \u2014 \u0394\u0394Ct",
  theme = app_theme,
  header = tags$head(
    tags$style(HTML(app_css)),
    tags$script(HTML("
      Shiny.addCustomMessageHandler('copy_to_clipboard', function(msg) {
        if (navigator.clipboard && navigator.clipboard.writeText) {
          navigator.clipboard.writeText(msg.text);
        } else {
          // Legacy fallback
          var ta = document.createElement('textarea');
          ta.value = msg.text; document.body.appendChild(ta);
          ta.select(); document.execCommand('copy');
          document.body.removeChild(ta);
        }
      });
    "))
  ),

  bslib::nav_panel(
    tagList(icon("chart-column"), " Analysis"),
    bslib::layout_sidebar(
      sidebar = bslib::sidebar(
        width  = 380,
        open   = TRUE,
        shared_sidebar()
      ),
      analysis_main_panel()
    )
  ),

  bslib::nav_panel(
    tagList(icon("link"), " Repeated Measures"),
    bslib::layout_sidebar(
      sidebar = bslib::sidebar(width = 320, open = TRUE, rm_sidebar()),
      rm_main_panel()
    )
  ),

  bslib::nav_panel(
    tagList(icon("circle-question"), " Help"),
    div(style = "padding: 20px; max-width: 960px; margin: auto;",
        help_content())
  ),

  # Right-side navbar items: data status badge, dark mode toggle, GitHub link.
  bslib::nav_spacer(),
  bslib::nav_item(uiOutput("data_status_badge", inline = TRUE)),
  bslib::nav_item(bslib::input_dark_mode(id = "dark_mode")),
  bslib::nav_item(
    tags$a(icon("github"), " GitHub",
           href = "https://github.com/karthi-ntu/qPCR_analysis",
           target = "_blank",
           style = "color: white; margin-right: 10px;")
  )
)

# ---- Server ---------------------------------------------------------------
server <- function(input, output, session) {

  rv <- reactiveValues(wide = empty_wide())

  long_data    <- reactive({ to_long_df(rv$wide) })
  has_subgroup <- reactive({ "Subgroup" %in% names(rv$wide) })

  # -- Top-bar data status badge -------------------------------------------
  output$data_status_badge <- renderUI({
    data_status_badge_html(rv$wide, has_subgroup())
  })

  # -- Data summary card (shows above the data table after paste) ---------
  output$data_summary_ui <- renderUI({
    wide <- rv$wide
    if (nrow(wide) == 0) {
      return(div(class = "data-summary empty",
                 icon("circle-info"), " No data yet. Paste from Excel or load an example above."))
    }
    gene_cols <- setdiff(names(wide), c("Sample", "Group", "Subgroup", "Ct_reference"))
    groups    <- unique(wide$Group[nzchar(wide$Group)])
    sub_txt <- if (has_subgroup()) {
      subs <- unique(wide$Subgroup[nzchar(wide$Subgroup)])
      paste0(" \u00b7 Subgroups: ", paste(subs, collapse = ", "))
    } else ""
    na_cells <- sum(sapply(gene_cols, function(g) sum(is.na(wide[[g]]))))
    na_msg <- if (na_cells > 0)
      span(class = "text-warning", " \u00b7 ", na_cells, " NA value(s)")
    cls <- if (na_cells > 0) "data-summary warn" else "data-summary"
    div(class = cls,
        icon("check-circle"), " ",
        strong(nrow(wide), "samples"), " \u00b7 ",
        strong(length(gene_cols)), if (length(gene_cols) == 1) " gene" else " genes",
        " (", paste(gene_cols, collapse = ", "), ") \u00b7 ",
        strong(length(groups)), " groups (", paste(groups, collapse = ", "), ")",
        sub_txt, na_msg)
  })

  # -- Guided vs Expert mode: open/close accordion panels ------------------
  # When Guided is on, only Step 1 stays open. When off, all 6 open at once.
  observeEvent(input$guided_mode, {
    all_steps <- c("data", "experiment", "stats", "plot", "styling", "export")
    if (isTRUE(input$guided_mode)) {
      # Close everything except Data.
      for (s in setdiff(all_steps, "data"))
        bslib::accordion_panel_close("sidebar_accordion", values = s)
      bslib::accordion_panel_open("sidebar_accordion", values = "data")
    } else {
      for (s in all_steps)
        bslib::accordion_panel_open("sidebar_accordion", values = s)
    }
  }, ignoreInit = TRUE)

  # -- Example data loaders --------------------------------------------------
  observeEvent(input$load_ex_1factor, {
    rv$wide <- example_data_1factor()
    showNotification("Loaded 1-factor example (Control vs Treated, 2 genes)",
                     type = "message", duration = 3)
  })
  observeEvent(input$load_ex_2factor, {
    rv$wide <- example_data_2factor()
    showNotification("Loaded 2-factor example (Young/Aged \u00d7 Control/Thaps)",
                     type = "message", duration = 3)
  })

  # -- Copy-to-clipboard for methods text -----------------------------------
  observeEvent(input$copy_methods, {
    txt <- methods_txt()
    req(nzchar(txt))
    session$sendCustomMessage("copy_to_clipboard", list(text = txt))
    showNotification("Methods copied to clipboard",
                     type = "message", duration = 2)
  })

  # -- Shared table ----------------------------------------------------------
  output$table <- renderDT({
    datatable(
      rv$wide,
      editable  = TRUE,
      rownames  = FALSE,
      selection = "none",
      options   = list(dom = "t", pageLength = 50, scrollY = "300px")
    )
  })

  observeEvent(input$table_cell_edit, {
    rv$wide <- editData(rv$wide, input$table_cell_edit, proxy = "table",
                        rownames = FALSE)
    num_cols <- setdiff(names(rv$wide), c("Sample", "Group", "Subgroup"))
    for (col in num_cols)
      rv$wide[[col]] <- suppressWarnings(as.numeric(rv$wide[[col]]))
  })

  observeEvent(input$add_row, {
    mcols <- meta_cols(rv$wide)
    gene_cols <- setdiff(names(rv$wide), mcols)
    new_row <- data.frame(Sample = "", Group = "", stringsAsFactors = FALSE)
    if ("Subgroup" %in% mcols) new_row$Subgroup <- ""
    new_row$Ct_reference <- NA_real_
    for (g in gene_cols) new_row[[g]] <- NA_real_
    new_row <- new_row[, names(rv$wide), drop = FALSE]
    rv$wide <- rbind(rv$wide, new_row)
  })

  observeEvent(input$clear_tbl, { rv$wide <- empty_wide() })

  observeEvent(input$pasted_clipboard, {
    txt <- input$pasted_clipboard
    req(!is.null(txt), nchar(trimws(txt)) > 0)
    txt   <- gsub("\r\n", "\n", txt)
    txt   <- gsub("\r",   "\n", txt)
    lines <- strsplit(txt, "\n")[[1]]
    lines <- lines[nchar(trimws(lines)) > 0]
    req(length(lines) >= 2)

    headers <- trimws(strsplit(lines[1], "\t")[[1]])
    req(length(headers) >= 4)

    h3_low  <- tolower(headers[3])
    has_sub <- h3_low %in% c("subgroup", "factor2", "factor 2")
    meta_n  <- if (has_sub) 4 else 3
    gene_names <- headers[(meta_n + 1):length(headers)]
    req(length(gene_names) >= 1)

    rows <- lapply(lines[-1], function(line) {
      fields <- trimws(strsplit(line, "\t")[[1]])
      if (length(fields) < meta_n + 1) return(NULL)
      ct_ref <- suppressWarnings(as.numeric(fields[meta_n]))
      lapply(seq_along(gene_names), function(i) {
        r <- data.frame(
          Sample       = fields[1],
          Group        = fields[2],
          Gene         = gene_names[i],
          Ct_target    = suppressWarnings(as.numeric(fields[meta_n + i])),
          Ct_reference = ct_ref,
          stringsAsFactors = FALSE
        )
        if (has_sub) r$Subgroup <- fields[3]
        r
      })
    })

    parsed <- do.call(rbind, Filter(Negate(is.null), unlist(rows, recursive = FALSE)))
    if (!is.null(parsed) && nrow(parsed) > 0)
      rv$wide <- to_wide_df(parsed)
  })

  # -- Per-group color pickers (dynamic based on detected groups) ------------
  # OKABE_ITO is sourced from R/plot.R; duplicate the hex list here so we can
  # populate the default hex values in the textInputs.
  OKABE_ITO_UI <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442",
                    "#0072B2", "#D55E00", "#CC79A7", "#000000")

  # The levels that the plot colors: Subgroup when has_subgroup, else Group.
  fill_levels <- reactive({
    df <- rv$wide
    if (nrow(df) == 0) return(character(0))
    if (has_subgroup() && "Subgroup" %in% names(df)) {
      lv <- unique(df$Subgroup); lv[nzchar(lv)]
    } else {
      lv <- unique(df$Group);    lv[nzchar(lv)]
    }
  })

  output$color_picker_ui <- renderUI({
    lvs <- fill_levels()
    if (length(lvs) == 0) return(p(em("(add data to set colors)"),
                                   style = "font-size:11px; color:#888;"))
    palette <- rep_len(OKABE_ITO_UI, length(lvs))
    tagList(
      tags$label("Group colors (hex)", style = "font-weight:bold;"),
      lapply(seq_along(lvs), function(i) {
        level <- lvs[i]
        current <- isolate(input[[paste0("color_", make.names(level))]])
        if (is.null(current) || !nzchar(current)) current <- palette[i]
        fluidRow(
          column(5, div(style = paste0(
            "height:28px; margin-top:18px; border:1px solid #888; border-radius:3px; ",
            "background:", current, ";"), ""),
            tags$small(level, style = "display:block; font-size:10px; overflow:hidden;")),
          column(7, textInput(paste0("color_", make.names(level)),
                              label = NULL, value = current,
                              placeholder = "#RRGGBB"))
        )
      })
    )
  })

  # Reset all pickers to the Okabe-Ito defaults.
  observeEvent(input$reset_colors, {
    lvs <- fill_levels()
    palette <- rep_len(OKABE_ITO_UI, length(lvs))
    for (i in seq_along(lvs)) {
      updateTextInput(session, paste0("color_", make.names(lvs[i])),
                      value = palette[i])
    }
  })

  # Build the named override vector passed to plot functions.
  color_override <- reactive({
    lvs <- fill_levels()
    if (length(lvs) == 0) return(NULL)
    vals <- vapply(lvs, function(lv) {
      v <- input[[paste0("color_", make.names(lv))]]
      if (is.null(v)) NA_character_ else v
    }, character(1))
    names(vals) <- lvs
    vals[!is.na(vals)]
  })

  text_sizes <- reactive({
    list(
      title      = input$ts_title      %||% 18,
      axis_title = input$ts_axis_title %||% 14,
      axis_text  = input$ts_axis_text  %||% 12,
      legend     = input$ts_legend     %||% 12,
      facet      = input$ts_facet      %||% 14,
      sig_bar    = input$ts_sig_bar    %||% 3.5
    )
  })

  # -- 2-factor plot layout question (Prism-wizard style) --------------------
  # Only shown when a Subgroup column is detected. Lets the user pick which
  # factor drives the x-axis and which is the dodged/colored grouping.
  output$layout_question_ui <- renderUI({
    if (!has_subgroup()) return(NULL)
    tagList(
      tags$div(style = "padding:6px 8px; margin:4px 0; background:#f4f8ff; border-left:3px solid #4a7bd4; border-radius:3px;",
        tags$b("Two-factor design detected."),
        tags$p(style = "font-size:12px; margin:4px 0 0 0;",
               "Prism-grouped layout is used: one factor on the x-axis, the other shown as dodged colored clusters with nested labels below.")
      ),
      radioButtons("x_axis_var", "Which factor on the x-axis?",
                   choices = c("Group (e.g. Young / Aged)"              = "Group",
                               "Subgroup (e.g. Carrier / Thapsigargin)" = "Subgroup"),
                   selected = "Group")
    )
  })

  # -- Control group dropdown (populated from detected Group column) --------
  # Falls back to a disabled placeholder when no data has been pasted yet.
  # Preserves the user's current selection on re-render via isolate().
  output$ctrl_group_ui <- renderUI({
    grps <- unique(rv$wide$Group)
    grps <- grps[nzchar(grps)]
    if (length(grps) == 0) {
      return(selectInput("ctrl_group", "Control group name",
                         choices = c("(paste data first)" = ""),
                         selected = ""))
    }
    current <- isolate(input$ctrl_group)
    sel <- if (!is.null(current) && nzchar(current) && current %in% grps)
             current else grps[1]
    selectInput("ctrl_group", "Control group name",
                choices = grps, selected = sel)
  })

  # -- Control subgroup selector (only appears when Subgroup column present) --
  output$ctrl_subgroup_ui <- renderUI({
    if (!has_subgroup()) return(NULL)
    subs <- unique(rv$wide$Subgroup)
    subs <- subs[nzchar(subs)]
    if (length(subs) == 0) return(NULL)
    current <- isolate(input$ctrl_subgroup)
    sel <- if (!is.null(current) && current %in% c("", subs)) current else ""
    selectInput(
      "ctrl_subgroup",
      "Control subgroup (optional)",
      choices = c("(use Group mean - mix all subgroups)" = "", subs),
      selected = sel
    )
  })

  # -- Validation ------------------------------------------------------------
  validation_msg <- reactive({
    df <- long_data()
    if (nrow(df) == 0) return("Add at least one gene column and enter data.")
    if (any(is.na(df$Ct_target) | is.na(df$Ct_reference)))
      return("All Ct values must be numeric.")
    # input$ctrl_group comes from a dynamic selectInput (renderUI) and is NULL
    # until the UI has rendered and Shiny has synced its value. Tolerate that
    # instead of throwing "argument is of length zero".
    ctrl <- input$ctrl_group %||% ""
    if (!nzchar(ctrl)) return("Waiting for control group selection\u2026")
    if (!ctrl %in% df$Group)
      return(paste0("Control group '", ctrl, "' not found in Group column."))
    grp_counts <- table(df$Gene, df$Group)
    if (any(grp_counts[grp_counts > 0] < 2))
      return("At least 2 samples per group per gene required.")
    NULL
  })

  output$error_msg <- renderUI({
    msg <- validation_msg()
    if (!is.null(msg))
      div(style = "color:red; margin-top:8px;", icon("exclamation-triangle"), " ", msg)
  })

  # -- Analysis (Tab 1) ------------------------------------------------------
  # The processed frame — applies technical-replicate averaging BEFORE delta-Ct
  # when the user has picked "technical replicates" mode. This matches MIQE
  # guidelines: tech reps collapse to one biological observation per sample.
  processed_long <- reactive({
    df <- long_data()
    if (identical(input$rep_mode, "technical")) df <- average_tech_reps(df)
    df
  })

  analyzed <- reactive({
    req(is.null(validation_msg()))
    df <- compute_delta_ct(processed_long())
    out <- compute_fold_change(df, control_group = input$ctrl_group,
                               control_subgroup = input$ctrl_subgroup)
    # Preserve tech_counts attribute through delta/fold-change for methods text.
    attr(out, "tech_counts") <- attr(processed_long(), "tech_counts")
    attr(out, "rep_mode")    <- input$rep_mode
    out
  })

  stats_result <- reactive({
    req(analyzed())
    run_stats(analyzed(), paired = input$paired,
              test_type     = input$test_type,
              has_subgroup  = has_subgroup(),
              control_group = input$ctrl_group)
  })

  genes_in_data <- reactive({ req(analyzed()); unique(analyzed()$Gene) })

  # -- Plot metadata chip (shown in the plot card header) ------------------
  output$plot_metadata_chip <- renderUI({
    if (is.null(validation_msg()) && !is.null(stats_result())) {
      test_name <- unique(stats_result()$Test)[1] %||% ""
      n_txt <- tryCatch(replicate_summary(analyzed(),
                                          rep_mode = input$rep_mode %||% "biological"),
                        error = function(e) "")
      tags$span(style = "font-size:12px; color:#666;",
                if (nzchar(n_txt)) paste(n_txt, " \u00b7 ", sep = "") else "",
                if (nzchar(test_name)) test_name else "")
    }
  })

  # -- Significance bar controls (with select/deselect all) -----------------
  output$hide_comps_ui <- renderUI({
    req(stats_result())
    sig <- stats_result()
    sig <- sig[sig$Significant == "Yes" & grepl(" vs ", sig$Comparison, fixed = TRUE), ]
    if (nrow(sig) == 0) return(NULL)
    choices <- unique(paste0(sig$Gene, ": ", sig$Comparison))
    tagList(
      tags$label("Significance bars to show", style = "font-weight:bold;"),
      fluidRow(
        column(6, actionButton("sig_select_all",   "\u2713 Select all",
                               style = "width:100%; font-size:11px;")),
        column(6, actionButton("sig_deselect_all", "\u2717 Deselect all",
                               style = "width:100%; font-size:11px;"))
      ),
      br(),
      checkboxGroupInput("visible_comps", label = NULL,
                         choices  = choices, selected = choices)
    )
  })

  observeEvent(input$sig_select_all, {
    sig <- stats_result()
    sig <- sig[sig$Significant == "Yes" & grepl(" vs ", sig$Comparison, fixed = TRUE), ]
    choices <- unique(paste0(sig$Gene, ": ", sig$Comparison))
    updateCheckboxGroupInput(session, "visible_comps", selected = choices)
  })

  observeEvent(input$sig_deselect_all, {
    updateCheckboxGroupInput(session, "visible_comps", selected = character(0))
  })

  visible_comps_for_gene <- function(gene) {
    sel <- input$visible_comps
    if (is.null(sel)) return(character(0))
    prefix <- paste0(gene, ": ")
    sub(prefix, "", sel[startsWith(sel, prefix)], fixed = TRUE)
  }

  # Label used on the y-axis and in methods text: composite when a control
  # subgroup is picked ("Young | Control"), otherwise just the group.
  baseline_label <- reactive({
    cs <- input$ctrl_subgroup
    if (!is.null(cs) && nzchar(cs) && has_subgroup()) {
      paste(input$ctrl_group, cs, sep = " | ")
    } else {
      input$ctrl_group
    }
  })

  # -- Plots (Tab 1) ---------------------------------------------------------
  output$plots_ui <- renderUI({
    req(genes_in_data())
    if (identical(input$plot_layout, "combined")) {
      n      <- length(genes_in_data())
      ncol   <- max(1, as.integer(input$facet_ncol %||% 3))
      nrow   <- ceiling(n / ncol)
      height <- paste0(max(350, 350 * nrow), "px")
      plotOutput("plot_combined", height = height)
    } else {
      do.call(tagList, lapply(genes_in_data(), function(gene)
        plotOutput(paste0("plot_", make.names(gene)), height = "450px")))
    }
  })

  output$plot_combined <- renderPlot({
    req(genes_in_data())
    make_combined_plot(analyzed(), stats_result(),
                       error_type    = input$error_type,
                       plot_type     = input$plot_type,
                       show_sig      = input$show_sig,
                       sig_comparisons_all = input$visible_comps,
                       control_group = baseline_label(),
                       rotate_x      = input$rotate_x,
                       aspect_ratio  = input$aspect_ratio,
                       ncol          = max(1, as.integer(input$facet_ncol %||% 3)),
                       has_subgroup  = has_subgroup(),
                       fill_override = color_override(),
                       text_sizes    = text_sizes(),
                       font_family   = input$font_family %||% "",
                       sig_format    = input$sig_format %||% "exact",
                       x_axis_var    = input$x_axis_var %||% "Group",
                       y_min = if (is.na(input$y_min)) NULL else input$y_min,
                       y_max = if (is.na(input$y_max)) NULL else input$y_max)
  })

  observe({
    req(genes_in_data())
    lapply(genes_in_data(), function(gene) {
      local({
        g <- gene
        output[[paste0("plot_", make.names(g))]] <- renderPlot({
          make_barplot(analyzed(), stats_result(), gene = g,
                       error_type    = input$error_type,
                       plot_type     = input$plot_type,
                       show_sig      = input$show_sig,
                       sig_comparisons = visible_comps_for_gene(g),
                       control_group = baseline_label(),
                       rotate_x      = input$rotate_x,
                       aspect_ratio  = input$aspect_ratio,
                       has_subgroup  = has_subgroup(),
                       fill_override = color_override(),
                       text_sizes    = text_sizes(),
                       font_family   = input$font_family %||% "",
                       sig_format    = input$sig_format %||% "exact",
                       x_axis_var    = input$x_axis_var %||% "Group",
                       y_min = if (is.na(input$y_min)) NULL else input$y_min,
                       y_max = if (is.na(input$y_max)) NULL else input$y_max)
        })
      })
    })
  })

  # -- Stats / Normality tables (Tab 1) --------------------------------------
  output$stats_table <- renderDT({
    req(stats_result())
    df <- stats_result()
    fmt <- pick_sig_formatter(input$sig_format %||% "exact")
    df$p_value <- fmt(df$p_value)
    datatable(df, rownames = FALSE, options = list(dom = "t")) |>
      formatStyle("Significant", target = "row",
                  backgroundColor = styleEqual("Yes", "#d4edda"))
  })

  methods_txt <- reactive({
    req(stats_result(), analyzed())
    generate_methods_text(
      df = analyzed(), stats_df = stats_result(),
      test_type = input$test_type, paired = input$paired,
      control_group = baseline_label(), has_subgroup = has_subgroup(),
      rep_mode = input$rep_mode %||% "biological"
    )
  })

  output$methods_text <- renderText({ methods_txt() })

  # -- Downloads (Tab 1) -----------------------------------------------------
  render_individual_plots <- function() {
    lapply(genes_in_data(), function(g)
      make_barplot(analyzed(), stats_result(), gene = g,
                   error_type    = input$error_type,
                   plot_type     = input$plot_type,
                   show_sig      = input$show_sig,
                   sig_comparisons = visible_comps_for_gene(g),
                   control_group = baseline_label(),
                   rotate_x      = input$rotate_x,
                   aspect_ratio  = input$aspect_ratio,
                   has_subgroup  = has_subgroup(),
                   fill_override = color_override(),
                   text_sizes    = text_sizes(),
                   font_family   = input$font_family %||% "",
                   sig_format    = input$sig_format %||% "exact",
                   x_axis_var        = input$x_axis_var        %||% "Group",
                   y_min = if (is.na(input$y_min)) NULL else input$y_min,
                   y_max = if (is.na(input$y_max)) NULL else input$y_max))
  }

  render_combined_plot <- function() {
    make_combined_plot(analyzed(), stats_result(),
                       error_type    = input$error_type,
                       plot_type     = input$plot_type,
                       show_sig      = input$show_sig,
                       sig_comparisons_all = input$visible_comps,
                       control_group = baseline_label(),
                       rotate_x      = input$rotate_x,
                       aspect_ratio  = input$aspect_ratio,
                       ncol          = max(1, as.integer(input$facet_ncol %||% 3)),
                       has_subgroup  = has_subgroup(),
                       fill_override = color_override(),
                       text_sizes    = text_sizes(),
                       font_family   = input$font_family %||% "",
                       sig_format    = input$sig_format %||% "exact",
                       x_axis_var        = input$x_axis_var        %||% "Group",
                       y_min = if (is.na(input$y_min)) NULL else input$y_min,
                       y_max = if (is.na(input$y_max)) NULL else input$y_max)
  }

  # Always render a single combined facet plot for download
  # (Shinylive-compatible: single ggplot + ggsave works reliably; raw png()/pdf()
  # devices can fail silently in WebAssembly.)
  download_plot <- function() {
    make_combined_plot(analyzed(), stats_result(),
                       error_type          = input$error_type,
                       plot_type           = input$plot_type,
                       show_sig            = input$show_sig,
                       sig_comparisons_all = input$visible_comps,
                       control_group       = baseline_label(),
                       rotate_x            = input$rotate_x,
                       aspect_ratio        = input$aspect_ratio,
                       ncol = max(1, as.integer(input$facet_ncol %||% 3)),
                       has_subgroup        = has_subgroup(),
                       fill_override       = color_override(),
                       text_sizes          = text_sizes(),
                       font_family         = input$font_family %||% "",
                       sig_format          = input$sig_format %||% "exact",
                       x_axis_var          = input$x_axis_var        %||% "Group",
                       y_min = if (is.na(input$y_min)) NULL else input$y_min,
                       y_max = if (is.na(input$y_max)) NULL else input$y_max)
  }

  # Download dimensions — manuscript-ready: use user-provided mm + DPI.
  download_dims <- function() {
    list(w_mm = input$export_w_mm %||% 174,
         h_mm = input$export_h_mm %||% 120,
         dpi  = as.integer(input$export_dpi %||% 300))
  }

  output$dl_png <- downloadHandler(
    filename    = function() paste0("qpcr_plot_", Sys.Date(), ".png"),
    contentType = "image/png",
    content     = function(file) {
      req(genes_in_data())
      d <- download_dims()
      ggsave(file, download_plot(),
             width = d$w_mm, height = d$h_mm, units = "mm",
             dpi = d$dpi, device = "png", bg = "white")
    }
  )

  output$dl_pdf <- downloadHandler(
    filename    = function() paste0("qpcr_plot_", Sys.Date(), ".pdf"),
    contentType = "application/pdf",
    content     = function(file) {
      req(genes_in_data())
      d <- download_dims()
      ggsave(file, download_plot(),
             width = d$w_mm, height = d$h_mm, units = "mm",
             device = "pdf", bg = "white")
    }
  )

  output$dl_svg <- downloadHandler(
    filename    = function() paste0("qpcr_plot_", Sys.Date(), ".svg"),
    contentType = "image/svg+xml",
    content     = function(file) {
      req(genes_in_data())
      d <- download_dims()
      ggsave(file, download_plot(),
             width = d$w_mm, height = d$h_mm, units = "mm",
             device = "svg", bg = "white")
    }
  )

  output$dl_tiff <- downloadHandler(
    filename    = function() paste0("qpcr_plot_", Sys.Date(), ".tiff"),
    contentType = "image/tiff",
    content     = function(file) {
      req(genes_in_data())
      d <- download_dims()
      ggsave(file, download_plot(),
             width = d$w_mm, height = d$h_mm, units = "mm",
             dpi = d$dpi, device = "tiff", bg = "white",
             compression = "lzw")
    }
  )

  output$dl_csv <- downloadHandler(
    filename    = function() paste0("qpcr_stats_", Sys.Date(), ".csv"),
    contentType = "text/csv",
    content     = function(file) {
      req(stats_result())
      write.csv(stats_result(), file, row.names = FALSE)
    }
  )

  output$dl_methods <- downloadHandler(
    filename    = function() paste0("qpcr_methods_", Sys.Date(), ".txt"),
    contentType = "text/plain",
    content     = function(file) {
      req(methods_txt())
      writeLines(methods_txt(), file)
    }
  )

  # -- Tab 2: Repeated Measures ----------------------------------------------
  rm_valid <- reactive({
    df <- long_data()
    if (nrow(df) == 0) return(list(ok = FALSE, msg = "No data. Paste your data in the Analysis tab first."))
    if (any(is.na(df$Ct_target) | is.na(df$Ct_reference)))
      return(list(ok = FALSE, msg = "Some Ct values are non-numeric."))
    if (!is_repeated_measures(df))
      return(list(ok = FALSE,
                  msg = "This data doesn't look like repeated measures. Each Sample appears in only one Group. Use the Analysis tab instead."))
    ctrl <- input$ctrl_group %||% ""
    if (!nzchar(ctrl))
      return(list(ok = FALSE, msg = "Waiting for control group selection\u2026"))
    if (!ctrl %in% df$Group)
      return(list(ok = FALSE,
                  msg = paste0("Control group '", ctrl, "' not found.")))
    list(ok = TRUE, msg = NULL)
  })

  output$rm_status <- renderUI({
    v <- rm_valid()
    if (!v$ok) div(style = "color:#b00; padding:10px; border:1px solid #f3c0c0; background:#fde8e8; border-radius:4px;",
                   icon("info-circle"), " ", v$msg)
  })

  rm_analyzed <- reactive({
    req(rm_valid()$ok)
    df <- compute_delta_ct(long_data())
    compute_fold_change(df, control_group = input$ctrl_group,
                        control_subgroup = input$ctrl_subgroup)
  })

  rm_stats <- reactive({
    req(rm_analyzed())
    run_paired_stats(rm_analyzed())
  })

  output$rm_plot_ui <- renderUI({
    req(rm_valid()$ok, rm_analyzed())
    n      <- length(unique(rm_analyzed()$Gene))
    ncol   <- max(1, as.integer(input$facet_ncol %||% 3))
    nrow   <- ceiling(n / ncol)
    height <- paste0(max(400, 400 * nrow), "px")
    plotOutput("rm_plot", height = height)
  })

  render_rm_plot <- function() {
    make_paired_plot(rm_analyzed(), rm_stats(),
                     error_type    = input$error_type,
                     show_sig      = input$show_sig,
                     sig_comparisons_all = input$visible_comps,
                     control_group = baseline_label(),
                     rotate_x      = input$rotate_x,
                     aspect_ratio  = input$aspect_ratio,
                     ncol          = max(1, as.integer(input$facet_ncol %||% 3)),
                     text_sizes    = text_sizes(),
                     font_family   = input$font_family %||% "",
                     sig_format    = input$sig_format %||% "exact",
                     y_min = if (is.na(input$y_min)) NULL else input$y_min,
                     y_max = if (is.na(input$y_max)) NULL else input$y_max)
  }

  output$rm_plot <- renderPlot({ render_rm_plot() })

  output$rm_stats_table <- renderDT({
    req(rm_stats())
    df <- rm_stats()
    df$p_value <- format_pvalue(df$p_value)
    datatable(df, rownames = FALSE, options = list(dom = "t")) |>
      formatStyle("Significant", target = "row",
                  backgroundColor = styleEqual("Yes", "#d4edda"))
  })

  rm_download_dims <- function() {
    list(w_mm = input$export_w_mm %||% 174,
         h_mm = input$export_h_mm %||% 120,
         dpi  = as.integer(input$export_dpi %||% 300))
  }

  output$rm_dl_png <- downloadHandler(
    filename    = function() paste0("qpcr_rm_plot_", Sys.Date(), ".png"),
    contentType = "image/png",
    content     = function(file) {
      req(rm_valid()$ok)
      d <- rm_download_dims()
      ggsave(file, render_rm_plot(),
             width = d$w_mm, height = d$h_mm, units = "mm",
             dpi = d$dpi, device = "png", bg = "white")
    }
  )

  output$rm_dl_pdf <- downloadHandler(
    filename    = function() paste0("qpcr_rm_plot_", Sys.Date(), ".pdf"),
    contentType = "application/pdf",
    content     = function(file) {
      req(rm_valid()$ok)
      d <- rm_download_dims()
      ggsave(file, render_rm_plot(),
             width = d$w_mm, height = d$h_mm, units = "mm",
             device = "pdf", bg = "white")
    }
  )

  output$rm_dl_svg <- downloadHandler(
    filename    = function() paste0("qpcr_rm_plot_", Sys.Date(), ".svg"),
    contentType = "image/svg+xml",
    content     = function(file) {
      req(rm_valid()$ok)
      d <- rm_download_dims()
      ggsave(file, render_rm_plot(),
             width = d$w_mm, height = d$h_mm, units = "mm",
             device = "svg", bg = "white")
    }
  )

  output$rm_dl_tiff <- downloadHandler(
    filename    = function() paste0("qpcr_rm_plot_", Sys.Date(), ".tiff"),
    contentType = "image/tiff",
    content     = function(file) {
      req(rm_valid()$ok)
      d <- rm_download_dims()
      ggsave(file, render_rm_plot(),
             width = d$w_mm, height = d$h_mm, units = "mm",
             dpi = d$dpi, device = "tiff", bg = "white",
             compression = "lzw")
    }
  )

  output$rm_dl_csv <- downloadHandler(
    filename    = function() paste0("qpcr_rm_stats_", Sys.Date(), ".csv"),
    contentType = "text/csv",
    content     = function(file) {
      req(rm_stats())
      write.csv(rm_stats(), file, row.names = FALSE)
    }
  )
}

shinyApp(ui, server)
