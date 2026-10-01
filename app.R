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

source("R/data.R")
source("R/analysis.R")
source("R/plot.R")

`%||%` <- function(a, b) if (is.null(a)) b else a

# Numeric input guard: NULL / NA / non-finite -> default.
num_or <- function(x, default, min = -Inf, max = Inf) {
  x <- suppressWarnings(as.numeric(x))
  if (length(x) != 1 || is.na(x) || !is.finite(x)) return(default)
  min(max(x, min), max)
}

# --- Example data (for the "Load example" buttons) -------------------------
example_data_1factor <- function() {
  data.frame(
    Sample       = paste0("mouse", 1:8),
    Group        = rep(c("Control", "Treated"), each = 4),
    Ct_reference = c(13.5, 13.6, 13.4, 13.5, 13.5, 13.4, 13.6, 13.5),
    HERPUD1      = c(25.1, 25.3, 25.2, 25.0, 22.1, 22.3, 22.0, 22.4),
    ATF4         = c(24.5, 24.7, 24.6, 24.4, 21.8, 22.0, 21.9, 22.1),
    stringsAsFactors = FALSE
  )
}

example_data_2factor <- function() {
  data.frame(
    Sample       = paste0("cell", 1:12),
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
app_css <- "
body, .form-control, .btn, .selectize-input, .dataTable { font-family: Arial, Helvetica, sans-serif; }
.paste-dropzone {
  border: 2px dashed #00A8A8;
  border-radius: 8px;
  padding: 12px;
  background: rgba(0, 168, 168, 0.05);
  text-align: center;
  margin-bottom: 10px;
}
.paste-dropzone .pz-title { font-weight: 600; color: #007373; margin-bottom: 4px; }
.paste-dropzone .pz-sub   { font-size: 12px; opacity: 0.8; }
.data-summary {
  padding: 8px 12px; border-radius: 6px; margin: 8px 0 12px; font-size: 13px;
  background: rgba(0, 168, 168, 0.08); border-left: 3px solid #00A8A8;
}
.data-summary.warn  { background: rgba(217, 119, 6, 0.10); border-left-color: #d97706; }
.data-summary.empty { background: rgba(120, 120, 120, 0.08); border-left-color: #aaa; opacity: 0.8; }
.data-badge {
  font-size: 12px; padding: 2px 8px; border-radius: 10px;
  background: rgba(0, 168, 168, 0.15); color: #007373; margin-right: 10px;
}
.data-badge.empty { background: rgba(120,120,120,0.12); color: #777; }
.status-box { padding: 10px 12px; border-radius: 6px; margin-bottom: 10px; font-size: 13px; }
.status-box.error { background: rgba(220, 53, 69, 0.10); border-left: 3px solid #dc3545; }
.status-box.warn  { background: rgba(217, 119, 6, 0.10); border-left: 3px solid #d97706; }
.status-box.info  { background: rgba(0, 168, 168, 0.08); border-left: 3px solid #00A8A8; }
.status-box ul { margin: 4px 0 0 0; padding-left: 18px; }
.methods-copy-wrap { position: relative; }
.methods-copy-btn { position: absolute; top: 6px; right: 6px; z-index: 5; font-size: 11px; padding: 2px 8px; }
.plot-card .card-header { font-size: 13px; padding: 6px 12px; background: transparent; border-bottom: 1px solid rgba(128,128,128,0.25); }
.color-row { display: flex; align-items: center; gap: 8px; margin-bottom: 4px; }
.color-row input[type=color] { width: 44px; height: 28px; padding: 0; border: 1px solid #888; border-radius: 4px; background: none; }
.color-row span { font-size: 12px; }
.small-note { font-size: 11px; opacity: 0.75; }
.shiny-plot-output img { display: block; margin: 0 auto; max-width: none; }
.plot-card .card-body { overflow-x: auto; }
"

# Navbar data-status badge.
data_status_badge_html <- function(wide, has_sub) {
  if (nrow(wide) == 0)
    return(span(class = "data-badge empty", icon("circle"), " no data"))
  n_genes  <- length(gene_cols(wide))
  n_groups <- length(unique(wide$Group[nzchar(wide$Group)]))
  txt <- paste(nrow(wide), "rows,", n_genes, if (n_genes == 1) "gene," else "genes,",
               n_groups, "groups")
  if (has_sub) txt <- paste(txt, "(2-factor)")
  span(class = "data-badge", icon("check-circle"), " ", txt)
}

sidebar_tip <- function(trigger, msg, placement = "top", ...) {
  bslib::tooltip(trigger, msg, placement = placement, ...)
}

# --- Sidebar ----------------------------------------------------------------
shared_sidebar <- function() {
  tagList(
    # Global paste listener: pastes anywhere on the page (outside text inputs)
    # are sent to the server.
    tags$script(HTML("
      document.addEventListener('paste', function(e) {
        var el = document.activeElement;
        var tag = el ? el.tagName.toLowerCase() : '';
        if (tag === 'input' || tag === 'textarea' || (el && el.isContentEditable)) return;
        var text = e.clipboardData.getData('text/plain');
        if (text && text.trim().length > 0) {
          Shiny.setInputValue('pasted_clipboard', text, {priority: 'event'});
          e.preventDefault();
        }
      });
      Shiny.addCustomMessageHandler('copy_to_clipboard', function(msg) {
        function done(ok) { Shiny.setInputValue('copy_result', {ok: ok, t: Date.now()}, {priority: 'event'}); }
        if (navigator.clipboard && navigator.clipboard.writeText) {
          navigator.clipboard.writeText(msg.text).then(function() { done(true); }, function() { done(false); });
        } else {
          try {
            var ta = document.createElement('textarea');
            ta.value = msg.text; document.body.appendChild(ta);
            ta.select(); var ok = document.execCommand('copy');
            document.body.removeChild(ta); done(!!ok);
          } catch (err) { done(false); }
        }
      });
    ")),

    sidebar_tip(
      checkboxInput("guided_mode",
                    tagList(icon("graduation-cap"), " Guided mode (open sections step by step)"),
                    value = TRUE),
      "On: only the current step is open. Off: all sections open at once."
    ),

    bslib::accordion(
      id = "sidebar_accordion",
      open = c("data"),
      multiple = TRUE,

      # --- Step 1: Data -------------------------------------------------------
      bslib::accordion_panel(
        value = "data",
        title = tagList(icon("database"), " 1. Data"),
        div(class = "paste-dropzone",
            div(class = "pz-title", icon("clipboard"), " Paste your Excel data"),
            div(class = "pz-sub",
                "Copy the cells in Excel, then press ", tags$kbd("Ctrl"), "+", tags$kbd("V"),
                " anywhere on this page (", tags$kbd("Cmd"), "+", tags$kbd("V"), " on Mac). ",
                "Header: ", code("Sample | Group | [Subgroup] | Ct_ref | Gene1 | ..."))
        ),
        textAreaInput("paste_box", NULL, rows = 3,
                      placeholder = "Or paste here and click Load (works on tablets and when Ctrl+V is blocked)"),
        actionButton("load_paste_box", "Load pasted text", icon = icon("upload"),
                     class = "btn-sm", style = "width:100%; margin-bottom:8px;"),
        fluidRow(
          column(6, actionButton("load_ex_1factor", "Load 1-factor example",
                                 icon = icon("flask"), class = "btn-sm", style = "width:100%;")),
          column(6, actionButton("load_ex_2factor", "Load 2-factor example",
                                 icon = icon("flask-vial"), class = "btn-sm", style = "width:100%;"))
        ),
        uiOutput("data_summary_ui"),
        DTOutput("table"),
        p(class = "small-note", style = "margin:4px 0;",
          "Double-click a cell to edit it. Click rows to select them for deletion."),
        fluidRow(
          column(4, actionButton("add_row", "Add row", icon = icon("plus"),
                                 class = "btn-sm", style = "width:100%;")),
          column(4, actionButton("delete_rows", "Delete selected", icon = icon("minus"),
                                 class = "btn-sm", style = "width:100%;")),
          column(4, actionButton("clear_tbl", "Clear table", icon = icon("trash"),
                                 class = "btn-sm", style = "width:100%;"))
        ),
        uiOutput("error_msg")
      ),

      # --- Step 2: Experiment -------------------------------------------------
      bslib::accordion_panel(
        value = "experiment",
        title = tagList(icon("vial"), " 2. Experiment"),
        sidebar_tip(
          uiOutput("ctrl_group_ui"),
          "The group used as the baseline (fold change = 1). Chosen from your Group column."
        ),
        sidebar_tip(
          uiOutput("ctrl_subgroup_ui"),
          "Two-factor designs: which Subgroup of the control group is the baseline (for example Young + Control). The pooled option mixes all subgroups of the control group and is rarely what you want."
        ),
        sidebar_tip(
          radioButtons("rep_mode", "Replicate type",
                       choices = c("Biological (one row per sample)" = "biological",
                                   "Technical (average rows that share a Sample name)" = "technical"),
                       selected = "biological"),
          "Biological: every row is an independent sample. Technical: rows with the same Sample name are averaged before delta Ct (MIQE)."
        ),
        uiOutput("layout_question_ui")
      ),

      # --- Step 3: Statistics -------------------------------------------------
      bslib::accordion_panel(
        value = "stats",
        title = tagList(icon("square-root-variable"), " 3. Statistics"),
        sidebar_tip(
          radioButtons("test_type", "Statistical test",
                       choices = c("Parametric (t-test / ANOVA + Tukey)"           = "parametric",
                                   "Nonparametric (Mann-Whitney / Kruskal-Wallis)" = "nonparametric",
                                   "vs control only: pooled t-tests (Holm)"        = "parametric_vs_control",
                                   "vs control only: Mann-Whitney (Holm)"          = "nonparametric_vs_control"),
                       selected = "parametric"),
          "All tests run on delta Ct values. The vs-control options compare every group with the control only (Dunnett-type), which is Prism's usual choice for dose-response designs."
        ),
        conditionalPanel(
          "input.test_type == 'parametric'",
          sidebar_tip(
            checkboxInput("var_equal", "Assume equal variances (Student's t) for 2 groups", value = FALSE),
            "Off = Welch's t-test (safer with small, unequal groups). On = Student's t-test, Prism's default."
          )
        ),
        sidebar_tip(
          checkboxInput("paired", "Paired samples (2 groups, matched by Sample name)", value = FALSE),
          "Only applies to exactly 2 groups without a Subgroup column. Samples are matched by their Sample name, so use the same name in both groups."
        ),
        sidebar_tip(
          radioButtons("sig_format", "P-value display",
                       choices = c("Exact (p = 0.0123)" = "exact",
                                   "Stars (*, **, ***, ****)" = "stars"),
                       selected = "exact"),
          "Stars follow Prism: * < 0.05, ** < 0.01, *** < 0.001, **** < 0.0001; ns otherwise."
        )
      ),

      # --- Step 4: Plot -------------------------------------------------------
      bslib::accordion_panel(
        value = "plot",
        title = tagList(icon("chart-column"), " 4. Plot"),
        sidebar_tip(
          sliderInput("plot_zoom", "Plot size on screen (%)",
                      min = 25, max = 300, value = 100, step = 5),
          "Zoom for viewing only. Downloads use the millimetre size set under Step 6."
        ),
        radioButtons("plot_type", "Plot type",
                     choices = c("Scatter + error bars" = "scatter",
                                 "Column + dots" = "column"),
                     inline = TRUE),
        radioButtons("y_scale", "Y-axis",
                     choices = c("log2 fold change" = "log2", "Fold change (linear)" = "linear"),
                     inline = TRUE),
        radioButtons("plot_layout", "Plot layout",
                     choices = c("Individual" = "individual", "Combined grid" = "combined"),
                     selected = "individual"),
        conditionalPanel(
          "input.plot_layout == 'combined'",
          numericInput("facet_ncol", "Columns in grid", value = 3, min = 1, max = 6, step = 1)
        ),
        sidebar_tip(
          radioButtons("error_type", "Error bars",
                       choices = c("SD" = "SD", "SEM" = "SEM", "95% CI" = "CI95"),
                       inline = TRUE),
          "95% CI uses the t-distribution (Prism default)."
        ),
        checkboxInput("show_sig", "Show significance brackets", value = TRUE),
        checkboxInput("rotate_x", "Rotate x-axis labels 45 degrees", value = FALSE),
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
                      choices = c("Arial"           = "Arial",
                                  "Helvetica"       = "Helvetica",
                                  "Times New Roman" = "Times New Roman",
                                  "Courier New"     = "Courier New",
                                  "Device default (sans)" = ""),
                      selected = "Arial"),
          "Named fonts need to be installed on the computer that renders the plot; when missing, the device default sans font is used."
        ),
        sliderInput("ts_title",      "Plot title size (pt)",      6, 28, 18, 1),
        sliderInput("ts_axis_title", "Axis title size (pt)",      6, 24, 14, 1),
        sliderInput("ts_axis_text",  "Axis tick text size (pt)",  6, 22, 12, 1),
        sliderInput("ts_legend",     "Legend text size (pt)",     6, 20, 12, 1),
        sliderInput("ts_facet",      "Facet (gene) title size (pt)", 6, 24, 14, 1),
        sliderInput("ts_sig_bar",    "Bracket text size (mm)",    2, 10, 3.5, 0.5),
        uiOutput("color_picker_ui"),
        actionButton("reset_colors", "Reset colors to Okabe-Ito",
                     icon = icon("rotate-left"), class = "btn-sm",
                     style = "width:100%; margin-top:6px;")
      ),

      # --- Step 6: Export -----------------------------------------------------
      bslib::accordion_panel(
        value = "export",
        title = tagList(icon("download"), " 6. Export size"),
        p(class = "small-note", style = "margin-bottom:4px;",
          "Exact dimensions for manuscript figures. Common widths: ",
          strong("85 mm"), " (single column), ", strong("174 mm"),
          " (double column), ", strong("89 / 183 mm"), " (Nature)."),
        fluidRow(
          column(6, numericInput("export_w_mm", "Width (mm)", value = 174, min = 30, max = 400, step = 1)),
          column(6, numericInput("export_h_mm", "Height (mm)", value = 120, min = 30, max = 400, step = 1))
        ),
        radioButtons("export_dpi", "DPI (PNG / TIFF)",
                     choices = c("300" = 300, "600" = 600, "1200" = 1200),
                     selected = 300, inline = TRUE),
        p(class = "small-note",
          "SVG and PDF are vector formats; DPI is ignored. Downloads always contain the combined grid of all genes.")
      )
    )
  )
}

# --- Help tab content ------------------------------------------------------
help_content <- function() {
  tagList(
    h2("qPCR Analysis: help"),
    h3("Which tab should I use?"),
    tags$ul(
      tags$li(strong("Analysis"), ": standard qPCR. Each sample belongs to one group (different animals or cultures per condition)."),
      tags$li(strong("Repeated Measures"), ": the same animal or culture measured under 2 or more conditions. Each Sample name appears in 2 or more Groups."),
      tags$li(strong("Help"), ": this page.")
    ),

    h3("Step 1: paste your data"),
    p("Copy the block from Excel (including the header row) and press Ctrl+V anywhere on the page, ",
      "or paste into the text box in the sidebar and click ", strong("Load pasted text"), ". ",
      "Tab, comma and semicolon separated text are all accepted. Decimal commas (25,3) are converted. ",
      "Cells that read Undetermined, N/A or are blank are treated as missing and reported in the data summary."),
    tags$pre(
      "Sample   Group    Ct_ref   Gene1    Gene2\n",
      "mouse1   Control  13.50    20.51    25.40\n",
      "mouse2   Control  13.56    20.62    25.23\n",
      "mouse3   Control  13.64    20.53    25.28\n",
      "mouse4   Treated  17.60    20.54    25.24\n",
      "mouse5   Treated  17.58    20.51    26.25\n",
      "mouse6   Treated  18.38    21.02    26.39"
    ),
    p("The third column is the reference (housekeeping) gene Ct. Its header can be Ct_ref, Reference gene, Reference or Housekeeping. ",
      "Every column after it is treated as a target gene."),
    p(strong("Two-factor design"), " (for example Young/Aged x Control/Treated): add a column named ",
      code("Subgroup"), " (also accepted: Factor2, Treatment, Condition) as the 3rd column, between Group and the reference Ct:"),
    tags$pre(
      "Sample   Group   Subgroup  Ct_ref   Gene1   Gene2\n",
      "cell1    Young   Control   13.50    20.51   25.40\n",
      "cell2    Young   Control   13.56    20.62   25.23\n",
      "cell3    Young   Control   13.64    20.53   25.28\n",
      "cell4    Young   Treated   17.60    17.20   22.10\n",
      "cell5    Young   Treated   17.58    17.15   22.30\n",
      "cell6    Young   Treated   18.38    17.40   22.40\n",
      "cell7    Aged    Control   13.55    20.35   25.15\n",
      "cell8    Aged    Control   13.60    20.40   25.20\n",
      "cell9    Aged    Control   13.65    20.30   25.10\n",
      "cell10   Aged    Treated   17.62    16.10   21.25\n",
      "cell11   Aged    Treated   17.55    16.25   21.40\n",
      "cell12   Aged    Treated   18.40    16.35   21.55"
    ),
    p(strong("Technical replicates"), ": paste one row per well and give the wells of one sample the same Sample name. ",
      "Then choose ", em("Technical"), " under Replicate type: the app averages Ct_target and Ct_ref of those rows before calculating delta Ct. ",
      "Leave the default ", em("Biological"), " when every row is an independent sample."),
    p(strong("Excluding a value"), ": select the row in the table and click ", strong("Delete selected"), ", ",
      "or double-click a Ct cell and clear it. A missing Ct removes that sample for that gene only; the other genes are unaffected."),

    h3("Step 2: baseline"),
    p(strong("Control group"), " is the reference for fold change (its mean fold change is 1). ",
      "In two-factor designs also pick the ", strong("Control subgroup"), " so that the baseline is one cell, for example Young + Control. ",
      "The pooled option averages all subgroups of the control group, which is only meaningful when the subgroups are expected to be equal at baseline."),

    h3("Step 3: statistical tests"),
    p("All tests are run on delta Ct values (log2 scale), which are closer to normally distributed than linear fold changes. ",
      "Effect sizes and 95% confidence intervals in the statistics table are given on the log2 fold-change scale."),
    tags$ul(
      tags$li(strong("delta-delta Ct"), " (Livak and Schmittgen, 2001): delta Ct = Ct(target) - Ct(reference); ",
              "delta-delta Ct = delta Ct - mean(delta Ct of the baseline); fold change = 2^(-delta-delta Ct)."),
      tags$li(strong("2 groups, parametric"), ": Welch's t-test (default) or Student's t-test (tick 'Assume equal variances'). ",
              "Tick 'Paired' for matched samples (paired t-test)."),
      tags$li(strong("2 groups, nonparametric"), ": Mann-Whitney U test, or Wilcoxon signed-rank when paired. ",
              "Exact p values are used for small samples without ties."),
      tags$li(strong("3+ groups, parametric"), ": one-way ANOVA followed by Tukey's HSD for all pairwise comparisons."),
      tags$li(strong("3+ groups, nonparametric"), ": Kruskal-Wallis test followed by pairwise Mann-Whitney U tests with Holm-Bonferroni correction. ",
              "This is not Dunn's test, which is what Prism runs; the two usually agree closely."),
      tags$li(strong("vs control only"), ": every other group is compared with the control only, using pooled-variance t-tests ",
              "(the ANOVA error term, as in Dunnett's test) or Mann-Whitney U tests, with Holm-Bonferroni correction across those comparisons. ",
              "This is a Dunnett-type procedure; Prism's Dunnett test uses a slightly different (usually less conservative) correction."),
      tags$li(strong("Two-factor design"), " (Subgroup column): two-way ANOVA with Type II sums of squares reporting the Group effect, the Subgroup effect ",
              "and their interaction, followed by Tukey's HSD between all Group x Subgroup cells. ",
              "With 'Nonparametric' the cells are compared with Kruskal-Wallis plus pairwise Mann-Whitney (Holm), because there is no nonparametric two-way ANOVA. ",
              "With 'vs control only' every cell is compared with the control cell (Control group + Control subgroup)."),
      tags$li(strong("Paired"), " only applies to exactly 2 groups without a Subgroup column. In other cases the app tells you it ran an unpaired analysis.")
    ),
    p("Whenever a setting cannot be applied, a note appears above the plot and the methods text describes what was actually run."),

    h3("Repeated Measures tab"),
    p("Same data format, but ", strong("each Sample name must appear in 2 or more Groups"), ":"),
    tags$pre(
      "Sample   Group      Ct_ref   Gene1    Gene2\n",
      "mouse1   Baseline   13.50    20.51    25.40\n",
      "mouse1   Treated    13.56    17.70    23.19\n",
      "mouse2   Baseline   13.64    20.53    25.28\n",
      "mouse2   Treated    14.22    17.67    23.12\n",
      "mouse3   Baseline   13.60    20.55    25.30\n",
      "mouse3   Treated    14.12    17.67    23.10"
    ),
    p("Each animal gets one colour and a line connects its values across conditions. Paired t-tests are run between every pair of conditions ",
      "and Holm-corrected within each gene. Because samples are matched by name, make sure numbering does not restart per group."),

    h3("Plot and export"),
    tags$ul(
      tags$li(strong("Y-axis"), ": log2 fold change (symmetric, recommended) or linear fold change."),
      tags$li(strong("Error bars"), ": SD, SEM or 95% CI (t-distribution)."),
      tags$li(strong("Significance brackets"), ": choose which comparisons to draw in the list under Step 4. ",
              "Significant comparisons are ticked by default; tick others to show them as ns or with their p value."),
      tags$li(strong("Two-factor layout"), ": one factor on the x-axis with nested labels, the other as coloured clusters (Prism grouped style). ",
              "Choose which factor goes on the x-axis under Step 2."),
      tags$li(strong("Colours"), ": Okabe-Ito colourblind-safe palette by default; change each level with the colour pickers under Step 5."),
      tags$li(strong("Font"), ": Arial by default. Fonts must be installed on the computer that draws the plot."),
      tags$li(strong("Downloads"), ": PNG and TIFF at the chosen DPI, PDF and SVG as vectors, in the width and height set under Step 6. ",
              "Downloads always contain the combined grid of all genes. ",
              strong("Results"), " exports every sample's delta Ct, delta-delta Ct and fold change; ",
              strong("Summary"), " exports the per-group means; ",
              strong("Stats"), " exports the statistics table; ",
              strong("Methods"), " exports the methods paragraph.")
    ),

    h3("References"),
    tags$ul(
      tags$li("Livak KJ, Schmittgen TD (2001). Analysis of relative gene expression data using real-time quantitative PCR and the 2^(-delta delta Ct) method. ", em("Methods"), " 25(4):402-408."),
      tags$li("Bustin SA et al. (2009). The MIQE guidelines: minimum information for publication of quantitative real-time PCR experiments. ", em("Clin Chem"), " 55(4):611-622.")
    ),
    hr(),
    p(class = "small-note",
      "Built with R Shiny and Shinylive. Source: ",
      a("github.com/karthi-ntu/qPCR_analysis", href = "https://github.com/karthi-ntu/qPCR_analysis", target = "_blank"))
  )
}

# ---- UI -------------------------------------------------------------------
app_theme <- bslib::bs_theme(
  version      = 5,
  primary      = "#00A8A8",
  "navbar-bg"  = "#00A8A8",
  base_font    = "Arial, Helvetica, sans-serif",
  heading_font = "Arial, Helvetica, sans-serif"
)

download_row <- function(prefix, extra = list()) {
  btns <- list(
    downloadButton(paste0(prefix, "dl_png"),  "PNG",  icon = icon("image"), class = "btn-sm"),
    downloadButton(paste0(prefix, "dl_pdf"),  "PDF",  icon = icon("file-pdf"), class = "btn-sm"),
    downloadButton(paste0(prefix, "dl_svg"),  "SVG",  icon = icon("bezier-curve"), class = "btn-sm"),
    downloadButton(paste0(prefix, "dl_tiff"), "TIFF", icon = icon("image"), class = "btn-sm")
  )
  div(style = "display:flex; flex-wrap:wrap; gap:6px;", btns, extra)
}

analysis_main_panel <- function() {
  tagList(
    uiOutput("main_status"),
    bslib::card(
      class = "plot-card",
      bslib::card_header(
        div(style = "display:flex; justify-content:space-between; align-items:center;",
            div(icon("chart-column"), strong(" Plot")),
            uiOutput("plot_metadata_chip", inline = TRUE))
      ),
      bslib::card_body(min_height = "420px", uiOutput("plots_ui"))
    ),
    bslib::card(
      bslib::card_header(icon("square-root-variable"), strong(" Statistics")),
      bslib::card_body(
        p(class = "small-note", "Differences and 95% CI are on the log2 fold-change scale (positive = higher in the first-named group)."),
        DTOutput("stats_table"))
    ),
    bslib::card(
      bslib::card_header(icon("file-lines"), strong(" Methods (copy into your manuscript)")),
      bslib::card_body(
        div(class = "methods-copy-wrap",
            actionButton("copy_methods", "Copy", icon = icon("copy"),
                         class = "btn-sm btn-outline-primary methods-copy-btn"),
            verbatimTextOutput("methods_text"))
      )
    ),
    bslib::card(
      bslib::card_header(icon("download"), strong(" Downloads")),
      bslib::card_body(
        download_row("", list(
          downloadButton("dl_results", "Results (per sample)", icon = icon("table"), class = "btn-sm"),
          downloadButton("dl_summary", "Summary (per group)",  icon = icon("table"), class = "btn-sm"),
          downloadButton("dl_csv",     "Stats",   icon = icon("table"), class = "btn-sm"),
          downloadButton("dl_methods", "Methods", icon = icon("file-lines"), class = "btn-sm")
        ))
      )
    )
  )
}

rm_main_panel <- function() {
  tagList(
    bslib::card(
      class = "plot-card",
      bslib::card_header(icon("link"), strong(" Paired plot (sample-matched)")),
      bslib::card_body(min_height = "420px", uiOutput("rm_status"), uiOutput("rm_plot_ui"))
    ),
    bslib::card(
      bslib::card_header(icon("square-root-variable"), strong(" Paired statistics")),
      bslib::card_body(DTOutput("rm_stats_table"))
    ),
    bslib::card(
      bslib::card_header(icon("download"), strong(" Downloads")),
      bslib::card_body(download_row("rm_", list(
        downloadButton("rm_dl_csv", "Stats", icon = icon("table"), class = "btn-sm"))))
    )
  )
}

rm_sidebar <- function() {
  tagList(
    p(style = "font-size:13px;",
      icon("info-circle"), " This tab uses the data pasted in the ", strong("Analysis"),
      " tab and treats ", strong("each Sample name as one animal or culture"),
      ". Each Sample must appear in 2 or more Groups."),
    p(style = "font-size:13px;",
      "Baseline, replicate type, error bars, y-axis, font, sizes and export settings come from the Analysis sidebar."),
    hr(),
    p(class = "small-note", "Dots are coloured by sample; lines connect the same sample across groups. Black bars show the group mean and error.")
  )
}

ui <- bslib::page_navbar(
  title = tagList(icon("dna"), " qPCR Analysis"),
  window_title = "qPCR Analysis (delta-delta Ct)",
  theme = app_theme,
  header = tags$head(tags$style(HTML(app_css))),

  bslib::nav_panel(
    tagList(icon("chart-column"), " Analysis"),
    bslib::layout_sidebar(
      sidebar = bslib::sidebar(width = 380, open = TRUE, shared_sidebar()),
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
    div(style = "padding: 20px; max-width: 960px; margin: auto;", help_content())
  ),

  bslib::nav_spacer(),
  bslib::nav_item(uiOutput("data_status_badge", inline = TRUE)),
  bslib::nav_item(bslib::input_dark_mode(id = "dark_mode")),
  bslib::nav_item(
    tags$a(icon("github"), " GitHub",
           href = "https://github.com/karthi-ntu/qPCR_analysis",
           target = "_blank", style = "color: white; margin-right: 10px;")
  )
)

# ---- Server ---------------------------------------------------------------
server <- function(input, output, session) {

  rv <- reactiveValues(wide = empty_wide(), prev_wide = NULL,
                       colors = list(), color_epoch = 0,
                       selected_comps = NULL)

  long_data    <- reactive({ to_long_df(rv$wide) })
  has_subgroup <- reactive({ "Subgroup" %in% names(rv$wide) })

  output$data_status_badge <- renderUI({ data_status_badge_html(rv$wide, has_subgroup()) })

  # -- Data summary (above the table) --------------------------------------
  output$data_summary_ui <- renderUI({
    wide <- rv$wide
    if (nrow(wide) == 0) {
      return(div(class = "data-summary empty", icon("circle-info"),
                 " No data yet. Paste from Excel or load an example."))
    }
    gcols  <- gene_cols(wide)
    groups <- unique(wide$Group[nzchar(wide$Group)])
    sub_txt <- if (has_subgroup()) {
      subs <- unique(wide$Subgroup[nzchar(wide$Subgroup)])
      paste0("; Subgroups: ", paste(subs, collapse = ", "))
    } else ""
    missing <- missing_ct_cells(long_data())
    dup_samples <- sum(duplicated(wide[, intersect(c("Sample", "Group", "Subgroup"), names(wide))]))
    cls <- if (length(missing) > 0) "data-summary warn" else "data-summary"
    div(class = cls,
        icon("check-circle"), " ",
        strong(nrow(wide), "rows"), "; ",
        strong(length(gcols)), if (length(gcols) == 1) " gene" else " genes",
        " (", paste(gcols, collapse = ", "), "); ",
        strong(length(groups)), " groups (", paste(groups, collapse = ", "), ")", sub_txt,
        if (dup_samples > 0) tags$div(class = "small-note",
          icon("info-circle"), " ", dup_samples, " row(s) share a Sample name with another row. ",
          "Choose 'Technical' under Replicate type if these are technical replicates."),
        if (length(missing) > 0) tags$div(class = "small-note", style = "color:#b45309;",
          icon("triangle-exclamation"), " Missing Ct (excluded for that gene): ",
          paste(missing, collapse = "; ")))
  })

  # -- Guided vs Expert mode ------------------------------------------------
  observeEvent(input$guided_mode, {
    all_steps <- c("data", "experiment", "stats", "plot", "styling", "export")
    if (isTRUE(input$guided_mode)) {
      for (s in setdiff(all_steps, "data"))
        bslib::accordion_panel_close("sidebar_accordion", values = s)
      bslib::accordion_panel_open("sidebar_accordion", values = "data")
    } else {
      for (s in all_steps) bslib::accordion_panel_open("sidebar_accordion", values = s)
    }
  }, ignoreInit = TRUE)

  # In guided mode, open the Experiment step once data is present.
  observeEvent(rv$wide, {
    if (isTRUE(input$guided_mode) && nrow(rv$wide) > 0)
      bslib::accordion_panel_open("sidebar_accordion", values = "experiment")
  }, ignoreInit = TRUE)

  # -- Loading data -----------------------------------------------------------
  set_data <- function(wide, source_label) {
    rv$prev_wide <- rv$wide
    rv$wide <- wide
    rv$colors <- list()
    showNotification(
      tagList(paste0(source_label, ": ", nrow(wide), " rows loaded. "),
              actionLink("undo_load", "Undo")),
      type = "message", duration = 6, id = "load_note")
  }

  observeEvent(input$undo_load, {
    if (!is.null(rv$prev_wide)) {
      rv$wide <- rv$prev_wide
      rv$prev_wide <- NULL
      removeNotification("load_note")
      showNotification("Previous data restored.", type = "message", duration = 3)
    }
  })

  observeEvent(input$load_ex_1factor, { set_data(example_data_1factor(), "1-factor example") })
  observeEvent(input$load_ex_2factor, { set_data(example_data_2factor(), "2-factor example") })

  handle_paste <- function(txt) {
    res <- parse_pasted_text(txt)
    if (!isTRUE(res$ok)) {
      showNotification(tagList(strong("Could not read the pasted data. "), res$error),
                       type = "error", duration = 10)
      return(invisible(FALSE))
    }
    set_data(res$wide, "Pasted data")
    for (w in res$warnings)
      showNotification(w, type = "warning", duration = 8)
    invisible(TRUE)
  }

  observeEvent(input$pasted_clipboard, { handle_paste(input$pasted_clipboard) })
  observeEvent(input$load_paste_box, {
    if (handle_paste(input$paste_box)) updateTextAreaInput(session, "paste_box", value = "")
  })

  # -- Copy methods -------------------------------------------------------------
  observeEvent(input$copy_methods, {
    txt <- methods_txt()
    req(nzchar(txt))
    session$sendCustomMessage("copy_to_clipboard", list(text = txt))
  })
  observeEvent(input$copy_result, {
    if (isTRUE(input$copy_result$ok))
      showNotification("Methods copied to clipboard.", type = "message", duration = 2)
    else
      showNotification("Copy blocked by the browser. Select the text and copy it manually.",
                       type = "warning", duration = 5)
  })

  # -- Data table -------------------------------------------------------------
  output$table <- renderDT({
    datatable(
      rv$wide,
      editable  = TRUE,
      rownames  = FALSE,
      selection = "multiple",
      options   = list(dom = "t", pageLength = 100, scrollY = "300px", scrollX = TRUE)
    )
  })

  observeEvent(input$table_cell_edit, {
    rv$wide <- editData(rv$wide, input$table_cell_edit, rownames = FALSE)
    for (col in setdiff(names(rv$wide), c("Sample", "Group", "Subgroup")))
      rv$wide[[col]] <- parse_ct(rv$wide[[col]])
  })

  observeEvent(input$add_row, {
    wide <- rv$wide
    if (length(gene_cols(wide)) == 0) wide$Gene1 <- numeric(nrow(wide))
    new_row <- as.data.frame(lapply(wide, function(col) if (is.numeric(col)) NA_real_ else ""),
                             stringsAsFactors = FALSE)
    if (nrow(wide) == 0) new_row <- new_row[1, , drop = FALSE]
    rv$wide <- rbind(wide, new_row)
  })

  observeEvent(input$delete_rows, {
    sel <- input$table_rows_selected
    if (length(sel) == 0) {
      showNotification("Click one or more rows in the table first.", type = "warning", duration = 4)
      return()
    }
    rv$prev_wide <- rv$wide
    rv$wide <- rv$wide[-sel, , drop = FALSE]
    showNotification(tagList(paste0(length(sel), " row(s) deleted. "), actionLink("undo_load", "Undo")),
                     type = "message", duration = 6, id = "load_note")
  })

  observeEvent(input$clear_tbl, {
    rv$prev_wide <- rv$wide
    rv$wide <- empty_wide()
  })

  # -- Control group / subgroup selectors ---------------------------------------
  output$ctrl_group_ui <- renderUI({
    grps <- unique(rv$wide$Group)
    grps <- grps[nzchar(grps)]
    if (length(grps) == 0) {
      return(selectInput("ctrl_group", "Control group",
                         choices = c("(paste data first)" = ""), selected = ""))
    }
    current <- isolate(input$ctrl_group)
    # Prefer a group whose name looks like a control when nothing was chosen.
    guess <- grps[grepl("^(control|ctrl|ctl|vehicle|veh|carrier|wt|sham|baseline|untreated|mock)", grps, ignore.case = TRUE)]
    sel <- if (!is.null(current) && nzchar(current) && current %in% grps) current
           else if (length(guess)) guess[1] else grps[1]
    selectInput("ctrl_group", "Control group", choices = grps, selected = sel)
  })

  output$ctrl_subgroup_ui <- renderUI({
    if (!has_subgroup()) return(NULL)
    subs <- unique(rv$wide$Subgroup)
    subs <- subs[nzchar(subs)]
    if (length(subs) == 0) return(NULL)
    current <- isolate(input$ctrl_subgroup)
    guess <- subs[grepl("^(control|ctrl|ctl|vehicle|veh|carrier|dmso|untreated|mock|baseline)", subs, ignore.case = TRUE)]
    sel <- if (!is.null(current) && current %in% c("", subs)) current
           else if (length(guess)) guess[1] else subs[1]
    selectInput("ctrl_subgroup", "Control subgroup (baseline cell)",
                choices = c(subs, "(pool all subgroups of the control group)" = ""),
                selected = sel)
  })

  output$layout_question_ui <- renderUI({
    if (!has_subgroup()) return(NULL)
    tagList(
      div(class = "status-box info",
          tags$b("Two-factor design detected. "),
          "One factor goes on the x-axis; the other is shown as coloured clusters with nested labels."),
      radioButtons("x_axis_var", "Factor on the x-axis",
                   choices = c("Group (e.g. Young / Aged)" = "Group",
                               "Subgroup (e.g. Control / Treated)" = "Subgroup"),
                   selected = isolate(input$x_axis_var) %||% "Group")
    )
  })

  # -- Processing -------------------------------------------------------------
  # Rows with a missing Ct are dropped for that gene only.
  clean_long <- reactive({
    df <- long_data()
    df[!is.na(df$Ct_target) & !is.na(df$Ct_reference), , drop = FALSE]
  })

  processed_long <- reactive({
    df <- clean_long()
    if (identical(input$rep_mode, "technical")) df <- average_tech_reps(df)
    df
  })

  cell_key <- function(df) {
    if ("Subgroup" %in% names(df)) paste(df$Group, df$Subgroup, sep = " | ") else df$Group
  }

  # Returns NULL when analysis can proceed, otherwise list(level, msg).
  validation_msg <- reactive({
    raw <- long_data()
    if (nrow(raw) == 0) return(list(level = "info", msg = "Paste data or load an example to start."))
    df <- processed_long()
    if (nrow(df) == 0) return(list(level = "error", msg = "No usable rows: every row has a missing Ct value."))
    ctrl <- input$ctrl_group %||% ""
    if (!nzchar(ctrl)) return(list(level = "info", msg = "Waiting for the control group selection."))
    if (!ctrl %in% df$Group)
      return(list(level = "error", msg = paste0("Control group '", ctrl, "' is not in the Group column.")))
    cs <- input$ctrl_subgroup
    if (has_subgroup() && !is.null(cs) && nzchar(cs) &&
        !any(df$Group == ctrl & df$Subgroup == cs))
      return(list(level = "error",
                  msg = paste0("There are no rows with Group '", ctrl, "' and Subgroup '", cs,
                               "'. Pick a different control subgroup.")))
    n_cells <- length(unique(cell_key(df)))
    if (n_cells < 2)
      return(list(level = "error", msg = "At least 2 groups are needed to compare anything."))
    for (g in unique(df$Gene)) {
      sub <- df[df$Gene == g, ]
      idx <- if (has_subgroup() && !is.null(cs) && nzchar(cs)) sub$Group == ctrl & sub$Subgroup == cs
             else sub$Group == ctrl
      if (!any(idx))
        return(list(level = "error", msg = paste0("Gene ", g, " has no usable control values.")))
    }
    NULL
  })

  # Soft warnings that do not block the analysis.
  data_warnings <- reactive({
    df <- processed_long()
    if (nrow(df) == 0) return(character(0))
    notes <- character(0)
    counts <- table(df$Gene, cell_key(df))
    small <- which(counts > 0 & counts < 2, arr.ind = TRUE)
    if (nrow(small) > 0) {
      cells <- unique(colnames(counts)[small[, 2]])
      notes <- c(notes, paste0("Only 1 sample in: ", paste(cells, collapse = ", "),
                               ". Tests involving these groups give no p value."))
    }
    missing <- missing_ct_cells(long_data(), max_items = 4)
    if (length(missing) > 0)
      notes <- c(notes, paste0("Missing Ct values were excluded: ", paste(missing, collapse = "; "), "."))
    notes
  })

  analyzed <- reactive({
    req(is.null(validation_msg()))
    df <- compute_delta_ct(processed_long())
    out <- tryCatch(
      compute_fold_change(df, control_group = input$ctrl_group,
                          control_subgroup = if (has_subgroup()) input$ctrl_subgroup else NULL),
      error = function(e) { validate(need(FALSE, conditionMessage(e))); NULL })
    attr(out, "tech_counts") <- attr(processed_long(), "tech_counts")
    out
  })

  stats_result <- reactive({
    req(analyzed())
    run_stats(analyzed(), paired = isTRUE(input$paired),
              test_type        = input$test_type %||% "parametric",
              has_subgroup     = has_subgroup(),
              control_group    = input$ctrl_group,
              control_subgroup = if (has_subgroup()) input$ctrl_subgroup else NULL,
              var_equal        = isTRUE(input$var_equal))
  })

  genes_in_data <- reactive({ req(analyzed()); unique(analyzed()$Gene) })

  baseline_label <- reactive({
    cs <- input$ctrl_subgroup
    if (has_subgroup() && !is.null(cs) && nzchar(cs)) paste(input$ctrl_group, cs, sep = " | ")
    else input$ctrl_group
  })

  # -- Status and notes shown in the main panel ----------------------------
  output$error_msg <- renderUI({
    v <- validation_msg()
    if (!is.null(v) && v$level == "error")
      div(class = "status-box error", icon("triangle-exclamation"), " ", v$msg)
  })

  output$main_status <- renderUI({
    v <- validation_msg()
    if (!is.null(v)) {
      cls <- if (v$level == "error") "status-box error" else "status-box info"
      ic  <- if (v$level == "error") icon("triangle-exclamation") else icon("circle-info")
      return(div(class = cls, ic, " ", v$msg))
    }
    notes <- c(attr(stats_result(), "notes"), data_warnings())
    if (length(notes) == 0) return(NULL)
    div(class = "status-box warn", icon("circle-info"), strong(" Please note"),
        tags$ul(lapply(notes, tags$li)))
  })

  output$plot_metadata_chip <- renderUI({
    if (is.null(validation_msg()) && !is.null(stats_result())) {
      tests <- setdiff(unique(stats_result()$Test), NA)
      tests <- tests[!grepl("ANOVA|Kruskal", tests)]
      test_name <- if (length(tests)) tests[1] else unique(stats_result()$Test)[1] %||% ""
      n_txt <- tryCatch(replicate_summary(analyzed(), rep_mode = input$rep_mode %||% "biological"),
                        error = function(e) "")
      tags$span(class = "small-note", paste(c(n_txt, test_name), collapse = "; "))
    }
  })

  # -- Bracket selection --------------------------------------------------------
  comparison_choices <- reactive({
    st <- stats_result()
    st <- st[grepl(" vs ", st$Comparison, fixed = TRUE), , drop = FALSE]
    list(all = unique(paste0(st$Gene, ": ", st$Comparison)),
         sig = unique(paste0(st$Gene, ": ", st$Comparison)[st$Significant == "Yes"]))
  })

  # The selection is mirrored server-side because a checkbox group with
  # nothing ticked reports NULL, which is indistinguishable from "not rendered
  # yet". Every time the statistics change the selection resets to the
  # significant comparisons.
  output$hide_comps_ui <- renderUI({
    req(stats_result())
    ch <- comparison_choices()
    rv$selected_comps <- ch$sig
    if (length(ch$all) == 0) return(NULL)
    tagList(
      tags$label("Significance brackets to draw", style = "font-weight:bold;"),
      div(style = "display:flex; gap:4px;",
          actionButton("sig_select_sig", "Significant", class = "btn-sm", style = "flex:1; font-size:11px;"),
          actionButton("sig_select_all", "All",         class = "btn-sm", style = "flex:1; font-size:11px;"),
          actionButton("sig_deselect_all", "None",      class = "btn-sm", style = "flex:1; font-size:11px;")),
      br(),
      checkboxGroupInput("visible_comps", label = NULL, choices = ch$all, selected = ch$sig)
    )
  })
  observeEvent(input$visible_comps, {
    rv$selected_comps <- input$visible_comps %||% character(0)
  }, ignoreNULL = FALSE, ignoreInit = TRUE)
  observeEvent(input$sig_select_all, {
    updateCheckboxGroupInput(session, "visible_comps", selected = comparison_choices()$all) })
  observeEvent(input$sig_select_sig, {
    updateCheckboxGroupInput(session, "visible_comps", selected = comparison_choices()$sig) })
  observeEvent(input$sig_deselect_all, {
    updateCheckboxGroupInput(session, "visible_comps", selected = character(0)) })

  visible_comps <- reactive({ rv$selected_comps })

  # -- Colours ------------------------------------------------------------------
  fill_levels <- reactive({
    df <- rv$wide
    if (nrow(df) == 0) return(character(0))
    if (has_subgroup() && !identical(input$x_axis_var, "Subgroup")) {
      lv <- unique(df$Subgroup)
    } else {
      lv <- unique(df$Group)
    }
    lv[nzchar(lv)]
  })

  output$color_picker_ui <- renderUI({
    lvs <- fill_levels()
    rv$color_epoch
    if (length(lvs) == 0) return(p(class = "small-note", "(add data to set colours)"))
    palette <- rep_len(OKABE_ITO, length(lvs))
    tagList(
      tags$label("Colours", style = "font-weight:bold;"),
      lapply(seq_along(lvs), function(i) {
        level <- lvs[i]
        current <- isolate(rv$colors[[level]]) %||% palette[i]
        div(class = "color-row",
            tags$input(type = "color", value = current,
                       oninput = sprintf("Shiny.setInputValue('color_pick', {level: %s, value: this.value}, {priority: 'event'})",
                                         jsonlite::toJSON(level, auto_unbox = TRUE))),
            span(level))
      })
    )
  })

  observeEvent(input$color_pick, {
    rv$colors[[input$color_pick$level]] <- input$color_pick$value
  })

  observeEvent(input$reset_colors, {
    rv$colors <- list()
    rv$color_epoch <- rv$color_epoch + 1
  })

  color_override <- reactive({
    lvs <- fill_levels()
    vals <- unlist(rv$colors)
    if (is.null(vals)) return(NULL)
    vals[names(vals) %in% lvs]
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

  # Arguments shared by every plot call.
  plot_args <- reactive({
    list(
      error_type    = input$error_type %||% "SD",
      plot_type     = input$plot_type %||% "scatter",
      show_sig      = isTRUE(input$show_sig),
      control_group = baseline_label(),
      rotate_x      = isTRUE(input$rotate_x),
      aspect_ratio  = num_or(input$aspect_ratio, 1, 0.2, 5),
      has_subgroup  = has_subgroup(),
      fill_override = color_override(),
      text_sizes    = text_sizes(),
      font_family   = input$font_family %||% DEFAULT_FONT,
      sig_format    = input$sig_format %||% "exact",
      x_axis_var    = input$x_axis_var %||% "Group",
      y_scale       = input$y_scale %||% "log2",
      y_min         = if (is.finite(num_or(input$y_min, NA))) input$y_min else NULL,
      y_max         = if (is.finite(num_or(input$y_max, NA))) input$y_max else NULL
    )
  })

  facet_ncol <- reactive({ as.integer(num_or(input$facet_ncol, 3, 1, 6)) })

  combined_plot <- reactive({
    do.call(make_combined_plot,
            c(list(analyzed(), stats_result(), ncol = facet_ncol(),
                   sig_comparisons_all = visible_comps()), plot_args()))
  })

  # -- Plots --------------------------------------------------------------------
  # Plot heights follow the container width and the aspect ratio, so a
  # square panel is never clipped when the browser window is wide.
  # (plotOutput(height = "auto") + renderPlot(height = fn) lets the server
  # decide the height from session$clientData.)
  # The on-screen size comes from the "Plot size on screen" slider: each panel
  # is drawn about 380 px wide at 100%, capped by the available width.
  plot_dims_px <- function(output_id, n_panels, ncol, per_panel_extra = 130) {
    # The reported width is 0 or missing while the output is hidden or has
    # not been laid out yet (e.g. right after switching layout); fall back to
    # a sensible width so the plot never renders at 0 px.
    avail  <- session$clientData[[paste0("output_", output_id, "_width")]]
    if (is.null(avail) || !is.finite(avail) || avail < 200) avail <- 900
    zoom   <- num_or(input$plot_zoom, 100, 25, 300) / 100
    ncol_e <- max(1, min(ncol, n_panels))
    nrow   <- ceiling(n_panels / ncol_e)
    # Above 100% the plot may grow past the card (it scrolls sideways).
    width  <- round(min(avail * max(1, zoom), ncol_e * 380 * zoom + 110))
    panel_w <- max(120, (width - 110) / ncol_e - 20)
    aspect  <- num_or(input$aspect_ratio, 1, 0.2, 5)
    list(w = width,
         h = max(240, round(nrow * (panel_w * aspect + per_panel_extra))))
  }

  output$plots_ui <- renderUI({
    req(genes_in_data())
    if (identical(input$plot_layout, "combined")) {
      plotOutput("plot_combined", height = "auto")
    } else {
      do.call(tagList, lapply(genes_in_data(), function(gene)
        plotOutput(paste0("plot_", make.names(gene)), height = "auto")))
    }
  })

  output$plot_combined <- renderPlot(
    { req(genes_in_data()); combined_plot() },
    width  = function() plot_dims_px("plot_combined", length(genes_in_data()), facet_ncol())$w,
    height = function() plot_dims_px("plot_combined", length(genes_in_data()), facet_ncol())$h
  )

  registered_gene_plots <- character(0)
  observe({
    req(genes_in_data())
    for (gene in setdiff(genes_in_data(), registered_gene_plots)) {
      local({
        g  <- gene
        id <- paste0("plot_", make.names(g))
        output[[id]] <- renderPlot({
          req(g %in% genes_in_data())
          do.call(make_barplot,
                  c(list(analyzed(), stats_result(), gene = g,
                         sig_comparisons = visible_comps()), plot_args()))
        }, width  = function() plot_dims_px(id, 1, 1)$w,
           height = function() plot_dims_px(id, 1, 1)$h)
      })
      registered_gene_plots <<- c(registered_gene_plots, gene)
    }
  })

  # -- Stats table, methods ---------------------------------------------------
  output$stats_table <- renderDT({
    req(stats_result())
    df  <- stats_result()
    fmt <- pick_sig_formatter(input$sig_format %||% "exact")
    df$p_display <- fmt(df$p_value)
    df <- df[, c("Gene", "Comparison", "Test", "Statistic", "df",
                 "Diff_log2FC", "CI95_low", "CI95_high", "p_display", "Significant")]
    names(df) <- c("Gene", "Comparison", "Test", "Statistic", "df",
                   "Diff (log2FC)", "CI low", "CI high", "p", "Significant")
    datatable(df, rownames = FALSE, options = list(dom = "t", pageLength = 200, scrollX = TRUE)) |>
      formatStyle("Significant", target = "row",
                  backgroundColor = styleEqual("Yes", "rgba(0, 168, 168, 0.18)"))
  })

  methods_txt <- reactive({
    req(stats_result(), analyzed())
    generate_methods_text(
      df = analyzed(), stats_df = stats_result(),
      paired = isTRUE(input$paired), control_group = baseline_label(),
      has_subgroup = has_subgroup(), rep_mode = input$rep_mode %||% "biological",
      sig_format = input$sig_format %||% "exact", var_equal = isTRUE(input$var_equal)
    )
  })
  output$methods_text <- renderText({ methods_txt() })

  # -- Downloads ------------------------------------------------------------------
  download_dims <- function() {
    list(w_mm = num_or(input$export_w_mm, 174, 30, 400),
         h_mm = num_or(input$export_h_mm, 120, 30, 400),
         dpi  = as.integer(num_or(input$export_dpi, 300, 72, 1200)))
  }

  plot_download <- function(id, fmt, plot_fun, prefix) {
    force(id); force(fmt); force(plot_fun); force(prefix)
    output[[id]] <- downloadHandler(
      filename    = function() paste0(prefix, Sys.Date(), ".", fmt),
      content     = function(file) {
        d <- download_dims()
        save_plot_file(file, plot_fun(), fmt, d$w_mm, d$h_mm, dpi = d$dpi,
                       font_family = input$font_family %||% DEFAULT_FONT)
      }
    )
  }
  for (fmt in c("png", "pdf", "svg", "tiff")) {
    plot_download(paste0("dl_", fmt),    fmt, function() { req(genes_in_data()); combined_plot() }, "qpcr_plot_")
    plot_download(paste0("rm_dl_", fmt), fmt, function() { req(rm_valid()$ok); rm_plot() },        "qpcr_rm_plot_")
  }

  results_table <- reactive({
    df <- analyzed()
    cols <- intersect(c("Sample", "Group", "Subgroup", "Gene", "Ct_target", "Ct_reference",
                        "delta_ct", "delta_delta_ct", "fold_change", "log2_fold_change"), names(df))
    df[, cols]
  })

  output$dl_results <- downloadHandler(
    filename = function() paste0("qpcr_results_per_sample_", Sys.Date(), ".csv"),
    content  = function(file) { req(analyzed()); write.csv(results_table(), file, row.names = FALSE) }
  )
  output$dl_summary <- downloadHandler(
    filename = function() paste0("qpcr_summary_per_group_", Sys.Date(), ".csv"),
    content  = function(file) { req(analyzed()); write.csv(summarize_groups(analyzed()), file, row.names = FALSE) }
  )
  output$dl_csv <- downloadHandler(
    filename = function() paste0("qpcr_stats_", Sys.Date(), ".csv"),
    content  = function(file) { req(stats_result()); write.csv(stats_result(), file, row.names = FALSE) }
  )
  output$dl_methods <- downloadHandler(
    filename = function() paste0("qpcr_methods_", Sys.Date(), ".txt"),
    content  = function(file) { req(methods_txt()); writeLines(methods_txt(), file) }
  )

  # -- Repeated Measures tab -------------------------------------------------------
  rm_valid <- reactive({
    df <- processed_long()
    if (nrow(long_data()) == 0)
      return(list(ok = FALSE, msg = "No data. Paste your data in the Analysis tab first."))
    if (nrow(df) == 0)
      return(list(ok = FALSE, msg = "No usable rows: every row has a missing Ct value."))
    if (!is_repeated_measures(df))
      return(list(ok = FALSE, msg = "This does not look like repeated measures: each Sample name appears in only one Group. Use the Analysis tab instead."))
    ctrl <- input$ctrl_group %||% ""
    if (!nzchar(ctrl)) return(list(ok = FALSE, msg = "Waiting for the control group selection."))
    if (!ctrl %in% df$Group)
      return(list(ok = FALSE, msg = paste0("Control group '", ctrl, "' not found.")))
    list(ok = TRUE, msg = NULL)
  })

  output$rm_status <- renderUI({
    v <- rm_valid()
    if (!v$ok) div(class = "status-box info", icon("info-circle"), " ", v$msg)
  })

  rm_analyzed <- reactive({
    req(rm_valid()$ok)
    df <- compute_delta_ct(processed_long())
    tryCatch(compute_fold_change(df, control_group = input$ctrl_group,
                                 control_subgroup = if (has_subgroup()) input$ctrl_subgroup else NULL),
             error = function(e) { validate(need(FALSE, conditionMessage(e))); NULL })
  })

  rm_stats <- reactive({ req(rm_analyzed()); run_paired_stats(rm_analyzed()) })

  output$rm_plot_ui <- renderUI({
    req(rm_valid()$ok, rm_analyzed())
    plotOutput("rm_plot", height = "auto")
  })

  rm_plot <- reactive({
    a <- plot_args()
    a$plot_type <- NULL; a$has_subgroup <- NULL; a$fill_override <- NULL; a$x_axis_var <- NULL
    do.call(make_paired_plot,
            c(list(rm_analyzed(), rm_stats(), ncol = facet_ncol(), sig_comparisons_all = NULL), a))
  })

  output$rm_plot <- renderPlot(
    { rm_plot() },
    width  = function() plot_dims_px("rm_plot", length(unique(rm_analyzed()$Gene)), facet_ncol(), 150)$w,
    height = function() plot_dims_px("rm_plot", length(unique(rm_analyzed()$Gene)), facet_ncol(), 150)$h
  )

  output$rm_stats_table <- renderDT({
    req(rm_stats())
    df  <- rm_stats()
    fmt <- pick_sig_formatter(input$sig_format %||% "exact")
    df$p_display <- fmt(df$p_value)
    df <- df[, c("Gene", "Comparison", "Test", "Statistic", "df",
                 "Diff_log2FC", "CI95_low", "CI95_high", "p_display", "Significant")]
    names(df) <- c("Gene", "Comparison", "Test", "Statistic", "df",
                   "Diff (log2FC)", "CI low", "CI high", "p", "Significant")
    datatable(df, rownames = FALSE, options = list(dom = "t", pageLength = 200, scrollX = TRUE)) |>
      formatStyle("Significant", target = "row",
                  backgroundColor = styleEqual("Yes", "rgba(0, 168, 168, 0.18)"))
  })

  output$rm_dl_csv <- downloadHandler(
    filename = function() paste0("qpcr_rm_stats_", Sys.Date(), ".csv"),
    content  = function(file) { req(rm_stats()); write.csv(rm_stats(), file, row.names = FALSE) }
  )
}

shinyApp(ui, server)
