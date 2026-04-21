library(shiny)
library(ggplot2)
library(DT)
library(ggsignif)

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

# --- Sidebar (shared controls) ---------------------------------------------
shared_sidebar <- function() {
  tagList(
    textInput("ctrl_group", "Control group name", value = "Control"),
    radioButtons("test_type", "Statistical test",
                 choices = c("Parametric (t-test / ANOVA)" = "parametric",
                             "Non-parametric (Wilcoxon / Kruskal)" = "nonparametric"),
                 selected = "parametric"),
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
    radioButtons("error_type", "Error bars",
                 choices = c("SD", "SEM"), inline = TRUE),
    checkboxInput("paired", "Paired samples (2 groups only)", value = FALSE),
    checkboxInput("show_sig", "Show significance bars", value = TRUE),
    checkboxInput("rotate_x", "Rotate x-axis labels 45\u00b0", value = FALSE),
    sliderInput("aspect_ratio", "Plot aspect ratio (height/width)",
                min = 0.5, max = 2.5, value = 1.0, step = 0.1),
    fluidRow(
      column(6, numericInput("y_min", "Y-axis min", value = NA, step = 0.5)),
      column(6, numericInput("y_max", "Y-axis max", value = NA, step = 0.5))
    ),
    uiOutput("hide_comps_ui"),
    hr(),
    h5("Data"),
    p(style = "font-size:12px; color:#888;",
      "Paste from Excel. Header: ",
      strong("Sample | Group | [Subgroup] | Reference gene | Gene1 | Gene2 | ..."),
      br(), em("Subgroup column is optional (triggers two-way ANOVA).")),
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
    DTOutput("table"),
    br(),
    fluidRow(
      column(6, actionButton("add_row",   "Add row",     icon = icon("plus"))),
      column(6, actionButton("clear_tbl", "Clear table", icon = icon("trash")))
    ),
    hr(),
    uiOutput("error_msg")
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
    p(em("Optional:"), " add a ", code("Subgroup"), " column (as the 3rd column, before ",
      code("Reference gene"), ") to run two-way ANOVA. Example: ",
      code("Sample | Group | Subgroup | Reference gene | Gene1 | Gene2")),
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
      tags$li(strong("Two-way ANOVA"), " \u2014 when Subgroup column present; tests Group, Subgroup, and interaction.")
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
ui <- navbarPage(
  title = "qPCR Analysis \u2014 \u0394\u0394Ct Method",
  tabPanel(
    "Analysis",
    sidebarLayout(
      sidebarPanel(width = 4, shared_sidebar()),
      mainPanel(
        width = 8,
        uiOutput("plots_ui"),
        br(),
        h4("Statistics"),
        DTOutput("stats_table"),
        br(),
        h4("Methods (copy-paste for papers)"),
        verbatimTextOutput("methods_text"),
        br(),
        fluidRow(
          column(3, downloadButton("dl_svg",     "Plot (SVG)")),
          column(2, downloadButton("dl_png",     "PNG")),
          column(2, downloadButton("dl_pdf",     "PDF")),
          column(2, downloadButton("dl_csv",     "Stats CSV")),
          column(3, downloadButton("dl_methods", "Methods TXT"))
        ),
        p(style = "font-size:12px; color:#888; margin-top:8px;",
          strong("SVG"), " is recommended for the online (GitHub Pages) version \u2014 ",
          "opens in any browser, Illustrator, Inkscape; infinite resolution. ",
          em("PNG/PDF only work reliably in the local version."))
      )
    )
  ),
  tabPanel(
    "Repeated Measures / Animal",
    sidebarLayout(
      sidebarPanel(
        width = 4,
        p(style = "font-size:13px;",
          "This tab uses the same data pasted in the ", strong("Analysis"), " tab, ",
          "but treats ", strong("each Sample name as one animal"), ". ",
          "Each Sample must appear in 2+ Groups."),
        hr(),
        h5("RM plot style"),
        p(style = "font-size:12px; color:#888;",
          "Dots are colored by animal; thin lines connect the same animal across groups.")
      ),
      mainPanel(
        width = 8,
        uiOutput("rm_status"),
        uiOutput("rm_plot_ui"),
        br(),
        h4("Paired statistics"),
        DTOutput("rm_stats_table"),
        br(),
        fluidRow(
          column(3, downloadButton("rm_dl_svg", "Plot (SVG)")),
          column(3, downloadButton("rm_dl_png", "PNG")),
          column(3, downloadButton("rm_dl_pdf", "PDF")),
          column(3, downloadButton("rm_dl_csv", "Stats CSV"))
        ),
        p(style = "font-size:12px; color:#888; margin-top:8px;",
          em("SVG works on GitHub Pages; PNG/PDF only work locally."))
      )
    )
  ),
  tabPanel("Help", fluidPage(help_content()))
)

# ---- Server ---------------------------------------------------------------
server <- function(input, output, session) {

  rv <- reactiveValues(wide = empty_wide())

  long_data    <- reactive({ to_long_df(rv$wide) })
  has_subgroup <- reactive({ "Subgroup" %in% names(rv$wide) })

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

  # -- Validation ------------------------------------------------------------
  validation_msg <- reactive({
    df <- long_data()
    if (nrow(df) == 0) return("Add at least one gene column and enter data.")
    if (any(is.na(df$Ct_target) | is.na(df$Ct_reference)))
      return("All Ct values must be numeric.")
    if (!input$ctrl_group %in% df$Group)
      return(paste0("Control group '", input$ctrl_group, "' not found in Group column."))
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
  analyzed <- reactive({
    req(is.null(validation_msg()))
    df <- compute_delta_ct(long_data())
    compute_fold_change(df, control_group = input$ctrl_group)
  })

  stats_result <- reactive({
    req(analyzed())
    run_stats(analyzed(), paired = input$paired,
              test_type = input$test_type,
              has_subgroup = has_subgroup())
  })

  genes_in_data <- reactive({ req(analyzed()); unique(analyzed()$Gene) })

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
                       control_group = input$ctrl_group,
                       rotate_x      = input$rotate_x,
                       aspect_ratio  = input$aspect_ratio,
                       ncol          = max(1, as.integer(input$facet_ncol %||% 3)),
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
                       control_group = input$ctrl_group,
                       rotate_x      = input$rotate_x,
                       aspect_ratio  = input$aspect_ratio,
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
    df$p_value <- format_pvalue(df$p_value)
    datatable(df, rownames = FALSE, options = list(dom = "t")) |>
      formatStyle("Significant", target = "row",
                  backgroundColor = styleEqual("Yes", "#d4edda"))
  })

  methods_txt <- reactive({
    req(stats_result(), analyzed())
    generate_methods_text(
      df = analyzed(), stats_df = stats_result(),
      test_type = input$test_type, paired = input$paired,
      control_group = input$ctrl_group, has_subgroup = has_subgroup()
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
                   control_group = input$ctrl_group,
                   rotate_x      = input$rotate_x,
                   aspect_ratio  = input$aspect_ratio,
                   y_min = if (is.na(input$y_min)) NULL else input$y_min,
                   y_max = if (is.na(input$y_max)) NULL else input$y_max))
  }

  render_combined_plot <- function() {
    make_combined_plot(analyzed(), stats_result(),
                       error_type    = input$error_type,
                       plot_type     = input$plot_type,
                       show_sig      = input$show_sig,
                       sig_comparisons_all = input$visible_comps,
                       control_group = input$ctrl_group,
                       rotate_x      = input$rotate_x,
                       aspect_ratio  = input$aspect_ratio,
                       ncol          = max(1, as.integer(input$facet_ncol %||% 3)),
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
                       control_group       = input$ctrl_group,
                       rotate_x            = input$rotate_x,
                       aspect_ratio        = input$aspect_ratio,
                       ncol = max(1, as.integer(input$facet_ncol %||% 3)),
                       y_min = if (is.na(input$y_min)) NULL else input$y_min,
                       y_max = if (is.na(input$y_max)) NULL else input$y_max)
  }

  download_dims <- function() {
    n    <- length(genes_in_data())
    ncol <- max(1, as.integer(input$facet_ncol %||% 3))
    nrow <- ceiling(n / ncol)
    list(w = 4.5 * ncol, h = 4.5 * nrow)
  }

  output$dl_svg <- downloadHandler(
    filename    = function() paste0("qpcr_plot_", Sys.Date(), ".svg"),
    contentType = "image/svg+xml",
    content     = function(file) {
      req(genes_in_data())
      d <- download_dims()
      ggsave(file, download_plot(), width = d$w, height = d$h,
             device = "svg", bg = "white")
    }
  )

  output$dl_png <- downloadHandler(
    filename    = function() paste0("qpcr_plot_", Sys.Date(), ".png"),
    contentType = "image/png",
    content     = function(file) {
      req(genes_in_data())
      d <- download_dims()
      ggsave(file, download_plot(), width = d$w, height = d$h,
             dpi = 300, device = "png", bg = "white")
    }
  )

  output$dl_pdf <- downloadHandler(
    filename    = function() paste0("qpcr_plot_", Sys.Date(), ".pdf"),
    contentType = "application/pdf",
    content     = function(file) {
      req(genes_in_data())
      d <- download_dims()
      ggsave(file, download_plot(), width = d$w, height = d$h,
             device = "pdf", bg = "white")
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
    if (!input$ctrl_group %in% df$Group)
      return(list(ok = FALSE,
                  msg = paste0("Control group '", input$ctrl_group, "' not found.")))
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
    compute_fold_change(df, control_group = input$ctrl_group)
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
                     control_group = input$ctrl_group,
                     rotate_x      = input$rotate_x,
                     aspect_ratio  = input$aspect_ratio,
                     ncol          = max(1, as.integer(input$facet_ncol %||% 3)),
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
    n    <- length(unique(rm_analyzed()$Gene))
    ncol <- max(1, as.integer(input$facet_ncol %||% 3))
    nrow <- ceiling(n / ncol)
    list(w = 5 * ncol, h = 5 * nrow)
  }

  output$rm_dl_svg <- downloadHandler(
    filename    = function() paste0("qpcr_rm_plot_", Sys.Date(), ".svg"),
    contentType = "image/svg+xml",
    content     = function(file) {
      req(rm_valid()$ok)
      d <- rm_download_dims()
      ggsave(file, render_rm_plot(), width = d$w, height = d$h,
             device = "svg", bg = "white")
    }
  )

  output$rm_dl_png <- downloadHandler(
    filename    = function() paste0("qpcr_rm_plot_", Sys.Date(), ".png"),
    contentType = "image/png",
    content     = function(file) {
      req(rm_valid()$ok)
      d <- rm_download_dims()
      ggsave(file, render_rm_plot(), width = d$w, height = d$h,
             dpi = 300, device = "png", bg = "white")
    }
  )

  output$rm_dl_pdf <- downloadHandler(
    filename    = function() paste0("qpcr_rm_plot_", Sys.Date(), ".pdf"),
    contentType = "application/pdf",
    content     = function(file) {
      req(rm_valid()$ok)
      d <- rm_download_dims()
      ggsave(file, render_rm_plot(), width = d$w, height = d$h,
             device = "pdf", bg = "white")
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
