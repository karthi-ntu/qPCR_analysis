library(shiny)
library(ggplot2)
library(DT)
library(ggsignif)

source("R/analysis.R")
source("R/plot.R")

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

ui <- fluidPage(
  titlePanel("qPCR Analysis \u2014 \u0394\u0394Ct Method"),
  sidebarLayout(
    sidebarPanel(
      width = 4,
      textInput("ctrl_group", "Control group name", value = "Control"),
      radioButtons("test_type", "Statistical test",
                   choices = c("Auto (by normality)" = "auto",
                               "Parametric (t-test / ANOVA)" = "parametric",
                               "Non-parametric (Wilcoxon / Kruskal)" = "nonparametric"),
                   selected = "parametric"),
      radioButtons("plot_type", "Plot type",
                   choices = c("Scatter + error bars" = "scatter",
                               "Column + dots" = "column"),
                   inline = TRUE),
      radioButtons("error_type", "Error bars",
                   choices = c("SD", "SEM"), inline = TRUE),
      checkboxInput("paired", "Paired samples (2 groups only)", value = FALSE),
      checkboxInput("show_sig", "Show significance bars", value = TRUE),
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
    ),
    mainPanel(
      width = 8,
      uiOutput("plots_ui"),
      br(),
      h4("Statistics"),
      DTOutput("stats_table"),
      br(),
      h4("Normality (Shapiro-Wilk)"),
      DTOutput("norm_table"),
      br(),
      h4("Methods (copy-paste for papers)"),
      verbatimTextOutput("methods_text"),
      br(),
      fluidRow(
        column(3, downloadButton("dl_png", "Plot (PNG)")),
        column(3, downloadButton("dl_pdf", "Plot (PDF)")),
        column(3, downloadButton("dl_csv", "Stats (CSV)")),
        column(3, downloadButton("dl_methods", "Methods (TXT)"))
      )
    )
  )
)

server <- function(input, output, session) {

  rv <- reactiveValues(wide = empty_wide())

  long_data <- reactive({ to_long_df(rv$wide) })

  has_subgroup <- reactive({ "Subgroup" %in% names(rv$wide) })

  # -- Table ------------------------------------------------------------------
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
    # Reorder columns to match rv$wide
    new_row <- new_row[, names(rv$wide), drop = FALSE]
    rv$wide <- rbind(rv$wide, new_row)
  })

  observeEvent(input$clear_tbl, {
    rv$wide <- empty_wide()
  })

  # -- Paste loader ----------------------------------------------------------
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

    # Detect Subgroup: if 3rd header (case-insensitive) is "subgroup" / "factor2"
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

  # -- Validation -------------------------------------------------------------
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

  # -- Analysis ---------------------------------------------------------------
  analyzed <- reactive({
    req(is.null(validation_msg()))
    df <- compute_delta_ct(long_data())
    compute_fold_change(df, control_group = input$ctrl_group)
  })

  norm_result <- reactive({
    req(analyzed())
    run_normality(analyzed())
  })

  stats_result <- reactive({
    req(analyzed())
    run_stats(analyzed(), paired = input$paired,
              test_type = input$test_type,
              norm_df = norm_result(),
              has_subgroup = has_subgroup())
  })

  genes_in_data <- reactive({
    req(analyzed())
    unique(analyzed()$Gene)
  })

  # -- Dynamic UI: hide comparisons ------------------------------------------
  output$hide_comps_ui <- renderUI({
    req(stats_result())
    sig <- stats_result()
    sig <- sig[sig$Significant == "Yes" & grepl(" vs ", sig$Comparison, fixed = TRUE), ]
    if (nrow(sig) == 0) return(NULL)
    choices <- unique(paste0(sig$Gene, ": ", sig$Comparison))
    checkboxGroupInput("visible_comps",
                       label = "Significance bars to show",
                       choices  = choices,
                       selected = choices)
  })

  # Filter comparisons per-gene based on visible_comps selection
  visible_comps_for_gene <- function(gene) {
    sel <- input$visible_comps
    if (is.null(sel)) return(NULL)
    prefix <- paste0(gene, ": ")
    comps  <- sub(prefix, "", sel[startsWith(sel, prefix)], fixed = TRUE)
    comps
  }

  # -- Plots ------------------------------------------------------------------
  output$plots_ui <- renderUI({
    req(genes_in_data())
    do.call(tagList, lapply(genes_in_data(), function(gene)
      plotOutput(paste0("plot_", make.names(gene)), height = "450px")))
  })

  observe({
    req(genes_in_data())
    lapply(genes_in_data(), function(gene) {
      local({
        g <- gene
        output[[paste0("plot_", make.names(g))]] <- renderPlot({
          make_barplot(analyzed(), stats_result(), gene = g,
                       error_type = input$error_type,
                       plot_type  = input$plot_type,
                       show_sig   = input$show_sig,
                       sig_comparisons = visible_comps_for_gene(g),
                       y_min = if (is.na(input$y_min)) NULL else input$y_min,
                       y_max = if (is.na(input$y_max)) NULL else input$y_max)
        })
      })
    })
  })

  # -- Stats / Normality tables ----------------------------------------------
  output$stats_table <- renderDT({
    req(stats_result())
    df <- stats_result()
    df$p_value <- format_pvalue(df$p_value)
    datatable(df, rownames = FALSE,
              options = list(dom = "t")) |>
      formatStyle("Significant", target = "row",
                  backgroundColor = styleEqual("Yes", "#d4edda"))
  })

  output$norm_table <- renderDT({
    req(norm_result())
    datatable(norm_result(), rownames = FALSE,
              options = list(dom = "t")) |>
      formatStyle("Normal", target = "row",
                  backgroundColor = styleEqual("No", "#f8d7da"))
  })

  # -- Methods text -----------------------------------------------------------
  methods_txt <- reactive({
    req(stats_result(), analyzed(), norm_result())
    generate_methods_text(
      df = analyzed(), stats_df = stats_result(),
      norm_df = norm_result(), test_type = input$test_type,
      paired = input$paired, control_group = input$ctrl_group,
      has_subgroup = has_subgroup()
    )
  })

  output$methods_text <- renderText({ methods_txt() })

  # -- Downloads --------------------------------------------------------------
  render_all_plots <- function() {
    lapply(genes_in_data(), function(g)
      make_barplot(analyzed(), stats_result(), gene = g,
                   error_type = input$error_type,
                   plot_type  = input$plot_type,
                   show_sig   = input$show_sig,
                   sig_comparisons = visible_comps_for_gene(g),
                   y_min = if (is.na(input$y_min)) NULL else input$y_min,
                   y_max = if (is.na(input$y_max)) NULL else input$y_max))
  }

  output$dl_png <- downloadHandler(
    filename = function() paste0("qpcr_plot_", Sys.Date(), ".png"),
    content  = function(file) {
      req(genes_in_data())
      n     <- length(genes_in_data())
      plots <- render_all_plots()
      png(file, width = 3200, height = 2000 * n, res = 600)
      on.exit(dev.off(), add = TRUE)
      for (p in plots) print(p)
    }
  )

  output$dl_pdf <- downloadHandler(
    filename = function() paste0("qpcr_plot_", Sys.Date(), ".pdf"),
    content  = function(file) {
      req(genes_in_data())
      n     <- length(genes_in_data())
      plots <- render_all_plots()
      pdf(file, width = 7, height = 5 * n)
      on.exit(dev.off(), add = TRUE)
      for (p in plots) print(p)
    }
  )

  output$dl_csv <- downloadHandler(
    filename = function() paste0("qpcr_stats_", Sys.Date(), ".csv"),
    content  = function(file) {
      req(stats_result())
      write.csv(stats_result(), file, row.names = FALSE)
    }
  )

  output$dl_methods <- downloadHandler(
    filename = function() paste0("qpcr_methods_", Sys.Date(), ".txt"),
    content  = function(file) {
      req(methods_txt())
      writeLines(methods_txt(), file)
    }
  )
}

shinyApp(ui, server)
