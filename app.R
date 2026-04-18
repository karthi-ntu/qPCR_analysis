library(shiny)
library(ggplot2)
library(DT)
library(ggsignif)

source("R/analysis.R")
source("R/plot.R")

empty_wide <- function() {
  data.frame(Sample = character(0), Group = character(0),
             Ct_reference = numeric(0), stringsAsFactors = FALSE)
}

empty_long <- function() {
  data.frame(Sample = character(0), Group = character(0),
             Gene = character(0), Ct_target = numeric(0),
             Ct_reference = numeric(0), stringsAsFactors = FALSE)
}

to_long_df <- function(wide) {
  gene_cols <- setdiff(names(wide), c("Sample", "Group", "Ct_reference"))
  if (length(gene_cols) == 0 || nrow(wide) == 0) return(empty_long())
  rows <- lapply(seq_len(nrow(wide)), function(i) {
    lapply(gene_cols, function(g) {
      data.frame(
        Sample       = wide$Sample[i],
        Group        = wide$Group[i],
        Gene         = g,
        Ct_target    = suppressWarnings(as.numeric(wide[[g]][i])),
        Ct_reference = suppressWarnings(as.numeric(wide$Ct_reference[i])),
        stringsAsFactors = FALSE
      )
    })
  })
  do.call(rbind, Filter(Negate(is.null), unlist(rows, recursive = FALSE)))
}

to_wide_df <- function(long) {
  genes <- unique(long$Gene[nchar(trimws(long$Gene)) > 0])
  if (length(genes) == 0 || nrow(long) == 0) return(empty_wide())
  unique_idx   <- !duplicated(paste(long$Sample, long$Group, sep = "\t"))
  wide <- long[unique_idx, c("Sample", "Group", "Ct_reference")]
  for (g in genes) {
    gene_sub <- long[long$Gene == g, c("Sample", "Group", "Ct_target")]
    names(gene_sub)[3] <- g
    wide <- merge(wide, gene_sub, by = c("Sample", "Group"), all.x = TRUE, sort = FALSE)
  }
  wide
}

ui <- fluidPage(
  titlePanel("qPCR Analysis \u2014 \u0394\u0394Ct Method"),
  sidebarLayout(
    sidebarPanel(
      width = 4,
      textInput("ctrl_group", "Control group name", value = "Control"),
      radioButtons("error_type", "Error bars",
                   choices = c("SD", "SEM"), inline = TRUE),
      checkboxInput("paired", "Paired samples (2 groups only)", value = FALSE),
      fluidRow(
        column(6, numericInput("y_min", "Y-axis min", value = NA, step = 0.5)),
        column(6, numericInput("y_max", "Y-axis max", value = NA, step = 0.5))
      ),
      hr(),
      h5("Data"),
      p(style = "font-size:12px; color:#888;",
        "Select data in Excel (include header row: ",
        strong("Sample | Group | Reference gene | Gene1 | Gene2 | ..."),
        "), press Ctrl+C, then press Ctrl+V anywhere on this page."),
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
      fluidRow(
        column(4, downloadButton("dl_png", "Download Plot (PNG)")),
        column(4, downloadButton("dl_pdf", "Download Plot (PDF)")),
        column(4, downloadButton("dl_csv", "Download Stats (CSV)"))
      )
    )
  )
)

server <- function(input, output, session) {

  rv <- reactiveValues(wide = empty_wide())

  long_data <- reactive({ to_long_df(rv$wide) })

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
    num_cols <- setdiff(names(rv$wide), c("Sample", "Group"))
    for (col in num_cols)
      rv$wide[[col]] <- suppressWarnings(as.numeric(rv$wide[[col]]))
  })

  # -- Add / clear rows -------------------------------------------------------
  observeEvent(input$add_row, {
    gene_cols <- setdiff(names(rv$wide), c("Sample", "Group", "Ct_reference"))
    new_row <- data.frame(Sample = "", Group = "", Ct_reference = NA_real_,
                          stringsAsFactors = FALSE)
    for (g in gene_cols) new_row[[g]] <- NA_real_
    rv$wide <- rbind(rv$wide, new_row)
  })

  observeEvent(input$clear_tbl, {
    rv$wide <- empty_wide()
  })

  # -- Paste loader (global Ctrl+V) ------------------------------------------
  observeEvent(input$pasted_clipboard, {
    txt <- input$pasted_clipboard
    req(!is.null(txt), nchar(trimws(txt)) > 0)
    txt   <- gsub("\r\n", "\n", txt)
    txt   <- gsub("\r",   "\n", txt)
    lines <- strsplit(txt, "\n")[[1]]
    lines <- lines[nchar(trimws(lines)) > 0]
    req(length(lines) >= 2)

    headers    <- trimws(strsplit(lines[1], "\t")[[1]])
    req(length(headers) >= 4)
    gene_names <- headers[4:length(headers)]

    rows <- lapply(lines[-1], function(line) {
      fields <- trimws(strsplit(line, "\t")[[1]])
      if (length(fields) < 4) return(NULL)
      ct_ref <- suppressWarnings(as.numeric(fields[3]))
      lapply(seq_along(gene_names), function(i) {
        data.frame(Sample = fields[1], Group = fields[2], Gene = gene_names[i],
                   Ct_target    = suppressWarnings(as.numeric(fields[3 + i])),
                   Ct_reference = ct_ref, stringsAsFactors = FALSE)
      })
    })

    parsed <- do.call(rbind, Filter(Negate(is.null), unlist(rows, recursive = FALSE)))
    if (!is.null(parsed) && nrow(parsed) > 0)
      rv$wide <- to_wide_df(parsed)
  })

  # -- Validation -------------------------------------------------------------
  validation_msg <- reactive({
    df <- long_data()
    if (nrow(df) == 0)
      return("Add at least one gene column and enter data.")
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

  stats_result <- reactive({
    req(analyzed())
    run_stats(analyzed(), paired = input$paired)
  })

  genes_in_data <- reactive({
    req(analyzed())
    unique(analyzed()$Gene)
  })

  # -- Plots ------------------------------------------------------------------
  output$plots_ui <- renderUI({
    req(genes_in_data())
    do.call(tagList, lapply(genes_in_data(), function(gene)
      plotOutput(paste0("plot_", make.names(gene)), height = "400px")))
  })

  observe({
    req(genes_in_data())
    lapply(genes_in_data(), function(gene) {
      local({
        g <- gene
        output[[paste0("plot_", make.names(g))]] <- renderPlot({
          make_barplot(analyzed(), stats_result(), gene = g,
                       error_type = input$error_type,
                       y_min = if (is.na(input$y_min)) NULL else input$y_min,
                       y_max = if (is.na(input$y_max)) NULL else input$y_max)
        })
      })
    })
  })

  # -- Stats table ------------------------------------------------------------
  output$stats_table <- renderDT({
    req(stats_result())
    datatable(stats_result(), rownames = FALSE,
              options = list(dom = "t")) |>
      formatStyle("Significant", target = "row",
                  backgroundColor = styleEqual("Yes", "#d4edda"))
  })

  # -- Downloads --------------------------------------------------------------
  output$dl_png <- downloadHandler(
    filename = function() paste0("qpcr_plot_", Sys.Date(), ".png"),
    content  = function(file) {
      req(genes_in_data())
      n     <- length(genes_in_data())
      plots <- lapply(genes_in_data(), function(g)
        make_barplot(analyzed(), stats_result(), gene = g,
                     error_type = input$error_type,
                     y_min = if (is.na(input$y_min)) NULL else input$y_min,
                     y_max = if (is.na(input$y_max)) NULL else input$y_max))
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
      plots <- lapply(genes_in_data(), function(g)
        make_barplot(analyzed(), stats_result(), gene = g,
                     error_type = input$error_type,
                     y_min = if (is.na(input$y_min)) NULL else input$y_min,
                     y_max = if (is.na(input$y_max)) NULL else input$y_max))
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
}

shinyApp(ui, server)
