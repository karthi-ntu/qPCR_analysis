library(shiny)
library(ggplot2)
library(DT)
library(ggsignif)

source("R/analysis.R")
source("R/plot.R")

empty_table <- function() {
  data.frame(
    Sample       = character(0),
    Group        = character(0),
    Gene         = character(0),
    Ct_target    = numeric(0),
    Ct_reference = numeric(0),
    stringsAsFactors = FALSE
  )
}

ui <- fluidPage(
  titlePanel("qPCR Analysis \u2014 \u0394\u0394Ct Method"),
  sidebarLayout(
    sidebarPanel(
      width = 4,
      textInput("ctrl_group",  "Control group name",            value = "Control"),
      radioButtons("error_type", "Error bars",
                   choices = c("SD", "SEM"), inline = TRUE),
      checkboxInput("paired", "Paired samples (2 groups only)", value = FALSE),
      fluidRow(
        column(6, numericInput("y_min", "Y-axis min", value = NA, step = 0.5)),
        column(6, numericInput("y_max", "Y-axis max", value = NA, step = 0.5))
      ),
      hr(),
      h5("Paste from Excel"),
      p(style = "font-size:12px; color:#888;",
        "Copy cells from Excel in this column order:",
        br(), strong("Sample | Group | Gene | Ct_target | Ct_reference")),
      textAreaInput("paste_data", label = NULL,
                    placeholder = "Paste here (Ctrl+V), then click Load",
                    rows = 4, width = "100%"),
      actionButton("load_paste", "Load pasted data", icon = icon("upload"),
                   style = "width:100%; margin-bottom:6px;"),
      hr(),
      h5("Or enter row by row"),
      DTOutput("table"),
      br(),
      actionButton("add_row",   "Add row",    icon = icon("plus")),
      actionButton("clear_tbl", "Clear table", icon = icon("trash")),
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

  rv <- reactiveValues(data = empty_table())

  # -- Editable data table ----------------------------------------------------
  output$table <- renderDT({
    datatable(
      rv$data,
      editable  = TRUE,
      rownames  = FALSE,
      selection = "none",
      options   = list(dom = "t", pageLength = 50, scrollY = "300px")
    )
  })

  observeEvent(input$table_cell_edit, {
    info     <- input$table_cell_edit
    rv$data  <- editData(rv$data, info, proxy = "table", rownames = FALSE)
    rv$data$Ct_target    <- suppressWarnings(as.numeric(rv$data$Ct_target))
    rv$data$Ct_reference <- suppressWarnings(as.numeric(rv$data$Ct_reference))
  })

  observeEvent(input$add_row, {
    rv$data <- rbind(rv$data, data.frame(
      Sample       = "",
      Group        = "",
      Gene         = "",
      Ct_target    = NA_real_,
      Ct_reference = NA_real_,
      stringsAsFactors = FALSE
    ))
  })

  observeEvent(input$clear_tbl, {
    rv$data <- empty_table()
  })

  observeEvent(input$load_paste, {
    txt <- trimws(input$paste_data)
    req(nchar(txt) > 0)
    lines <- strsplit(txt, "\n")[[1]]
    rows <- lapply(lines, function(line) {
      fields <- strsplit(trimws(line), "\t")[[1]]
      if (length(fields) < 5) return(NULL)
      data.frame(
        Sample       = trimws(fields[1]),
        Group        = trimws(fields[2]),
        Gene         = trimws(fields[3]),
        Ct_target    = suppressWarnings(as.numeric(fields[4])),
        Ct_reference = suppressWarnings(as.numeric(fields[5])),
        stringsAsFactors = FALSE
      )
    })
    parsed <- do.call(rbind, Filter(Negate(is.null), rows))
    if (!is.null(parsed) && nrow(parsed) > 0) {
      rv$data <- parsed
      updateTextAreaInput(session, "paste_data", value = "")
    }
  })

  # -- Validation -------------------------------------------------------------
  validation_msg <- reactive({
    df <- rv$data
    if (nrow(df) == 0 || all(df$Gene == ""))
      return("Enter at least one row of data.")
    df <- df[df$Gene != "", ]
    if (any(is.na(df$Ct_target) | is.na(df$Ct_reference)))
      return("All Ct_target and Ct_reference values must be numeric.")
    if (!input$ctrl_group %in% df$Group)
      return(paste0("Control group '", input$ctrl_group, "' not found in Group column."))
    grp_counts <- table(df$Gene, df$Group)
    if (any(grp_counts[grp_counts > 0] < 2))
      return("At least 2 samples per group per gene required for statistical testing.")
    NULL
  })

  output$error_msg <- renderUI({
    msg <- validation_msg()
    if (!is.null(msg))
      div(style = "color:red; margin-top:8px;", icon("exclamation-triangle"), " ", msg)
  })

  # -- Core analysis ----------------------------------------------------------
  analyzed <- reactive({
    req(is.null(validation_msg()))
    df <- rv$data[rv$data$Gene != "", ]
    df <- compute_delta_ct(df)
    df <- compute_fold_change(df, control_group = input$ctrl_group)
    df
  })

  stats_result <- reactive({
    req(analyzed())
    run_stats(analyzed(), paired = input$paired)
  })

  summary_result <- reactive({
    req(analyzed())
    # Preserve control-first order for each gene
    df <- analyzed()
    gene_group_order <- unique(df[, c("Gene", "Group")])
    sm <- summarize_groups(df)
    sm$Group <- factor(sm$Group,
                       levels = unique(gene_group_order$Group[
                         gene_group_order$Group %in% sm$Group]))
    sm$Group <- as.character(sm$Group)
    sm
  })

  genes_in_data <- reactive({
    req(analyzed())
    unique(analyzed()$Gene)
  })

  # -- Per-gene plots (dynamic UI) --------------------------------------------
  output$plots_ui <- renderUI({
    req(genes_in_data())
    plot_output_list <- lapply(genes_in_data(), function(gene) {
      plotOutput(paste0("plot_", make.names(gene)), height = "400px")
    })
    do.call(tagList, plot_output_list)
  })

  observe({
    req(genes_in_data())
    lapply(genes_in_data(), function(gene) {
      local({
        g <- gene
        output[[paste0("plot_", make.names(g))]] <- renderPlot({
          make_barplot(summary_result(), stats_result(),
                       gene = g, error_type = input$error_type,
                       y_min = if (is.na(input$y_min)) NULL else input$y_min,
                       y_max = if (is.na(input$y_max)) NULL else input$y_max)
        })
      })
    })
  })

  # -- Statistics table -------------------------------------------------------
  output$stats_table <- renderDT({
    req(stats_result())
    datatable(stats_result(), rownames = FALSE,
              options = list(dom = "t")) |>
      formatStyle("Significant",
                  target          = "row",
                  backgroundColor = styleEqual("Yes", "#d4edda"))
  })

  # -- Downloads --------------------------------------------------------------
  output$dl_png <- downloadHandler(
    filename = function() paste0("qpcr_plot_", Sys.Date(), ".png"),
    content  = function(file) {
      req(genes_in_data())
      n <- length(genes_in_data())
      plots <- lapply(genes_in_data(), function(g)
        make_barplot(summary_result(), stats_result(),
                     gene = g, error_type = input$error_type,
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
      n <- length(genes_in_data())
      plots <- lapply(genes_in_data(), function(g)
        make_barplot(summary_result(), stats_result(),
                     gene = g, error_type = input$error_type,
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
