# set up -----------------------------
# load libraries
library(shiny)
library(tidyverse)
library(data.table)

# source volcano plot script
source("/apps/volcano/volcano.R")

# look for data file
tsv <- list.files("/work", pattern = ".tsv", full.names = TRUE)
csv <- list.files("/work", pattern = ".csv", full.names = TRUE)
txt <- list.files("/work", pattern = ".txt", full.names = TRUE)

if(length(tsv) > 0) {
  data <- read.table(tsv, header = TRUE, sep = "\t")
}

if(length(csv) > 0) {
  data <- read.table(csv, header = TRUE, sep = ",")
}

if(length(txt) > 0) {
  data <- fread(txt)
}

data <- as.data.frame(data)

# functions
pvalue_candidate_f <- function(x) {
  if (class(data[[x]]) == "numeric") {
    if (max(data[[x]], na.rm = TRUE) <= 1) {
      if (min(data[[x]], na.rm = TRUE) >= 0) {
        return(TRUE)
      }
    }
  }
  return(FALSE)
}

logfc_candidate_f <- function(x) {
  if (class(data[[x]]) == "numeric") {
        return(TRUE)}
  return(FALSE)
}

gene_candidate_f <- function(x) {
  if (class(data[[x]]) == "character") {
    return(TRUE)}
  return(FALSE)
}

# check data columns to ID most likely candidates for pval, logfc, and gene ID
pval_cols <- names(data)[sapply(names(data), pvalue_candidate_f)]
logfc_cols <- names(data)[sapply(names(data), logfc_candidate_f)]
gene_cols <- names(data)[sapply(names(data), gene_candidate_f)]

# ui -----------------------------
ui <- fluidPage(
    # tab bar on main panel
    tabsetPanel(
      # VOLCANO PLOT PANEL -----
      tabPanel("Volcano Plot",
               h2("Interactive Volcano Plot"),
               sidebarLayout(
                 
                 # VOLCANO PLOT SIDE PANEL ------
                 sidebarPanel(width = 4,
                              
                              # SELECT AXES LABELS -----
                              h4("Select volcano plot axes:"),
                             
                              # select column for pval
                              selectInput("pvalue_col",
                                          "Input column for significance (y axis)",
                                          pval_cols,
                                          multiple = FALSE),
                                        
                              # select column for fold change
                              selectInput("logfc_col",
                                          "Input column for effect size (x axis)",
                                          logfc_cols,
                                          multiple = FALSE),
                              # SET PVAL AND LOGFC THRESHOLDS ----- 
                              h4("Set significance and effect size thresholds:"),
                                        
                              # set pvalue threshold 
                              sliderInput("pvalue_threshold",
                                          "Set significance threshold",
                                          min = 0,
                                          max = 1,
                                          value = .05),
                              
                              # set logfc threshold
                              uiOutput("logfc_slider"),
                              
                              # HIGHLIGHT GENES -----
                              h4("Highlight features of interest:"),
                              
                              # select column for gene ID input
                              selectInput("gene_col",
                                          "Select input column for feature label",
                                          gene_cols,
                                          multiple = FALSE),
                              
                              # gene selector menu
                              uiOutput("gene_selector"),
                              
                              # CUSTOMIZE PLOT -----
                              h4("Customize plot:"),
                              
                              # show/hide logfc and pval line
                              checkboxInput("show_pvalue_threshold",
                                            "Show significance threshold line",
                                            value = TRUE),
                              
                              # show/hide logfc lines
                              checkboxInput("show_logfc_threshold",
                                            "Show effect size threshold line",
                                            value = TRUE),
                              
                              # color differentially expressed genes
                              checkboxInput("color_by_de",
                                            "Color significantly different features",
                                            TRUE),
                              
                              # output ui for axis label inputs
                              uiOutput("y_axis_labeler"),
                              uiOutput("x_axis_labeler"),
                              
                              # label legend
                              textInput("legend_title",
                                        "Specify legend title",
                                        value = "Differentially Expressed")),
                 
                 # VOLCANO PLOT MAIN PANEL -----
                 mainPanel(
                   # output info from click
                   p("Hover over points to view features label, effect size, and significance."),
                   verbatimTextOutput("click_info",
                                      placeholder = TRUE),
                   
                   # output ggplot volcano
                   plotOutput("volcano_plot",
                                      width = "100%",
                                      height = "600px",
                                      hover = "volcano_hover"),
                           
                           # Download button for plot
                           downloadButton('download_volcano', 'Download volcano plot as PDF'),
                           
                           br(),
                           br(),
                           
                           # HIGHLIGHTED GENES TABLE -----
                           dataTableOutput("gene_highlight_tbl"))
                 
                 ) # end sidebarLayout
               ), # end volcano plot tabPanel
      
      # DATA PANEL -----
      tabPanel("Data",
               sidebarLayout(
                 
                 # DATA PANEL SIDEBAR
                 sidebarPanel(width = 3,
                              
                              # some text explanation
                              em("Threshold for what is considered differentially expressed is set in Volcano Plot tab by using 
                                 the significance and effect size sliders"),
                              
                              # Show differentiall expressed genes only
                              checkboxInput("show_de",
                                            "Show only significantly different features",
                                            FALSE)),
                 
                 # DATA PANEL MAIN PANEL
                 mainPanel(dataTableOutput("gene_data")))
               ) # end data tab panel
                ) # end tabsetPanel
) # end fluidPage

# server -------------------------
server <- function(input, output) {
  
  # IDENTIFY DIFFERENTIALLY EXPRESSED GENES -----
  
  # render UI for logfc slider
  # min and max set reactively with logfc based on selected logfc input col
  output$logfc_slider <- renderUI({
    sliderInput("logfc_threshold",
                "Select effect size threshold",
                min = 0,
                max = round(max(data[[input$logfc_col]])),
                value = 2,
                step = .5)
  })
  
  # use columns and thresholds selected in UI
  is_de <- reactive({
    abs(data[[input$logfc_col]]) >= input$logfc_threshold & data[[input$pvalue_col]] <= input$pvalue_threshold
  })
  
  
  # FILTERABLE DATAFRAME BY DE GENE -----
  
  # reactively filter data based on checkbox
  de_gene_data <- reactive({
    if (input$show_de) {
      filter(data, is_de())
    } else {
      data
    }
  })
  
  # render data frame of gene data
  output$gene_data <- renderDataTable(
    de_gene_data()
  )
  
  # X AND Y AXES LABELER -----
  
  # capture pvalue column selected and default value with it
  reactive_pvalue_value <- reactive({
    paste0("-log10(", input$pvalue_col, ")")
  })
  
  # enter custom x (logfc) axis label
  output$x_axis_labeler <- renderUI({
    textInput("x_axis_lab",
              "Specify X axis label",
              value = input$logfc_col,
              placeholder = "ex: Log Fold Change")
  })
    
    # enter custom x (logfc) axis label
    output$y_axis_labeler <- renderUI({
      textInput("y_axis_lab",
                "Specify Y axis label",
                value = reactive_pvalue_value(),
                placeholder = "ex: -log10(FDR)")
    })
  
  # HIGHLIGHTED GENE TABLE -----

  # select genes to highlight
  output$gene_selector <- renderUI({
    selectInput("highlight_genes",
                "Select feature(s) to highlight",
                sort(data[[input$gene_col]]),
                multiple = TRUE,
                selectize= TRUE)
  })
  
  # reactive function that subsets data by highlighted_gene vector
  highlight_gene_data <- reactive({
    if (length(input$highlight_genes) > 0) {
      highlight_gene_data <- data[data[[input$gene_col]] %in% input$highlight_genes, c(input$gene_col, input$logfc_col, input$pvalue_col)]
    } else {
      highlight_gene_data <- data.frame(NA, NA, NA)
      names(highlight_gene_data) <- c(input$gene_col, input$logfc_col, input$pvalue_col)
    }
  })

  # render a data table of highlighted genes info
  output$gene_highlight_tbl <- renderDataTable({
    highlight_gene_data()
  })
  
  # PLOT AND RENDER VOLCANO -----
  
  # volcano plot in reactive function (is this necessary?? can't be sure.)
  reactive_volcano <- reactive({
    plotVolcano(data = data, 
                logfc_col = input$logfc_col, 
                pval_col = input$pvalue_col,
                gene_col = input$gene_col,
                pval_thresh = input$pvalue_threshold,
                logfc_thresh = input$logfc_threshold,
                de_vec = is_de(),
                color_by_de = input$color_by_de,
                show_logfc_thresh = input$show_logfc_threshold,
                show_pvalue_thresh = input$show_pvalue_threshold,
                highlight_genes = input$highlight_genes,
                x_label = input$x_axis_lab,
                y_label = input$y_axis_lab,
                legend_title = input$legend_title)
  })

  # output volcano plot
  output$volcano_plot <- renderPlot({
    reactive_volcano()
  })
  
  # DISPLAY GENE INFO ON HOVER OVER -----
  
  # Create -log10 pvalue column
  data_w_log_pval <- reactive({
    # make new cols and select
    reduced_data <- data %>%
      mutate(log_pval = -log10(data[[input$pvalue_col]]))
  })
  
  # Collect nearpoint info and reduce to only gene_col, logfc_col and pvalue_col
  point_info <- reactive({
     nearpoint_out <- nearPoints(data_w_log_pval(), input$volcano_hover, xvar = input$logfc_col, yvar = data$log_pval)
     nearpoint_out %>%
       select(input$gene_col, input$logfc_col, input$pvalue_col)
   })
  
  # render printed text
  output$click_info <- renderPrint({
    point_info()
  })

  # DOWNLOAD HANDLER -----
  
  output$download_volcano <- downloadHandler(
    filename = function() {
      paste0("volcano-plot-", Sys.Date(), ".pdf")
    },

    content = function(file) {
      ggsave(file, reactive_volcano(), device = "pdf", width = 10, height = 5, units = "in")
    })
}

# build app ----------------------
shinyApp(ui, server)
