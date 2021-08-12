# set up -----------------------------
# load libraries
library(shiny)
library(tidyverse)
library(plotly)

# source volcano plot script
source("volcano.R")

# read in data
data <- read.table("210709_KO_vs_WT.txt", sep = "\t", header = T)

# functions
pvalue_candidate_f <- function(x) {
  if (class(data[[x]]) == "numeric") {
    if (max(data[[x]]) <= 1) {
      if (min(data[[x]]) >= 0) {
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
                                          "Input column for P-value (y axis)",
                                          pval_cols,
                                          multiple = FALSE),
                                        
                              # select column for fold change
                              selectInput("logfc_col",
                                          "Input column for effect size (x axis)",
                                          logfc_cols,
                                          multiple = FALSE),
                              # SET PVAL AND LOGFC THRESHOLDS ----- 
                              h4("Set differential gene thresholds:"),
                                        
                              # set pvalue threshold 
                              sliderInput("pvalue_threshold",
                                          "Set P value / FDR threshold",
                                          min = 0,
                                          max = 1,
                                          value = .05),
                              
                              # set logfc threshold
                              uiOutput("logfc_slider"),
                              
                              # CUSTOMIZE PLOT -----
                              h4("Customize plot:"),
                              
                              # show/hide logfc and pval line
                              checkboxInput("show_pvalue_threshold",
                                            "Show P value threshold line?",
                                            value = TRUE),
                              
                              # show/hide logfc lines
                              checkboxInput("show_logfc_threshold",
                                            "Show effect size threshold line?",
                                            value = TRUE),
                              
                              # color differentially expressed genes
                              checkboxInput("color_by_de",
                                            "Color differentially expressed genes",
                                            TRUE),
                              
                              # enter custom x (logfc) axis label
                              textInput("x_axis_lab",
                                        "Specify X axis label",
                                        placeholder = "ex: Log Fold Change"),
                              
                              # enter cutsom y (pval) axis label
                              textInput("y_axis_lab",
                                        "Specify Y axis label",
                                        placeholder = "ex: FDR"),
                              
                              # HIGHLIGHT GENES -----
                              h4("Highlight genes of interest:"),
                              
                              # select column for gene ID input
                              selectInput("gene_col",
                                          "Select input column for gene ID",
                                          gene_cols,
                                          multiple = FALSE),
                              
                              # gene selector menu
                              uiOutput("gene_selector")),
                 
                 # VOLCANO PLOT MAIN PANEL -----
                 mainPanel(plotOutput("volcano_plot",
                                      width = "100%",
                                      height = "600px"),
                           
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
               h1("Provided dataset"),
               sidebarLayout(
                 
                 # DATA PANEL SIDEBAR
                 sidebarPanel(width = 3,
                              
                              # Show differentiall expressed genes only
                              checkboxInput("show_de",
                                            "Show only differentially expressed genes",
                                            FALSE)),
                 
                 # DATA PANEL MAIN PANEL
                 mainPanel(dataTableOutput("gene_data")))
               ) # end data tab panel
                ) # end tabsetPanel
) # end fluidPage

# server -------------------------
server <- function(input, output) {
  
  # reactively mark DE genes from using pval and logfc
  is_de <- reactive({
    abs(data[[input$logfc_col]]) >= input$logfc_threshold & data[[input$pvalue_col]] <= input$pvalue_threshold
  })
  
  # reactively filter dataframe based on checkbox
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
  
  # select genes to highlight
  output$gene_selector <- renderUI({
    selectInput("highlight_genes",
                "Highlight gene(s)",
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
                y_label = input$y_axis_lab)
  })

  # output volcano plot
  output$volcano_plot <- renderPlot({
    reactive_volcano()
  })
  
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