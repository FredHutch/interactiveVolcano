# load libraries
library(shiny)
library(tidyverse)
library(plotly)

# source scripts
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
  # main panel
  mainPanel(
    # tab bar on main panel
    tabsetPanel(
                tabPanel("Splash page",
                         h3("welcome to the carousel volcano plot app, 
                           info about how this visualizaiton is created, where to find code, etc")), # end tabPanel
                
                # Volcano plot panel
                tabPanel("Volcano Plot",
                         h1("Interactive Volcano Plot"),
                         sidebarLayout(
                           sidebarPanel(width = 3,
                                        # select column for pval
                                        selectInput("pvalue_col",
                                                    "Select input column for P value",
                                                    pval_cols,
                                                    multiple = FALSE),
                                        # select column for fold change
                                        selectInput("logfc_col",
                                                    "Select input column for fold change",
                                                    logfc_cols,
                                                    multiple = FALSE),
                                        # select column for gene ID
                                        selectInput("gene_col",
                                                    "Select input column for gene ID",
                                                    gene_cols,
                                                    multiple = FALSE),
                                        # set pvalue threshold 
                                        sliderInput("pvalue_threshold",
                                                    "Select P value / FDR threshold",
                                                    min = 0,
                                                    max = 1,
                                                    value = .05),
                                        # set logfc threshold
                                        uiOutput("logfc_slider"),
                                        # show/hide logfc and pval lines
                                        checkboxInput("show_pvalue_threshold",
                                                      "Show P value threshold?",
                                                      value = TRUE),
                                        checkboxInput("show_logfc_threshold",
                                                      "Show log fold change threshold?",
                                                      value = TRUE)),
                                        # select genes to highlight
                                        #uiOutput("gene_selector")),
                         mainPanel(plotOutput("volcano_plot",
                                      click = "volcano_click"),
                                   verbatimTextOutput("info"))
                         ) # end sidebarLayout
                         ), # end tabPanel
                # Data panel
                tabPanel("Data",
                         h1("Provided dataset"),
                         sidebarLayout(
                           sidebarPanel(width = 3,
                                        checkboxInput("show_de",
                                                       "Show only differentially expressed genes",
                                                       FALSE)),
                           mainPanel(dataTableOutput("gene_data"))))
                ) # end tabsetPanel
    ) # end mainPanel
) # end fluidPage

# server -------------------------
server <- function(input, output) {
  
  # reactively mark DE genes from using pval and logfc
  is_de <- reactive({
    abs(data[[input$logfc_col]]) >= input$logfc_threshold & data[[input$pvalue_col]] <= input$pvalue_threshold
  })
  
  # reactively filter dataframe based on checkbox
  filtered_gene_data <- reactive({
    if (input$show_de) {
      filter(data, is_de())
    } else {
      data
    }
  })
  
  # render data frame of gene data
  output$gene_data <- renderDataTable(
    filtered_gene_data()
  )
  
  # render UI for logfc slider
  # min and max set reactively with logfc based on selected logfc input col
  output$logfc_slider <- renderUI({
    sliderInput("logfc_threshold",
                "Select log fold change threshold",
                min = 0,
                max = round(max(data[[input$logfc_col]])),
                value = 2)
  })
  
  # # gene
  # output$gene_selector <- renderUI({
  #   selectInput("selected_genes",
  #               "Highlight gene(s)",
  #               sort(data[[input$gene_col]]),
  #               multiple = TRUE,
  #               selectize= TRUE)
  # })
  
  # output volcano plot
  output$volcano_plot <- renderPlot({
    plotVolcano(data = data, 
                logfc_col = input$logfc_col, 
                pval_col = input$pvalue_col,
                gene_col = input$gene_col,
                pval_thresh = input$pvalue_threshold,
                logfc_thresh = input$logfc_threshold,
                de_vec = is_de(),
                show_logfc_thresh = input$show_logfc_threshold,
                show_pvalue_thresh = input$show_pvalue_threshold)
                
  })
  
  # output volcano_click info as text  
  output$info <- renderText({
    paste0(input$volcano_click)
    })
}

# build app ----------------------
shinyApp(ui, server)