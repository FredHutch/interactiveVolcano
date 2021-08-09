# load libraries
library(shiny)
library(tidyverse)
library(plotly)

# source scripts
source("volcano.R")

# read in data
data <- read.table("./data/210709_KO_vs_WT.txt", sep = "\t", header = T)

# ui -----------------------------
ui <- fluidPage(
  # title
  #titlePanel("Data Core Interactive Volcano App"),
  # main panel
  mainPanel(
    # tab bar on main panel
    tabsetPanel(
                tabPanel("Splash page",
                         p("welcome to the carousel volcano plot app, 
                           info about how this visualizaiton is created, where to find code, etc")),
                # Data panel
                tabPanel("Data",
                         h1("Provided dataset"),
                         dataTableOutput("gene_data")),
                # Volcano plot panel
                tabPanel("Volcano Plot",
                         h1("Interactive Volcano Plot"),
                         sidebarLayout(
                           sidebarPanel(width = 3,
                                        selectInput("gene_list",
                                                    "Highlight a gene",
                                                    sort(data[[gene_col]]),
                                                    multiple = TRUE,
                                                    selectize= TRUE)),
                         mainPanel(plotOutput("volcano_plot",
                                      click = "volcano_click"),
                                   verbatimTextOutput("info"))))
                ) # end tabsetPanel
    ) # end mainPanel
) # end fluidPage

# server -------------------------
server <- function(input, output) {
  # output uploaded dataframe
  output$gene_data <- renderDataTable(
    data
  )
  
  # volcano plot output
  output$volcano_plot <- renderPlot({
    plotVolcano(data, "logFC", "FDR", "Genename")
  })
  
  # output nearPoint to volcano plot click as datatable  
  output$info <- renderText(
    paste0(input$volcano_click))
}

# build app ----------------------
shinyApp(ui, server)