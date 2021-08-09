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
                         h3("welcome to the carousel volcano plot app, 
                           info about how this visualizaiton is created, where to find code, etc")),
                # Data panel
                tabPanel("Data",
                         h1("Provided dataset"),
                         sidebarLayout(
                           sidebarPanel( width = 3,
                                         checkboxInput("show_de",
                                                       "Show only differentially expressed genes",
                                                       FALSE)),
                           mainPanel(dataTableOutput("gene_data")))),
                         
                # Volcano plot panel
                tabPanel("Volcano Plot",
                         h1("Interactive Volcano Plot"),
                         sidebarLayout(
                           sidebarPanel(width = 3,
                                        selectInput("gene_list",
                                                    "Highlight gene(s)",
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
  
  # filter dataframe based on checkbox
  filtered_gene_data <- reactive({
    if (input$show_de) {
      filter(data, is.DE != 0)
    } else {
      data
    }
  })
  
  # output dataframe
  output$gene_data <- renderDataTable(
    filtered_gene_data()
  )

  # output volcano plot
  output$volcano_plot <- renderPlot({
    plotVolcano(data, "logFC", "FDR", "Genename")
  })
  
  # output volcano_click info as text  
  output$info <- renderText({
    paste0(input$volcano_click)
    })
}

# build app ----------------------
shinyApp(ui, server)