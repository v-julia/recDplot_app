library(shiny)
library(DT)
library(shinyjs)
library(shinycssloaders)
library(shinythemes)
source("scripts/rec_plots.R")

options(shiny.maxRequestSize = 40*1024^2)


# Run the application 
#shinyApp(ui = ui, server = server)

#shinyApp(
#  ui = fluidPage(DT::dataTableOutput('tbl')),
#  server = function(input, output) {
#    output$tbl = DT::renderDataTable(
#      {iris}, options = list(pageLength = 20)
#    )
#  }
#)#