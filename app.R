#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)

# Define UI for application that draws a histogram
ui <- fluidPage(
   
   # Application title
   titlePanel("Matrix with recombination distance plots"),
   
   sidebarLayout(
     sidebarPanel(
       file_alignment <- fileInput(
         "upload", "Upload an alignment in fasta-format", accept = c(".fasta", ".fas")
       ),
       window <- numericInput(
         "window", 
         "Size of sliding window", 
         min = 50, max = 10000, value = 200
       ),
       step <- numericInput(
         "step", 
         "Size of step", 
         min = 10, max = 10000, value = 50
       )
     ),
     mainPanel(
       p("Short description of this app"),
       plotOutput("rmse_matrix")
       )
   )
)
   
   

# Define server logic required to draw a histogram
server <- function(input, output) {
   
   output$distPlot <- renderPlot({
      # generate bins based on input$bins from ui.R
      x    <- faithful[, 2] 
      bins <- seq(min(x), max(x), length.out = input$bins + 1)
      
      # draw the histogram with the specified number of bins
      hist(x, breaks = bins, col = 'darkgray', border = 'white')
   })
}

# Run the application 
shinyApp(ui = ui, server = server)

