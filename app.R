library(shiny)
library(DT)
source("scripts/rec_plots.R")


# Define UI for application that draws a histogram
ui <- navbarPage(
        "Exploring recombination in genetic sequences using distance plots",
         tabPanel("RMSE matrix plot",
                  sidebarLayout(
                    sidebarPanel(
                      
                      fileInput(
                        "file_alignment", "Upload alignment in fasta-format", accept = c(".fasta", ".fas")
                      ),
                      numericInput(
                        "window", 
                        "Size of sliding window", 
                        min = 50, max = 10000, value = 200
                      ),
                      numericInput(
                        "step", 
                        "Size of step", 
                        min = 10, max = 10000, value = 50
                      ),
                      actionButton("goButton", "Run")
                    ),
                    mainPanel(
                      plotOutput("rmse_matrix_plot")
                    )
                  )
         ),
         
         tabPanel("Distance plot",
                  sidebarLayout(
                    sidebarPanel(
                      h4("Positions of fragment 1"),
                      fluidRow(
                        column(5,
                               numericInput(
                                 "start1", 
                                 "Start", 
                                 min = 1, max = 1000000, value = 1
                               )
                        ),
                        column(5,
                               numericInput(
                                 "end1", 
                                 "End", 
                                 min = 50, max = 10000, value = 100
                               )
                        )
                      ),
                      h4("Positions of fragment 2"),
                      fluidRow(
                        column(5,
                               numericInput(
                                 "start2", 
                                 "Start", 
                                 min = 1, max = 1000000, value = 100
                               )
                        ),
                         column(5,
                                numericInput(
                                  "end2", 
                                  "End", 
                                  min = 10, max = 1000000, value = 200
                                )
                         )
                      ),
                      actionButton("goButton2", "Run")
                    ),
                    mainPanel(
                      fluidRow(
                        splitLayout(
                          cellWidths = 500, 
                          tagList(h4("Control distance plot"), plotOutput("control")
                          ),
                          tagList(h4("Distance plot"), plotOutput("dist_plot",
                                                                  brush = brushOpts(id = "plot1_brush"),
                                                                  click = clickOpts(id = "plot1_click")
                                                                  )
                          )
                        )
                    ),
                    fluidRow(
                      #verbatimTextOutput("mm"),
                      #DT::dataTableOutput("min_max"),
                      DT::dataTableOutput("brush_info")
                    )
                  )
         )
        ),
        tabPanel("About",
                   mainPanel(
                     p("The description of this app")
                   )
        )
)


   


   
   
   #sidebarLayout(
  #   sidebarPanel(
  #     fileInput(
  #       "file_alignment", "Upload alignment in fasta-format", accept = c(".fasta")
  #     ),
  #     numericInput(
  #       "window", 
  #       "Size of sliding window", 
  #       min = 50, max = 10000, value = 200
  #     ),
  #     numericInput(
  #       "step", 
  #       "Size of step", 
  #       min = 10, max = 10000, value = 50
  #     ),
  #     actionButton("goButton", "Run")
  #
  #     
  #   ),
  #   
  #   
  #   mainPanel(
  #     p("Short description of this app"),
  #     textOutput("selected_var"),
  #     plotOutput("rmse_matrix_plot", click = "plot_click"),
  #     verbatimTextOutput("info")
  #     
  #     
  #     )
      
     
     
  # )

   


# Define server logic required to draw a histogram
server <- function(input, output) {

  observeEvent(input$goButton, {
    aln <- read.dna(as.character(input$file_alignment$datapath), format="fasta")
    output$rmse_matrix_plot <- renderPlot({
      
      #draws rmse matrix
      matrix_rmse = plot_rmse(aln, input$step, input$window, "pdist", "pairwise")
      heatmap.2(as.matrix(matrix_rmse[nrow(matrix_rmse):1,]), Rowv = FALSE, Colv = "Rowv", dendrogram = 'none', col=matlab.like, tracecol=NA)
    })
    
    #np <- nearPoints(matrix_rmse, input$plot_click)
    
    #output$info <- renderText({
    #  paste0("x=", input$plot_click$row, "\ny=", input$plot_click$col)
    #})
    
  }
  )
  
  
  observeEvent(
    input$goButton2, {
      aln = read.dna(as.character(input$file_alignment$datapath), format="fasta")
      output$control <- renderPlot(
                                  plot_control(aln)
                        )
      
      l = plot_dist_test(aln, input$start1,input$end1,input$start2,input$end2)
      
      
      output$dist_plot <- renderPlot(l[[1]])
      df = l[[2]]
      
      observeEvent(input$plot1_brush, {
        brushed_points <- brushedPoints(df, input$plot1_brush)
        min_1 = input$plot1_brush$xmin
        max_1 = input$plot1_brush$xmax
        
        min_2 = input$plot1_brush$ymin
        max_2 = input$plot1_brush$ymax
        
        #output$mm <- renderPrint(min_1)
        
        #output$mm <- renderText(paste(toString(min(brushed_points[2,])), toString(max(brushed_points[2,])), sep=","))
        
        output$min_max <- DT::renderDataTable(brushed_points)
        
        output$brush_info <- DT::renderDataTable({
          #brushedPoints(df, input$plot1_brush)
          #typeof(l[[3]])

          #brushed_points[[2]]
          #c(min_1, max_1, min_2, max_2)
          find_recomb_names(l[[3]], min_1, max_1, l[[4]], min_2, max_2)
          
        })
        
      })

      
    }
  )
}

# Run the application 
shinyApp(ui = ui, server = server)