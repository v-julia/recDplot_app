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
                        min = 50, max = 10000, value = 500
                      ),
                      numericInput(
                        "step", 
                        "Size of step", 
                        min = 10, max = 10000, value = 200
                      ),
                      actionButton("goButton", "Run")
                    ),
                    mainPanel(
                      plotOutput("rmse_matrix_plot"),
                      width = 6
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

   
# plots heatmap with RMSE in pairwise distance comparison plot for each pair of genomic regions
# dna_object -  list of DNA sequences (class DNAbin)
# step
# window - length of genomic regions to compare
# method - method of calculation distances ("pdist", "JC", "Kimura", "TN")
# modification - pairwise deletion of positions with gaps or not
# updateProgress - function 

#returns matrix with rmse values for each pair f=of genomic regions

create_rmse = function(dna_object, step,window, method, modification=NA, updateProgress = NULL){
  
  length_aln = length(dna_object[1,]) #length of alignment
  num_seq = length(dna_object[,1]) # number of sequences in alignment
  
  starts = seq(from=0, to=length_aln-window, by = step) # start positions of genomic regions
  starts[1]=1
  ends = seq(from=window, to = length_aln, by = step) # end positions of genomic regions
  if (length_aln%%step>step){ends=c(ends,length_aln)}
  
  df_intervals = cbind(starts,ends) #intervals
  
  #names = apply(df_intervals, 1, function(x){paste(toString(x[1]),toString(x[2]),sep="_")})
  
  #dataframe to store RMSE values of each comparison
  rmse_df = data.frame(matrix(ncol=length(starts), nrow = length(starts)))
  colnames(rmse_df)=starts
  rownames(rmse_df)=starts
  
  
  
  #list of distance matrices for each pair of genomic regions
  n = nrow(df_intervals)
  
  dist_matrices = list()
  for (i in 1:n){
    slice = dna_object[1:num_seq, seq(from = df_intervals[i,"starts"], to = df_intervals[i,"ends"], by=1)]
    
    #dist_matrices[[i]] = dist.dna(slice,  as.matrix = TRUE,  model = "JC69")
    if (method == "pdist"){
      if (modification=="pairwise"){
        dist_matrices[[i]] = dist.gene(slice, method = "percentage",  pairwise.deletion = TRUE)}
      else {
        dist_matrices[[i]] = dist.gene(slice, method = "percentage",  pairwise.deletion = FALSE)}
    }
    
    else {
      if (method == "JC"){dist_matrices[[i]] = dist.dna(slice,  as.matrix = TRUE,  model = "JC69")}
      if (method == "Kimura"){dist_matrices[[i]] = dist.dna(slice,  as.matrix = TRUE,  model = "K80")}
      if (method == "TN"){dist_matrices[[i]] = dist.dna(slice,  as.matrix = TRUE,  model = "K80")}
      #else{print("Unknown method")}
      
    }
    
    if (is.function(updateProgress)){
      incProgress(0.5*(1/n), detail = "Calculating genetic distances in windows")

    }
  }
  n = nrow(df_intervals)
  print(length(dist_matrices))
  for (i in 1:n){
    for (j in 1:(n - i + 1)){
      #for (j in 1:(n)){
      #print(paste(toString(i), toString(j), sep=","))
      #fits pairwise distance comparison plots linear model, calculates rmse
      rmse_i_j = (rmse(lm(dist_matrices[[j]]~dist_matrices[[i]])) + rmse(lm(dist_matrices[[i]]~dist_matrices[[j]]))) /2.0
      #rmse_i_j = rmse(lm(dist_matrices[[j]]~dist_matrices[[i]]))
      rmse_df[i,j] = rmse_i_j
      rmse_df[n-j+1,n-i+1] = rmse_i_j
      
      if (is.function(updateProgress)){
        incProgress(0.5*(2/(n*(n-1))), detail = "Calculating rmse for windows pairs")

      }
      
    }
  }
  #print(rmse_df)
  #colnames(rmse_df)
  return(rmse_df)
  
}




# Define server logic required to draw a histogram
server <- function(input, output) {

  observeEvent(input$goButton, {
    
    aln <- read.dna(as.character(input$file_alignment$datapath), format="fasta", as.character=TRUE)
    aln[aln=='-'] <- NA
    output$rmse_matrix_plot <- renderPlot({
      
      
      # Create a Progress object
      #progress <- shiny::Progress$new()
      
      # Make sure it closes when we exit this reactive, even if there's an error
      #on.exit(progress$close())
      
      #progress$set(message = "Making rmse matrix plot", value = 0)
      withProgress(message = 'Creating rmse distance matrix', value = 0, {
        
        # Create a closure to update progress.
        # Each time this is called:
        # - If `value` is NULL, it will move the progress bar 1/5 of the remaining
        #   distance. If non-NULL, it will set the progress to that value.
        # - It also accepts optional detail text.
        updateProgress <- function(value = NULL, detail = NULL) {
          if (is.null(value)) {
            value <- progress$getValue()
            value <- value + (progress$getMax() - value) / 5
          }
          progress$set(value = value, detail = detail)
        }  
        
        
        #draws rmse matrix
        matrix_rmse = create_rmse(aln, input$step, input$window, "pdist", "pairwise", updateProgress)
        

        # updating progress
        #progress$inc(0.5)
        
        heatmap.2(as.matrix(matrix_rmse[nrow(matrix_rmse):1,]), Rowv = FALSE, Colv = "Rowv", dendrogram = 'none', col=matlab.like, tracecol=NA)
        
        incProgress(0.1, detail = "heatmap finished")
      
      }
      )
      #width = 300
      #height = 900
    })
    
    #np <- nearPoints(matrix_rmse, input$plot_click)
    
    #output$info <- renderText({
    #  paste0("x=", input$plot_click$row, "\ny=", input$plot_click$col)
    #})
    
    
    
    
  }
  )
  
  
  observeEvent(
    input$goButton2, {

      
      withProgress(message = 'Creating distance plots', value = 0, {
      #progress1 <- shiny::Progress$new()
      
      
      #progress1$set(message = "Making control distance plot", value = 0)
      
      aln = read.dna(as.character(input$file_alignment$datapath), format="fasta", as.character=TRUE)
      aln[aln=='-'] <- NA
      output$control <- renderPlot({
                                  # Create a Progress object
                                  #progress1 <- shiny::Progress$new()
                                  
                                  # Make sure it closes when we exit this reactive, even if there's an error
                                  #on.exit(progress1$close())
                                  #progress1$set(message = "Making control distance plot", value = 0)
                                  
                                  plot_control(aln)
                                  }
                                  
                                  
                        )
      incProgress(0.4, detail = "Control plot finished")
      Sys.sleep(0.5)
      #progress2 <- shiny::Progress$new()
      
      # Make sure it closes when we exit this reactive, even if there's an error
      #on.exit(progress2$close())
      #progress2$set(message = "Making distance plot", value = 0)
      
      #progress1$set(message = "Making distance plot", value = 50)
      
      l = plot_dist_test(aln, input$start1,input$end1,input$start2,input$end2)
      
      incProgress(0.4, detail = "Distance plot finished")
      Sys.sleep(0.5)
      
      
      
      output$dist_plot <- renderPlot(l[[1]])
      df = l[[2]]
      
      incProgress(0.2, detail = "Creating plots")
      #Sys.sleep(1)
      }
      )
      
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