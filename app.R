library(shiny)
library(DT)
library(shinyjs)
source("scripts/rec_plots.R")

options(shiny.maxRequestSize = 40*1024^2)





# Define UI for application that draws a histogram
ui <- navbarPage(
        "Exploring recombination in genetic sequences using distance plots",
         tabPanel("RMSE matrix plot",
                  sidebarLayout(
                    sidebarPanel(
                      
                      fileInput(
                        "file_alignment", "Upload alignment in fasta-format", accept = c(".fasta", ".fas")
                      ),
                      radioButtons("default_alignment", "Or reproduce the results for Coronavirus genera",
                                   choices = c("Alphacoronavirus" = "alphac",
                                               "Betacoronavirus" = "betac",
                                               "Gammacoronavirus" = "gammac",
                                               "Deltacoronavirus" = "deltac"),
                                   selected = "alphac"),
                      
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
        tabPanel("Figures from paper",
                 
                 tags$head(tags$script(HTML("$(document).on('click', '.fig_button', function (e) {
                                e.stopPropagation()
                                if (typeof BUTTON_CLICK_COUNT == 'undefined') {
                                  BUTTON_CLICK_COUNT = 1; 
                                } else {
                                  BUTTON_CLICK_COUNT ++;
                                }
                                Shiny.onInputChange('js.button_clicked', 
                                  e.target.id + '_' + BUTTON_CLICK_COUNT);
                             });"))),
                 sidebarLayout(
                   sidebarPanel(
                     h4("Figure 3. Correspondence of pairwise nucleotide distances (PDC plots) between ORF1a and ORF1b (a), ORF1ab and spike (b), S1 and S2 regions of spike (c)."),
                     fluidRow(
                       h5("a"),
                       actionButton(inputId = "A-1ab-button", label = NULL, class="fig_button", style = "width: 150px; height: 120px;
                            background: url('AlphaORF1avsb.png');  background-size: cover; background-position: center;"),
                       actionButton(inputId = "B-1ab-button", label = NULL, class="fig_button", style = "width: 150px; height: 120px;
                            background: url('betaORF1avsb.png');  background-size: cover; background-position: center;"),
                       actionButton(inputId = "G-1ab-button", label = NULL, class="fig_button", style = "width: 150px; height: 120px;
                            background: url('gammaORF1avsb.png');  background-size: cover; background-position: center;"),
                       actionButton(inputId = "D-1ab-button", label = NULL, class="fig_button", style = "width: 150px; height: 120px;
                            background: url('deltaORF1avsb.png');  background-size: cover; background-position: center;")
                       
                       ),
                     fluidRow(
                       h5("b"),
                       actionButton(inputId = "A-1abs-button", label = NULL, class="fig_button", style = "width: 150px; height: 120px;
                            background: url('Alpha_ORF1abvsS.png');  background-size: cover; background-position: center;"),
                       actionButton(inputId = "B-1abs-button", label = NULL, class="fig_button", style = "width: 150px; height: 120px;
                            background: url('beta_ORF1abvsS.png');  background-size: cover; background-position: center;"),
                       actionButton(inputId = "G-1abs-button", label = NULL, class="fig_button", style = "width: 150px; height: 120px;
                            background: url('gamma_ORF1abvsS.png');  background-size: cover; background-position: center;"),
                       actionButton(inputId = "D-1abs-button", label = NULL, class="fig_button", style = "width: 150px; height: 120px;
                            background: url('delta_ORF1abvsS.png');  background-size: cover; background-position: center;")
                       
                     ),
                     
                     fluidRow(
                       h5("c"),
                       actionButton(inputId = "A-ss-button", label = NULL, class="fig_button", style = "width: 150px; height: 120px;
                            background: url('AlphaS1vsS2.png');  background-size: cover; background-position: center;"),
                       actionButton(inputId = "B-ss-button", label = NULL, class="fig_button", style = "width: 150px; height: 120px;
                            background: url('betaS1vsS2.png');  background-size: cover; background-position: center;"),
                       actionButton(inputId = "G-ss-button", label = NULL, class="fig_button", style = "width: 150px; height: 120px;
                            background: url('gammaS1vsS2.png');  background-size: cover; background-position: center;"),
                       actionButton(inputId = "D-ss-button", label = NULL, class="fig_button", style = "width: 150px; height: 120px;
                            background: url('deltaS1vsS2.png');  background-size: cover; background-position: center;")
                       
                     ),

                     
                   ),
                   mainPanel(
                     width = 6,
                     fluidRow(
                       
                       splitLayout(
                         cellWidths = 450, 
                         tagList(h4("Control plot"), plotOutput("fig3_control",
                         )
                         ),

                         tagList(h4("Distance plot"), plotOutput("fig3_plot",
                                                                 brush = brushOpts(id = "fig3_brush"),
                                                                 click = clickOpts(id = "fig3_click")
                         )
                         )
                       )
                     ),
                     fluidRow(
                       #verbatimTextOutput("mm"),
                       #DT::dataTableOutput("min_max"),
                       DT::dataTableOutput("fig3_brush_info")
                     )
                   )
                   
                   
                   
                 ),
                 
           sidebarLayout(

             sidebarPanel(
               h4("Figure 4. Correspondence of pairwise nucleotide distances in different regions of spike protein: S1-NTD and S1-CTD domains (a), two halves of  S1-NTD domain (b),two halves of  S1-CTD domain (c), S1-NTD domain and S2 region (d), S1-CTD and S2 region (e), two halves of S2 region (f)."),
               fluidRow(
                 h5("a"),
                 actionButton(inputId = "A-nc-button", label = NULL, class="fig_button", style = "width: 150px; height: 120px;
                            background: url('Alpha_NTDvsCTD.png');  background-size: cover; background-position: center;"),
                 actionButton(inputId = "B-nc-button", label = NULL, class="fig_button", style = "width: 150px; height: 120px;
                            background: url('beta_NTDvsCTD.png');  background-size: cover; background-position: center;"),
                 actionButton(inputId = "G-nc-button", label = NULL, class="fig_button", style = "width: 150px; height: 120px;
                            background: url('gamma_NTDvsCTD.png');  background-size: cover; background-position: center;"),
                 actionButton(inputId = "D-nc-button", label = NULL, class="fig_button", style = "width: 150px; height: 120px;
                            background: url('delta_NTDvsCTD.png');  background-size: cover; background-position: center;")
                 
               ),
               fluidRow(
                 h5("b"),
                 actionButton(inputId = "A-ntd-button", label = NULL, class="fig_button", style = "width: 150px; height: 120px;
                            background: url('Alpha_NTD.png');  background-size: cover; background-position: center;"),
                 actionButton(inputId = "B-ntd-button", label = NULL, class="fig_button", style = "width: 150px; height: 120px;
                            background: url('beta_NTD.png');  background-size: cover; background-position: center;"),
                 actionButton(inputId = "G-ntd-button", label = NULL, class="fig_button", style = "width: 150px; height: 120px;
                            background: url('gamma_NTD.png');  background-size: cover; background-position: center;"),
                 actionButton(inputId = "D-ntd-button", label = NULL, class="fig_button", style = "width: 150px; height: 120px;
                            background: url('delta_NTD.png');  background-size: cover; background-position: center;")
                 
               ),
               
               fluidRow(
                 h5("c"),
                 actionButton(inputId = "A-ctd-button", label = NULL, class="fig_button", style = "width: 150px; height: 120px;
                            background: url('Alpha_CTD.png');  background-size: cover; background-position: center;"),
                 actionButton(inputId = "B-ctd-button", label = NULL, class="fig_button", style = "width: 150px; height: 120px;
                            background: url('beta_CTD.png');  background-size: cover; background-position: center;"),
                 actionButton(inputId = "G-ctd-button", label = NULL, class="fig_button", style = "width: 150px; height: 120px;
                            background: url('gamma_CTD.png');  background-size: cover; background-position: center;"),
                 actionButton(inputId = "D-ctd-button", label = NULL, class="fig_button", style = "width: 150px; height: 120px;
                            background: url('delta_CTD.png');  background-size: cover; background-position: center;")
                 
               ),
               fluidRow(
                 h5("d"),
                 actionButton(inputId = "A-ns2-button", label = NULL, class="fig_button", style = "width: 150px; height: 120px;
                            background: url('Alpha_NTDvsS2.png');  background-size: cover; background-position: center;"),
                 actionButton(inputId = "B-ns2-button", label = NULL, class="fig_button", style = "width: 150px; height: 120px;
                            background: url('beta_NTDvsS2.png');  background-size: cover; background-position: center;"),
                 actionButton(inputId = "G-ns2-button", label = NULL, class="fig_button", style = "width: 150px; height: 120px;
                            background: url('gamma_NTDvsS2.png');  background-size: cover; background-position: center;"),
                 actionButton(inputId = "D-ns2-button", label = NULL, class="fig_button", style = "width: 150px; height: 120px;
                            background: url('delta_NTDvsS2.png');  background-size: cover; background-position: center;")
                 
               ),
               fluidRow(
                 h5("e"),
                 actionButton(inputId = "A-cs2-button", label = NULL, class="fig_button", style = "width: 150px; height: 120px;
                            background: url('Alpha_CTDvsS2.png');  background-size: cover; background-position: center;"),
                 actionButton(inputId = "B-cs2-button", label = NULL, class="fig_button", style = "width: 150px; height: 120px;
                            background: url('beta_CTDvsS2.png');  background-size: cover; background-position: center;"),
                 actionButton(inputId = "G-cs2-button", label = NULL, class="fig_button", style = "width: 150px; height: 120px;
                            background: url('gamma_CTDvsS2.png');  background-size: cover; background-position: center;"),
                 actionButton(inputId = "D-cs2-button", label = NULL, class="fig_button", style = "width: 150px; height: 120px;
                            background: url('delta_CTDvsS2.png');  background-size: cover; background-position: center;")
                 
               ),
               fluidRow(
                 h5("f"),
                 actionButton(inputId = "A-s2-button", label = NULL, class="fig_button", style = "width: 150px; height: 120px;
                            background: url('Alpha_S2.png');  background-size: cover; background-position: center;"),
                 actionButton(inputId = "B-s2-button", label = NULL, class="fig_button", style = "width: 150px; height: 120px;
                            background: url('beta_S2.png');  background-size: cover; background-position: center;"),
                 actionButton(inputId = "G-s2-button", label = NULL, class="fig_button", style = "width: 150px; height: 120px;
                            background: url('gamma_S2.png');  background-size: cover; background-position: center;"),
                 actionButton(inputId = "D-s2-button", label = NULL, class="fig_button", style = "width: 150px; height: 120px;
                            background: url('delta_S2.png');  background-size: cover; background-position: center;")
                 
               ),
             ),
             mainPanel(
               width = 6,
               fluidRow(
                 splitLayout(
                     #cellWidths = 450, 
                     tagList(h4("Control plot"), plotOutput("fig4_control",
                     )
                     ),
                     
                     tagList(h4("Distance plot"), plotOutput("fig4_plot",
                                                             brush = brushOpts(id = "fig4_brush"),
                                                             click = clickOpts(id = "fig4_click")
                     )
                     )
                 )  
                 
               ),
               fluidRow(
                 #verbatimTextOutput("mm"),
                 #DT::dataTableOutput("min_max"),
                 DT::dataTableOutput("fig4_brush_info")
               )
                 
               ),
               
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
    for (j in (i):(n)){
    #for (j in 1:(n - i + 1)){
      #for (j in 1:(n)){
      #print(paste(toString(i), toString(j), sep=","))
      #fits pairwise distance comparison plots linear model, calculates rmse
      rmse_i_j = (rmse(lm(dist_matrices[[j]]~dist_matrices[[i]])) + rmse(lm(dist_matrices[[i]]~dist_matrices[[j]]))) /2.0
      #rmse_i_j = rmse(lm(dist_matrices[[j]]~dist_matrices[[i]]))
      rmse_df[i,j] = rmse_i_j
      rmse_df[j,i] = rmse_i_j
      #rmse_df[n-j+1,n-i+1] = rmse_i_j
      
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


    if (! is.null(input$file_alignment$datapath)){
      file = input$file_alignment$datapath
    } else if (input$default_alignment == "alphac"){
      
      file = 'data/alpha.fasta'
    } else if (input$default_alignment == "betac") {
      file = 'data/beta.fasta'
    } else if (input$default_alignment == "gammac") {
      file = 'data/gamma.fasta'
    } else if (input$default_alignment == "deltac") {
      file = 'data/delta.fasta'
    }

    print(file)
    aln <- read.dna(as.character(file), format="fasta", as.character=TRUE)
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

      if (! is.null(input$file_alignment$datapath)){
        file = input$file_alignment$datapath
      } else if (input$default_alignment == "alphac"){
        file = 'data/alpha.fasta'
      } else if (input$default_alignment == "betac") {
        file = 'data/beta.fasta'
      } else if (input$default_alignment == "gammac") {
        file = 'data/gamma.fasta'
      } else if (input$default_alignment == "deltac") {
        file = 'data/delta.fasta'
      }
      
      print(file)
      withProgress(message = 'Creating distance plots', value = 0, {
      #progress1 <- shiny::Progress$new()
      
      
      #progress1$set(message = "Making control distance plot", value = 0)
      
      aln = read.dna(as.character(file), format="fasta", as.character=TRUE)
      aln[aln=='-'] <- NA
      
      aln_control = cbind(aln[,input$start1:input$end1], aln[,input$start2:input$end2])
      output$control <- renderPlot({
                                  # Create a Progress object
                                  #progress1 <- shiny::Progress$new()
                                  
                                  # Make sure it closes when we exit this reactive, even if there's an error
                                  #on.exit(progress1$close())
                                  #progress1$set(message = "Making control distance plot", value = 0)
                                  
                                  plot_control(aln_control)
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
          #find_recomb_names(l[[4]], min_2, max_2, l[[3]], min_1, max_1)
          
        })
        
      })
      
      
      
    }
  )
  

  
  observeEvent(input$js.button_clicked, {
    
    uid = strsplit(input$js.button_clicked, "_")
    #print(uid)
    button = uid[[1]][1]
    chs = strsplit(button, "-")
    print(chs)
    genus = chs[[1]][1]
    region = chs[[1]][2]
    print(genus)
    if (genus == "A"){
      aln_fig = read.dna(as.character("data/alpha.fasta"), format="fasta", as.character=TRUE)
      if (region == "1ab"){
        f = 3
        min1 = 1
        max1 = 11242
        min2 = 11243
        max2 = 19255
      }
      else if (region == "1abs"){
        f = 3
        min1 = 1
        max1 = 19255
        min2 = 19256
        max2 = 22710
      }
      else if (region == "ss"){
        f = 3
        min1 = 19255
        max1 = 20871
        min2 = 20872
        max2 = 22710
      }
      else if (region == "nc"){
        f = 4
        min1 = 19612
        max1 = 20190
        min2 = 20272
        max2 = 22710
      }
      else if (region == "ntd"){
        f = 4
        min1 = 19612
        max1 = 19901
        min2 = 19902
        max2 = 20190
      }
      else if (region == "ctd"){
        f = 4
        min1 = 20272
        max1 = 20426
        min2 = 20427
        max2 = 20580
      }
      else if (region == "ns2"){
        f = 4
        min1 = 19612
        max1 = 20190
        min2 = 20872
        max2 = 22710
      }
      else if (region == "cs2"){
        f = 4
        min1 = 20272
        max1 = 20580
        min2 = 20872
        max2 = 22710
      }
      else if (region == "s2"){
        f = 4
        min1 = 20872
        max1 = 21800
        min2 = 21800
        max2 = 22710
      }
    }
    else if (genus == "B"){
      aln_fig = read.dna(as.character("data/beta.fasta"), format="fasta", as.character=TRUE)
      if (region == "1ab"){
        f = 3
        min1 = 1
        max1 = 11484
        min2 = 11485
        max2 = 19398
      }
      else if (region == "1abs"){
        f = 3
        min1 = 1
        max1 = 19398
        min2 = 19399
        max2 = 22746
      }
      else if (region == "ss"){
        f = 3
        min1 = 19468
        max1 = 21060
        min2 = 21061
        max2 = 22743
      }
      else if (region == "nc"){
        f = 4
        min1 = 19429
        max1 = 20124
        min2 = 20212
        max2 = 20646
      }
      else if (region == "ntd"){
        f = 4
        min1 = 19429
        max1 = 19750
        min2 = 19751
        max2 = 20124
      }
      else if (region == "ctd"){
        f = 4
        min1 = 20212
        max1 = 20400
        min2 = 20401
        max2 = 20646
      }
      else if (region == "ns2"){
        f = 4
        min1 = 19429
        max1 = 20124
        min2 = 21061
        max2 = 22743
      }
      else if (region == "cs2"){
        f = 4
        min1 = 20212
        max1 = 20646
        min2 = 21061
        max2 = 22743
      }
      else if (region == "s2"){
        f = 4
        min1 = 21061
        max1 = 21900
        min2 = 21901
        max2 = 22743
      }
    }
    else if (genus == "G"){
      aln_fig = read.dna(as.character("data/gamma.fasta"), format="fasta", as.character=TRUE)
      if (region == "1ab"){
        f = 3
        min1 = 1
        max1 = 11754
        min2 = 11755
        max2 =19629
      }
      else if (region == "1abs"){
        f = 3
        min1 = 1
        max1 = 19629
        min2 = 19630
        max2 = 22689
      }
      else if (region == "ss"){
        f = 3
        min1 = 19714
        max1 = 20862
        min2 = 20863
        max2 = 22689
      }
      else if (region == "nc"){
        f = 4
        min1 = 19684
        max1 = 20142
        min2 = 20182
        max2 = 20553
      }
      else if (region == "ntd"){
        f = 4
        min1 = 19684
        max1 = 19900
        min2 = 19901
        max2 = 20142
      }
      else if (region == "ctd"){
        f = 4
        min1 = 20182
        max1 = 20350
        min2 = 20351
        max2 = 20553
      }
      else if (region == "ns2"){
        f = 4
        min1 = 19684
        max1 = 20142
        min2 = 20863
        max2 = 22689
      }
      else if (region == "cs2"){
        f = 4
        min1 = 20182
        max1 = 20553
        min2 = 20863
        max2 = 22689
      }
      else if (region == "s2"){
        f = 4
        min1 = 20863
        max1 = 21800
        min2 = 21800
        max2 = 22689
      }
    }
    else if (genus == "D"){
      aln_fig = read.dna(as.character("data/delta.fasta"), format="fasta", as.character=TRUE)
      if (region == "1ab"){
        f = 3
        min1 = 1
        max1 = 10608
        min2 = 10609
        max2 =18531
      }
      else if (region == "1abs"){
        f = 3
        min1 = 1
        max1 = 18531
        min2 = 18532
        max2 = 21915
      }
      else if (region == "ss"){
        f = 3
        min1 = 18697
        max1 = 20007
        min2 = 20029
        max2 = 21876
      }
      else if (region == "nc"){
        f = 4
        min1 = 18697
        max1 = 19338
        min2 = 19408
        max2 = 19791
      }
      else if (region == "ntd"){
        f = 4
        min1 = 18697
        max1 = 19400
        min2 = 19401
        max2 = 19791
      }
      else if (region == "ctd"){
        f = 4
        min1 = 19408
        max1 = 19600
        min2 = 19601
        max2 = 19791
      }
      else if (region == "ns2"){
        f = 4
        min1 = 18697
        max1 = 19791
        min2 = 20029
        max2 = 21876
      }
      else if (region == "cs2"){
        f = 4
        min1 = 19408
        max1 = 19791
        min2 = 20029
        max2 = 21876
      }
      else if (region == "s2"){
        f = 4
        min1 = 20029
        max1 = 20900
        min2 = 20901
        max2 = 21876
      }
    }
    aln_fig[aln_fig=='-'] <- NA

    #l_fig = plot_dist_test(aln_fig, min1,max1,min2,max2)
    #df = l_fig[[2]]
    
    if (f==3){
      withProgress(message = 'Creating distance plots', value = 0, {
        
        aln_control =  cbind(aln_fig[,min1:max1], aln_fig[,min2:max2])
        output$fig3_control <- renderPlot({
          plot_control(aln_control)
        })
        incProgress(0.4, detail = "Control plot finished")
        Sys.sleep(0.5)
        
        l_fig = plot_dist_test(aln_fig, min1,max1,min2,max2)
        df = l_fig[[2]]
        incProgress(0.4, detail = "Distance plot finished")
        Sys.sleep(0.5)
        
        output$fig3_plot <- renderPlot(l_fig[[1]])
        
        
        observeEvent(input$fig3_brush, {
          brushed_points3 <- brushedPoints(df, input$fig3_brush)
          min_1 = input$fig3_brush$xmin
          max_1 = input$fig3_brush$xmax
          
          min_2 = input$fig3_brush$ymin
          max_2 = input$fig3_brush$ymax
          output$min_max <- DT::renderDataTable(brushed_points3)
          
          output$fig3_brush_info <- DT::renderDataTable({
            find_recomb_names(l_fig[[3]], min_1, max_1, l_fig[[4]], min_2, max_2)
          })
        })
        incProgress(0.2, detail = "Creating plots")
      })
      
    }
    if (f==4){
      withProgress(message = 'Creating distance plots', value = 0, {
      
        aln_control =  cbind(aln_fig[,min1:max1], aln_fig[,min2:max2])
        output$fig4_control <- renderPlot({
          plot_control(aln_control)
        })
        incProgress(0.4, detail = "Control plot finished")
        Sys.sleep(0.5)
        
        l_fig = plot_dist_test(aln_fig, min1,max1,min2,max2)
        df = l_fig[[2]]
        incProgress(0.4, detail = "Distance plot finished")
        Sys.sleep(0.5)
        
        output$fig4_plot <- renderPlot(l_fig[[1]])
        observeEvent(input$fig4_brush, {
          brushed_points4 <- brushedPoints(df, input$fig4_brush)
          min_1 = input$fig4_brush$xmin
          max_1 = input$fig4_brush$xmax
          
          min_2 = input$fig4_brush$ymin
          max_2 = input$fig4_brush$ymax
          output$min_max <- DT::renderDataTable(brushed_points4)
          
          output$fig4_brush_info <- DT::renderDataTable({
            find_recomb_names(l_fig[[3]], min_1, max_1, l_fig[[4]], min_2, max_2)
          })
        })
        incProgress(0.2, detail = "Creating plots")
      })
    }
  
  }
  )
  
}

# Run the application 
shinyApp(ui = ui, server = server)