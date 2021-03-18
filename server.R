library(shiny)
library(DT)
library(shinyjs)
library(shinycssloaders)
library(shinythemes)
source("scripts/rec_plots.R")

options(shiny.maxRequestSize = 40*1024^2)


server <- function(input, output) {
# ---------RMSE PLOT GO BUTTON
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
      
    },
    width = 800,
    height = 500)
    
    #np <- nearPoints(matrix_rmse, input$plot_click)
    
    #output$info <- renderText({
    #  paste0("x=", input$plot_click$row, "\ny=", input$plot_click$col)
    #})
    
    
    
    
  }
  )
# ---------DISTANCE PLOT GO BUTTON  
  
  observeEvent(
    input$goButton2, {

      if (! is.null(input$file_alignment2$datapath)){
        file = input$file_alignment2$datapath
      } else if (input$default_alignment2 == "alphac"){
        file = 'data/alpha.fasta'
      } else if (input$default_alignment2 == "betac") {
        file = 'data/beta.fasta'
      } else if (input$default_alignment2 == "gammac") {
        file = 'data/gamma.fasta'
      } else if (input$default_alignment2 == "deltac") {
        file = 'data/delta.fasta'
      }
      
      print(file)
      withProgress(message = 'Creating distance plots', value = 0, {

      aln = read.dna(as.character(file), format="fasta", as.character=TRUE)
      aln[aln=='-'] <- NA
      
      aln_control = cbind(aln[,input$start1:input$end1], aln[,input$start2:input$end2])
      output$control <- renderPlot({
                                  
                                  plot_control(aln_control)
                                  }
                                  
                                  
                        )
      incProgress(0.4, detail = "Control plot finished")
      Sys.sleep(0.5)

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
          

          
        },extensions = 'Buttons',
        options = list(pageLength = 20, dom = 'Bfrtip', buttons = c('copy', 'csv'))
        )
        
      })
      
      
      
    }
  )
  
# GO BUTTONS FOR FIGURES FROM PAPER
  
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
        min1 = 19256
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
      else if (region == "1abS1"){
        f = 5
        min1 = 1
        max1 = 19255
        min2 = 19256
        max2 = 20871
      }
      else if (region == "1abS2"){
        f = 5
        min1 = 1
        max1 = 19255
        min2 = 20872
        max2 = 22710
      }
      else if (region == "1abntd"){
        f = 5
        min1 = 1
        max1 = 19255
        min2 = 19612
        max2 = 20190
      }
      else if (region == "1abctd"){
        f = 5
        min1 = 1
        max1 = 19255
        min2 = 20272
        max2 = 20580
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
      else if (region == "1abS1"){
        f = 5
        min1 = 1
        max1 = 19398
        min2 = 19468
        max2 = 21060
      }
      else if (region == "1abS2"){
        f = 5
        min1 = 1
        max1 = 19398
        min2 = 21061
        max2 = 22743
      }
      else if (region == "1abntd"){
        f = 5
        min1 = 1
        max1 = 19398
        min2 = 19429
        max2 = 20124
      }
      else if (region == "1abctd"){
        f = 5
        min1 = 1
        max1 = 19398
        min2 = 20212
        max2 = 20646
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
      else if (region == "1abS1"){
        f = 5
        min1 = 1
        max1 = 19629
        min2 = 19714
        max2 = 20862
      }
      else if (region == "1abS2"){
        f = 5
        min1 = 1
        max1 = 19629
        min2 = 20863
        max2 = 22689
      }
      else if (region == "1abntd"){
        f = 5
        min1 = 1
        max1 = 19629
        min2 = 19684
        max2 = 20142
      }
      else if (region == "1abctd"){
        f = 5
        min1 = 1
        max1 = 19629
        min2 = 20182
        max2 = 20553
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
        max1 = 19017
        min2 = 19018
        max2 = 19338
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
        max1 = 19338
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
      else if (region == "1abS1"){
        f = 5
        min1 = 1
        max1 = 18531
        min2 = 18697
        max2 = 20007
      }
      else if (region == "1abS2"){
        f = 5
        min1 = 1
        max1 = 18531
        min2 = 20029
        max2 = 21876
      }
      else if (region == "1abntd"){
        f = 5
        min1 = 1
        max1 = 18531
        min2 = 18697
        max2 = 19338
      }
      else if (region == "1abctd"){
        f = 5
        min1 = 1
        max1 = 18531
        min2 = 19408
        max2 = 19791
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
          }, options = list(pageLength = 20))
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
          }, options = list(pageLength = 20))
        })
        incProgress(0.2, detail = "Creating plots")
      })
    }
    if (f==5){
      withProgress(message = 'Creating distance plots', value = 0, {
        
        aln_control =  cbind(aln_fig[,min1:max1], aln_fig[,min2:max2])
        output$fig5_control <- renderPlot({
          plot_control(aln_control)
        })
        incProgress(0.4, detail = "Control plot finished")
        Sys.sleep(0.5)
        
        l_fig = plot_dist_test(aln_fig, min1,max1,min2,max2)
        df = l_fig[[2]]
        incProgress(0.4, detail = "Distance plot finished")
        Sys.sleep(0.5)
        
        output$fig5_plot <- renderPlot(l_fig[[1]])
        observeEvent(input$fig5_brush, {
          brushed_points5 <- brushedPoints(df, input$fig5_brush)
          min_1 = input$fig5_brush$xmin
          max_1 = input$fig5_brush$xmax
          
          min_2 = input$fig5_brush$ymin
          max_2 = input$fig5_brush$ymax
          output$min_max <- DT::renderDataTable(brushed_points5)
          
          output$fig5_brush_info <- DT::renderDataTable({
            find_recomb_names(l_fig[[3]], min_1, max_1, l_fig[[4]], min_2, max_2)
          }, options = list(pageLength = 20))
        })
        incProgress(0.2, detail = "Creating plots")
      })
    }
  }
  )
  
}
