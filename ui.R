# Define UI for application
ui <- navbarPage(
        "Exploring recombination in genetic sequences using distance plots",
         theme = "styles.css",
#------------------RMSE PLOT TAB
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

                    )
                  )
         ),
#------------------DISTANCE PLOT TAB
         tabPanel("Distance plot",
                  sidebarLayout(
                    sidebarPanel(
                      h4("Upload alignment in fasta-format"),
                      fileInput(
                        "file_alignment2", "", accept = c(".fasta", ".fas")
                      ),
                      h4("Or reproduce the results for Coronavirus genera"),
                      radioButtons("default_alignment2", "",
                                   choices = c("Alphacoronavirus" = "alphac",
                                               "Betacoronavirus" = "betac",
                                               "Gammacoronavirus" = "gammac",
                                               "Deltacoronavirus" = "deltac"),
                                   selected = "alphac"),
                      
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
#------------------FIGURES FROM PAPER TAB
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
                     h4("Click the figures below to reproduce the plots. Select the dots on the distance plot to show the virus pairs that correspond to them."),
                     h5("Figure 3. Correspondence of pairwise nucleotide distances (PDC plots) between ORF1a and ORF1b (a), ORF1ab and spike (b), S1 and S2 regions of spike (c)."),
                     fluidRow(
                       h5("a"),
                       actionButton(inputId = "A-1ab-button", label = NULL, class="fig_button", style = "width: 140px; height: 110px;
                            background: url('img/AlphaORF1avsb.png');  background-size: cover; background-position: center;"),
                       actionButton(inputId = "B-1ab-button", label = NULL, class="fig_button", style = "width: 140px; height: 110px;
                            background: url('img/BetaORF1avsb.png');  background-size: cover; background-position: center;"),
                       actionButton(inputId = "G-1ab-button", label = NULL, class="fig_button", style = "width: 140px; height: 110px;
                            background: url('img/GammaORF1avsb.png');  background-size: cover; background-position: center;"),
                       actionButton(inputId = "D-1ab-button", label = NULL, class="fig_button", style = "width: 140px; height: 110px;
                            background: url('img/DeltaORF1avsb.png');  background-size: cover; background-position: center;")
                       
                       ),
                     fluidRow(
                       h5("b"),
                       actionButton(inputId = "A-1abs-button", label = NULL, class="fig_button", style = "width: 140px; height: 110px;
                            background: url('img/Alpha_ORF1abvsS.png');  background-size: cover; background-position: center;"),
                       actionButton(inputId = "B-1abs-button", label = NULL, class="fig_button", style = "width: 140px; height: 110px;
                            background: url('img/Beta_ORF1abvsS.png');  background-size: cover; background-position: center;"),
                       actionButton(inputId = "G-1abs-button", label = NULL, class="fig_button", style = "width: 140px; height: 110px;
                            background: url('img/Gamma_ORF1abvsS.png');  background-size: cover; background-position: center;"),
                       actionButton(inputId = "D-1abs-button", label = NULL, class="fig_button", style = "width: 140px; height: 110px;
                            background: url('img/delta_ORF1abvsS.png');  background-size: cover; background-position: center;")
                       
                     ),
                     
                     fluidRow(
                       h5("c"),
                       actionButton(inputId = "A-ss-button", label = NULL, class="fig_button", style = "width: 140px; height: 110px;
                            background: url('img/AlphaS1vsS2.png');  background-size: cover; background-position: center;"),
                       actionButton(inputId = "B-ss-button", label = NULL, class="fig_button", style = "width: 140px; height: 110px;
                            background: url('img/BetaS1vsS2.png');  background-size: cover; background-position: center;"),
                       actionButton(inputId = "G-ss-button", label = NULL, class="fig_button", style = "width: 140px; height: 110px;
                            background: url('img/GammaS1vsS2.png');  background-size: cover; background-position: center;"),
                       actionButton(inputId = "D-ss-button", label = NULL, class="fig_button", style = "width: 140px; height: 110px;
                            background: url('img/deltaS1vsS2.png');  background-size: cover; background-position: center;")
                       
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
               h5("Figure 4. Correspondence of pairwise nucleotide distances in different regions of spike protein: S1-NTD and S1-CTD domains (a), two halves of  S1-NTD domain (b),two halves of  S1-CTD domain (c), S1-NTD domain and S2 region (d), S1-CTD and S2 region (e), two halves of S2 region (f)."),
               fluidRow(
                 h5("a"),
                 actionButton(inputId = "A-nc-button", label = NULL, class="fig_button", style = "width: 140px; height: 110px;
                            background: url('img/Alpha_NTDvsCTD.png');  background-size: cover; background-position: center;"),
                 actionButton(inputId = "B-nc-button", label = NULL, class="fig_button", style = "width: 140px; height: 110px;
                            background: url('img/beta_NTDvsCTD.png');  background-size: cover; background-position: center;"),
                 actionButton(inputId = "G-nc-button", label = NULL, class="fig_button", style = "width: 140px; height: 110px;
                            background: url('img/gamma_NTDvsCTD.png');  background-size: cover; background-position: center;"),
                 actionButton(inputId = "D-nc-button", label = NULL, class="fig_button", style = "width: 140px; height: 110px;
                            background: url('img/delta_NTDvsCTD.png');  background-size: cover; background-position: center;")
                 
               ),
               fluidRow(
                 h5("b"),
                 actionButton(inputId = "A-ntd-button", label = NULL, class="fig_button", style = "width: 140px; height: 110px;
                            background: url('img/Alpha_NTD.png');  background-size: cover; background-position: center;"),
                 actionButton(inputId = "B-ntd-button", label = NULL, class="fig_button", style = "width: 140px; height: 110px;
                            background: url('img/beta_NTD.png');  background-size: cover; background-position: center;"),
                 actionButton(inputId = "G-ntd-button", label = NULL, class="fig_button", style = "width: 140px; height: 110px;
                            background: url('img/gamma_NTD.png');  background-size: cover; background-position: center;"),
                 actionButton(inputId = "D-ntd-button", label = NULL, class="fig_button", style = "width: 140px; height: 110px;
                            background: url('img/delta_NTD.png');  background-size: cover; background-position: center;")
                 
               ),
               
               fluidRow(
                 h5("c"),
                 actionButton(inputId = "A-ctd-button", label = NULL, class="fig_button", style = "width: 140px; height: 110px;
                            background: url('img/Alpha_CTD.png');  background-size: cover; background-position: center;"),
                 actionButton(inputId = "B-ctd-button", label = NULL, class="fig_button", style = "width: 140px; height: 110px;
                            background: url('img/beta_CTD.png');  background-size: cover; background-position: center;"),
                 actionButton(inputId = "G-ctd-button", label = NULL, class="fig_button", style = "width: 140px; height: 110px;
                            background: url('img/gamma_CTD.png');  background-size: cover; background-position: center;"),
                 actionButton(inputId = "D-ctd-button", label = NULL, class="fig_button", style = "width: 140px; height: 110px;
                            background: url('img/delta_CTD.png');  background-size: cover; background-position: center;")
                 
               ),
               fluidRow(
                 h5("d"),
                 actionButton(inputId = "A-ns2-button", label = NULL, class="fig_button", style = "width: 140px; height: 110px;
                            background: url('img/Alpha_NTDvsS2.png');  background-size: cover; background-position: center;"),
                 actionButton(inputId = "B-ns2-button", label = NULL, class="fig_button", style = "width: 140px; height: 110px;
                            background: url('img/beta_NTDvsS2.png');  background-size: cover; background-position: center;"),
                 actionButton(inputId = "G-ns2-button", label = NULL, class="fig_button", style = "width: 140px; height: 110px;
                            background: url('img/gamma_NTDvsS2.png');  background-size: cover; background-position: center;"),
                 actionButton(inputId = "D-ns2-button", label = NULL, class="fig_button", style = "width: 140px; height: 110px;
                            background: url('img/delta_NTDvsS2.png');  background-size: cover; background-position: center;")
                 
               ),
               fluidRow(
                 h5("e"),
                 actionButton(inputId = "A-cs2-button", label = NULL, class="fig_button", style = "width: 140px; height: 110px;
                            background: url('img/Alpha_CTDvsS2.png');  background-size: cover; background-position: center;"),
                 actionButton(inputId = "B-cs2-button", label = NULL, class="fig_button", style = "width: 140px; height: 110px;
                            background: url('img/beta_CTDvsS2.png');  background-size: cover; background-position: center;"),
                 actionButton(inputId = "G-cs2-button", label = NULL, class="fig_button", style = "width: 140px; height: 110px;
                            background: url('img/gamma_CTDvsS2.png');  background-size: cover; background-position: center;"),
                 actionButton(inputId = "D-cs2-button", label = NULL, class="fig_button", style = "width: 140px; height: 110px;
                            background: url('img/delta_CTDvsS2.png');  background-size: cover; background-position: center;")
                 
               ),
               fluidRow(
                 h5("f"),
                 actionButton(inputId = "A-s2-button", label = NULL, class="fig_button", style = "width: 140px; height: 110px;
                            background: url('img/Alpha_S2.png');  background-size: cover; background-position: center;"),
                 actionButton(inputId = "B-s2-button", label = NULL, class="fig_button", style = "width: 140px; height: 110px;
                            background: url('img/beta_S2.png');  background-size: cover; background-position: center;"),
                 actionButton(inputId = "G-s2-button", label = NULL, class="fig_button", style = "width: 140px; height: 110px;
                            background: url('img/gamma_S2.png');  background-size: cover; background-position: center;"),
                 actionButton(inputId = "D-s2-button", label = NULL, class="fig_button", style = "width: 140px; height: 110px;
                            background: url('img/delta_S2.png');  background-size: cover; background-position: center;")
                 
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
               
             ),
         sidebarLayout(
           sidebarPanel(
             h5("Figure S5. Correspondence of pairwise nucleotide distances (PDC plots) between ORF1ab and domains of spike protein: S1 (a), S2 (b), S1-NTD (c), S1-CTD (d)."),
             fluidRow(
               h5("a"),
               actionButton(inputId = "A-1abS1-button", label = NULL, class="fig_button", style = "width: 140px; height: 110px;
                          background: url('img/AlphaORF1abvsS1.png');  background-size: cover; background-position: center;"),
               actionButton(inputId = "B-1abS1-button", label = NULL, class="fig_button", style = "width: 140px; height: 110px;
                          background: url('img/BetaORF1abvsS1.png');  background-size: cover; background-position: center;"),
               actionButton(inputId = "G-1abS1-button", label = NULL, class="fig_button", style = "width: 140px; height: 110px;
                          background: url('img/GammaORF1abvsS1.png');  background-size: cover; background-position: center;"),
               actionButton(inputId = "D-1abS1-button", label = NULL, class="fig_button", style = "width: 140px; height: 110px;
                          background: url('img/delta_ORF1abvsS1.png');  background-size: cover; background-position: center;")
               
             ),
             fluidRow(
               h5("b"),
               actionButton(inputId = "A-1abS2-button", label = NULL, class="fig_button", style = "width: 140px; height: 110px;
                          background: url('img/AlphaORF1abvsS2.png');  background-size: cover; background-position: center;"),
               actionButton(inputId = "B-1abS2-button", label = NULL, class="fig_button", style = "width: 140px; height: 110px;
                          background: url('img/BetaORF1abvsS2.png');  background-size: cover; background-position: center;"),
               actionButton(inputId = "G-1abS2-button", label = NULL, class="fig_button", style = "width: 140px; height: 110px;
                          background: url('img/GammaORF1abvsS2.png');  background-size: cover; background-position: center;"),
               actionButton(inputId = "D-1abS2-button", label = NULL, class="fig_button", style = "width: 140px; height: 110px;
                          background: url('img/delta_ORF1abvsS2.png');  background-size: cover; background-position: center;")
               
             ),
             
             fluidRow(
               h5("c"),
               actionButton(inputId = "A-1abntd-button", label = NULL, class="fig_button", style = "width: 140px; height: 110px;
                          background: url('img/AlphaORF1abvsNTD.png');  background-size: cover; background-position: center;"),
               actionButton(inputId = "B-1abntd-button", label = NULL, class="fig_button", style = "width: 140px; height: 110px;
                          background: url('img/BetaORF1abvsNTD.png');  background-size: cover; background-position: center;"),
               actionButton(inputId = "G-1abntd-button", label = NULL, class="fig_button", style = "width: 140px; height: 110px;
                          background: url('img/GammaORF1abvsNTD.png');  background-size: cover; background-position: center;"),
               actionButton(inputId = "D-1abntd-button", label = NULL, class="fig_button", style = "width: 140px; height: 110px;
                          background: url('img/delta_ORF1abvsNTD.png');  background-size: cover; background-position: center;")
               
             ),
             fluidRow(
               h5("c"),
               actionButton(inputId = "A-1abctd-button", label = NULL, class="fig_button", style = "width: 140px; height: 110px;
                          background: url('img/AlphaORF1abvsCTD.png');  background-size: cover; background-position: center;"),
               actionButton(inputId = "B-1abctd-button", label = NULL, class="fig_button", style = "width: 140px; height: 110px;
                          background: url('img/BetaORF1abvsCTD.png');  background-size: cover; background-position: center;"),
               actionButton(inputId = "G-1abctd-button", label = NULL, class="fig_button", style = "width: 140px; height: 110px;
                          background: url('img/GammaORF1abvsCTD.png');  background-size: cover; background-position: center;"),
               actionButton(inputId = "D-1abctd-button", label = NULL, class="fig_button", style = "width: 140px; height: 110px;
                          background: url('img/delta_ORF1abvsCTD.png');  background-size: cover; background-position: center;")
               
             ),
             
             
           ),
           mainPanel(
             width = 6,
             fluidRow(
               
               splitLayout(
                 cellWidths = 450, 
                 tagList(h4("Control plot"), plotOutput("fig5_control",
                 )
                 ),
                 
                 tagList(h4("Distance plot"), plotOutput("fig5_plot",
                                                         brush = brushOpts(id = "fig5_brush"),
                                                         click = clickOpts(id = "fig5_click")
                 )
                 )
               )
             ),
             fluidRow(
               #verbatimTextOutput("mm"),
               #DT::dataTableOutput("min_max"),
               DT::dataTableOutput("fig5_brush_info")
             )
           )
           
           
           
         )
           
             
           
                 
                 
        ),
#------------------ABOUT APP TAB
        tabPanel("About",
                   mainPanel(
                     p("The description of this app")
                   )
        )

)
