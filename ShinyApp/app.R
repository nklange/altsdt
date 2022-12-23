source("R/helper.R")


SB2020_e1 <- readRDS("Data/Data_SB2020_e1.rds") %>%
  mutate(oldnew = isold)
SB2020_e2 <- readRDS("Data/Data_SB2020_e2.rds")%>%
  mutate(oldnew = isold)
SB2020_e3 <- readRDS("Data/Data_SB2020_e3.rds")%>%
  mutate(oldnew = isold)

SB2020_e123_fits <- readRDS("Fits/SB2020_e123_fits.rds")

SB2021_e1 <- readRDS("Data/Data_SB2022_e1.rds")
SB2021_e1_fits <- readRDS("Fits/SB2021_e12_fits.rds")
SB2021_e2 <- readRDS("Data/Data_SB2022_e2.rds")
SB2021_e3 <- readRDS("Data/Data_SB2022_e3.rds") %>%
  mutate(id = as.character(id)) %>%
  mutate(exp = "SB2021_e3")
SB2021_e3_fits <- bind_rows(readRDS("Fits/SB2021_e3_modelfits_baselineUVSDT.rds"),
                            readRDS("Fits/SB2021_e3_modelfits_ExGaussNormEVSDT.rds"),
                            readRDS("Fits/SB2021_e3_modelfits_GumbelEVSDT.rds"))%>%
  mutate(id = as.character(id)) %>%
  mutate(exp = "SB2021_e3")


MWW2007 <- readRDS("Data/MWW2007_sixpoint_preprocesseddata.rds") %>%
  mutate(rating = response) %>%
  mutate(oldnew = ifelse(oldnew == "Old",1,0)) %>%
  filter(rating %in% c(1:6))
MWW2007_fits <- readRDS("Fits/MWW2007_fits.rds") %>%
  filter(model != "GaussianUVSDT_equi") %>%
  group_by(id,model) %>%
  filter(objective == min(objective)) %>%
  ungroup()
Gsqs <- readRDS("Fits/Gsquared_allexps.rds") %>%
  mutate(model = ifelse(model == "baselineUVSDT","GaussianUVSDT",model))
Preds <- readRDS("Fits/Preds_allexps.rds")

ui <-
  navbarPage(title = "SDT PARAMETERIZATIONS",

             #useWaitress(),
             theme = "bootstrap.css",
             shinyjs::useShinyjs(),
             # tabPanel(title = "Aggregate",
             #          sidebarLayout(
             #
             #            sidebarPanel( width=2,
             #              checkboxGroupInput("selectmodel", h3("Models:"),
             #                                 choices = list("UVSDT" = "UVSDT",
             #                                   "EVSDT" = "EVSDT",
             #                                   "Gumbel" = "Gumbel"),
             #                                 selected = "UVSDT")
             #            ),
             #            mainPanel(width=10,
             #              fluidRow(
             #                div(style = "padding:20px",
             #                    h2("Deviance"),
             #                    htmlOutput("aggDev_info"),
             #                column(6,
             #                       plotOutput("aggDev_exps1", height = 400)
             #                       ),
             #                column(6,
             #                       plotOutput("aggDev_exps2", height = 400)
             #                )
             #                )
             #              ),
             #              fluidRow(
             #                div(style = "padding:20px",
             #                    h2("K-Fold Crossvalidation"),
             #                    htmlOutput("aggKFCV_info"),
             #                    column(6,
             #                           plotOutput("aggKFCV_exps1", height = 400)
             #                    ),
             #                    column(6,
             #                           plotOutput("aggKFCV_exps2", height = 400)
             #                    )
             #                )
             #              ),
             #              fluidRow(
             #                div(style = "padding:20px",
             #                    h2("Leave-Participant-Out Crossvalidation"),
             #                    htmlOutput("aggLOP_info"),
             #                    column(6,
             #                           plotOutput("aggLOP_exps1", height = 400)
             #                    ),
             #                    column(6,
             #                           plotOutput("aggLOP_exps2", height = 400)
             #                    )
             #                )
             #              )
             #            )
             #          )
             # ),

             tabPanel(title = "Individual Datasets",
                      sidebarLayout(
                        sidebarPanel(
                          width = 2,

                          selectInput("exp", h3("Experiments"),
                                       choices = list("SB2020_e1", "SB2020_e2",
                                                      "SB2020_e3",
                                                      "SB2021_e1", "SB2021_e2",
                                                      "SB2021_e3", "MWW2007_e1",
                                                      "MWW2007_e2","MWW2007_e3",
                                                      "D2007_e1a","D2007_e1b"),
                                       selected = "SB2020_e1"),
                          conditionalPanel(
                            condition = "input.exp == 'SB2020_e1'",
                            title = "Individual",
                            radioButtons("SB2020_e1", "Individual:",
                                         choices = as.list(unique(SB2020_e1$id)),selected="SB2020_e1_1")
                          ),
                          conditionalPanel(
                            condition = "input.exp == 'SB2020_e2'",
                            title = "Individual",
                            radioButtons("SB2020_e2", "Individual:",
                                         choices = as.list(unique(SB2021_e2$id)),selected="SB2020_e2_1")
                          ),
                          conditionalPanel(
                            condition = "input.exp == 'SB2020_e3'",
                            title = "Individual",
                            radioButtons("SB2020_e3", "Individual:",
                                         choices = as.list(unique(SB2020_e3$id)),selected="SB2020_e3_1")
                          ),
                          conditionalPanel(
                            condition = "input.exp == 'SB2021_e1'",
                            title = "Individual",
                            radioButtons("SB2021_e1", "Individual:",
                                        choices = as.list(unique(SB2021_e1$id)),selected="SB2021_e1_1")
                          ),
                          conditionalPanel(
                            condition = "input.exp == 'SB2021_e2'",
                            title = "Individual",
                            radioButtons("SB2021_e2", "Individual:",
                                         choices = as.list(unique(SB2021_e2$id)),selected="SB2021_e2_1")
                          ),
                          conditionalPanel(
                            condition = "input.exp == 'SB2021_e3'",
                            title = "Individual",
                            radioButtons("SB2021_e3", "Individual:",
                                         choices = as.list(unique(SB2021_e3$id)),selected="10597")
                          ),
                          conditionalPanel(
                            condition = "input.exp == 'MWW2007_e1'",
                            title = "Individual",
                            radioButtons("MWW2007_e1", "Individual:",
                                         choices = as.list(unique(MWW2007 %>%
                                                                    filter(exp == "MWW2007_e1") %>%
                                                                    .$id)),selected="MWW2007_e1_1")
                          ),
                          conditionalPanel(
                            condition = "input.exp == 'MWW2007_e2'",
                            title = "Individual",
                            radioButtons("MWW2007_e2", "Individual:",
                                         choices = as.list(unique(MWW2007 %>%
                                                                    filter(exp == "MWW2007_e2") %>%
                                                                    .$id)),selected="MWW2007_e2_15")
                          ),
                          conditionalPanel(
                            condition = "input.exp == 'MWW2007_e3'",
                            title = "Individual",
                            radioButtons("MWW2007_e3", "Individual:",
                                         choices = as.list(unique(MWW2007 %>%
                                                                    filter(exp == "MWW2007_e3") %>%
                                                                    .$id)),selected="MWW2007_e3_31")
                          )


                        ),
                        mainPanel(width=10,
                                  fluidRow(
                                    column(12,
                                           div(style = "padding:20px",
                                               h2("AIC"),
                                               div(style = "text-align:right",
                                                   a(id = "toggleDev", "Show explanation", href = "#")
                                               ),
                                               shinyjs::hidden(
                                                 div(id = "devianceinfo",
                                                     htmlOutput("devianceinfo")
                                                 )

                                               ),

                                               plotOutput("devianceplot")


                                           )



                                    )
                                  ),
                                  fluidRow(
                                    column(12,
                                           div(style = "padding:20px",
                                               h2("ROC and zROC"),
                                               div(style = "text-align:right",
                                                   a(id = "toggleROC", "Show explanation", href = "#")
                                               ),
                                               shinyjs::hidden(
                                                 div(id = "ROCinfo",
                                                     htmlOutput("ROCinfo")
                                                 )

                                               ),

                                               uiOutput("ui")


                                           )



                                    )
                                  ),

                                  fluidRow(
                                    column(12,
                                           div(style = "padding:20px",
                                               h2("Observed and predicted responding & residuals"),
                                               div(style = "text-align:right",
                                                   a(id = "toggleRes", "Show explanation", href = "#")
                                               ),
                                               shinyjs::hidden(
                                                 div(id = "residualsinfo",
                                                     htmlOutput("residualsinfo")
                                                 )

                                               ),

                                               plotOutput("residualsplot",height=800)


                                           )



                                    )
                                  ),

                                  fluidRow(
                                    column(12,
                                           div(style = "padding:20px",
                                               h2("Parameter estimates"),
                                               dataTableOutput("parestimates")


                                           )



                                    )
                                  )

                        )
                      )
             )


  )

server <- function(input, output, session) {




  shinyjs::onclick("toggleROC",
                   shinyjs::toggle(id = "ROCinfo", anim = TRUE))

  output$ROCinfo <- renderUI({


    str_e1 <- "ROC and zROC curves show the observed data (black points), and model predictions (lines) based on the parameters of the maximum likelihood estimate of each model."
    str_e2 <- "ROC curves are created by constrasting responses to old items to responses to new items. Where the experimental condition varies by study-test block, models are fitted to each condition, with each condition resulting in an independent set of parameter estimates. Where the experimental condition varies within a study-test block, and multiple sets of new items exist, one set of new items was designated the reference set."

    HTML(paste(str_e1, str_e2, sep = '<p/><p/>'))

  })


  output$ui <- renderUI({

    output$rocplot <- renderPlot({

      if(input$exp == "SB2020_e1"){

        #inp<-which(exps==input$SB2021_e1)
        makeROC(SB2020_e1 %>% filter(id==input$SB2020_e1),
                SB2020_e123_fits %>% filter(id == input$SB2020_e1))
      } else if(input$exp == "SB2020_e2"){

        #inp<-which(exps==input$SB2020_e1)
        makeROC(SB2020_e2 %>% filter(id==input$SB2020_e2),
                SB2020_e123_fits %>% filter(id == input$SB2020_e2))
      } else if(input$exp == "SB2020_e3"){

        #inp<-which(exps==input$SB2020_e1)
        makeROC(SB2020_e3 %>% filter(id==input$SB2020_e3),
                SB2020_e233_fits %>% filter(id == input$SB2020_e3))

      } else if(input$exp == "SB2021_e1"){

        #inp<-which(exps==input$SB2021_e1)
        makeROC(SB2021_e1 %>% filter(id==input$SB2021_e1),
              SB2021_e1_fits %>% filter(id == input$SB2021_e1))
      } else if(input$exp == "SB2021_e2"){

        #inp<-which(exps==input$SB2021_e1)
        makeROC(SB2021_e2 %>% filter(id==input$SB2021_e2),
                SB2021_e1_fits %>% filter(id == input$SB2021_e2))
      } else if(input$exp == "SB2021_e3"){

        #inp<-which(exps==input$SB2021_e1)
        makeROC(SB2021_e3 %>% filter(id==input$SB2021_e3),
                SB2021_e3_fits %>% filter(id == input$SB2021_e3))

      } else if(input$exp == "MWW2007_e1"){

        #inp<-which(exps==input$SB2021_e1)
        makeROC(MWW2007 %>% filter(id==input$MWW2007_e1),
                MWW2007_fits %>% filter(id == input$MWW2007_e1))

      } else if(input$exp == "MWW2007_e2"){

        #inp<-which(exps==input$SB2021_e1)
        makeROC(MWW2007 %>% filter(id==input$MWW2007_e2),
                MWW2007_fits %>% filter(id == input$MWW2007_e2))

      } else if(input$exp == "MWW2007_e3"){

        #inp<-which(exps==input$SB2021_e1)
        makeROC(MWW2007 %>% filter(id==input$MWW2007_e3),
                MWW2007_fits %>% filter(id == input$MWW2007_e3))

      } else if(input$exp == "D2007_e1a"){

        #inp<-which(exps==input$SB2021_e1)
        makeROC(D2007 %>% filter(id==input$D2007_e1a),
                D2007_fits %>% filter(id == input$D2007_e1a))

      } else if(input$exp == "D2007_e1b"){

        #inp<-which(exps==input$SB2021_e1)
        makeROC(D2007 %>% filter(id==input$D2007_e1b),
                D2007_fits %>% filter(id == input$D2007_e1b))

      }



    })

    # if(input$radio == 2) {
    #   exps <- exps2
    #   dexp <- d2
    #   selectee <- input$exp2
    # } else {
    #   exps <- exps1
    #   dexp <- d1
    #   selectee <- input$exp1
    # }
    # inp<-which(exps==selectee)
    #
    # ln <-  length(unique(dexp[[inp]] %>% .$condition))
    # nn <- ifelse((ln > 3 & ln <= 6) , 2, ifelse(ln > 6,3,1))

    plotOutput("rocplot", height = 800)

  })

  shinyjs::onclick("toggleDev",
                   shinyjs::toggle(id = "devianceinfo", anim = TRUE))







  output$devianceinfo <- renderUI({


    str_e1 <- "Quantitative fit, expressed as Delta AIC where Delta AIC = 0 shows the winning model. For each model and condition, the p-value associated with the G^2 statistics is additionally shown where p > .05 suggests that the observed data is adequately captured by the model."


    HTML(paste(str_e1, sep = '<p/><p/>'))

  })

  output$density <- renderPlot({

    plotdensity



  })

  output$devianceplot <- renderPlot({


    if(input$exp == "SB2020_e1"){

      #inp<-which(exps==input$SB2021_e1)
      makeDeviance(SB2020_e123_fits %>% filter(id == input$SB2020_e1),
                   Gsqs %>% filter(id == input$SB2020_e1))
    } else if(input$exp == "SB2020_e2"){

      #inp<-which(exps==input$SB2021_e1)
      makeDeviance(SB2020_e123_fits %>% filter(id == input$SB2020_e2),
                   Gsqs %>% filter(id == input$SB2020_e2))
    }else if(input$exp == "SB2020_e3"){

      #inp<-which(exps==input$SB2021_e1
      makeDeviance(SB2020_e123_fits %>% filter(id == input$SB2020_e3),
                   Gsqs %>% filter(id == input$SB2020_e3))
    }else if(input$exp == "SB2021_e1"){

      #inp<-which(exps==input$SB2021_e1)
      makeDeviance(SB2021_e1_fits %>% filter(id == input$SB2021_e1),
                   Gsqs %>% filter(id == input$SB2021_e1))
    } else if(input$exp == "SB2021_e2"){

      #inp<-which(exps==input$SB2021_e1)
      makeDeviance(SB2021_e1_fits %>% filter(id == input$SB2021_e2),
                   Gsqs %>% filter(id == input$SB2021_e2))
    }else if(input$exp == "SB2021_e3"){

      #inp<-which(exps==input$SB2021_e1)
      makeDeviance(SB2021_e3_fits %>% filter(id == input$SB2021_e3),
                   Gsqs %>% filter(id == input$SB2021_e3))
    }else if(input$exp == "MWW2007_e1"){

      #inp<-which(exps==input$SB2021_e1)
      makeDeviance(MWW2007_fits %>% filter(id == input$MWW2007_e1),
                   Gsqs %>% filter(id == input$MWW2007_e1))
    } else if(input$exp == "MWW2007_e2"){

      #inp<-which(exps==input$SB2021_e2)
      makeDeviance(MWW2007_fits %>% filter(id == input$MWW2007_e2),
                   Gsqs %>% filter(id == input$MWW2007_e2))
    } else if(input$exp == "MWW2007_e3"){

      #inp<-which(exps==input$SB2021_e)
      makeDeviance(MWW2007_fits %>% filter(id == input$MWW2007_e3),
                   Gsqs %>% filter(id == input$MWW2007_e3))
    } else if(input$exp == "D2007_e1a"){

      #inp<-which(exps==input$SB2021_e)
      makeDeviance(D2007_fits %>% filter(id == input$D2007_e1a),
                   Gsqs %>% filter(id == input$D2007_e1a))
    }else if(input$exp == "D2007_e1b"){

      #inp<-which(exps==input$SB2021_e)
      makeDeviance(D2007_fits %>% filter(id == input$D2007_e1b),
                   Gsqs %>% filter(id == input$D2007_e1b))
    }


  })

  shinyjs::onclick("toggleRes",
                   shinyjs::toggle(id = "residualsinfo", anim = TRUE))




  output$residualsinfo <- renderUI({


    str_e1 <- "Difference in predicted confidence bin probabilities and observed proportions of rating responses."


    HTML(paste(str_e1, sep = '<p/><p/>'))

  })

  output$residualsplot <- renderPlot({


    if(input$exp == "SB2020_e1"){

      #inp<-which(exps==input$SB2021_e1)
      makeFreqprediction(SB2020_e1 %>% filter(id == input$SB2020_e1),
                         Preds$SB2020_e123 %>% filter(id ==  input$SB2020_e1))
    } else if(input$exp == "SB2020_e2"){

      #inp<-which(exps==input$SB2021_e1)
      makeFreqprediction(SB2021_e2 %>% filter(id == input$SB2020_e2),
                         Preds$SB2020_e123 %>% filter(id ==  input$SB2020_e2))
    }else if(input$exp == "SB2020_e3"){

      #inp<-which(exps==input$SB2021_e1)
      makeFreqprediction(SB2020_e3 %>% filter(id == input$SB2020_e3),
                         Preds$SB2020_e123 %>% filter(id ==  input$SB2020_e3))
    }else if(input$exp == "SB2021_e1"){

      #inp<-which(exps==input$SB2021_e1)
      makeFreqprediction(SB2021_e1 %>% filter(id == input$SB2021_e1),
                         Preds$SB2021_e12 %>% filter(id ==  input$SB2021_e1))
    } else if(input$exp == "SB2021_e2"){

      #inp<-which(exps==input$SB2021_e1)
      makeFreqprediction(SB2021_e2 %>% filter(id == input$SB2021_e2),
                         Preds$SB2021_e12 %>% filter(id ==  input$SB2021_e2))
    }else if(input$exp == "SB2021_e3"){

      #inp<-which(exps==input$SB2021_e1)
      makeFreqprediction(SB2021_e3 %>% filter(id == input$SB2021_e3),
                         Preds$SB2021_e3 %>% filter(id ==  input$SB2021_e3))
    }else if(input$exp == "MWW2007_e1"){

      #inp<-which(exps==input$SB2021_e1)
      makeFreqprediction(MWW2007 %>% filter(id == input$MWW2007_e1),
                         Preds$MWW2007 %>% filter(id ==  input$MWW2007_e1))
    } else if(input$exp == "MWW2007_e2"){

      #inp<-which(exps==input$SB2021_e1)
      makeFreqprediction(MWW2007 %>% filter(id == input$MWW2007_e2),
                         Preds$MWW2007 %>% filter(id ==  input$MWW2007_e2))
    } else if(input$exp == "MWW2007_e3"){

      #inp<-which(exps==input$SB2021_e1)
      makeFreqprediction(MWW2007 %>% filter(id == input$MWW2007_e3),
                         Preds$MWW2007 %>% filter(id ==  input$MWW2007_e3))
    } else if(input$exp == "D2007_e1a"){

      #inp<-which(exps==input$SB2021_e1)
      makeFreqprediction(D2007 %>% filter(id == input$D2007_e1a),
                         Preds$D2007 %>% filter(id ==  input$D2007_e1a))
    } else if(input$exp == "D2007_e1b"){

      #inp<-which(exps==input$SB2021_e1)
      makeFreqprediction(D2007 %>% filter(id == input$D2007_e1b),
                         Preds$D2007 %>% filter(id ==  input$D2007_e1b))
    }





  })
  tabs <- reactive({
    #actual reactive
    if(input$exp == "SB2020_e1"){
      makePar(SB2020_e123_fits %>% filter(id == input$SB2020_e1))
    } else if(input$exp == "SB2020_e2"){
      makePar(SB2020_e123_fits %>% filter(id == input$SB2020_e2))
    } else if(input$exp == "SB2020_e3"){
      makePar(SB2020_e123_fits %>% filter(id == input$SB2020_e3))
    } else if(input$exp == "SB2021_e1"){
      makePar(SB2021_e1_fits %>% filter(id == input$SB2021_e1))
    } else if(input$exp == "SB2021_e2"){
      makePar(SB2021_e1_fits %>% filter(id == input$SB2021_e2))
    } else if(input$exp == "MWW2007_e1"){
    makePar(MWW2007_fits %>% filter(id == input$MWW2007_e1))
    } else if(input$exp == "MWW2007_e2"){
      makePar(MWW2007_fits %>% filter(id == input$MWW2007_e2))
    } else if(input$exp == "MWW2007_e3"){
      makePar(MWW2007_fits %>% filter(id == input$MWW2007_e3))
    } else if(input$exp == "SB2021_e3"){
      makePar(SB2021_e3_fits %>% filter(id == input$SB2021_e3))
    } else if(input$exp == "D2007_e1a"){
      makePar(D2007_fits %>% filter(id == input$D2007_e1a))
    } else if(input$exp == "D2007_e1b"){
      makePar(D2007_fits %>% filter(id == input$D2007_e1b))
    }
  })

  output$parestimates <- renderDataTable(


    tabs(),

    rownames= FALSE,
    options = list(searching = FALSE,
                   paging = FALSE,

                   initComplete = JS(
                     "function(settings, json) {",
                     "$(this.api().table().header()).css({'background-color': '#000', 'color': '#fff'});",
                     "}")
    )

#}
)











}

shinyApp(ui, server)

# add in: KFCV, LOP into same plot. empty if necessary
# add number of data sets included of all possible ones
# add explanation at the top