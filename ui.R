library(shiny)
library(shinyjs)
library(visNetwork)
library(shinythemes)
library(shinycssloaders)
library(shinyWidgets)
library(DT)
library(ggplot2)
library(car)
library(nortest)
library(tseries)
library(RcmdrMisc)
library(lmtest)
library(slickR)
library(igraph)

library(shiny)

source("Components/Comp_Network.R")
shinyUI(fluidPage(title = "BIO114",theme = shinytheme("cyborg"),setBackgroundImage(src = "wallpaper.jpg", shinydashboard = F),tags$head( includeScript("https://d3js.org/d3.v3.min.js") ),
                  tags$style(HTML("td {
  border: 1px solid black;
}")),
                  includeScript(path ="http://d3js.org/d3.v3.min.js"),
                  navbarPage(title = div(),
                             #"Biologia Molecular y Biotecnologia (BIO-114)",style="color:white;text-align:center"
                             tabPanel(div(img(src="LogoGrupo.svg",width="200px",height="86px")),
                                      sidebarPanel(br(),
                                                   fileInput("fileNetRoots", "Choose csv File", accept = ".csv"),
                                                   textInput("seqID_R_N",p("ID gene database:",style="color:white; text-align:center")),
                                                   textInput("geneName_R_N",p("Node name:",style="color:white; text-align:center")),
                                                   actionButton("buttonAddNetwork", "Add Gene"),
                                                   br(),br(),
                                                   actionButton("buttonDeleteNetwork", "Delete Gene"),
                                                   br(),
                                                   hr(),
                                                   br(),
                                                   numericInput("CorteCorNetR",value = 0,max = 1, min = -1, step = 0.1,p("Cutoff correlation value:",style="color:white; text-align:center")),
                                                   numericInput("CortePvalueNetR",value = 0.05,min = 0, max = 0.06,step = 0.01,p("Cutoff p-value value:",style="color:white; text-align:center")),
                                                   actionButton("buttonFilterNetwork", "Filter"),
                                                   br(),
                                                   hr(),
                                                   br(),
                                                   actionButton("buttonEmptyNetwork", "Delete Network"),

                                                   # Only show this panel if the plot type is a histogram
                                                   ),
                                      mainPanel(
                                        visNetworkOutput("networkRoots"),
                                        style = "height:700px;background-color:#FCFBF5;"
                                      ),
                                      br(),br(),
                                      ),
                             
                             tabPanel(div("Root Apex",style="text-align: center;"),
                                      tabsetPanel(
                                        tabPanel("Expression View",br(),
                                      sidebarLayout(
                                        sidebarPanel(
                                          useShinyjs(),
                                          textInput("seqID_R",p("ID gene database:",style="color:white; text-align:center")),
                                          imageOutput("legendtHeatR",height = "150px"),br(),
                                          div(tags$img(src="Root_legend.png",width="310px",height="320px")),
                                          
                                          width = 3
                                        ),
                                        mainPanel(
                                          fluidRow(
                                            uiOutput("descriptionSeq_R"),br()
                                          ),
                                          fluidRow(
                                            
                                            column(fluidRow(
                                              h5(p("Ammonium condition (3mM)",style="color:white;text-align:center"),style="background-color:black;border-radius: 8px"),
                                              imageOutput("heatRoots_N")),width = 4,style="text-align: center;"),
  
                                            column(fluidRow(
                                              useShinyjs(),
                                              h5(p("Differential expression",style="color:white;text-align:center"),style="background-color:black;border-radius: 8px"),
                                              tableOutput("diffRootTable"),br(),
                                              tableOutput("tableGOR"),br(),
                                              plotOutput("barplotRoots"), align="center"),br(),br(),
                                              width = 4,style="text-align: center;"),
                                            column(fluidRow(
                                              h5(p("Control condition",style="color:white;text-align:center"),style="background-color:black;border-radius: 8px"),
                                              imageOutput("heatRoots_C")),width = 4,style="text-align: center;")
                                          )
                                        ))),tabPanel("Correlation",br(),
                                                     sidebarPanel(
                                                       actionButton("buttonCorR", "Generate correlation table"),
                                                       width = 3
                                                     ),mainPanel(
                                                       DT::dataTableOutput("tableCorR") %>% withSpinner(color="#0dc5c1")
                                                     ))
                                      ,tabPanel("GO search",br(),
                                                sidebarPanel(
                                                  textInput("goR",p("Gene Onthology ID:",style="color:white; text-align:center")),
                                                ),mainPanel(
                                                  DT::dataTableOutput("searchGO_R") %>% withSpinner(color="#0dc5c1")
                                                )
                                      ),
                                      tabPanel("General",br(),sidebarLayout(
                                        sidebarPanel(
                                          textInput("goRTotal",p("Gene Onthology ID:",style="color:white; text-align:center")),
                                          selectInput("moduleRTotal", h5("Module selection"), 
                                                      choices = moduleList(), selected = 1),
                                          checkboxGroupInput("checkGroupR", 
                                                             h5("Select diferential condition"), 
                                                             choices = list("RC Ammonium > RC Control" = "RC,Ammonium condition", 
                                                                            "RC Control > RC Ammonium" = "RC,Control condition",
                                                                            "RDC Ammonium > RDC Control" = "RDC,Ammonium condition", 
                                                                            "RDC Control > RDC Ammonium" = "RDC,Control condition",
                                                                            "RM Ammonium > RM Control" = "RM,Ammonium condition", 
                                                                            "RM Control > RM Ammonium" = "RM,Control condition",
                                                                            "RDV Ammonium > RDV Control" = "RDV,Ammonium condition", 
                                                                            "RDV Control > RDV Ammonium" = "RDV,Control condition", 
                                                                            "RC > RDC" = "RC-RDC,RC",
                                                                            "RC < RDC" = "RC-RDC,RDC",
                                                                            "RC > RDV" = "RC-RDV,RC",
                                                                            "RC < RDV" = "RC-RDV,RDV",
                                                                            "RC > RM" = "RC-RM,RC",
                                                                            "RC < RM" = "RC-RM,RM",
                                                                            "RDC > RDV" = "RDC-RDV,RDC",
                                                                            "RDC < RDV" = "RDC-RDV,RDV",
                                                                            "RDC > RM" = "RDC-RM,RDC",
                                                                            "RDC < RM" = "RDC-RM,RM",
                                                                            "RDV > RM" = "RDV-RM,RDV",
                                                                            "RDV < RM" = "RDV-RM,RM")),
                                          br(),
                                          actionButton("buttonGeneralR", "Generate table"),
                                          br(),
                                          hr(),
                                          numericInput("clusterHeatR","HCluster cutree",value = 5),
                                          br(),
                                          fileInput("fileHeatmapRoots", "Choose csv File", accept = ".csv"),
                                          downloadButton("downloadHeatRoots", "Download heatmap"),
                                          width = 3
                                        ),mainPanel(
                                          DT::dataTableOutput("tableGeneralRoots") %>% withSpinner(color="#0dc5c1")
                                        )))
                                      )
                                      
                             )), class='flex-center'
                  
                  
))