library(shiny)
library(shinydashboard)


ui <- dashboardPage( # header

  header <- dashboardHeader(),


  sidebar <- dashboardSidebar(

    sidebarUserPanel("GenomeTornadoPlot2",
                     subtitle = a(href = "#", icon("circle", class = "text-success"), "Online"),
                     image = "userimage.png"
    ),

    sidebarMenu(
      id = "tabs",
      menuItem("Tornado Plot", icon = icon("th"), tabName = "To_pl")
    )

  ),

  body <- dashboardBody(

    tabItems(
      tabItem("headers",
              h1("GenomeTornadoPlot Easy2Use")
      ),


      # tornado plot page
      tabItem("To_pl",

              h1("GenomeTornadoPlot Easy2Use") ,
              box(title = "Genes", status = "primary",
                  br(),
                  br(),
                  fileInput("file1", "Choose CSV File",
                            accept = c(
                              "text/csv",
                              "text/comma-separated-values,text/plain",
                              ".csv")
                  ),
                  selectInput("gene1", label = h3("Gene A"),
                              choices = list("FOXP1" = "FOXP1", "MBD2" = "MBD2", "POLI" = "POLI"),
                              selected = "FOXP1"),
                  br(),
                  br(),
                  selectInput("gene2", label = h3("Gene B"),
                              choices = list("RYBP" = "RYBP", "MBD2" = "MBD2", "POLI" = "POLI"),
                              selected = "RYBP"),
                  br(),
                  br(),
                  br(),

              ),
              box(title = "Options", status = "primary",
                  radioButtons("cnv.type", label = h3("CNV types"),
                              choices = list("deletion" = "del", "duplication" = "dup", "both" = "both"),
                              selected = "del"),
                  radioButtons("sorting_method", label = h3("Sorting method"),
                               choices = list("by cohort" = "cohort", "by ploidy" = "ploidy"),
                               selected = "cohort"),
                  radioButtons("color_method", label = h3("Color method"),
                               choices = list("by cohort" = "cohort", "by ploidy" = "ploidy"),
                               selected = "cohort"),
                  radioButtons("Plot_type", label = h3("Plot type"),
                               choices = list("Twin" = 1, "Normal" = 2),
                               selected = 1),
                  radioButtons("max.length", label = h3("definition of focal"),
                               choices = list("1e6" = 1e6, "1e7" = 1e7),
                               selected = 1e7),

              ),
              br(),
              br(),
              actionButton("Draw_tornado","Tornado!"),
              plotOutput("pplot1"),
              downloadButton('downloadPlot', 'Download Plot')


      )
    )
  )
) # header
