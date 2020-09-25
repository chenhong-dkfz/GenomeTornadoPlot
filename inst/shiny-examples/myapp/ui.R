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
            selectInput("gene1", label = h3("Gene A"),
                        choices = list("MLLT3" = "MLLT3", "MBD2" = "MBD2", "POLI" = "POLI"),
                        selected = "MBD2"),
            br(),
            br(),
            selectInput("gene2", label = h3("Gene B"),
                        choices = list("MLLT3" = "MLLT3", "MBD2" = "MBD2", "POLI" = "POLI"),
                        selected = "POLI"),
            br(),
            br(),
            br()
            ),
            box(title = "Options", status = "primary",
            radioButtons("Score", label = h3("Score type"),
                         choices = list("Edge score" = 1, "FS/local rank" = 2),
                         selected = 1),
            radioButtons("Focal_threshold", label = h3("Focal threshold"),
                         choices = list("1e-6" = 1, "1e-7" = 2),
                         selected = 1),
            radioButtons("Plot_type", label = h3("Plot type"),
                         choices = list("Twin" = 1, "Normal" = 2),
                         selected = 1)
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
