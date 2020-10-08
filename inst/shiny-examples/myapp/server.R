server <- function(input, output, session) {


  # /Users/chen/Tornado_manuscript/gene_candidates/FOXP1_RYBP.txt


  observeEvent(input$file1, {
    inFile <- input$file1
    dty <- fread(inFile$datapath,data.table = F)
    colnames(dty) <- c("Chromosome","Start","End","Score","Gene","Cohort","PID")
    updateSelectInput(session, "gene1", label = "Gene A", choices = dty$Gene)
    updateSelectInput(session, "gene2", label = "Gene B", choices = dty$Gene)

  })


  observeEvent(input$Draw_tornado,{
    inFile <- input$file1
    print(inFile)
    dty <- fread(inFile$datapath,data.table = F)
    colnames(dty) <- c("Chromosome","Start","End","Score","Gene","Cohort","PID")
    dty$length <- dty$End-dty$Start
    dty <- dty[dty$length<=1000000,]
    dty$Chromosome <- 3 #NEAT1 chr11.
    input_gene_1 <- input$gene1
    input_gene_2 <- input$gene2
    plot_type <- as.numeric(input$Plot_type)
    print(plot_type)
    print(input_gene_1)
    print(input_gene_2)
    print("x")
    sdt <- MakeData(CNV=dty,gene_name_1 = input_gene_1,gene_name_2 = input_gene_2,score.type = "del")
    plotlist1s <- TornadoPlots(sdt,color.method="cohort",sort.method="length",SaveAsObject = T)
    out1 <<- plotlist1s[[plot_type]]
    output$pplot1 <- renderPlot({grid.arrange(out1)})
  })

  output$downloadPlot <- downloadHandler(
    filename = function() {
      paste('out1', '.pdf', sep='')
    },
    content = function(file) {
      pdf(file)

      grid.arrange(out1)

      dev.off()

    }
  )

}
