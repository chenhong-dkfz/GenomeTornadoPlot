server <- function(input, output, session) {


  # /Users/chen/Tornado_manuscript/gene_candidates/FOXP1_RYBP.txt


  observeEvent(input$file1, {
    inFile <- input$file1
    # if(input$pcawg_chr=="user"){
    #   inFile <- input$file1
    #   print(inFile)
    #   dty <- fread(inFile$datapath,data.table = F)
    # }else{
    #   chr.nr <- input$pcawg_chr
    #   data.name <- paste0("data/",chr.nr,".Rdata")
    #   print(data.name)
    #   data("chr3",package = "tornado.test.1")
    #   dty <- cnv_chr
    #   print("loading finished")
    # }
    dty <- fread(inFile$datapath,data.table = F)
    print("loading finished!")
    colnames(dty) <- c("Chromosome","Start","End","Score","Gene","Cohort","PID")
    print("loading finished!")
    updateSelectInput(session, "gene1", label = "Gene A", choices = dty$Gene)
    print("loading finished!")
    updateSelectInput(session, "gene2", label = "Gene B", choices = dty$Gene)
    print("loading finished!")

  })



  observeEvent(input$Draw_tornado,{
    # if(input$pcawg_chr=="user"){
    # inFile <- input$file1
    # print(inFile)
    # dty <- fread(inFile$datapath,data.table = F)
    # }else{
    #   chr.nr <- input$pcawg_chr
    #   data.name <- paste0("data/",chr.nr,".Rdata")
    #   load(data.name)
    #   dty <- cnv_chr
    # }
    inFile <- input$file1
    #print(inFile)
    dty <- fread(inFile$datapath,data.table = F)
    colnames(dty) <- c("Chromosome","Start","End","Score","Gene","Cohort","PID")
    dty$length <- dty$End-dty$Start
    max.length <- input$max.length
    dty <- dty[dty$length<=10000000,]
    dty <- dty[dty$length<=max.length,]
    input_gene_1 <- input$gene1
    input_gene_2 <- input$gene2
    cnv.type <- input$cnv.type
    plot_type <- as.numeric(input$Plot_type)
    color_method <- as.character(input$color_method)
    sort_method <- as.character(input$sorting_method)
    print(plot_type)
    print(input_gene_1)
    print(input_gene_2)
    print("x")
    sdt <- MakeData(CNV=dty,gene_name_1 = input_gene_1,gene_name_2 = input_gene_2,score.type = "del")
    plotlist1s <- TornadoPlots(sdt,cnv.type_1 = cnv.type,cnv.type_2 = cnv.type, color.method=color_method,sort.method=sort_method,SaveAsObject = T)
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
