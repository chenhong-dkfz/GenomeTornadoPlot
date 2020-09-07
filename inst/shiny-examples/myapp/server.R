server <- function(input, output) {

  observeEvent(input$Draw_tornado,{
    input_gene_1 <- input$gene1
    input_gene_2 <- input$gene2
    plot_type <- as.numeric(input$Plot_type)
    print(plot_type)
    print(input_gene_2)
    print("x")
    sdt <- MakeData(CNV=dt,gene_name_1 = input_gene_1,gene_name_2 = input_gene_2,score.type = "del")
    plotlist1s <- TornadoPlots(sdt,color.method="cohort",sort.method="length",SaveAsObject = T)
    out1 <- plotlist1s[[plot_type]]
    output$pplot1 <- renderPlot({grid.arrange(out1)})
  })

  output$downloadPlot <- downloadHandler(
    filename = function() {
      paste('out1', '.pdf', sep='')
    },
    content = function(file) {
      pdf(file)

      plot(grid.arrange(out1))

      dev.off()
    }
  )

}
