#'Plot Tornado Twin Plots
#'
#' @param CNV.input Grange object
#' @param gene.name string
#' @param pid string
#' @param legend int
#' @param legend.names string
#' @param out.dir string
#' @param file.type string
#' @param pixel.per.cnv int
#' @param color list
#' @param display string
#' @param gene.anno string
#' @param start.gene string
#' @param end.gene string
#' @param sort.method string
#' @param color.method string
#'

setMethod("TornadoPlots",signature("CNV_twin"),function(object,gene.name,pids,title,legend.type,legend.names,
                                                        out.dir,file.type,pixel.per.cnv,color,display,
                                                        gene.anno,start.gene,end.gene,color.method,sort.method,SaveAsObject){
  paralist0 <- PlotTwinsInit(twin.cnv=object,
                             title=title,legend.type=legend.type,legend.names=legend.names,
                             out.dir=out.dir,color=color,
                             gene.anno=gene.anno,start.gene,end.gene,
                             color.method=color.method,sort.method=sort.method)
  if(SaveAsObject==TRUE){
    plotlist0 <- PlotTwins(paralist=paralist0,SaveAsObject=SaveAsObject)
    #plot0 <- plotlist0[[1]]
    #plot1 <- plotlist0[[2]]
  }else{
    print("Output image is saved!!")
  }
})

