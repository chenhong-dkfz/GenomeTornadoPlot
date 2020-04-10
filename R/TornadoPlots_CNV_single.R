#'Plot Tornado Plots
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
#' @export



setMethod("TornadoPlots",signature("CNV_single"),function(object,gene.name,pids,title,legend,legend.names,
                                                          out.dir,file.type,pixel.per.cnv,color,display,
                                                          gene.anno,start.gene,end.gene,color.method,sort.method,SaveAsObject){
  if(missing(SaveAsObject)){SaveAsObject = FALSE}
  if(missing(color.method)){color.method = "cohort"}
  if(missing(sort.method)){sort.method = "length"}
  paralist0 <- CNV.by.method(object,gene.name,pids,title,legend,legend.names,
                             out.dir,file.type,pixel.per.cnv,color,display,
                             gene.anno,start.gene,end.gene,color.method,sort.method)
  SaveAsObject = TRUE
  if(SaveAsObject==TRUE){
    plotlist0 <- plotCnvs.cohort(paralist=paralist0,SaveAsObject = SaveAsObject)
    #plot0 <- plotlist0[[1]]
    #plot1 <- plotlist0[[2]]
  }else{
    print("Output image is saved!!")
  }
})
