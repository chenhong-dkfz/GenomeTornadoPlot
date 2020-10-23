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
                                                          gene.anno,start.gene,end.gene,color.method,sort.method,SaveAsObject,
                                                          multi_panel){
  if(missing(SaveAsObject)){SaveAsObject = TRUE}
  if(missing(color.method)){color.method = "cohort"}
  if(missing(sort.method)){sort.method = "length"}
  if(missing(multi_panel)){multi_panel = FALSE}
  paralist0 <- CNV.by.method(object,gene.name,pids,title,legend,legend.names,
                             out.dir,file.type,pixel.per.cnv,color,display,
                             gene.anno,start.gene,end.gene,color.method,sort.method)
  if(SaveAsObject==TRUE){
    if(multi_panel==FALSE){
      plotlist0 <- PlotTwins(paralist=paralist0,SaveAsObject=SaveAsObject)
    }else{
      plot_multipanel_single(paralist=paralist0)
    }
  }else{
    print("Output image is saved!!")
  }
})
