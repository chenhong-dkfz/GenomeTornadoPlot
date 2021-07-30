#'Plot Tornado Plots
#'
#' @param CNV.input Output of MakeData function, Grange object
#' @param gene.name gene name, string
#' @param pids the given patient id, string
#' @param legend type of legend display, int
#' @param legend.names string
#' @param out.dir path of output plot, string
#' @param file.type output plot file type, string
#' @param pixel.per.cnv int
#' @param color colors of tornado, list
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
                                                          multi_panel,zoomed){
  if(missing(SaveAsObject)){SaveAsObject = TRUE}
  if(missing(color.method)){color.method = "cohort"}
  if(missing(sort.method)){sort.method = "length"}
  if(missing(multi_panel)){multi_panel = FALSE}
  if(missing(zoomed)){multi_panel = FALSE}
  paralist0 <- CNV.by.method(object,gene.name,pids,title,legend,legend.names,
                             out.dir,file.type,pixel.per.cnv,color,display,
                             gene.anno,start.gene,end.gene,color.method,sort.method,zoomed)
  if(SaveAsObject==TRUE){
    if(multi_panel==FALSE){
      plotlist0 <- plotCnvs.cohort(paralist=paralist0,SaveAsObject=SaveAsObject)
    }else{
      plot_multipanel_single(paralist=paralist0)
    }
  }else{
    print("Output image is saved!!")
  }
})
