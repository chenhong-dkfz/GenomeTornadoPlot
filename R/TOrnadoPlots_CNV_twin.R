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
                                                        out.dir,file.type,pixel.per.cnv,color,display,cnv.type_1,cnv.type_2,
                                                        gene.anno,start.gene,end.gene,color.method,sort.method,SaveAsObject,
                                                        multi_panel,font.size.factor,file,format){
  if(missing(SaveAsObject)){SaveAsObject = TRUE}
  if(missing(format)){format = "tiff"}
  if(missing(file)){file = paste0("output_plot.",format)}
  if(missing(color.method)){color.method = "cohort"}
  if(missing(sort.method)){sort.method = "length"}
  if(missing(multi_panel)){multi_panel = FALSE}
  if(missing(font.size.factor)){font.size.factor = 1}
  paralist0 <- PlotTwinsInit(twin.cnv=object,
                             title=title,legend.type=legend.type,legend.names=legend.names,
                             out.dir=out.dir,color=color,cnv.type_1=cnv.type_1,cnv.type_2=cnv.type_2,
                             gene.anno=gene.anno,start.gene,end.gene,
                             color.method=color.method,sort.method=sort.method,
                             SaveAsObject=SaveAsObject,format=format,file=file)

  if(SaveAsObject==TRUE){
    if(multi_panel==FALSE){
    plotlist0 <- PlotTwins(paralist=paralist0,SaveAsObject=SaveAsObject)
    }else{
    plot_multipanel_twin(paralist=paralist0)
    }
  }else{
    print("Output image is saved!!")
  }
})

