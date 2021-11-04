#'Plot Tornado Twin Plots
#'
#' @param CNV.input Grange object
#' @param gene.name gene name. Character.
#' @param legend type of legend display. Character.
#' @param legend.names string
#' @param path path of output plot. character.
#' @param formath output plot file type. character.
#' @param SaveAsObject if true, the plot will be saved in R object. If false, it will be saved in a file. Boolean.
#' @param color colors of tornado.List.
#' @param zoomed if true, the plot will be shown in zoomed-in sight. Boolean.
#' @param multi_panel If true, a multiple panel plot will be shown.
#' @param orient vertical or horizontal arrange CNVs and chromosome. It should be "v" or "h".
#' @param font.size.factor multiply for rescale font size. Double.
#' @param sort.method the method for coloring the CNVs. Character.
#' @param color.method the method for coloring the CNVs. Character.
#'

setMethod("TornadoPlots",signature("CNV_twin"),function(object,gene.name,title,legend.type,legend.names,
                                                        pixel.per.cnv,color,cnv.type_1,cnv.type_2,
                                                        gene.anno,start.gene,end.gene,color.method,sort.method,SaveAsObject,
                                                        multi_panel,font.size.factor,file,format,path,orient,zoomed){
  if(missing(SaveAsObject)){SaveAsObject = TRUE}
  if(missing(format)){format = "tiff"}
  if(missing(path)){path = ""}
  if(missing(color.method)){color.method = "cohort"}
  if(missing(sort.method)){sort.method = "length"}
  if(missing(multi_panel)){multi_panel = FALSE}
  if(missing(font.size.factor)){font.size.factor = 1}
  if(missing(zoomed)){zoomed = FALSE}
  if(missing(orient)){orient = "v"}

  paralist0 <- PlotTwinsInit(twin.cnv=object,
                             title=title,legend.type=legend.type,legend.names=legend.names,
                             color=color,cnv.type_1=cnv.type_1,cnv.type_2=cnv.type_2,
                             gene.anno=gene.anno,start.gene=start.gene,end.gene=end.gene,
                             color.method=color.method,sort.method=sort.method,
                             SaveAsObject=SaveAsObject,format=format,path=path,
                             zoomed=zoomed, orient=orient)

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

