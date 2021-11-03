#'Plot Tornado Plots
#'
#' @param CNV.input Output of MakeData function, Grange object
#' @param gene.name gene name. Character.
#' @param pids the given patient id. Character.
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
#' @export



setMethod("TornadoPlots",signature("CNV_single"),function(object,gene.name,pids,title,legend,legend.names,
                                                          out.dir,file.type,pixel.per.cnv,color,display,
                                                          gene.anno,start.gene,end.gene,color.method,sort.method,SaveAsObject,
                                                          multi_panel,zoomed,font.size.factor,path,format,orient){
  if(missing(SaveAsObject)){SaveAsObject = TRUE}
  if(missing(format)){format = "tiff"}
  if(missing(path)){path = ""}
  if(missing(color.method)){color.method = "cohort"}
  if(missing(sort.method)){sort.method = "length"}
  if(missing(multi_panel)){multi_panel = FALSE}
  if(missing(zoomed)){zoomed = FALSE}
  if(missing(font.size.factor)){font.size.factor = 1}
  if(missing(orient)){orient = "v"}


  paralist0 <- CNV.by.method(CNV.input=object,gene.name=gene.name,pids=pids,title=title,legend=legend,
                             legend.names=legend.names,
                             out.dir=out.dir,file.type=file.type,pixel.per.cnv=pixel.per.cnv,
                             color=color,display=display,
                             gene.anno=gene.anno,start.gene=start.gene,end.gene=end.gene,
                             color.method=color.method,sort.method=sort.method,zoomed=zoomed,
                             SaveAsObject=SaveAsObject,format=format,path=path,orient=orient)
  if(SaveAsObject==TRUE){
    if(multi_panel==FALSE){
      plotlist0 <- plotCnvs.cohort(paralist=paralist0,SaveAsObject=SaveAsObject,font.size.factor=font.size.factor)
    }else{
      plot_multipanel_single(paralist=paralist0,font.size.factor=font.size.factor,orient=orient)
    }
  }else{
    print("Output image is saved!!")
  }
})
