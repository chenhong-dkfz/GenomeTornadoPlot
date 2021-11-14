#'Plot Tornado Plots
#'
#' @param CNV.input Output of MakeData function, Grange object
#' @param gene.name gene name. Character.
#' @param title the title shown in the plots. Character.
#' @param legend type of legend display. Character.
#' @param path path of output plot. character.
#' @param formath output plot file type. character.
#' @param SaveAsObject if true, the plot will be saved in R object. If false, it will be saved in a file. Boolean.
#' @param color colors of tornado.List.
#' @param pixel.per.cnv thickness of CNV lines in pixel. Int.
#' @param zoomed indicates if plot will be shown in zoomed-in sight. It should be "global","region" or "gene".Character
#' @param multi_panel If true, a multiple panel plot will be shown.
#' @param orient vertical or horizontal arrange CNVs and chromosome. It should be "v" or "h".
#' @param font.size.factor multiply for rescale font size. Double.
#' @param sort.method the method for coloring the CNVs. Character.
#' @param color.method the method for coloring the CNVs. Character.
#'
#' @export



setMethod("TornadoPlots",signature("CNV_single"),function(object,gene.name,title,legend,
                                                          pixel.per.cnv,color,
                                                          gene.anno,
                                                          start.gene,end.gene,
                                                          color.method,sort.method,SaveAsObject,drop.low.amp,
                                                          multi_panel,zoomed,font.size.factor,path,format,orient){
  if(missing(SaveAsObject)){SaveAsObject = TRUE}
  if(missing(format)){format = "tiff"}
  if(missing(path)){path = ""}
  if(missing(color.method)){color.method = "cohort"}
  if(missing(sort.method)){sort.method = "length"}
  if(missing(multi_panel)){multi_panel = FALSE}
  if(missing(zoomed)){zoomed = "global"}
  if(missing(font.size.factor)){font.size.factor = 1}
  if(missing(orient)){orient = "v"}
  if(missing(drop.low.amp)){drop.low.amp = FALSE}


  paralist0 <- CNV.by.method(CNV.input=object,gene.name=gene.name,title=title,legend=legend,
                             out.dir=out.dir,file.type=file.type,pixel.per.cnv=pixel.per.cnv,
                             color=color,
                             gene.anno=gene.anno,start.gene=start.gene,end.gene=end.gene,
                             color.method=color.method,sort.method=sort.method,zoomed=zoomed,
                             SaveAsObject=SaveAsObject,format=format,path=path,orient=orient,
                             drop.low.amp = drop.low.amp)
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
