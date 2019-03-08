
setMethod("TornadoPlots",signature("CNV_twin"),function(object,gene.name,pids,title,legend,legend.names,
                                                        out.dir,file.type,pixel.per.cnv,color,display,
                                                        gene.anno,start.gene,end.gene,color.method,sort.method,SaveAsObject){
  paralist0 <- PlotTwinsInit()
  if(SaveAsObject==TRUE){
    plot0 <- PlotTwins(paralist=paralist0,SaveAsObject=SaveAsObject)
  }else{
    print("Output image is saved!!")
  }
})

