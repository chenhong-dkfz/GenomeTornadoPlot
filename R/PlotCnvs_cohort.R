
plotCnvs.cohort <- function(paralist,SaveAsObject){
  chrom = unlist(paralist["chrom"])
  startPos = unlist(paralist["startPos"])
  endPos = unlist(paralist["endPos"])
  rescore = unlist(paralist["rescore"])
  score = unlist(paralist["score"])
  sorting = unlist(paralist["sorting"])
  title = unlist(paralist["title"])
  pixel.per.cnv = unlist(paralist["pixel.per.cnv"])
  plot.type = unlist(paralist["plot.type"])
  # getcolor
  legend = unlist(paralist["legend"])
  legend.names = unlist(paralist["legend.names"])
  color = unlist(paralist["color"])
  score.values = unlist(paralist["score.values"])
  n = unlist(paralist["n"])
  display = unlist(paralist["display"])
  gene.anno = unlist(paralist["gene.anno"])
  cnv.type = unlist(paralist["cnv.type"])
  start.gene = unlist(paralist["start.gene"])
  end.gene = unlist(paralist["end.gene"])
  sort.method = unlist(paralist["sort.method"])
  color.method = unlist(paralist["color.method"])
  score = unlist(paralist["score"])
  pids = unlist(paralist["pids"])
  cohort = unlist(paralist["cohort"])



  chroms <- chrom[sorting]
  starts <- startPos[sorting]
  ends <- endPos[sorting]
  cohorts <- cohort[sorting]
  #cohorts <- rescore[sorting] #?
  cohorts <- droplevels.factor(cohorts, exclude = if(anyNA(levels(cohorts)))NULL else NA)  ## erase factor levels = 0 (turns out very important for color plotting)
  cnv.number <-  length(chroms) # number of lines in input
  chromWidth <- round((pixel.per.cnv * cnv.number) * 0.1)
  f.score <- focallity.score(m=length(starts),starts = starts,ends = ends)


  if (length(unique(chroms)) > 1){
    print(unique(chroms))
    print("More than one chromosome id - use other function")
    return()
  }

  y <- lengthChromosome(chroms[1],"bases") + 10000000



  plot.new()
  #png("t1.png",width = 1024,height=768,units = "px")
  tiff(file="t1.tiff", width=12, height=8,units="in", compression="lzw", res=150)
  par(c(5,3,4,4))
  pixelPerChrom <- chromWidth + (pixel.per.cnv)*(cnv.number+1)+10 # determines space between chromsomes
  x.size <- pixelPerChrom
  y.size <- y+100
  plot(c(0,x.size),c(0,y.size),type="n",xaxt="n",yaxt="n",xlab="SVs",ylab="Chromosomal location",main=title)
  chrStr <- paste("chr",toString(chroms[1]))
  text(c((chromWidth/2)),c(0),labels=c(chrStr))
  if(gene.anno == TRUE)    ###added RT gene.anno arg
  {
    m <- mean(start.gene,end.gene)
    text(c(-1),c(y-m+(y*0.035)),labels=c(cnv.type),cex=0.5)
    paintCytobands(chroms[1],pos=c(chromWidth,y),units="bases",width=chromWidth-7,orientation="v",legend=FALSE)
    rect(7,y-start.gene,chromWidth-1,y-end.gene,col="gray50", border = "gray50")
    lines(c(0.6,2.25),c(y-m+(y*0.02),y-m),col="gray50")
    lines(c(2.25,5),c(y-m,y-m),col="gray50")
  }else{
    paintCytobands(chroms[1],pos=c(chromWidth,y),units="bases",width=chromWidth,orientation="v",legend=FALSE)
  }

  plotCnv.cohort(chroms,starts,ends,y,
                 chromWidth=chromWidth,pixel.per.cnv=pixel.per.cnv,score=score,
                 cohorts=cohort,startPoint=chromWidth,method=color.method,color=color)


  # legend position decision (top or bottom)
  centro <- c(125,93.3,91,50.4,48.4,61,59.9,45.6,49,40.2,53.7,35.8,17.9,17.6,19,36.6,24,17.2,26.5,27.5,13.2,14.7,60.6,12.5)*1000000
  length <- lengthChromosome(c(1:22,"X","Y"),"bases")/2
  genome <- data.frame(chromosome=c(1:22,"X","Y"), centromere=centro, length=length) # dataframe containing chromosome and centromere position info

  mean.pos <- mean(c(starts,ends)) # mean position of all CNV´s
  #centroo <- genome[genome$chromosome %in% chroms,] # centromere position in the current chromosome
  half.length <- genome[genome$chromosome %in% chroms,] # half of the length of the current chromosome

  # mean CNV is over the centromere -> legend is plotted bottomright
  if(mean.pos < half.length$length){
    xtr <- "bottomright"
    xtf <- c(4,24,20.5,3)
    text(c(pixelPerChrom/2),c(y-10),labels = paste("score: ",f.score),cex=1.2) # score on opposite
  }

  if(mean.pos > half.length$length){
    xtr <- "topright"
    xtf <- c(21.5,24,4,3)
    text(c(pixelPerChrom/2),c(10),labels = paste("score: ",f.score),cex=0.75)
  }



  # legend type decision ----------------------------------------------------------------------------
  if(color.method=="cohort" | color.method=="length"){
    legend.color <- GetColor(method="cohort",color=color,cohorts=cohort)
  }
  if(color.method=="ploidy"){
    legend.color <- GetColor(method="ploidy",color=color,cohorts=cohorts)
  }

  if(legend=="missing" || legend==1){
    legend(xtr,legend=unique(cohorts),col=GetColor(color=color,cohorts=cohorts,q=F,method="by.cohort"),cex=0.75,pch=16) # normal legend
  }

  print(color.method)
  print(table(cohorts))

  if(legend==2){
    par(new=T,mar=xtf )
    #par(new=T,mar=c(2,12,10,1))
    if(color.method=="cohort" | color.method=="length"){
      pie(table(cohorts),col=legend.color,cex=1) # piechart legend
    }else if(color.method=="ploidy"){
      tb <- table(score)
      dp.list <-c("bi-del","mo-del","diploidy","gain-low","gain-mid","gain-high","n/a")
      for(i in 1:length(names(tb))){names(tb)[i] <- dp.list[as.integer(names(tb)[i])]}
      pie(tb,col=legend.color,cex=1)
    }
  }

  # ploidy print pass
  #
  dev.off()

  if(SaveAsObject==TRUE){
    img <- readTIFF("t1.tiff")
    g <- rasterGrob(img, interpolate=TRUE)
    return(g)
  }

}
