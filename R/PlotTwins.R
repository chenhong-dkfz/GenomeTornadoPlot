
PlotTwins <- function(paralist,SaveAsObject = SaveAsObject){

  chroms_1 <- unlist(paralist["chrom_1"])
  starts_1 <- unlist(paralist["startPos_1"])
  ends_1 <- unlist(paralist["endPos_1"])
  scores_1 <- unlist(paralist["score_1"])
  f.score_1 <- focallity.score(m=length(starts_1),starts = starts_1,ends = ends_1)

  chroms_2 <- unlist(paralist["chrom_2"])
  starts_2 <- unlist(paralist["startPos_2"])
  ends_2 <- unlist(paralist["endPos_2"])
  scores_2 <- unlist(paralist["score_2"])
  f.score_2 <- focallity.score(m=length(starts_2),starts = starts_2,ends = ends_2)

  cnv.number <- (length(chroms_1)+length(chroms_2)) # number of lines in input
  chromWidth <- round((pixel.per.cnv * cnv.number) * 0.1)

  if (length(unique(chroms_1)) > 1){
    print(unique(chroms_1))
    print("More than one chromosome id - use other function")
    return()
  }


  if (length(unique(chroms_2)) > 1){
    print(unique(chroms_2))
    print("More than one chromosome id - use other function")
    return()
  }


  y <- lengthChromosome(chroms_1[1],"bases") + 10000000  ## are you sure??

  # plot parameters -----------------------------------------------------------------------------------------------------------------
  plot.new()
  png("t2.png",width = 1024,height=768,units = "px")

  par(c(5,3,4,4))
  pixelPerChrom_1 <-  (pixel.per.cnv)*(length(chroms_1)+1)
  pixelPerChrom_2 <-  (pixel.per.cnv)*(length(chroms_2)+1)
  pixelPerChrom <- chromWidth+pixelPerChrom_1+pixelPerChrom_2+10 # determines space between chromsomes

  x.size <- pixelPerChrom
  y.size <- y+100
  plot(c(0,x.size),c(0,y.size),type="n",xaxt="n",yaxt="n",xlab="CNVs",ylab="Chromosomal location",main=title)
  chrStr <- paste("chr",toString(chroms_1[1]))
  text(c(pixelPerChrom_1+(chromWidth/2)),c(0),labels=c(chrStr))

  if(gene.anno == TRUE){
    paintCytobands(chroms_1[1],pos=c(pixelPerChrom_1+chromWidth,y),units="bases",width=chromWidth,orientation="v",legend=FALSE)
    m_1 <- mean(start.gene_1,end.gene_1)
    m_2 <- mean(start.gene_2,end.gene_2)

    text(c(pixelPerChrom_1+chromWidth+15),c(y-m_1+(y*0.045)),labels=c(cnv.type_1),cex=0.7)
    rect(pixelPerChrom_1+1,y-m_1,pixelPerChrom_1+chromWidth-1,y-m_1,col="gray50", border = "gray50")
    lines(c(pixelPerChrom_1+chromWidth+7,pixelPerChrom_1+chromWidth+4),c(y-m_1+(y*0.03),y-m_1),col="gray50")
    lines(c(pixelPerChrom_1+chromWidth+4,pixelPerChrom_1+chromWidth+1),c(y-m_1,y-m_1),col="gray50")

    text(c(pixelPerChrom_1-15),c(y-m_2-(y*0.045)),labels=c(cnv.type_2),cex=0.7)
    rect(pixelPerChrom_1+1,y-m_2,pixelPerChrom_1+chromWidth-1,y-m_2,col="gray50", border = "gray50")
    lines(c(pixelPerChrom_1-7,pixelPerChrom_1-4),c(y-m_2-(y*0.03),y-m_2),col="gray50")
    lines(c(pixelPerChrom_1-4,pixelPerChrom_1-1),c(y-m_2,y-m_2),col="gray50")

  }else{
    paintCytobands(chroms_1[1],pos=c(pixelPerChrom_1+chromWidth,y),units="bases",width=chromWidth,orientation="v",legend=FALSE)
  }

  plotCnv(chroms_1,starts_1,ends_1,y,scores_1,pixel.per.cnv=pixel.per.cnv,method="by.length",color=color,score.values = score.values_1,n=n_1,startPoint=(pixelPerChrom_1),direction = "left")
  plotCnv(chroms_2,starts_2,ends_2,y,scores_2,pixel.per.cnv=pixel.per.cnv,method="by.length",color=color,score.values = score.values_2,n=n_2,startPoint=(pixelPerChrom_1+chromWidth),direction = "right")

  #plotCnv.cohort <- function(chroms,starts,ends,y,chromWidth,pixel.per.cnv,cohorts,startPoint,color,method)

  # legend parameters ------------------------------------------------------------------------------------------------------
  df.color.ploidy <- data.frame(color=color[1:n], # colors according to getColor.ploidy/2
                                score=score.values[1:n], # score according to
                                names=legend.names[1:n])

  data.score <- data.frame(score=c(sort(unique(scores)))) # unique and present scores of inout data
  dtt <- df.color.ploidy[df.color.ploidy$score %in% data.score$score,] # subset only present scores of input data
  color <- as.vector(dtt$color)
  labs <- as.vector(dtt$names)

  # legend position decision (top or bottom)
  centro <- c(125,93.3,91,50.4,48.4,61,59.9,45.6,49,40.2,53.7,35.8,17.9,17.6,19,36.6,24,17.2,26.5,27.5,13.2,14.7,60.6,12.5)*1000000
  length <- lengthChromosome(c(1:22,"X","Y"),"bases")/2
  genome <- data.frame(chromosome=c(1:22,"X","Y"), centromere=centro, length=length) # dataframe containing chromosome and centromere position info

  mean.pos <- mean(c(starts,ends)) # mean position of all CNV´s
  #centroo <- genome[genome$chromosome %in% chroms,] # centromere position in the current chromosome
  half.length <- genome[genome$chromosome %in% chroms,] # half of the length of the current chromosome

  if(mean.pos < half.length$length) {
    xtr <- "bottomleft"
    xtr2 <- "bottomright"
    xtf <- c(4,5,20.5,22)
    xtf2 <- c(4,24,20.5,3)
    text(c(pixelPerChrom_1/2),c(y-10),labels = paste("score: ",f.score_1),cex=0.75)
    text(c(pixelPerChrom_1+chromWidth+(pixelPerChrom_2/2)),c(y-10),labels = paste("score: ",f.score_2),cex=0.75)

  }    # mean start end smaller than subset chrom centromer
  if(mean.pos > half.length$length){
    xtr <- "topleft"
    xtr2 <- "topright"
    xtf <- c(21.5,5,4,22)
    xtf2 <- c(21.5,24,4,3)
    text(c(pixelPerChrom_1/2),c(10),labels = paste("score: ",f.score_1),cex=0.75)
    text(c(pixelPerChrom_1+chromWidth+(pixelPerChrom_2/2)),c(10),labels = paste("score: ",f.score_2),cex=0.75)

  }

  # legend type decision ----------------------------------------------------------------------------
  if(legend=="missing" || legend==1){
    legend(xtr,legend=labs,col=color,cex=0.75,pch=16) # normal legend
    legend(xtr2,legend=labs,col=color,cex=0.75,pch=16)
  }
  if(legend==2) {
    par(new=T,mar=xtf )
    pie(table(score_1),labels=labs,col=color,cex=0.52) # piechart legend
    par(new=T,mar=xtf2)
    pie(table(score_2),labels=labs,col=color,cex=0.52)
  } else{} # no legend

  dev.off()
  if(SaveAsObject==TRUE){
    img <- readPNG("t2.png")
    g <- rasterGrob(img, interpolate=TRUE)
    return(g)
  }


}
