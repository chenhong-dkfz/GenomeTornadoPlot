PlotTwins <- function(paralist,SaveAsObject = SaveAsObject){

  chroms_1 <- unlist(paralist["chrom_1"])
  starts_1 <- unlist(paralist["startPos_1"])
  ends_1 <- unlist(paralist["endPos_1"])
  rescores_1 <- unlist(paralist["rescore_1"])
  f.score_1 <- unlist(paralist["f.score_1"])

  chroms_2 <- unlist(paralist["chrom_2"])
  starts_2 <- unlist(paralist["startPos_2"])
  ends_2 <- unlist(paralist["endPos_2"])
  rescores_2 <- unlist(paralist["rescore_2"])
  f.score_2 <- unlist(paralist["f.score_2"])

  chroms <- chroms_1

  sorting_1 <- unlist(paralist["sorting_1"])
  sorting_2 <- unlist(paralist["sorting_2"])

  cohort_1 <- unlist(paralist["cohort_1"])
  cohort_2 <- unlist(paralist["cohort_2"])
  cohort_1 <- factor(cohort_1)
  cohort_2 <- factor(cohort_2)
  print("cohort size:")
  print(cohort_1)
  print(cohort_2)

  repeat_1 <- unlist(paralist["repeat_1"])
  repeat_2 <- unlist(paralist["repeat_2"])


  gene.name_1 <- unlist(paralist["gene.name_1"])
  gene.name_2 <- unlist(paralist["gene.name_2"])

  pixel.per.cnv <- as.numeric(paralist["pixel.per.cnv"])
  cnv.number <- (length(chroms_1)+length(chroms_2)) # number of lines in input
  chromWidth <- round((pixel.per.cnv * cnv.number) * 0.1)
  gene.anno <- paralist$gene.anno
  if(gene.anno!=""){gene.anno<-FALSE}
  title <- paralist$title
  color.method <- paralist$color.method

  legend.type <- paralist$legend.type

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
  #png("t2.png",width = 1024,height=768,units = "px")
  tiff(file="t2.tiff", width=12, height=8,units="in", compression="lzw", res=150)

  par(c(5,3,4,4))
  pixelPerChrom_1 <-  (pixel.per.cnv)*(length(chroms_1)+1)
  pixelPerChrom_2 <-  (pixel.per.cnv)*(length(chroms_2)+1)
  pixelPerChrom <- chromWidth+pixelPerChrom_1+pixelPerChrom_2+10 # determines space between chromsomes

  x.size <- pixelPerChrom
  y.size <- y+100
  plot(c(0,x.size),c(0,y.size),type="n",xaxt="n",yaxt="n",xlab="CNVs",
       ylab="Chromosomal location",main=title)
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

  cohort_max <- sort(unique(c(levels(cohort_1),levels(cohort_2))))
  color.value <- GetColor(method=color.method,cohorts=cohort_max)

  plotCnv(chroms_1,starts_1,ends_1,y,rescores_1,pixel.per.cnv=pixel.per.cnv,
          sorting = sorting_1, cohort = cohort_1,cohort_max = cohort_max,
          color.value = color.value,
          color.method="cohort",score.values = score.values_1,n=n_1,
          startPoint=(pixelPerChrom_1),direction = "left")
  plotCnv(chroms_2,starts_2,ends_2,y,rescores_2,pixel.per.cnv=pixel.per.cnv,
          sorting = sorting_2, cohort = cohort_2,cohort_max = cohort_max,
          color.value = color.value,
          color.method="cohort",score.values = score.values_2,n=n_2,
          startPoint=(pixelPerChrom_1+chromWidth),direction = "right")

  # legend parameters ------------------------------------------------------------------------------------------------------


  if(color.method == "ploidy"){
    color.base <- colorRampPalette(c("red2","indianred4","royalblue4","steelblue1","chartreuse3","darkgreen","grey"))(7)
    if(length(color)<6){color <- color.base}
    if(n_1>=n_2){nmax=n_1}else{nmax=n_2}
    #legend.names = unlist(paralist["legend.names"])
    legend.names = c("deepdel","monodel","normal","gainlow","gainmid","gainhigh","unknown")
    df.color.ploidy <- data.frame(color=color[1:7], # colors according to getColor.ploidy/2
                                  score=c(1:7), # score according to ??
                                  names=legend.names)
    data.score <- data.frame(score=c(sort(unique(c(rescores_1,rescores_2))))) # unique and present scores of inout data
    dtt <- df.color.ploidy[df.color.ploidy$score %in% data.score$score,] # subset only present scores of input data
    color <- as.vector(dtt$color)
    labs <- as.vector(dtt$names)
    factor_1 <- rescores_1
    factor_2 <- rescores_2

  }else{
    #cohort.max <- unique(c(levels(cohort_1),levels(cohort_2)))
    cohort.dim <- length(cohort_max)
    df.color.cohort <- data.frame(color=color.value, # colors according to getColor.ploidy/2
                                  score=cohort_max, # score according to ??
                                  names=cohort_max)
    dtt <- df.color.cohort
    color <- as.vector(dtt$color)
    labs <- as.vector(dtt$names)
    factor_1 <- cohort_1
    factor_2 <- cohort_2
  }






  # legend position decision (top or bottom)
  centro <- c(125,93.3,91,50.4,48.4,61,59.9,45.6,49,40.2,53.7,35.8,17.9,17.6,19,36.6,24,17.2,26.5,27.5,13.2,14.7,60.6,12.5)*1000000
  length <- lengthChromosome(c(1:22,"X","Y"),"bases")/2
  genome <- data.frame(chromosome=c(1:22,"X","Y"), centromere=centro, length=length) # dataframe containing chromosome and centromere position info

  mean.pos <- mean(c(starts_1,ends_1)) # mean position of all CNV´s
  #centroo <- genome[genome$chromosome %in% chroms,] # centromere position in the current chromosome
  half.length <- genome[genome$chromosome %in% chroms_1,] # half of the length of the current chromosome

  if(mean.pos < half.length$length) {
    xtr <- "bottomleft"
    xtr2 <- "bottomright"
    xtf <- c(4,5,20.5,22)
    xtf2 <- c(4,24,20.5,0)
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

  print(legend.type)

  # legend type decision ----------------------------------------------------------------------------
  if(legend.type=="normal"){
    legend(xtr,legend=labs,col=color,cex=0.75,pch=16) # normal legend
    legend(xtr2,legend=labs,col=color,cex=0.75,pch=16)
    print("normal legend.")
  }
  if(legend.type =="pie") {
    par(new=T,mar=xtf )

    factor_1 <- droplevels.factor(factor_1, exclude = if(anyNA(levels(factor_1)))NULL else NA)
    factor_2 <- droplevels.factor(factor_2, exclude = if(anyNA(levels(factor_2)))NULL else NA)

    dtt_1 <- dtt[dtt$names%in%levels(factor_1),]
    color_1 <- as.vector(dtt_1$color)
    labs_1 <- dtt_1$score
    pie(table(factor_1),labels=labs_1,col=color_1,cex=0.45,radius = 0.6) # piechart legend

    par(new=T,mar=xtf2)
    dtt_2 <- dtt[dtt$names%in%levels(factor_2),]
    color_2 <- as.vector(dtt_2$color)
    labs_2 <- dtt_2$score
    pie(table(factor_2),labels=labs_2,col=color_2,cex=0.45,radius = 0.6)
    print("pie plot legend！")
  } else{} # no legend

  dev.off()
  if(SaveAsObject==TRUE){
    img <- readTIFF("t2.tiff")
    g <- rasterGrob(img, interpolate=TRUE)
  }

  # plot parameters for unique plots



  plot.new()
  #png("t1.png",width = 1024,height=768,units = "px")
  tiff(file="t3.tiff", width=12, height=8,units="in", compression="lzw", res=150)
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
    paintCytobands(chroms_1[1],pos=c(chromWidth,y),units="bases",width=chromWidth,orientation="v",legend=FALSE)
  }
  chroms <- c(chroms_1,chroms_2)
  starts <- c(starts_1,starts_2)
  ends <- c(ends_1,ends_2)
  score <- c(rescores_1,rescores_2)
  cohorts <- c(cohort_1,cohort_2)
  startPoint <- chromWidth
  method <- "repeat"
  repeats <- c(repeat_1,repeat_2)

  plotCnv.cohort(chroms,starts,ends,y,
                 chromWidth=chromWidth,pixel.per.cnv=pixel.per.cnv,score=score,
                 cohorts=cohorts,startPoint=chromWidth,method="repeat",color=color,rep=repeats)


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
    text(c(pixelPerChrom/2),c(y-10),labels = paste("score: ",f.score_1,",",f.score_2),cex=1.2) # score on opposite
  }

  if(mean.pos > half.length$length){
    xtr <- "topright"
    xtf <- c(21.5,24,4,3)
    text(c(pixelPerChrom/2),c(10),labels = paste("score: ",f.score_1,",",f.score_2),cex=0.75)
  }


  # legend type decision ----------------------------------------------------------------------------

  if(2==2 || legend=="pie"){
    legend.color <- c("red","grey","blue")
    par(new=T,mar=xtf )
    #par(new=T,mar=c(2,12,10,1))
    repu <- repeats!="C2"
    repeats_common <- repeats[repu]
    repeats_common <- factor(repeats_common,levels = c("U1","C1","U2"))
    pie(table(repeats_common),col=legend.color,cex=1,labels = c(gene.name_1,"both",gene.name_2))

  }

  # ploidy print pass
  #
  dev.off()

  if(SaveAsObject==TRUE){
    img <- readTIFF("t3.tiff")
    h <- rasterGrob(img, interpolate=TRUE)
  }

  results <- list(g,h)
  return(results)

}
