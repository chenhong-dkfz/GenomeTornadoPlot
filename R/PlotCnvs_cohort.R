plotCnvs.cohort <- function(paralist,SaveAsObject){
  gene_name = unlist(paralist["gene.name"])
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
  #score.values = unlist(paralist["score.values"])
  n = unlist(paralist["n"])
  display = unlist(paralist["display"])
  gene.anno = unlist(paralist["gene.anno"])
  cnv.type = unlist(paralist["cnv.type"])
  start.gene = unlist(paralist["start.gene"])
  end.gene = unlist(paralist["end.gene"])
  sort.method = unlist(paralist["sort.method"])
  color.method = unlist(paralist["color.method"])
  #score = unlist(paralist["score"])
  pids = unlist(paralist["pids"])
  cohort = unlist(paralist["cohort"])
  rescore = unlist(paralist["rescore"])
  sorting.color <- unlist(paralist["sorting.color"])
  sorting.plot.color <- cohort[sorting.color]
  f.score <- unlist(paralist["f.score"])

  # original

  chroms <- chrom[sorting]
  starts <- startPos[sorting]
  ends <- endPos[sorting]
  cohorts <- cohort[sorting]
  rescore <- rescore[sorting]
  cohorts <- droplevels.factor(cohorts, exclude = if(anyNA(levels(cohorts)))NULL else NA)  ## erase factor levels = 0 (turns out very important for color plotting)




  # single side #
  # plot parameters ---------------------------------------------------------------------

  if(cnv.type=="dup"){
  del_dup.index <- rescore>3
  }else{
  del_dup.index <- rescore<3
  }

  chroms_0 <- chroms[del_dup.index]
  cohort_0 <- cohorts[del_dup.index]
  starts_0 <- starts[del_dup.index]
  ends_0 <- ends[del_dup.index]
  rescore_0 <- rescore[del_dup.index]
  #sorting_1 <- sorting[del.index]
  sorting_0 <- order(ends_0 - starts_0,cohort_0)
  #score.values_0 <- score.values[del_dup.index]

  cnv.number_0 <-  length(chroms_0) # number of lines in input
  chromWidth_0 <- round((pixel.per.cnv * cnv.number_0) * 0.1)

  cohort_0 <- droplevels.factor(cohort_0, exclude = if(anyNA(levels(cohort_0)))NULL else NA)

  if (length(unique(chroms_0)) > 1){
    print(unique(chroms_0))
    print("More than one chromosome id - use other function")
    return()
  }

  y <- lengthChromosome(chroms_0[1],"bases") + 10000000


  legend.type = legend

  plot.new()
  tiff(file="t1.tiff", width=12, height=8,units="in", compression="lzw", res=150)
  par(c(5,3,4,4))
  pixelPerChrom <- chromWidth_0 + (pixel.per.cnv)*(cnv.number_0+1)+10 # determines space between chromsomes
  x.size <- pixelPerChrom
  y.size <- y+100
  nsample <- length(chroms_0)
  if(cnv.type=="dup"){
    tcnv = "duplication"
  }else{
    tcnv = "deletion"
  }
  ncohort <- length(unique(cohort_0))
  print(cohort_0)

  print(cohorts)
  title <- paste0(gene_name,": ",nsample," ",tcnv," events from ",ncohort," cohorts")
  plot(c(0,x.size),c(0,y.size),type="n",xaxt="n",yaxt="n",xlab="CNVs",ylab="Chromosomal location",main=title)
  chrStr <- paste("chr",toString(chroms_0[1]))
  text(c((chromWidth_0/2)),c(0),labels=c(chrStr))
  if(gene.anno == TRUE)    ###added RT gene.anno arg
  {
    m <- mean(start.gene,end.gene)
    text(c(-1),c(y-m+(y*0.035)),labels=c(cnv.type),cex=0.5)
    paintCytobands(chroms_0[1],pos=c(chromWidth_0,y),units="bases",width=chromWidth_0-7,orientation="v",legend=FALSE)
    rect(7,y-start.gene,chromWidth_0-1,y-end.gene,col="gray50", border = "gray50")
    lines(c(0.6,2.25),c(y-m+(y*0.02),y-m),col="gray50")
    lines(c(2.25,5),c(y-m,y-m),col="gray50")
  }else{
    paintCytobands(chroms_0[1],pos=c(chromWidth_0,y),units="bases",width=chromWidth_0,orientation="v",legend=FALSE)
  }

  plotCnv.cohort(chroms_0,starts_0,ends_0,y,
                 chromWidth=chromWidth_0,pixel.per.cnv=pixel.per.cnv,score=rescore_0,
                 cohorts=cohort_0,startPoint=chromWidth_0,method=color.method)


  # legend position decision (top or bottom)
  centro <- c(125,93.3,91,50.4,48.4,61,59.9,45.6,49,40.2,53.7,35.8,17.9,17.6,19,36.6,24,17.2,26.5,27.5,13.2,14.7,60.6,12.5)*1000000
  length <- lengthChromosome(c(1:22,"X","Y"),"bases")/2
  genome <- data.frame(chromosome=c(1:22,"X","Y"), centromere=centro, length=length) # dataframe containing chromosome and centromere position info

  mean.pos <- mean(c(starts_0,ends_0)) # mean position of all CNV´s
  #centroo <- genome[genome$chromosome %in% chroms,] # centromere position in the current chromosome
  half.length <- genome[genome$chromosome %in% chroms_0,] # half of the length of the current chromosome

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
    legend.color <- GetColor(method="cohort",color=color,cohorts=cohort_0)
  }
  if(color.method=="ploidy"){
    legend.color <- GetColor(method="ploidy",color=color,cohorts=cohort_0)
    print(legend.color)
    print("it is")
  }

  if(mean.pos < half.length$length) {
    sub.position <- c(.75, 1, .3, .7)
  }    # mean start end smaller than subset chrom centromer
  if(mean.pos > half.length$length){
    sub.position <- c(.75, 1, .6, 1)
  }

  if(color.method=="cohort" | color.method=="length"){
    if(legend==2 || legend=="pie"){
      par(new=T,mar=xtf)
      pie(table(cohort_0),col=legend.color,cex=1)
    }else{
      legend(xtr,legend=unique(cohort_0),col=legend.color,cex=0.75,pch=16) # normal legend
    }
  }else if(color.method=="ploidy"){

    tb <- table(rescore_0)
    print(rescore_0)
    dp.list <-c("bi-del","mo-del","diploidy","gain-low","gain-mid","gain-high","n/a")
    for(i in 1:length(names(tb))){names(tb)[i] <- dp.list[as.integer(names(tb)[i])]}
    legend.color.subset <- legend.color[sort(unique(rescore_0))]
    if(legend==2 || legend=="pie"){
      par(new=T,mar=xtf )
      pie(tb,col=legend.color.subset,cex=1)
    }else{

      legend(xtr,legend=unique(dp.list),col=legend.color,cex=0.75,pch=16)
    }
  }


  dev.off()

  if(SaveAsObject==TRUE){
    img <- readTIFF("t1.tiff")
    g <- rasterGrob(img, interpolate=TRUE)

  }




  # both side #
  # plot parameters -----------------------------------------------------------------------------------------------------------------


  del.index <- rescore<3
  dup.index <- rescore>3

  chroms_1 <- chroms[del.index]
  cnv.type_1 <- cnv.type
  cohort_1 <- cohorts[del.index]
  starts_1 <- starts[del.index]
  ends_1 <- ends[del.index]
  rescores_1 <- rescore[del.index]
  #sorting_1 <- sorting[del.index]
  sorting_1 <- order(ends_1 - starts_1,cohort_1)
  #score.values_1 <- score.values[del.index]

  chroms_2 <- chroms[dup.index]
  cnv.type_2 <- cnv.type
  cohort_2 <- cohorts[dup.index]
  starts_2 <- starts[dup.index]
  ends_2 <- ends[dup.index]
  rescores_2 <- rescore[dup.index]
  sorting_2 <- order(ends_2 - starts_2,cohort_2)
  #score.values_2 <- score.values[dup.index]

  cnv.number_d <-  length(chroms_1)+length(chroms_2) # number of lines in input
  chromWidth_d <- round((pixel.per.cnv * cnv.number_d) * 0.1)


  plot.new()
  #png("t5.png",width = 1024,height=768,units = "px")
  tiff(file="t5.tiff", width=12, height=8,units="in", compression="lzw", res=150)

  par(c(5,3,4,4))
  pixelPerChrom_1 <-  (pixel.per.cnv)*(length(chroms_1)+1)
  pixelPerChrom_2 <-  (pixel.per.cnv)*(length(chroms_2)+1)
  pixelPerChrom <- chromWidth_d+pixelPerChrom_1+pixelPerChrom_2+10 # determines space between chromsomes

  x.size <- pixelPerChrom
  y.size <- y+100
  ndel <- length(starts_1)
  ndup <- length(starts_2)
  title <- paste0(gene_name,": ",ndel," deletions and ",ndup," duplications")
  plot(c(0,x.size),c(0,y.size),type="n",xaxt="n",yaxt="n",xlab="CNVs",
       ylab="Chromosomal location",main=title)
  chrStr <- paste("chr",toString(chroms_1[1]))
  text(c(pixelPerChrom_1+(chromWidth_d/2)),c(0),labels=c(chrStr))

  if(gene.anno == TRUE){
    paintCytobands(chroms_1[1],pos=c(pixelPerChrom_1+chromWidth_d,y),units="bases",width=chromWidth_d,orientation="v",legend=FALSE)
    m_1 <- mean(start.gene_1,end.gene_1)
    m_2 <- mean(start.gene_2,end.gene_2)

    text(c(pixelPerChrom_1+chromWidth_d+15),c(y-m_1+(y*0.045)),labels=c(cnv.type_1),cex=0.7)
    rect(pixelPerChrom_1+1,y-m_1,pixelPerChrom_1+chromWidth_d-1,y-m_1,col="gray50", border = "gray50")
    lines(c(pixelPerChrom_1+chromWidth_d+7,pixelPerChrom_1+chromWidth_d+4),c(y-m_1+(y*0.03),y-m_1),col="gray50")
    lines(c(pixelPerChrom_1+chromWidth_d+4,pixelPerChrom_1+chromWidth_d+1),c(y-m_1,y-m_1),col="gray50")

    text(c(pixelPerChrom_1-15),c(y-m_2-(y*0.045)),labels=c(cnv.type_2),cex=0.7)
    rect(pixelPerChrom_1+1,y-m_2,pixelPerChrom_1+chromWidth_d-1,y-m_2,col="gray50", border = "gray50")
    lines(c(pixelPerChrom_1-7,pixelPerChrom_1-4),c(y-m_2-(y*0.03),y-m_2),col="gray50")
    lines(c(pixelPerChrom_1-4,pixelPerChrom_1-1),c(y-m_2,y-m_2),col="gray50")

  }else{
    paintCytobands(chroms_1[1],pos=c(pixelPerChrom_1+chromWidth_d,y),units="bases",width=chromWidth_d,orientation="v",legend=FALSE)
  }

  cohort_max <- sort(unique(c(levels(cohort_1),levels(cohort_2))))
  color.value <- GetColor(method=color.method,cohorts=cohort_max)


  sorting_ploidy_1 <- order(rescores_1,ends_1 - starts_1)
  sorting_ploidy_2 <- order(rescores_2,ends_2 - starts_2)

  plotCnv(chroms_1,starts_1,ends_1,y,rescores_1,pixel.per.cnv=pixel.per.cnv,
          sorting = sorting_ploidy_1, cohort = cohort_1,cohort_max = cohort_max,
          color.value = color.value,
          color.method="ploidy",n=n_1,
          startPoint=(pixelPerChrom_1),direction = "left")

  plotCnv(chroms_2,starts_2,ends_2,y,rescores_2,pixel.per.cnv=pixel.per.cnv,
          sorting = sorting_ploidy_2, cohort = cohort_2,cohort_max = cohort_max,
          color.value = color.value,
          color.method="ploidy",n=n_2,
          startPoint=(pixelPerChrom_1+chromWidth_d),direction = "right")



  # legend parameters ------------------------------------------------------------------------------------------------------



  #cohort.max <- unique(c(levels(cohort_1),levels(cohort_2)))
  cohort.dim <- length(cohort_max)






  if(color.method == "cohort"){
    df.color.cohort <- data.frame(color=color.value, # colors according to getColor.ploidy/2
                                  score=cohort_max, # score according to ??
                                  names=cohort_max)
    dtt <- df.color.cohort
    factor_1 <- cohort_1
    factor_2 <- cohort_2
  }else if(color.method == "ploidy"){
    df.color.ploidy <- data.frame(color=color.value, # colors according to getColor.ploidy/2
                                  score=c(1,2,3,4,5,6,7), # score according to ??
                                  names=c("bi-del","mo-del","diploidy","gain-low","gain-mid","gain-high","n/a")
    )
    dtt <- df.color.ploidy
    ploidy_levels <- c("bi-del","mo-del","diploidy","gain-low","gain-mid","gain-high","n/a")
    factor_1 <- ploidy_levels[rescores_1]
    factor_2 <- ploidy_levels[rescores_2]
  }



  color <- as.vector(dtt$color)
  labs <- as.vector(dtt$names)






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
    #xtf <- c(4,5,20.5,22)
    #xtf2 <- c(4,24,20.5,0)
    text(c(pixelPerChrom_1/2),c(y-10),labels = "deletions",cex=1)
    text(c(pixelPerChrom_1+chromWidth_d+(pixelPerChrom_2/2)),c(y-10),labels = "duplications",cex=1)

  }    # mean start end smaller than subset chrom centromer
  if(mean.pos > half.length$length){
    xtr <- "topleft"
    xtr2 <- "topright"
    xtf <- c(21.5,5,4,22)
    xtf2 <- c(21.5,24,4,3)
    text(c(pixelPerChrom_1/2),c(10),labels = "deletions",cex=1)
    text(c(pixelPerChrom_1+chromWidth_d+(pixelPerChrom_2/2)),c(10),labels = "duplications",cex=1)

  }

  print(legend.type)

  # legend type decision ----------------------------------------------------------------------------
  if(legend.type=="normal"){
    legend(xtr,legend=labs,col=color,cex=0.75,pch=16) # normal legend
    legend(xtr2,legend=labs,col=color,cex=0.75,pch=16)
    print("normal legend.")
  }
  if(legend.type =="pie") {
    xtf2 <- c(21.5,24,4,3)

    info1 <- factor_1
    info2 <- factor_2
    t1<-table(info1)
    t2<-table(info2)
    df01 <- data.frame(t1)
    rownames(df01) <- df01$info1
    df02 <- data.frame(t2)
    rownames(df02) <- df02$info2

    cohort_total <- unique(c(names(t1),names(t2)))
    df_cohort_total <- data.frame(matrix(rep(0,2*length(cohort_total)),ncol = 2))
    df_cohort_total$cohort <- cohort_total
    rownames(df_cohort_total)<-cohort_total
    df_cohort_total <- df_cohort_total[ order(row.names(df_cohort_total)), ]

    df_cohort_total1<-merge(df_cohort_total,df01,by="row.names",all.x=T)
    df_cohort_total2<-merge(df_cohort_total,df02,by="row.names",all.x=T)
    freqs <- cbind(df_cohort_total1$Freq,df_cohort_total2$Freq)
    rownames(freqs) <- rownames(df_cohort_total)
    freqs[is.na(freqs)] <- 0

    freqs <- freqs[order(-freqs[,1],freqs[,2]),]
    test0 <- t(freqs)

    if(mean.pos < half.length$length) {
      sub.position <- c(.75, 1, .3, .7)
    }    # mean start end smaller than subset chrom centromer
    if(mean.pos > half.length$length){
      sub.position <- c(.75, 1, .6, 1)
    }


    par(fig = sub.position , mar=c(0,0,5,5), new=TRUE)
    barplot(test0,
            #main="Deletions and Duplications",
            horiz=TRUE,
            xlab="cohorts", col=c("darkblue","red"),las=1,
            cex.main=0.5,cex.axis = 0.5,cex.names = 0.5,
            #legend = c("deletion","duplication"),
            beside=TRUE)
    legend("right",legend=c("DEL","DUP"),col=c("darkblue","red"),
           cex=0.5,pch=16,bty="n",text.width=2)

    print("pie plot legend！")
  } else{} # no legend

  dev.off()

  if(SaveAsObject==TRUE){
    img <- readTIFF("t5.tiff")
    h <- rasterGrob(img, interpolate=TRUE)
  }
  results <- list(g,h)
  return(results)

}
