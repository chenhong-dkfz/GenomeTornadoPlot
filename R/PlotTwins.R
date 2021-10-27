PlotTwins <- function(paralist,SaveAsObject,font.size.factor){
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

  repeat_1 <- unlist(paralist["repeat_1"])
  repeat_2 <- unlist(paralist["repeat_2"])


  sort.method = unlist(paralist["sort.method"])
  color.method = unlist(paralist["color.method"])

  gene.name_1 <- unlist(paralist["gene.name_1"])
  gene.name_2 <- unlist(paralist["gene.name_2"])

  t_gene_start_1 <- unlist(paralist["t_gene_start_1"])
  t_gene_end_1 <- unlist(paralist["t_gene_end_1"])
  t_gene_start_2 <- unlist(paralist["t_gene_start_2"])
  t_gene_end_2 <- unlist(paralist["t_gene_end_2"])

  zoomed <- unlist(paralist["zoomed"])
  orient <- unlist(paralist["orient"])
  path <- unlist(paralist["path"])
  format <- unlist(paralist["format"])
  SaveAsObject <- unlist(paralist["SaveAsObject"])



  n_1 <- unlist(paralist["n_1"])
  n_2 <- unlist(paralist["n_2"])

  pixel.per.cnv <- as.numeric(paralist["pixel.per.cnv"])
  cnv.number <- (length(chroms_1)+length(chroms_2)) # number of lines in input
  chromWidth <- round((pixel.per.cnv * cnv.number) * 0.1)
  gene.anno <- paralist$gene.anno
  if(gene.anno!=""){gene.anno<-FALSE}
  title <- paralist$title


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



  #######add on#######

  cnv.type_1 = unlist(paralist["cnv.type_1"])
  cnv.type_2 = unlist(paralist["cnv.type_2"])


  if(cnv.type_1=="dup"){
    del_dup.index_1 <- rescores_1>3
  }else if(cnv.type_1=="del") {
    del_dup.index_1 <- rescores_1<3
  }else{
    del_dup.index_1 <- rescores_1>=0
  }

  chroms_0_1 <- chroms_1[del_dup.index_1]
  cohort_0_1 <- cohort_1[del_dup.index_1]
  cohort_entropy_0_1 <- round(entropy(table(cohort_0_1),unit = "log2"),digits = 2)
  starts_0_1 <- starts_1[del_dup.index_1]
  ends_0_1 <- ends_1[del_dup.index_1]
  rescore_0_1 <- rescores_1[del_dup.index_1]
  sorting_0_1 <- sorting_1[del_dup.index_1]
  cnv.number_0_1 <-  length(chroms_0_1) # number of lines in input
  chromWidth_0_1 <- round((pixel.per.cnv * cnv.number_0_1) * 0.1)
  cohort_0_1 <- droplevels.factor(cohort_0_1, exclude = if(anyNA(levels(cohort_0_1)))NULL else NA)
  repeat_0_1 <- repeat_1[del_dup.index_1]
  if (length(unique(chroms_0_1)) > 1){
    print(unique(chroms_0_1))
    print("More than one chromosome id - use other function")
    return()
  }


  if(cnv.type_2=="dup"){
    del_dup.index_2 <- rescores_2>3
  }else if(cnv.type_2=="del") {
    del_dup.index_2 <- rescores_2<3
  }else{
    del_dup.index_2 <- rescores_2>=0
  }

  chroms_0_2 <- chroms_2[del_dup.index_2]
  cohort_0_2 <- cohort_2[del_dup.index_2]
  cohort_entropy_0_2 <- round(entropy(table(cohort_0_2),unit = "log2"),digits = 2)
  starts_0_2 <- starts_2[del_dup.index_2]
  ends_0_2 <- ends_2[del_dup.index_2]
  rescore_0_2 <- rescores_2[del_dup.index_2]
  sorting_0_2 <- sorting_2[del_dup.index_1]
  cnv.number_0_2 <-  length(chroms_0_2) # number of lines in input
  chromWidth_0_2 <- round((pixel.per.cnv * cnv.number_0_2) * 0.1)
  cohort_0_2 <- droplevels.factor(cohort_0_2, exclude = if(anyNA(levels(cohort_0_2)))NULL else NA)
  repeat_0_2 <- repeat_2[del_dup.index_2]

  if (length(unique(chroms_0_2)) > 1){
    print(unique(chroms_0_2))
    print("More than one chromosome id - use other function")
    return()
  }


  y <- lengthChromosome(chroms_0_1[1],"bases") + 10000000

  ######end addon############



  # plot parameters -----------------------------------------------------------------------------------------------------------------
  plot.new()

  if(SaveAsObject==TRUE){
    if(orient=="v"){
      tiff(file="t2.tiff", width=12, height=8,units="in", compression="lzw", res=150)}else{
        tiff(file="t2.tiff", width=8, height=12,units="in", compression="lzw", res=150)}
  }else{ # SaveAsObject==FALSE
    if(format=="tiff"){
      if(orient=="v"){
        tiff(file=paste0(path,"/bi-plot.tiff"), width=12, height=8,units="in", compression="lzw", res=150)}else{
          tiff(file=paste0(path,"/bi-plot.tiff"), width=8, height=12,units="in", compression="lzw", res=150)}
    }else{ # format = EPS
      if(orient=="v"){
        setEPS()
        postscript(file=paste0(path,"/bi-plot.tiff"),width=12,height=8)}else{
          setEPS()
          postscript(file=paste0(path,"/bi-plot.tiff"),width=8,height=12)
        }
    } # end else format = EPS
  }


  par(c(5,3,4,4))
  cnv.number_0 <- cnv.number_0_1 + cnv.number_0_2
  chromWidth_0 <- round((pixel.per.cnv * cnv.number_0) * 0.1)
  chromWidth_0 <- 10
  pixelPerChrom_0_1 <-  (pixel.per.cnv)*(length(chroms_0_1)+1)
  pixelPerChrom_0_2 <-  (pixel.per.cnv)*(length(chroms_0_2)+1)
  pixelPerChrom_0 <- chromWidth_0+pixelPerChrom_0_1+pixelPerChrom_0_2+10 # determines space between chromsomes

  #x.size <- pixelPerChrom_0
  #y.size <- y+100

  if(orient=="v"){
    x.size <- pixelPerChrom
    y.size <- y+100}else{
      y.size <- pixelPerChrom
      x.size <- y+100
    }

  #plot(c(0,x.size),c(0,y.size),type="n",xaxt="n",yaxt="n",xlab="CNVs",
  #     ylab="Chromosomal location",main=title)

  t_gene_start_1 = as.numeric(as.character(t_gene_start_1))
  t_gene_end_1 = as.numeric(as.character(t_gene_end_1))
  t_gene_start_2 = as.numeric(as.character(t_gene_start_2))
  t_gene_end_2 = as.numeric(as.character(t_gene_end_2))

  if(zoomed==TRUE){
    if(orient=="v"){
      y1 = min(y - t_gene_end_1,y-t_gene_end_2)
      y2 = max(y - t_gene_start_1,y-t_gene_start_2)
      plot(c(0,x.size),c(y1-1e7,y2+1e7),type="n",xaxt="n",yaxt="n",
           xlab="CNVs",ylab="Chromosomal location",main=title)
      lines(c(0,x.size), c(y-t_gene_start_1,y-t_gene_start_1), lty=3, lwd=1)
      lines(c(0,x.size), c(y-t_gene_end_1,y-t_gene_end_1), lty=3, lwd=1)
      lines(c(0,x.size), c(y-t_gene_start_2,y-t_gene_start_2), lty=3, lwd=1)
      lines(c(0,x.size), c(y-t_gene_end_2,y-t_gene_end_2), lty=3, lwd=1)
    }else{ # if orient is h
      y1 = max(t_gene_end_1,t_gene_end_2)
      y2 = min(t_gene_start_1,t_gene_start_2)
      plot(c(y2-1e7,y1+1e7),c(0,y.size),type="n",xaxt="n",yaxt="n",xlab="CNVs",ylab="Chromosomal location",main=title)
      lines(c(t_gene_start_1,t_gene_start_1), c(0,y.size),lty=3, lwd=1)
      lines(c(t_gene_end_2,t_gene_end_2), c(0,y.size),lty=3, lwd=1)
      lines(c(t_gene_start_2,t_gene_start_2), c(0,y.size),lty=3, lwd=1)
      lines(c(t_gene_end_1,t_gene_end_1), c(0,y.size),lty=3, lwd=1)
    }

  }else{
    if(orient=="v"){
      plot(c(0,x.size),c(0,y.size),type="n",xaxt="n",yaxt="n",xlab="CNVs",
           ylab="Chromosomal location",main=title)}else{
             plot(c(0,x.size),c(0,y.size),type="n",xaxt="n",yaxt="n",ylab="CNVs",
                  xlab="Chromosomal location",main=title)
           }
  }

  chrStr <- paste("chr",toString(chroms_0_1[1]))
  text(c(pixelPerChrom_0_1+(chromWidth_0/2)),c(0),labels=c(chrStr))



  if(gene.anno == TRUE){
    paintCytobands(chroms_0_1[1],pos=c(pixelPerChrom_0_1+chromWidth_0,y),units="bases",width=chromWidth_0,orientation="v",legend=FALSE)
    m_1 <- mean(c(t_gene_start_1,t_gene_end_1))
    m_2 <- mean(c(t_gene_start_2,t_gene_end_2))

    # annotation right side
    text(c(pixelPerChrom_0_1+chromWidth_0+15),c(y-m_1+(y*0.045)),labels=c(gene.name_2),cex=0.7)
    lines(c(pixelPerChrom_0_1+chromWidth_0+7,pixelPerChrom_0_1+chromWidth_0+4),c(y-m_1+(y*0.03),y-m_1),col="gray50")

    # annotation left side
    text(c(pixelPerChrom_0_1-15),c(y-m_2-(y*0.045)),labels=c(gene.name_1),cex=0.7)
    lines(c(pixelPerChrom_0_1-7,pixelPerChrom_0_1-4),c(y-m_2-(y*0.03),y-m_2),col="gray50")



  }else{ # acutally this is what really use
    #paintCytobands(chroms_0_1[1],pos=c(pixelPerChrom_0_1+chromWidth_0,y),units="bases",bands="major",width=chromWidth_0,orientation="v",legend=FALSE)
    if(orient=="v"){
      paintCytobands(chroms_1[1],pos=c(pixelPerChrom_0_1+chromWidth_0,y),units="bases",bands="major",width=chromWidth_0,orientation="v",legend=FALSE)}else{
        paintCytobands(chroms_1[1],pos=c(0,pixelPerChrom_0_1+1*chromWidth_0),units="bases",bands="major",width=chromWidth_0,orientation="h",legend=FALSE)
      }
  }


  cohort_max <- sort(unique(c(levels(cohort_0_1),levels(cohort_0_2))))
  color.value <- GetColor(method=color.method,cohorts=cohort_max)

  if(sort.method=="cohort"){
    sorting_x_1 <- order(cohort_0_1,ends_0_1 - starts_0_1)
    sorting_x_2 <- order(cohort_0_2,ends_0_2 - starts_0_2)
  }else if(sort.method=="ploidy"){
    sorting_x_1 <- order(rescore_0_1,ends_0_1 - starts_0_1)
    sorting_x_2 <- order(rescore_0_2,ends_0_2 - starts_0_2)
  }else{
    sorting_x_1 <- order(ends_0_1 - starts_0_1)
    sorting_x_2 <- order(ends_0_2 - starts_0_2)
  }

  # plotCnv(chroms_0_1,starts_0_1,ends_0_1,y,rescore_0_1,pixel.per.cnv=pixel.per.cnv,
  #         sorting = sorting_x_1, cohort = cohort_0_1,cohort_max = cohort_max,
  #         color.value = color.value,
  #         color.method=color.method,n=n_1,
  #         startPoint=(pixelPerChrom_0_1),direction = "left")
  # plotCnv(chroms_0_2,starts_0_2,ends_0_2,y,rescore_0_2,pixel.per.cnv=pixel.per.cnv,
  #         sorting = sorting_x_2, cohort = cohort_0_2,cohort_max = cohort_max,
  #         color.value = color.value,
  #         color.method=color.method,n=n_2,
  #         startPoint=(pixelPerChrom_0_1+chromWidth_0),direction = "right")

  if(orient=="v"){

    plotCnv(chroms_0_1,starts_0_1,ends_0_1,y,rescore_0_1,pixel.per.cnv=pixel.per.cnv,
                sorting = sorting_x_1, cohort = cohort_0_1,cohort_max = cohort_max,
                color.value = color.value,
                color.method=color.method,n=n_1,
                startPoint=(pixelPerChrom_0_1),direction = "left")

    plotCnv(chroms_0_2,starts_0_2,ends_0_2,y,rescore_0_2,pixel.per.cnv=pixel.per.cnv,
                sorting = sorting_x_2, cohort = cohort_0_2,cohort_max = cohort_max,
                color.value = color.value,
                color.method=color.method,n=n_2,
                startPoint=(pixelPerChrom_0_1+chromWidth_0),direction = "right")

  }else{
    plotCnv(chroms_0_1,starts_0_1,ends_0_1,y,rescore_0_1,pixel.per.cnv=pixel.per.cnv,
                sorting = sorting_x_1, cohort = cohort_0_1,cohort_max = cohort_max,
                color.value = color.value,
                color.method=color.method,n=n_1,
                startPoint=(pixelPerChrom_0_1),direction = "bottom")

    plotCnv(chroms_0_2,starts_0_2,ends_0_2,y,rescore_0_2,pixel.per.cnv=pixel.per.cnv,
                sorting = sorting_x_2, cohort = cohort_0_2,cohort_max = cohort_max,
                color.value = color.value,
                color.method=color.method,n=n_2,
                startPoint=(pixelPerChrom_0_1+chromWidth_0),direction = "top")
  }

  # legend parameters ------------------------------------------------------------------------------------------------------


  # dashline

  if(color.method == "ploidy"){
    color.base <- colorRampPalette(c("red2","indianred4","royalblue4","steelblue1","chartreuse3","darkgreen","grey"))(7)
    #if(length(color)<6){color <- color.base}
    color <- color.base
    if(n_1>=n_2){nmax=n_1}else{nmax=n_2}
    legend.names = c("bi-del","mo-del","CN<5","4<CN<9","CN>8")
    df.color.ploidy <- data.frame(color=color.value, # colors according to getColor.ploidy/2
                                  score=c(1,2,3,4,5), # score according to ??
                                  names=c("bi-del","mo-del","CN<5","4<CN<9","CN>8")
    )
    dtt <- df.color.ploidy
    ploidy_levels <- c("bi-del","mo-del","CN<5","4<CN<9","CN>8")
    factor_1 <- ploidy_levels[rescore_0_1]
    factor_2 <- ploidy_levels[rescore_0_2]


  }else{
    cohort.dim <- length(cohort_max)
    df.color.cohort <- data.frame(color=color.value, # colors according to getColor.ploidy/2
                                  score=cohort_max, # score according to ??
                                  names=cohort_max)
    dtt <- df.color.cohort
    color <- as.vector(dtt$color)
    labs <- as.vector(dtt$names)
    factor_1 <- cohort_0_1
    factor_2 <- cohort_0_2
  }


  #plot.x <- recordPlot(load=NULL, attach=NULL)

  if(zoomed==FALSE&orient=="v"){

    # legend position decision (top or bottom)
    centro <- c(125,93.3,91,50.4,48.4,61,59.9,45.6,49,40.2,53.7,35.8,17.9,17.6,19,36.6,24,17.2,26.5,27.5,13.2,14.7,60.6,12.5)*1000000
    length <- lengthChromosome(c(1:22,"X","Y"),"bases")/2
    genome <- data.frame(chromosome=c(1:22,"X","Y"), centromere=centro, length=length) # dataframe containing chromosome and centromere position info

    mean.pos <- mean(c(starts_0_1,ends_0_1)) # mean position of all CNV´s
    half.length <- genome[genome$chromosome %in% chroms_0_1,] # half of the length of the current chromosome

    if(mean.pos < half.length$length) {
      xtr <- "bottomleft"
      xtr2 <- "bottomright"
      xtf <- c(4,5,20.5,22)
      xtf2 <- c(4,24,20.5,0)
      text(c(pixelPerChrom_0_1/2),c(y-10),labels = paste("score: ",f.score_1," entropy:",cohort_entropy_0_1),cex=1)
      text(c(pixelPerChrom_0_1+chromWidth_0+(pixelPerChrom_0_2/2)),c(y-10),labels = paste("score: ",f.score_2," entropy:",cohort_entropy_0_2),cex=1)

    }    # mean start end smaller than subset chrom centromer
    if(mean.pos > half.length$length){
      xtr <- "topleft"
      xtr2 <- "topright"
      xtf <- c(21.5,5,4,22)
      xtf2 <- c(21.5,24,4,3)
      text(c(pixelPerChrom_0_1/2),c(10),labels = paste("score: ",f.score_1," entropy:",cohort_entropy_0_1),cex=1)
      text(c(pixelPerChrom_0_1+chromWidth_0+(pixelPerChrom_0_2/2)),c(10),labels = paste("score: ",f.score_2," entropy:",cohort_entropy_0_2),cex=1)

    }



    # may change here
    # legend type decision ----------------------------------------------------------------------------

    pixelPerChrom <- chromWidth_0+pixelPerChrom_0_1+pixelPerChrom_0_2+10

    x1 = 0.05
    x2 = pixelPerChrom_0_1/pixelPerChrom
    x3 = (pixelPerChrom_0_1+20+10)/pixelPerChrom
    x4 = 0.95

    d1 = x2-x1
    d2 = x4-x3
    if(d2>d1){
      x4 = x4 - (d2-d1)
    }else{
      x1 = x1 + (d1-d2)
    }


    if(mean.pos < half.length$length) {
      y2 = 0.4
      y1 = 0.2
    }    # mean start end smaller than subset chrom centromer
    if(mean.pos > half.length$length){
      y2 = 0.8
      y1 = 0.6
    }

    sub.position1 <- c(0.1,0.3,y1,y2)
    sub.position2 <- c(0.75,0.95, y1, y2)

    if(pixelPerChrom_0_1/pixelPerChrom<0.25){sub.position1 <- c(0.5,0.7,y1,y2)}
    if(pixelPerChrom_0_2/pixelPerChrom<0.25){sub.position2 <- c(0.3,0.5,y1,y2)}
    if(pixelPerChrom_0_1/pixelPerChrom<0.25 & pixelPerChrom_0_2/pixelPerChrom<0.25){
      sub.position1 <- c(0.1,0.3,y1,y2)
      sub.position2 <- c(0.75,0.95, y1, y2)
    }


    if(legend.type=="normal"){
      legend(xtr,legend=labs,col=color,cex=0.75,pch=16) # normal legend
      legend(xtr2,legend=labs,col=color,cex=0.75,pch=16)
      print("normal legend.")
    }
    if(legend.type =="pie") {
      #par(new=T,mar=xtf )
      par(fig = sub.position1 , mar=c(0,0,1,1), new=TRUE)
      factor_1 <- droplevels.factor(factor_1, exclude = if(anyNA(levels(factor_1)))NULL else NA)
      factor_2 <- droplevels.factor(factor_2, exclude = if(anyNA(levels(factor_2)))NULL else NA)

      dtt <-dtt[order(dtt$names),]

      dtt_1 <- dtt[dtt$names%in%levels(factor_1),]
      color_1 <- as.vector(dtt_1$color)
      labs_1 <- dtt_1$score
      freq <- data.frame(table(factor_1)/sum(table(factor_1)))
      if(nrow(freq[freq$Freq<0.05,])!=0){
        freq[freq$Freq<0.05,]$factor_1 <- NA
      }
      labs_1s <- freq$factor_1
      showtop = TRUE # here we only show the pie chart with more than 5% events
      if(showtop==TRUE){
        pie(table(factor_1),labels=labs_1s,col=color_1,cex=0.85,radius = 0.9) # piechart legend
      }else{
        pie(table(factor_1),labels=labs_1,col=color_1,cex=0.85,radius = 0.9) # piechart legend
      }

      #par(new=T,mar=xtf2)
      par(fig = sub.position2 , mar=c(0,0,1,1), new=TRUE)
      dtt_2 <- dtt[dtt$names%in%levels(factor_2),]
      color_2 <- as.vector(dtt_2$color)
      labs_2 <- dtt_2$score

      freq <- data.frame(table(factor_2)/sum(table(factor_2)))
      if(nrow(freq[freq$Freq<0.05,])!=0){
        freq[freq$Freq<0.05,]$factor_2 <- NA
      }
      labs_2s <- freq$factor_2
      if(showtop==TRUE){
        pie(table(factor_2),labels=labs_2s,col=color_2,cex=0.85,radius = 0.9) # piechart legend
      }else{
        pie(table(factor_2),labels=labs_2,col=color_2,cex=0.85,radius = 0.9) # piechart legend
      }


      print("pie plot legend！")
    } else{} # no legend
  }

  dev.off()
  if(SaveAsObject==TRUE){
    img <- readTIFF("t2.tiff")
    g <- rasterGrob(img, interpolate=TRUE)
    file.remove("t2.tiff")
  }



  # plot parameters for unique plots
  plot.new()
  #png("t1.png",width = 1024,height=768,units = "px")

  # save to object or to a file
  if(SaveAsObject==TRUE){
    if(orient=="v"){
      tiff(file="t3.tiff", width=12, height=8,units="in", compression="lzw", res=150)}else{
        tiff(file="t3.tiff", width=8, height=12,units="in", compression="lzw", res=150)}
  }else{ # SaveAsObject==FALSE
    if(format=="tiff"){
      if(orient=="v"){
        tiff(file=paste0(path,"/mixedplot.tiff"), width=12, height=8,units="in", compression="lzw", res=150)}else{
          tiff(file=paste0(path,"/mixedplot.tiff"), width=8, height=12,units="in", compression="lzw", res=150)}
    }else{ # format = EPS
      if(orient=="v"){
        setEPS()
        postscript(file=paste0(path,"/mixedplot.tiff"),width=12,height=8)}else{
          setEPS()
          postscript(file=paste0(path,"/mixedplot.tiff"),width=8,height=12)
        }
    } # end else format = EPS
  }


  par(c(5,3,4,4))

  chroms_0_1 <- chroms_1[del_dup.index_1]
  cohort_0_1 <- cohort_1[del_dup.index_1]
  cohort_entropy_0_1 <- round(entropy(table(cohort_0_1),unit = "log2"),digits = 2)
  starts_0_1 <- starts_1[del_dup.index_1]
  ends_0_1 <- ends_1[del_dup.index_1]
  rescore_0_1 <- rescores_1[del_dup.index_1]
  sorting_0_1 <- order(ends_0_1 - starts_0_1,cohort_0_1)
  cnv.number_0_1 <-  length(chroms_0_1) # number of lines in input
  chromWidth_0_1 <- round((pixel.per.cnv * cnv.number_0_1) * 0.1)
  cohort_0_1 <- droplevels.factor(cohort_0_1, exclude = if(anyNA(levels(cohort_0_1)))NULL else NA)

  chroms_0_2 <- chroms_2[del_dup.index_2]
  cohort_0_2 <- cohort_2[del_dup.index_2]
  cohort_entropy_0_2 <- round(entropy(table(cohort_0_2),unit = "log2"),digits = 2)
  starts_0_2 <- starts_2[del_dup.index_2]
  ends_0_2 <- ends_2[del_dup.index_2]
  rescore_0_2 <- rescores_2[del_dup.index_2]
  sorting_0_2 <- order(ends_0_2 - starts_0_2,cohort_0_2)
  cnv.number_0_2 <-  length(chroms_0_2) # number of lines in input
  chromWidth_0_2 <- round((pixel.per.cnv * cnv.number_0_2) * 0.1)
  cohort_0_2 <- droplevels.factor(cohort_0_2, exclude = if(anyNA(levels(cohort_0_2)))NULL else NA)

  cnv.number_0 <- cnv.number_0_1 + cnv.number_0_2
  chromWidth_0 <- round((pixel.per.cnv * cnv.number_0) * 0.1)

  chromWidth_0 <- 10



  pixelPerChrom <- chromWidth_0 + (pixel.per.cnv)*(cnv.number_0+1)+10 # determines space between chromsomes
  # x.size <- pixelPerChrom
  # y.size <- y+100

  if(orient=="v"){
    x.size <- pixelPerChrom
    y.size <- y+100}else{
      y.size <- pixelPerChrom
      x.size <- y+100
    }


  #plot(c(0,x.size),c(0,y.size),type="n",xaxt="n",yaxt="n",xlab="CNVs",ylab="Chromosomal location",main=title)


  if(zoomed==TRUE){
    if(orient=="v"){
      y1 = min(y - t_gene_end_1,y-t_gene_end_2)
      y2 = max(y - t_gene_start_1,y-t_gene_start_2)
      plot(c(0,x.size),c(y1-1e7,y2+1e7),type="n",xaxt="n",yaxt="n",xlab="CNVs",ylab="Chromosomal location",main=title)
      lines(c(0,x.size), c(y-t_gene_start_1,y-t_gene_start_1), lty=3, lwd=1)
      lines(c(0,x.size), c(y-t_gene_end_1,y-t_gene_end_1), lty=3, lwd=1)
      lines(c(0,x.size), c(y-t_gene_start_2,y-t_gene_start_2), lty=3, lwd=1)
      lines(c(0,x.size), c(y-t_gene_end_2,y-t_gene_end_2), lty=3, lwd=1)
    }else{ # if orient is h
      y1 = max(t_gene_end_1,t_gene_end_2)
      y2 = min(t_gene_start_1,t_gene_start_2)
      plot(c(y2-1e7,y1+1e7),c(0,y.size),type="n",xaxt="n",yaxt="n",xlab="CNVs",ylab="Chromosomal location",main=title)
      lines(c(t_gene_start_1,t_gene_start_1), c(0,y.size),lty=3, lwd=1)
      lines(c(t_gene_end_2,t_gene_end_2), c(0,y.size),lty=3, lwd=1)
      lines(c(t_gene_start_2,t_gene_start_2), c(0,y.size),lty=3, lwd=1)
      lines(c(t_gene_end_1,t_gene_end_1), c(0,y.size),lty=3, lwd=1)
    }

  }else{
    if(orient=="v"){
      plot(c(0,x.size),c(0,y.size),type="n",xaxt="n",yaxt="n",xlab="CNVs",
           ylab="Chromosomal location",main=title)}else{
             plot(c(0,x.size),c(0,y.size),type="n",xaxt="n",yaxt="n",ylab="CNVs",
                  xlab="Chromosomal location",main=title)
           }
  }







  chrStr <- paste("chr",toString(chroms[1]))
  text(c((chromWidth_0/2)),c(0),labels=c(chrStr))
  if(gene.anno == TRUE)    ###added RT gene.anno arg
  {
    m <- mean(start.gene,end.gene)
    text(c(-1),c(y-m+(y*0.035)),labels=c(cnv.type_1),cex=0.5)
    paintCytobands(chroms[1],pos=c(chromWidth_0,y),units="bases",width=chromWidth_0-7,orientation="v",legend=FALSE)
    rect(7,y-start.gene,chromWidth_0-1,y-end.gene,col="gray50", border = "gray50")
    lines(c(0.6,2.25),c(y-m+(y*0.02),y-m),col="gray50")
    lines(c(2.25,5),c(y-m,y-m),col="gray50")
  }else{
    #paintCytobands(chroms_1[1],pos=c(chromWidth_0,y),units="bases",width=chromWidth_0,orientation="v",legend=FALSE)
    if(orient=="v"){
      paintCytobands(chroms_1[1],pos=c(chromWidth_0,y),units="bases",bands="major",width=chromWidth_0,orientation="v",legend=FALSE)}else{
        paintCytobands(chroms_1[1],pos=c(0,chromWidth_0),units="bases",bands="major",width=chromWidth_0,orientation="h",legend=FALSE)
      }

  }
  chroms <- c(chroms_0_1,chroms_0_2)
  starts <- c(starts_0_1,starts_0_2)
  ends <- c(ends_0_1,ends_0_2)
  score <- c(rescore_0_1,rescore_0_2)
  cohorts <- c(cohort_0_1,cohort_0_2)
  startPoint <- chromWidth_0
  method <- "repeat"
  repeats <- c(repeat_0_1,repeat_0_2)

  plotCnv.cohort(chroms,starts,ends,y,
                     chromWidth=chromWidth_0,pixel.per.cnv=pixel.per.cnv,score=score,
                     cohorts=cohorts,startPoint=chromWidth_0,method="repeat",color=color,rep=repeats,
                     orient = orient)



  if(zoomed==FALSE&orient=="v"){

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
      text(c(pixelPerChrom/2),c(y-10),labels = paste("score: ",f.score_1,",",f.score_2,"; entropy: ",cohort_entropy_0_1,",",cohort_entropy_0_2),cex=1.2) # score on opposite
    }

    if(mean.pos > half.length$length){
      xtr <- "topright"
      xtf <- c(21.5,24,4,3)
      text(c(pixelPerChrom/2),c(10),labels = paste("score: ",f.score_1,",",f.score_2,"; entropy: ",cohort_entropy_0_1,",",cohort_entropy_0_2),cex=1.2)
    }


    # legend type decision ----------------------------------------------------------------------------

    if(2==2 || legend=="pie"){
      legend.color <- c("red","black","blue")
      par(new=T,mar=xtf )
      repu <- repeats!="C2"
      repeats_common <- repeats[repu]
      repeats_common <- factor(repeats_common,levels = c("U1","C1","U2"))
      pie(table(repeats_common),col=legend.color,cex=1,labels = c(gene.name_1,"both",gene.name_2))

    }
  }

  dev.off()

  if(SaveAsObject==TRUE){
    img <- readTIFF("t3.tiff")
    h <- rasterGrob(img, interpolate=TRUE)
    file.remove("t3.tiff")
  }

  if(SaveAsObject==TRUE){
    results <- list(g,h)}else{
      result <- "file saved!"
    }

  return(results)

}
