plot_multipanel_single <- function(paralist,font.size.factor,orient){


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

  t_gene_start <- unlist(paralist["t_gene_start"])
  t_gene_end <- unlist(paralist["t_gene_end"])
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
    del_dup.index <- rescore>=3
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


  #######
  del.index <- rescore<3
  dup.index <- rescore>=3
  cohort_1 <- cohorts[del.index]
  cohort_2 <- cohorts[dup.index]

  cohort_entropy_del <- round(entropy(table(cohort_1),unit = "log2"),digits=2)
  cohort_entropy_dup <- round(entropy(table(cohort_2),unit = "log2"),digits = 2)


  #####

  # tiff(file="multipanel_single.tiff", width=20, height=8,units="in", compression="lzw", res=150)


  ########################################################

  layout(mat = matrix(c(1,2,3,4,4,5), ncol = 3, byrow=T),heights = c(3, 4))
  #layout(mat = matrix(c(1,2,3,4,4,4), ncol = 3, byrow=T),heights = c(3, 4))

  ######################################################## data ready

  par(c(5,3,4,4))
  plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
  text(x = 0.5, y = 0.5, paste0("gene name: ",gene_name," \n",
                                "location: "," chr",chrom[1],":",t_gene_start,"-",t_gene_end,"\n",
                                "focality score of deletions: ",f.score,"\n",
                                #"focality score of deletions: ",f.score," (quantile)\n",
                                #"Casino score: NA (quantile)\n",
                                "duplication entropy:",cohort_entropy_dup,"\n",
                                "deleltion entropy: ",cohort_entropy_del,"\n",
                                "focal defininition threshold: 10Mbp"),
       cex = 1.6, col = "black")


  #############################################################################


  # both side #
  # plot parameters -----------------------------------------------------------------------------------------------------------------


  del.index <- rescore<3
  dup.index <- rescore>=3

  chroms_1 <- chroms[del.index]
  cnv.type_1 <- cnv.type
  cohort_1 <- cohorts[del.index]
  starts_1 <- starts[del.index]
  ends_1 <- ends[del.index]
  rescores_1 <- rescore[del.index]
  sorting_1 <- order(ends_1 - starts_1,cohort_1)
  chroms_2 <- chroms[dup.index]
  cnv.type_2 <- cnv.type
  cohort_2 <- cohorts[dup.index]
  starts_2 <- starts[dup.index]
  ends_2 <- ends[dup.index]
  rescores_2 <- rescore[dup.index]
  sorting_2 <- order(ends_2 - starts_2,cohort_2)
  cnv.number_d <-  length(chroms_1)+length(chroms_2) # number of lines in input
  chromWidth_d <- round((pixel.per.cnv * cnv.number_d) * 0.1)

  cohort_max <- sort(unique(c(levels(cohort_1),levels(cohort_2))))
  color.value <- GetColor(method="ploidy",cohorts=cohort_max)







  ######################################################## data ready
  # legend parameters ------------------------------------------------------------------------------------------------------


  cohort.dim <- length(cohort_max)
  color.value <- GetColor(method=color.method,cohorts=cohort_max)
  if(color.method == "cohort"){
    df.color.cohort <- data.frame(color=color.value, # colors according to getColor.ploidy/2
                                  score=cohort_max, # score according to ??
                                  names=cohort_max)
    dtt <- df.color.cohort
    ploidy_levels <- c("bi-del","mo-del","CN<5","4<CN<9","CN>8")
    factor_1 <- cohort_1
    factor_2 <- cohort_2
  }else if(color.method == "ploidy"){
    df.color.ploidy <- data.frame(color=color.value, # colors according to getColor.ploidy/2
                                  score=c(1,2,3,4,5), # score according to ??
                                  names=c("bi-del","mo-del","CN<5","4<CN<9","CN>8")
    )
    dtt <- df.color.ploidy
    ploidy_levels <- c("bi-del","mo-del","CN<5","4<CN<9","CN>8")
    factor_1 <- ploidy_levels[rescores_1]
    factor_2 <- ploidy_levels[rescores_2]
  }




  #############################################################################



  par(c(5,3,4,4))
  pixelPerChrom <- chromWidth_0 + (pixel.per.cnv)*(cnv.number_0+1)+10 # determines space between chromsomes

  #x.size <- pixelPerChrom
  #y.size <- y+100

  if(orient=="v"){
    x.size <- pixelPerChrom
    y.size <- y+100}else{
      y.size <- pixelPerChrom
      x.size <- y+100
    }



  nsample <- length(chroms_0)
  if(cnv.type=="dup"){
    tcnv = "duplication"
  }else{
    tcnv = "deletion"
  }
  ncohort <- length(unique(cohort_0))
  #print(cohort_0)

  #print(cohorts)
  title <- paste0(gene_name,": ",nsample," ",tcnv," events from ",ncohort," cohorts")


  if(orient=="v"){plot(c(0,x.size),c(0,y.size),type="n",xaxt="n",yaxt="n",xlab="CNVs",ylab="Chromosomal location",main=title,cex.lab=1.5)}else{
    plot(c(0,x.size),c(0,y.size),type="n",xaxt="n",yaxt="n",xlab="Chromosomal location",ylab="CNVs",main=title,cex.lab=1.5)
  }

  chrStr <- paste("chr",toString(chroms_0[1]))

  if(gene.anno == TRUE)    ###added RT gene.anno arg
  {
    m <- mean(start.gene,end.gene)
    text(c(-1),c(y-m+(y*0.035)),labels=c(cnv.type),cex=0.5)
    paintCytobands(chroms_0[1],pos=c(chromWidth_0,y),units="bases",width=chromWidth_0-7,orientation="v",legend=FALSE)
    rect(7,y-start.gene,chromWidth_0-1,y-end.gene,col="gray50", border = "gray50")
    lines(c(0.6,2.25),c(y-m+(y*0.02),y-m),col="gray50")
    lines(c(2.25,5),c(y-m,y-m),col="gray50")
  }else{

    if(orient=="v"){
      paintCytobands(chroms_0[1],pos=c(chromWidth_0,y),units="bases",width=chromWidth_0,orientation="v",legend=FALSE)}else{
        paintCytobands(chroms_0[1],pos=c(0,chromWidth_0),units="bases",width=chromWidth_0,orientation="h",legend=FALSE)
      }
  }

  plotCnv.cohort(chroms_0,starts_0,ends_0,y,
                     chromWidth=chromWidth_0,pixel.per.cnv=pixel.per.cnv,score=rescore_0,
                     cohorts=cohort_0,startPoint=chromWidth_0,method=color.method,orient = orient)


  # legend position decision (top or bottom)
  centro <- c(125,93.3,91,50.4,48.4,61,59.9,45.6,49,40.2,53.7,35.8,17.9,17.6,19,36.6,24,17.2,26.5,27.5,13.2,14.7,60.6,12.5)*1000000
  length <- lengthChromosome(c(1:22,"X","Y"),"bases")/2
  genome <- data.frame(chromosome=c(1:22,"X","Y"), centromere=centro, length=length) # dataframe containing chromosome and centromere position info

  mean.pos <- mean(c(starts_0,ends_0)) # mean position of all CNV´s
  #centroo <- genome[genome$chromosome %in% chroms,] # centromere position in the current chromosome
  half.length <- genome[genome$chromosome %in% chroms_0,] # half of the length of the current chromosome


  # legend type decision ----------------------------------------------------------------------------
  if(color.method=="cohort" | color.method=="length"){
    legend.color <- GetColor(method="cohort",color=color,cohorts=cohort_0)
  }
  if(color.method=="ploidy"){
    legend.color <- GetColor(method="ploidy",color=color,cohorts=cohort_0)
  }

  if(mean.pos < half.length$length) {
    sub.position <- c(.75, 1, .3, .7)
  }    # mean start end smaller than subset chrom centromer
  if(mean.pos > half.length$length){
    sub.position <- c(.75, 1, .6, 1)
  }

  if(color.method=="cohort" | color.method=="length"){
    if(legend==2 || legend=="pie"){

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
      showtop = TRUE
      if(showtop==TRUE){
        pie(table(factor_1),labels=labs_1s,col=color_1,cex=0.85*font.size.factor,radius = 0.9) # piechart legend
      }else{
        pie(table(factor_1),labels=labs_1,col=color_1,cex=0.85*font.size.factor,radius = 0.9) # piechart legend
      }


    }else{
      legend(xtr,legend=unique(cohort_0),col=legend.color,cex=0.75*font.size.factor,pch=16) # normal legend
    }
  }else if(color.method=="ploidy"){

    tb <- table(rescore_0)
    dp.list <-c("bi-del","mo-del","CN<5","4<CN<9","CN>8")
    for(i in 1:length(names(tb))){names(tb)[i] <- dp.list[as.integer(names(tb)[i])]}
    legend.color.subset <- legend.color[sort(unique(rescore_0))]
    if(legend==2 || legend=="pie"){
      pie(tb,col=legend.color.subset,cex=1*font.size.factor)
    }else{

      legend(xtr,legend=unique(dp.list),col=legend.color,cex=0.75*font.size.factor,pch=16)
    }
  }





  ############################################################################# second place here



  par(c(5,3,4,4))
  pixelPerChrom_1 <-  (pixel.per.cnv)*(length(chroms_1)+1)
  pixelPerChrom_2 <-  (pixel.per.cnv)*(length(chroms_2)+1)
  pixelPerChrom <- chromWidth_d+pixelPerChrom_1+pixelPerChrom_2+10 # determines space between chromsomes

  x.size <- pixelPerChrom
  y.size <- y+100
  ndel <- length(starts_1)
  ndup <- length(starts_2)
  title <- paste0(gene_name,": ",ndel," deletions and ",ndup," duplications")


  # modifies
  y1 = min(y - t_gene_end,y-max(ends_1),y-max(ends_2))
  y2 = max(y - t_gene_start,y-min(starts_1),y-min(starts_2))
  plot(c(0,x.size),c(y1*0.95,y2*1.05),type="n",xaxt="n",yaxt="n",xlab="CNVs",ylab="Chromosomal location",main=title,cex.lab=1.5)

  lines(c(0,x.size), c(y-t_gene_start,y-t_gene_start), lty=3, lwd=1)
  lines(c(0,x.size), c(y-t_gene_end,y-t_gene_end), lty=3, lwd=1)


  chrStr <- paste("chr",toString(chroms_1[1]))
  text(c(pixelPerChrom_1+(chromWidth_d/2)),c(0),labels=c(chrStr))


  paintCytobands(chroms_1[1],pos=c(pixelPerChrom_1+chromWidth_d,y),units="bases",width=chromWidth_d,orientation="v",legend=FALSE)
  m <- mean(t_gene_start,t_gene_end)

  if(sort.method=="cohort"){
    sorting_x_1 <- order(cohort_1,ends_1 - starts_1)
    sorting_x_2 <- order(cohort_2,ends_2 - starts_2)
  }else if(sort.method=="ploidy"){
    sorting_x_1 <- order(rescores_1,ends_1 - starts_1)
    sorting_x_2 <- order(rescores_2,ends_2 - starts_2)
  }else{
    sorting_x_1 <- order(ends_1 - starts_1)
    sorting_x_2 <- order(ends_2 - starts_2)
  }

  plotCnv(chroms_1,starts_1,ends_1,y,rescores_1,pixel.per.cnv=pixel.per.cnv,
          sorting = sorting_x_1, cohort = cohort_1,cohort_max = cohort_max,
          color.value = color.value,
          color.method=color.method,n=n_1,
          startPoint=(pixelPerChrom_1),direction = "left")

  plotCnv(chroms_2,starts_2,ends_2,y,rescores_2,pixel.per.cnv=pixel.per.cnv,
          sorting = sorting_x_2, cohort = cohort_2,cohort_max = cohort_max,
          color.value = color.value,
          color.method=color.method,n=n_2,
          startPoint=(pixelPerChrom_1+chromWidth_d),direction = "right")

  # legend("bottom",legend=ploidy_levels,col=color.value,cex=1*font.size.factor,
  #         pch=16,horiz=TRUE,box.lty=0,xjust=0.5,x.intersp = 0.3,y.intersp=0, bg="white")

  # here should be rewritten with windows

  #############################################################################

  color <- as.vector(dtt$color)
  labs <- as.vector(dtt$names)

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

  # barplot(test0,
  #         horiz=TRUE,
  #         xlab="cohorts", col=c("red","darkblue"),las=1,
  #         cex.main=0.5*font.size.factor,cex.axis = 0.5*font.size.factor,cex.names = 1*font.size.factor,
  #         beside=TRUE)
  # legend("center",legend=c("DEL","DUP"),col=c("red","darkblue"),
  #        cex=2*font.size.factor,pch=16,bty="n",text.width=3)


  test1 <- data.frame(t(test0))
  test1$sum <- test1$X1+test1$X2
  test2 <- data.frame(matrix(ncol = 1, nrow = 5))
  x <- c("bi-del","mo-del","CN<5","4<CN<9","CN>8")
  rownames(test2) <- x
  colnames(test2) <- "rank"
  test2$rank <- 1:nrow(test2)

  test1$rn <- rownames(test1)
  test2$rn <- rownames(test2)
  test3 <- merge(test1,test2,by="rn",all=TRUE)
  #test3$sum.x <- test3$sum.x + test3$sum.y
  test3[is.na(test3$sum)==TRUE,"sum"] <- 0
  test3<-test3[order(test3$rank),]
  test4 <- test3$sum
  names(test4) <- test3$rn

  color.value <- c("red2","indianred4","lightskyblue2","skyblue3","skyblue4")
  barplot(test4,
          #main="Deletions and Duplications",
          horiz=TRUE,
          xlab="cohorts",
          #col=c("red","darkblue"),
          col=color.value,
          las=1,
          cex.main=0.5*font.size.factor,cex.axis = 0.5*font.size.factor,cex.names = 0.5*font.size.factor,
          #legend = c("deletion","duplication"),
          beside=TRUE)




  print("pie plot legend！")

  #dev.off()

}
