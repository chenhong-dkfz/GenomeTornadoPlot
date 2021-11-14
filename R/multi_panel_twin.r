# it is also a half-finished version
plot_multipanel_twin <- function(paralist,font.size.factor,orient){

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
  #print("cohort size:")
  #print(cohort_1)
  #print(cohort_2)

  repeat_1 <- unlist(paralist["repeat_1"])
  repeat_2 <- unlist(paralist["repeat_2"])


  gene.name_1 <- unlist(paralist["gene.name_1"])
  gene.name_2 <- unlist(paralist["gene.name_2"])

  t_gene_start_1 <- unlist(paralist["t_gene_start_1"])
  t_gene_end_1 <- unlist(paralist["t_gene_end_1"])
  t_gene_start_2 <- unlist(paralist["t_gene_start_2"])
  t_gene_end_2 <- unlist(paralist["t_gene_end_2"])

  n_1 <- unlist(paralist["n_1"])
  n_2 <- unlist(paralist["n_2"])



  sort.method = unlist(paralist["sort.method"])
  color.method = unlist(paralist["color.method"])
  drop.low.amp <- unlist(paralist["drop.low.amp"])

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





  ######end addon############

#  tiff(file="multipanel_twin.tiff", width=25, height=8,units="in", compression="lzw", res=150)

  # the starts and text

  layout(mat = matrix(c(1,1,2,2,3,3,4,5,5,5,5,6), ncol = 6, byrow=T),heights = c(3, 4))
  #layout(mat = matrix(c(1,2,3,4,4,4), ncol = 3, byrow=T),heights = c(3, 4))



  ######## entropy #########

if(drop.low.amp==FALSE){
  del_dup.index_1 <- rescores_1>=3
  cohort_0_1 <- cohort_1[del_dup.index_1]
  cohort_entropy_dup_1 <- round(entropy(table(cohort_0_1),unit = "log2"),digits = 2)


  del_dup.index_2 <- rescores_2>=3
  cohort_0_2 <- cohort_2[del_dup.index_2]
  cohort_entropy_dup_2 <- round(entropy(table(cohort_0_2),unit = "log2"),digits = 2)


  del_dup.index_1 <- rescores_1<3
  cohort_0_1 <- cohort_1[del_dup.index_1]
  cohort_entropy_del_1 <- round(entropy(table(cohort_0_1),unit = "log2"),digits = 2)


  del_dup.index_2 <- rescores_2<3
  cohort_0_2 <- cohort_2[del_dup.index_2]
  cohort_entropy_del_2 <- round(entropy(table(cohort_0_2),unit = "log2"),digits = 2)
}else{
  del_dup.index_1 <- rescores_1>3
  cohort_0_1 <- cohort_1[del_dup.index_1]
  cohort_entropy_dup_1 <- round(entropy(table(cohort_0_1),unit = "log2"),digits = 2)


  del_dup.index_2 <- rescores_2>3
  cohort_0_2 <- cohort_2[del_dup.index_2]
  cohort_entropy_dup_2 <- round(entropy(table(cohort_0_2),unit = "log2"),digits = 2)


  del_dup.index_1 <- rescores_1<3
  cohort_0_1 <- cohort_1[del_dup.index_1]
  cohort_entropy_del_1 <- round(entropy(table(cohort_0_1),unit = "log2"),digits = 2)


  del_dup.index_2 <- rescores_2<3
  cohort_0_2 <- cohort_2[del_dup.index_2]
  cohort_entropy_del_2 <- round(entropy(table(cohort_0_2),unit = "log2"),digits = 2)
}



  ######################################################## data ready

  par(c(5,3,4,4))
  plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
  #text(x = 0.5, y = 0.5, paste0("what do you say?"),
  #     cex = 1.6, col = "black")
  text(x = 0.5, y = 0.5, paste0(gene.name_1," (chr",chroms_1[1],":",t_gene_start_1,"-",t_gene_end_1,")\n",
                                "focality score of deletions: ",f.score_1,"\n",
                                #"Casino score: NA (quantile)\n",
                                "duplication entropy:",cohort_entropy_dup_1,"\n",
                                "deleltion entropy: ",cohort_entropy_del_1,"\n",
                                gene.name_2," (chr",chroms_1[1],":",t_gene_start_2,"-",t_gene_end_2,")\n",
                                "focality score of deletions: ",f.score_2,"\n",
                                #"Casino score: NA (quantile)\n",
                                "duplication entropy:",cohort_entropy_dup_2,"\n",
                                "deleltion entropy: ",cohort_entropy_del_2,"\n",
                                "focal defininition threshold: 10Mb"),
       cex = 1.6, col = "black")



  ##################### the first plot
  #######add on#######

  cnv.type_1 = "dup"
  cnv.type_2 = "dup"


  # if(cnv.type_1=="dup"){
  #   del_dup.index_1 <- rescores_1>=3
  # }else if(cnv.type_1=="del") {
  #   del_dup.index_1 <- rescores_1<3
  # }else{
  #   del_dup.index_1 <- rescores_1>=0
  # }

  if(drop.low.amp==FALSE){
    if(cnv.type_1=="dup"){
      del_dup.index_1 <- rescores_1>=3
    }else if(cnv.type_1=="del") {
      del_dup.index_1 <- rescores_1<3
    }else{
      del_dup.index_1 <- rescores_1>=0
    }
  }else{
    if(cnv.type_1=="dup"){
      del_dup.index_1 <- rescores_1>3
    }else if(cnv.type_1=="del") {
      del_dup.index_1 <- rescores_1<3
    }else{
      del_dup.index_1 <- rescores_1>=0
    }
  }


  chroms_0_1 <- chroms_1[del_dup.index_1]
  cohort_0_1 <- cohort_1[del_dup.index_1]
  cohort_entropy_dup_1 <- round(entropy(table(cohort_0_1),unit = "log2"),digits = 2)
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

#
#   if(cnv.type_2=="dup"){
#     del_dup.index_2 <- rescores_2>=3
#   }else if(cnv.type_2=="del") {
#     del_dup.index_2 <- rescores_2<3
#   }else{
#     del_dup.index_2 <- rescores_2>=0
#   }
#
  if(drop.low.amp==FALSE){
    if(cnv.type_2=="dup"){
      del_dup.index_2 <- rescores_2>=3
    }else if(cnv.type_2=="del") {
      del_dup.index_2 <- rescores_2<3
    }else{
      del_dup.index_2 <- rescores_2>=0
    }
  }else{
    if(cnv.type_2=="dup"){
      del_dup.index_2 <- rescores_2>3
    }else if(cnv.type_2=="del") {
      del_dup.index_2 <- rescores_2<3
    }else{
      del_dup.index_2 <- rescores_2>=0
    }
  }



  chroms_0_2 <- chroms_2[del_dup.index_2]
  cohort_0_2 <- cohort_2[del_dup.index_2]
  cohort_entropy_dup_2 <- round(entropy(table(cohort_0_2),unit = "log2"),digits = 2)
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

  par(c(5,3,4,4))

  chroms_0_1 <- chroms_1[del_dup.index_1]
  cohort_0_1 <- cohort_1[del_dup.index_1]
  starts_0_1 <- starts_1[del_dup.index_1]
  ends_0_1 <- ends_1[del_dup.index_1]
  rescore_0_1 <- rescores_1[del_dup.index_1]
  sorting_0_1 <- order(ends_0_1 - starts_0_1,cohort_0_1)
  cnv.number_0_1 <-  length(chroms_0_1) # number of lines in input
  chromWidth_0_1 <- round((pixel.per.cnv * cnv.number_0_1) * 0.1)
  cohort_0_1 <- droplevels.factor(cohort_0_1, exclude = if(anyNA(levels(cohort_0_1)))NULL else NA)

  chroms_0_2 <- chroms_2[del_dup.index_2]
  cohort_0_2 <- cohort_2[del_dup.index_2]
  starts_0_2 <- starts_2[del_dup.index_2]
  ends_0_2 <- ends_2[del_dup.index_2]
  rescore_0_2 <- rescores_2[del_dup.index_2]
  sorting_0_2 <- order(ends_0_2 - starts_0_2,cohort_0_2)
  cnv.number_0_2 <-  length(chroms_0_2) # number of lines in input
  chromWidth_0_2 <- round((pixel.per.cnv * cnv.number_0_2) * 0.1)
  cohort_0_2 <- droplevels.factor(cohort_0_2, exclude = if(anyNA(levels(cohort_0_2)))NULL else NA)

  cnv.number_0 <- cnv.number_0_1 + cnv.number_0_2 - table(repeat_0_1)[1]
  chromWidth_0 <- round((pixel.per.cnv * cnv.number_0) * 0.1)


  pixelPerChrom <- chromWidth_0 + (pixel.per.cnv)*(cnv.number_0+1)+10 # determines space between chromsomes
  x.size <- pixelPerChrom
  y.size <- y+100
  plot(c(0,x.size),c(0,y.size),type="n",xaxt="n",yaxt="n",xlab="CNVs",ylab="Chromosomal location",main=paste0(gene.name_1," & ",gene.name_2,": duplications"))
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
    paintCytobands(chroms_1[1],pos=c(chromWidth_0,y),units="bases",width=chromWidth_0,orientation="v",legend=FALSE)
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
                 orient=orient)


  legend.color <- c("red","grey","blue")
  legend("bottom",legend=c(gene.name_1,"both",gene.name_2),col=legend.color,cex=0.75,pch=16,horiz=TRUE,box.lty=0)


  ################################# the second plot
  #######add on#######

  cnv.type_1 = "del"
  cnv.type_2 = "del"

if(drop.low.amp==FALSE){
  if(cnv.type_1=="dup"){
    del_dup.index_1 <- rescores_1>=3
  }else if(cnv.type_1=="del") {
    del_dup.index_1 <- rescores_1<3
  }else{
    del_dup.index_1 <- rescores_1>=0
  }
}else{
  if(cnv.type_1=="dup"){
    del_dup.index_1 <- rescores_1>3
  }else if(cnv.type_1=="del") {
    del_dup.index_1 <- rescores_1<3
  }else{
    del_dup.index_1 <- rescores_1>=0
  }
}

  chroms_0_1 <- chroms_1[del_dup.index_1]
  cohort_0_1 <- cohort_1[del_dup.index_1]
  cohort_entropy_del_1 <- round(entropy(table(cohort_0_1),unit = "log2"),digits = 2)
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

if(drop.low.amp==FALSE){
  if(cnv.type_2=="dup"){
    del_dup.index_2 <- rescores_2>=3
  }else if(cnv.type_2=="del") {
    del_dup.index_2 <- rescores_2<3
  }else{
    del_dup.index_2 <- rescores_2>=0
  }
}else{
  if(cnv.type_2=="dup"){
    del_dup.index_2 <- rescores_2>3
  }else if(cnv.type_2=="del") {
    del_dup.index_2 <- rescores_2<3
  }else{
    del_dup.index_2 <- rescores_2>=0
  }
}

  chroms_0_2 <- chroms_2[del_dup.index_2]
  cohort_0_2 <- cohort_2[del_dup.index_2]
  cohort_entropy_del_2 <- round(entropy(table(cohort_0_2),unit = "log2"),digits = 2)
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

  par(c(5,3,4,4))

  chroms_0_1 <- chroms_1[del_dup.index_1]
  cohort_0_1 <- cohort_1[del_dup.index_1]
  starts_0_1 <- starts_1[del_dup.index_1]
  ends_0_1 <- ends_1[del_dup.index_1]
  rescore_0_1 <- rescores_1[del_dup.index_1]
  sorting_0_1 <- order(ends_0_1 - starts_0_1,cohort_0_1)
  cnv.number_0_1 <-  length(chroms_0_1) # number of lines in input
  chromWidth_0_1 <- round((pixel.per.cnv * cnv.number_0_1) * 0.1)
  cohort_0_1 <- droplevels.factor(cohort_0_1, exclude = if(anyNA(levels(cohort_0_1)))NULL else NA)

  chroms_0_2 <- chroms_2[del_dup.index_2]
  cohort_0_2 <- cohort_2[del_dup.index_2]
  starts_0_2 <- starts_2[del_dup.index_2]
  ends_0_2 <- ends_2[del_dup.index_2]
  rescore_0_2 <- rescores_2[del_dup.index_2]
  sorting_0_2 <- order(ends_0_2 - starts_0_2,cohort_0_2)
  cnv.number_0_2 <-  length(chroms_0_2) # number of lines in input
  chromWidth_0_2 <- round((pixel.per.cnv * cnv.number_0_2) * 0.1)
  cohort_0_2 <- droplevels.factor(cohort_0_2, exclude = if(anyNA(levels(cohort_0_2)))NULL else NA)

  cnv.number_0 <- cnv.number_0_1 + cnv.number_0_2- table(repeat_0_1)[1]
  chromWidth_0 <- round((pixel.per.cnv * cnv.number_0) * 0.1)


  pixelPerChrom <- chromWidth_0 + (pixel.per.cnv)*(cnv.number_0+1)+10 # determines space between chromsomes
  x.size <- pixelPerChrom
  y.size <- y+100
  plot(c(0,x.size),c(0,y.size),type="n",xaxt="n",yaxt="n",xlab="CNVs",ylab="Chromosomal location",main=title)
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
    paintCytobands(chroms_1[1],pos=c(chromWidth_0,y),units="bases",width=chromWidth_0,orientation="v",legend=FALSE)
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
                 orient=orient)


  legend.color <- c("red","grey","blue")
  legend("bottom",legend=c(gene.name_1,"both",gene.name_2),col=legend.color,cex=0.75,pch=16,horiz=TRUE,box.lty=0)


  #############
  # pie plot 1


  cohort_max <- sort(unique(c(levels(cohort_0_1),levels(cohort_0_2))))
  color.value <- GetColor(method=color.method,cohorts=cohort_max)


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
    pie(table(factor_1),labels=labs_1s,col=color_1,cex=1,radius = 1.5) # piechart legend
  }else{
    pie(table(factor_1),labels=labs_1,col=color_1,cex=1,radius = 1.5) # piechart legend
  }







  ############################################################
  # the main plot

  par(c(5,3,4,4))
  cnv.number_0 <- cnv.number_0_1 + cnv.number_0_2
  chromWidth_0 <- round((pixel.per.cnv * cnv.number_0) * 0.1)
  pixelPerChrom_0_1 <-  (pixel.per.cnv)*(length(chroms_0_1)+1)
  pixelPerChrom_0_2 <-  (pixel.per.cnv)*(length(chroms_0_2)+1)
  pixelPerChrom_0 <- chromWidth_0+pixelPerChrom_0_1+pixelPerChrom_0_2+10 # determines space between chromsomes

  x.size <- pixelPerChrom_0
  y.size <- y+100


  t_gene_start_1 = as.numeric(as.character(t_gene_start_1))
  t_gene_end_1 = as.numeric(as.character(t_gene_end_1))
  t_gene_start_2 = as.numeric(as.character(t_gene_start_2))
  t_gene_end_2 = as.numeric(as.character(t_gene_end_2))

  # modifies
  y1 = min(y - t_gene_end_1,y - t_gene_end_2,y-max(ends_0_1),y-max(ends_0_2))
  y2 = max(y - t_gene_start_1,y - t_gene_start_2,y-min(starts_0_1),y-min(starts_0_2))
  plot(c(0,x.size),c(y1*0.95,y2*1.05),type="n",xaxt="n",yaxt="n",xlab="CNVs",ylab="Chromosomal location",main=title)

  lines(c(0,pixelPerChrom_0_1), c(y-t_gene_start_1,y-t_gene_start_1), lty=3, lwd=1)
  lines(c(0,pixelPerChrom_0_1), c(y-t_gene_end_1,y-t_gene_end_1), lty=3, lwd=1)

  lines(c(pixelPerChrom_0_1+chromWidth_0,x.size), c(y-t_gene_start_2,y-t_gene_start_2), lty=3, lwd=1)
  lines(c(pixelPerChrom_0_1+chromWidth_0,x.size), c(y-t_gene_end_2,y-t_gene_end_2), lty=3, lwd=1)


  #plot(c(0,x.size),c(0,y.size),type="n",xaxt="n",yaxt="n",xlab="CNVs",
  #     ylab="Chromosomal location",main=title)



  chrStr <- paste("chr",toString(chroms_0_1[1]))
  text(c(pixelPerChrom_0_1+(chromWidth_0/2)),c(0),labels=c(chrStr))


  if(gene.anno == TRUE){
    paintCytobands(chroms_0_1[1],pos=c(pixelPerChrom_0_1+chromWidth_0,y),units="bases",width=chromWidth_0,orientation="v",legend=FALSE)
    m_1 <- mean(c(t_gene_start_1,t_gene_end_1))
    m_2 <- mean(c(t_gene_start_2,t_gene_end_2))

    # annotation right side
    text(c(pixelPerChrom_0_1+chromWidth_0+15),c(y-m_1+(y*0.045)),labels=c(gene.name_2),cex=0.7)
    lines(c(pixelPerChrom_0_1+chromWidth_0+7,pixelPerChrom_0_1+chromWidth_0+4),c(y-m_1+(y*0.03),y-m_1),col="gray50")
    #lines(c(pixelPerChrom_1+chromWidth+4,pixelPerChrom_1+chromWidth+1),c(y-m_1,y-m_1),col="gray50")

    # annotation left side
    text(c(pixelPerChrom_0_1-15),c(y-m_2-(y*0.045)),labels=c(gene.name_1),cex=0.7)
    lines(c(pixelPerChrom_0_1-7,pixelPerChrom_0_1-4),c(y-m_2-(y*0.03),y-m_2),col="gray50")
    #lines(c(pixelPerChrom_1-4,pixelPerChrom_1-1),c(y-m_2,y-m_2),col="gray50")



  }else{ # acutally this is what really use
    paintCytobands(chroms_0_1[1],pos=c(pixelPerChrom_0_1+chromWidth_0,y),units="bases",bands="major",width=chromWidth_0,orientation="v",legend=FALSE)
  }




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

  # legend parameters ------------------------------------------------------------------------------------------------------


  # dashline




  ######### pie plot 2

  dtt_2 <- dtt[dtt$names%in%levels(factor_2),]
  color_2 <- as.vector(dtt_2$color)
  labs_2 <- dtt_2$score

  freq <- data.frame(table(factor_2)/sum(table(factor_2)))
  if(nrow(freq[freq$Freq<0.05,])!=0){
    freq[freq$Freq<0.05,]$factor_2 <- NA
  }
  labs_2s <- freq$factor_2
  if(showtop==TRUE){
    pie(table(factor_2),labels=labs_2s,col=color_2,cex=1,radius = 1.5) # piechart legend
  }else{
    pie(table(factor_2),labels=labs_2,col=color_2,cex=1,radius = 1.5) # piechart legend
  }

#dev.off()

}
