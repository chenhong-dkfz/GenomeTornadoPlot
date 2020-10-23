PlotTwinsInit <- function(twin.cnv,sort.method,color.method,
                          out.dir,title,
                          legend.type,legend.names,color,score.values_1,score.values_2,
                          gene.anno,cnv.type_1,cnv.type_2
){
  CNV_1 <- twin.cnv@matrix_1
  CNV_2 <- twin.cnv@matrix_2


  t_gene_start_1 = twin.cnv@t_gene_start_1
  t_gene_end_1 = twin.cnv@t_gene_end_1
  t_gene_start_2 = twin.cnv@t_gene_start_2
  t_gene_end_2 = twin.cnv@t_gene_end_2


  f.score_1 <- twin.cnv@gene_score_1
  f.score_2 <- twin.cnv@gene_score_2

  print(sort.method)

  chrom_1 <- as.vector(seqnames(CNV_1))
  start.CNV_1 <- start(CNV_1)
  end.CNV_1 <- end(CNV_1)
  score_1 <- CNV_1$Score
  cohort_1 <- CNV_1$Cohort
  pids_1 <- CNV_1$PID
  rep_1 <- CNV_1$rep


  chrom_2 <- as.vector(seqnames(CNV_2))
  start.CNV_2 <- start(CNV_2)
  end.CNV_2 <- end(CNV_2)
  score_2 <- CNV_2$Score
  cohort_2 <- CNV_2$Cohort
  pids_2 <- CNV_2$PID
  rep_2 <- CNV_2$rep

  max.length = twin.cnv@max.length





  gene.name_1 <- as.character(twin.cnv@gene_name_1)
  gene.name_2 <- as.character(twin.cnv@gene_name_2)

  if(missing(sort.method) & missing(color.method)){
    sort.method = "length"
    color.method = "length"
  }else{
    if(missing(sort.method)){
      sort.method = color.method
    }else if(missing(color.method)){
      color.method = sort.method
    }
  }



  index_1 <- (end.CNV_1 - start.CNV_1) < max.length # only events shorter than 10 M
  m_1 <- sum(index_1)
  startPos_1 <- start.CNV_1[index_1]
  endPos_1 <- end.CNV_1[index_1]
  score_1 <- score_1[index_1]
  rescore_1 <- unlist(lapply(score_1,MapPloidyClasses))
  score.values_1 <- as.character(sort(unique(rescore_1)))
  n_1 <- length(unique(rescore_1))

  index_2 <- (end.CNV_2 - start.CNV_2) < max.length # only events shorter than 10 M
  m_2 <- sum(index_2)
  startPos_2 <- start.CNV_2[index_2]
  endPos_2 <- end.CNV_2[index_2]
  score_2 <- score_2[index_2]
  rescore_2 <- unlist(lapply(score_2,MapPloidyClasses))
  score.values_2 <- as.character(sort(unique(rescore_2)))
  n_2 <- length(unique(rescore_2))

  repeat_1 <- rep_1[index_1]
  repeat_2 <- rep_2[index_2]

  cohort_1 <- cohort_1[index_1]
  cohort_2 <- cohort_2[index_2]
  ## default/optional parameter for gene name (default = "geneX")-------------------------------------------------------------------------------------------------------

  if(missing(cnv.type_1)){cnv.type_1 <- "del"}
  if(missing(cnv.type_2)){cnv.type_2 <- "del"}



  if(missing(title)){
    if(cnv.type_1 == "del"){
    title <- paste(gene.name_1,"&",gene.name_2,": deletions")
    }else if(cnv.type_1 == "dup"){
      title <- paste(gene.name_1,"&",gene.name_2,": duplications")
    }else{
      title <- paste(gene.name_1,"&",gene.name_2,": all CNVs")
    }
  }else{
    title <- title
  }



  ## default/optional parameter for legend (default = "pie" ,otherwise processed to "normal legend")
  if(missing(legend.type)){
    legend.type <- "pie"
  }

  ## default/optional parameter for legend.names
  if(missing(legend.names)){
    legend.names_1 <- c(score.values_1[1:n_1])
    legend.names_2 <- c(score.values_2[1:n_2])
    legend.names <- unique(c(legend.names_1,legend.names_2))
  }else{
    nmax <- max(n_1,n_2)
    legend.names <- c(legend.names)[1:nmax]
  }

  ## default/optional parameter for out.dir (defaultdirectory = "/package/TornadoCNV") --------------------------------------------------------------------------------
  if(missing(out.dir)){
    out.dir <- paste("twins.by.length","_",gene.name_1,"&",gene.name_2,".","png",sep = "")
  }else{
    out.dir <- paste(out.dir,"_",gene.name_1,"&",gene.name_2,".","png",sep="")
  }


  ## default/optional parameter for pixel.per.cnv (default = 5)---------------------------------------------------------------------------------------------------------
  #if(missing(pixel.per.cnv)){pixel.per.cnv <- 200/(m_1+m_2)}  ## better a equation dependened on the number of CNVs (index!)
  pixel.per.cnv <- 200/(m_1+m_2)
  ## sorting ----------------------------------------------------------------------------------------------------------------------------------------------------------
  if(sort.method=="length"){
    sorting_1 <- order(endPos_1 - startPos_1) # sort by length #### add ??
    sorting_2 <- order(endPos_2 - startPos_2)
  }
  if(sort.method=="ploidy"){
    sorting_1 <- order(rescore_1,endPos_1 - startPos_1)
    sorting_2 <- order(rescore_2,endPos_2 - startPos_2)
  }
  if(sort.method=="cohort"){
    sorting_1 <- order(cohort_1,endPos_1 - startPos_1)
    sorting_2 <- order(cohort_2,endPos_2 - startPos_2)
  }


  ## color ------------------------------------------------------------------------------------------------------------------------------------------------------------
  if(missing(color)){
    color <- "steelblue3"
  }else{
    color <- color
  }

  ## gene anno?----------------------------------------------------------------------------------------------------------------------------------------------------
  if(missing(gene.anno)){
    gene.anno <- ""
  }


  paralist <- list("gene.name_1"=gene.name_1,"gene.name_2"=gene.name_2,
                   "cnv.type_1"=cnv.type_1,"cnv.type_2"=cnv.type_2,
                   "title"=title,"legend.type"=legend.type,
                   "cohort_1"=cohort_1,"cohort_2"=cohort_2,
                   "legend.names"=legend.names,
                   "out.dir"=out.dir,"pixel.per.cnv"=pixel.per.cnv,
                   "color"=color,"sorting_1"=sorting_1,"sorting_2"=sorting_2,
                   #"start.gene_1"=start.gene_1,"end.gene_1"=end.gene_1,"start.gene_2"=start.gene_2,"end.gene_2"=end.gene_2,
                   "gene.anno"=gene.anno,"n_1"=n_1,"n_2"=n_2,
                   "chrom_1"=chrom_1,"start.CNV_1"=start.CNV_1,"end.CNV_1"=end.CNV_1,"rescore_1"=rescore_1,
                   "chrom_2"=chrom_2,"start.CNV_2"=start.CNV_2,"end.CNV_2"=end.CNV_2,"rescore_2"=rescore_2,
                   "index_1"=index_1,"index_2"=index_2,"m_1"=m_1,"m_2"=m_2,
                   "sort.method"=sort.method,"color.method"=color.method,
                   "startPos_1"=startPos_1,"endPos_1"=endPos_1,
                   "startPos_2"=startPos_2,"endPos_2"=endPos_2,
                   "repeat_1"=repeat_1,"repeat_2"=repeat_2,
                   "f.score_1"=f.score_1,"f.score_2"=f.score_2,
                   "t_gene_start_1" = t_gene_start_1,"t_gene_end_1" = t_gene_end_1,
                   "t_gene_start_2" = t_gene_start_2,"t_gene_end_2" = t_gene_end_2

  )


}
