plotCnv.cohort <- function(chroms,starts,ends,y,score,chromWidth,
                           pixel.per.cnv,cohorts,startPoint,color,method,rep){
  indX <- chroms == 'X'
  indY <- chroms == 'Y'
  len <- length(starts)
  if(missing(rep)){
    rep<-""
  }else if(is.null(rep) == TRUE){
    rep<-""
  }else{
    rep0 <- factor(rep,levels = c("U1","C1","C2","U2"))
    sorting_rep <- order(rep0,ends - starts)
    rep <- rep[sorting_rep]
    reindex <- rep!="C2"
    chroms <- chroms[sorting_rep][reindex]
    starts <- starts[sorting_rep][reindex]
    ends <- ends[sorting_rep][reindex]
    score <- score[sorting_rep][reindex]
    cohorts <- cohorts[sorting_rep][reindex]
    rep <- rep[reindex]

  }
  color.value <- GetColor(method=method,color=color,cohorts=cohorts)
  cohort.list <- sort(unique(cohorts))
  ploidy.list <- c("bi-del","mo-del","diploidy","gain-low","gain-mid","gain-high","n/a")
  ploidy.list <- 1:7
  repeat.list <- c("U1","C1","C2","U2")
  print(color.value)
  startPoint <- chromWidth

  if(method=="cohort"){
    class.list <- cohort.list
    class <- cohorts}
  if(method=="ploidy"){
    class.list <- ploidy.list
    class <- score}
  if(method=="repeat"){
    class.list <- repeat.list
    class <- rep
  }

  # Autosomes
  for(index in 1:len){
    if(method!="length"){
      #print(index)
      #print(ends[index]-starts[index])
      class.index <- match(class[index],class.list)
      #class.index <- class[sorting.color[index]]
      print(class.index)
      x <- startPoint + pixel.per.cnv*index
      lines(c(x,x),c(y-starts[index],y-ends[index]),
            col=color.value[class.index],lwd=pixel.per.cnv)
    }else{
      x <- startPoint + pixel.per.cnv*index
      lines(c(x,x),c(y-starts[index],y-ends[index]),
            col="black",lwd=pixel.per.cnv)
    }
  }

  # X Chromosome
  for(index in 1:len){
    if(indX[index] == TRUE){
      x <- startPoint +pixel.per.cnv*index
      lines(c(x,x),c(y-starts[index],y-ends[index]),col=cohorts[index],lwd=pixel.per.cnv)
    }
  }

  # Y Chromosome
  for(index in 1:len){
    if(indY[index] == TRUE){
      x <- startPoint + pixel.per.cnv*index
      lines(c(x,x),c(y-starts[index],y-ends[index]),col=cohorts[index],lwd=pixel.per.cnv)  ## RT cohorts for scores
    }
  }
}
