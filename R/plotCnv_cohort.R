
plotCnv.cohort <- function(chroms,starts,ends,y,score,chromWidth,pixel.per.cnv,cohorts,startPoint,color,method){
  indX <- chroms == 'X'
  indY <- chroms == 'Y'
  len <- length(starts)

  color.value <- GetColor(method=method,color=color,cohorts=cohorts)
  cohort.list <- sort(unique(cohorts))
  ploidy.list <- c("bi-del","mo-del","diploidy","gain-low","gain-mid","gain-high","n/a")
  ploidy.list <- 1:7


  if(method=="cohort"){
    class.list <- cohort.list
    class <- cohorts}
  if(method=="ploidy"){
    class.list <- ploidy.list
    class <- score
  }

  # Autosomes
  for(index in 1:len){
    if(method!="length"){
      class.index <- match(class[index],class.list)
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
