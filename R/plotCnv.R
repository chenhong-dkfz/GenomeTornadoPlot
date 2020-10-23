plotCnv <- function(chroms,starts,ends,y,scores,pixel.per.cnv,sorting,cohort,
                    cohort_max, color.value,
                    color.method,color,n,startPoint,direction){



  #if(direction=="left"){sorting=rev(sorting)}

  chroms <- chroms[sorting]
  starts <- starts[sorting]
  ends <- ends[sorting]
  cohorts <- cohort[sorting]
  scores <- scores[sorting]
  cnv.number <-  length(chroms) # number of lines in input



  indX <- chroms == 'X'
  indY <- chroms == 'Y'

  cohort.list <- cohort_max
  ploidy.list <- c("bi-del","mo-del","diploidy","gain-low","gain-mid","gain-high","n/a")
  ploidy.list <- 1:7

  if(color.method=="cohort"){
    class.list <- cohort.list
    class.cont <- cohorts}
  if(color.method=="ploidy"){
    class.list <- ploidy.list
    class.cont <- scores
  }

  len <- length(starts)
  if(missing(direction)){direction="right"}

  if(direction=="right"){
    for(index in 1:len){
      # Autosomes
      if(color.method!="length"){
        class.index <- match(class.cont[index],class.list)
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
   # for(index in 1:len) {
  #    if(indX[index] == TRUE){
  #      x <- startPoint + pixel.per.cnv*index
  #      lines(c(x,x),c(y-starts[index],y-ends[index]),col="black",lwd=pixel.per.cnv)
  #    }
  #  }

    # Y Chromosome
  #  for(index in 1:len){
  #    if(indY[index] == TRUE){
  #      x <- startPoint + pixel.per.cnv*index
  #      lines(c(x,x),c(y-starts[index],y-ends[index]),col="black",lwd=pixel.per.cnv)
  #    }
  #  }
  }else{ # if the lines are in the left side
    for(index in 1:len){
      # Autosomes
      if(color.method!="length"){
        class.index <- match(class.cont[index],class.list)
        x <- startPoint - pixel.per.cnv*index
        lines(c(x,x),c(y-starts[index],y-ends[index]),
              col=color.value[class.index],lwd=pixel.per.cnv)
      }else{
        x <- startPoint - pixel.per.cnv*index
        lines(c(x,x),c(y-starts[index],y-ends[index]),
              col="black",lwd=pixel.per.cnv)
      }
    }
  #  # X Chromosome
  #  for(index in 1:len) {
  #    if(indX[index] == TRUE){
  #      x <- startPoint - pixel.per.cnv*index
  #      lines(c(x,x),c(y-starts[index],y-ends[index]),col="black",lwd=pixel.per.cnv)
  #    }
  #  }

    # Y Chromosome
   # for(index in 1:len){
  #    if(indY[index] == TRUE){
  #      x <- startPoint - pixel.per.cnv*index
  #      lines(c(x,x),c(y-starts[index],y-ends[index]),col="black",lwd=pixel.per.cnv)
  #    }
  #  }
  }
}
