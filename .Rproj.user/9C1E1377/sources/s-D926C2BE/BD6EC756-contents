
plotCnv <- function(chroms,starts,ends,y,scores,pixel.per.cnv,method,color,score.values,n,startPoint,direction){

  indX <- chroms == 'X'
  indY <- chroms == 'Y'

  len <- length(starts)
  if(missing(direction)){direction="right"}
  if(direction=="right"){
    # Autosomes
    for(index in 1:len){
      x <- startPoint + pixel.per.cnv*index
      lines(c(x,x),c(y-starts[index],y-ends[index]),col="black",lwd=pixel.per.cnv)
    }

    # X Chromosome
    for(index in 1:len) {
      if(indX[index] == TRUE){
        x <- startPoint + pixel.per.cnv*index
        lines(c(x,x),c(y-starts[index],y-ends[index]),col="black",lwd=pixel.per.cnv)
      }
    }

    # Y Chromosome
    for(index in 1:len){
      if(indY[index] == TRUE){
        x <- startPoint + pixel.per.cnv*index
        lines(c(x,x),c(y-starts[index],y-ends[index]),col="black",lwd=pixel.per.cnv)
      }
    }
  }else{
    for(index in 1:len){
      x <- startPoint - pixel.per.cnv*index
      lines(c(x,x),c(y-starts[index],y-ends[index]),col="black",lwd=pixel.per.cnv)
    }

    # X Chromosome
    for(index in 1:len) {
      if(indX[index] == TRUE){
        x <- startPoint - pixel.per.cnv*index
        lines(c(x,x),c(y-starts[index],y-ends[index]),col="black",lwd=pixel.per.cnv)
      }
    }

    # Y Chromosome
    for(index in 1:len){
      if(indY[index] == TRUE){
        x <- startPoint - pixel.per.cnv*index
        lines(c(x,x),c(y-starts[index],y-ends[index]),col="black",lwd=pixel.per.cnv)
      }
    }
  }
}
