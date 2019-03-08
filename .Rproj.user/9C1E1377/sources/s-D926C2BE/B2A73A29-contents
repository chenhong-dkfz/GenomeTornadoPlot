focallity.score <- function(m,ends,starts){
  mean.length <-  (sum(ends - starts))/m
  range <- range(ends - starts)
  range.length <- range[2]-range[1]
  #sd.score <- round(sd(range.length/(ends - starts)),2)
  #sd.score <- round(sd(mean.length/(ends - starts)),2)
  mean.score <- round(mean(range.length/(ends - starts)),2)
  #mean.score <- round(mean(mean.length/(ends - starts)),2)
  #f.score <- paste(mean.score,"±",sd.score)
  f.score <- mean.score
}
