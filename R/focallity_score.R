# add new log method

focallity.score <- function(m,ends,starts,method){
  if(missing(method)){method="normal"}
  if(method=="normal"){
    mean.length <-  (sum(ends - starts))/m
    range <- range(ends - starts)
    range.length <- range[2]-range[1]
    mean.score <- round(mean(range.length/(ends - starts)),2)
    f.score <- mean.score
  }
  if(method=="log"){
    max.length <- max(abs(ends-starts))
    scores <- log10(max.length)-log10(ends-starts)
    f.score <- sum(scores)
  }
  return(f.score)
}
