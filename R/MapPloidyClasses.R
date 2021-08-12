# convert copy number to scores

MapPloidyClasses <- function(v){
  if(v == 0){class = 1 # homozygous deletion
  }else if(v<5){
    switch(v,
           "1" = {class = 2}, # heterozygous deletion
           "2" = {class = 3}, # diploidy
           "3" = {class = 4}, # copy number gain low
           "4" = {class = 4}
    )
  }else if(v>=5 & v<=8){
    class = 5 # copy number gain middle
  }else if(v>=9 & v<=99999999){
    class = 6 # copy number gain high
  }else{class = "unknown"}
  return(class)
}
