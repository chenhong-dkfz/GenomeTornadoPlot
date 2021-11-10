# convert copy number to scores

MapPloidyClasses <- function(v){
  if(v==0){
    cnv=1 # bi-del
  }else if(v==1){
    cnv=2 # mono-del
  }else if(v==3|v==4){
    cnv=3 # CN < 5
  }else if(v>=5 & v<=8){
    cnv=4 # CN 4~9
  }else if(v>=9 & v<=99999999){
    cnv=5 # CN > 8
  }else{cnv = 1} # not exist
  return(cnv)
}
