# GetColor: return a list of colors for CNVs
# method: the method of color arrangement. "confidence","ploidy","cohort","length","factor"
# score: ploidy of variation
# color: tbd
# score.values: tbd
# n: tbd
# greyscale: return color or only greyscale
# cohorts: cohort list for the study


GetColor <- function(method,score,color,score.values,n,greyscale,cohorts){
  if(missing(score)){score=0}
  if(missing(score.values)){score.values=0}
  if(missing(n)){n=0}
  if(missing(greyscale)){greyscale=FALSE}
  if(missing(cohorts)){cohort=c("all_patients")}
  switch(method,

         "confidence" = {
           color.value <- "black"
           switch (score,
                   "1" = {color.value = "grey"},
                   "2" = {color.value = "royalblue1"},
                   "3" = {color.value = "royalblue3"},
                   "4" = {color.value = "royalblue4"}
           )
         },

         "repeat" = {
           color.value <- "black"
           switch (score,
                   "U1" = {color.value = "red"},
                   "C1" = {color.value = "black"},
                   "U2" = {color.value = "blue"}
           )
           color.value <- c("red","black","black","blue")
         },

         "ploidy" = {
           color.base <- colorRampPalette(c("red2","indianred4","lightskyblue2","skyblue3","midnightblue"))(5)
           if(missing(color)){
             color1 <- color.base
           }else if(length(unique(color))!=1){
             color1 <- colorRampPalette(color)(5)
             # six types of ploidy, 2deletion+dip+3gain
           }else{color1 <- color.base}
           color.value <- color1
         },

         "ploidy2" = {
           # need redo
           if(missing(score)){
             color.value <- color[length(score.values)+1]
           }else{
             color.value <- color[score]
           }
         },

         "cohort" = {
           cohort.size <- length(unique(cohorts))
           if(missing(color)){
             color1 <- colorRampPalette(c("red2","indianred4","royalblue4","steelblue1","chartreuse3","darkgreen"))(cohort.size)
           }else if(color>0 & color<8 & is.integer(color)==TRUE){
             switch (color,
                     "1" = {color1 <- rainbow(cohort.size)},
                     "2" = {color1 <- colorRampPalette(c("seashell2","seagreen2","turquoise2","palevioletred2"))(cohort.size)},
                     "3" = {color1 <- colorRampPalette(c("gray7","mediumblue","deeppink4","sienna3"))(cohort.size)},
                     "4" = {color1 <- colorRampPalette(c("darkslateblue","darkslategray4","deeppink4","tan4","gray10"))(cohort.size)},
                     "5" = {color1 <- colorRampPalette(c("navajowhite3","orange3","olivedrab3"))(cohort.size)},
                     "6" = {color1 <- colorRampPalette(c("royalblue2","yellow1"))(cohort.size)},
                     "7" = {color1 <- gray.colors(cohort.size)}
             )
           }else if(greyscale==TRUE){
             color1 <- gray.colors(cohort.size)
             color1 <- palette(color1)
           }else{
             color1 <- colorRampPalette(c("red2","indianred4","royalblue4","steelblue1","chartreuse3","darkgreen"))(cohort.size)
           }
           color.value <- color1
         },

         "length" = {
           if(missing(color)){color="black"}
           color.value <- color
         },

         "factor" = {
           if(missing(color)){color="black"}
           color.value <- "black"
         }
  )
  return(color.value)
}

