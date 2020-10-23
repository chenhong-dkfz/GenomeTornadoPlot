# this function calculates entropy of cohort distributions

entropy.calculator <- function(cohorts){
  x = entropy(table(cohorts),unit = "log2")
  return(x)
}

