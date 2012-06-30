get.val <- function(x,y){
  file <- paste("~/proteinevoutk20/pkg/scratch/newton/RData/bg.",x,".",y,sep="")
  load(file)
##  par <- sapply(1:106, function(x) res[[x]]$par)
  val <- sapply(1:106, function(x) res[[x]]$objective)
  sum(val)
}

vget <- Vectorize(get.val,c("x","y")) #vectorized version of get.val, both arguments are vectorized
likelihood.array <- outer(1:20,1:20,vget) #apply vget to recombination of all values across 1:20, return 20*20 array
