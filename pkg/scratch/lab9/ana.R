get.val <- function(x,y){
  file <- paste("~/proteinevoutk20/pkd/scratch/newton/bg.",x,".",y,sep="")
  load(file)
##  par <- sapply(1:106, function(x) res[[x]]$par)
  val <- sapply(1:106, function(x) res[[x]]$objective)
  sum(val)
}
