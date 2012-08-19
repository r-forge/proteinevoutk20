source("~/proteinevoutk20/pkg/R/pphyproevo.R")
tips <- c(4,8,12)
BrLen <- c(100,50,20)
treegen <- function(tip,brlen){
  rtree(tip,rooted=TRUE,br=brlen*runif(2*tip-2))
}
trees <- mapply(treegen,tips,BrLen,SIMPLIFY=FALSE)
save.image(file="TreesForSim.RData",compress=TRUE)
