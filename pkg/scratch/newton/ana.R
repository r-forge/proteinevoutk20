filedir <- "~/proteinevoutk20/pkg/scratch/newton/grid_log/"
prefix <- "bglog"
suffix <- ".RData"
separator <- "_"

filedir <- "~/proteinevoutk20/pkg/scratch/newton/RData/FirstRun/"
prefix <- "bg"
suffix <- ""
separator <- "."

get.val <- function(x,y){
  file <- paste(filedir,prefix,separator,x,separator,y,suffix,sep="")
  #file <- paste("~/proteinevoutk20/pkg/scratch/newton/RData/FirstRun/bg.",x,".",y,sep="")
  load(file)
##  par <- sapply(1:106, function(x) res[[x]]$par)
  val <- sapply(1:106, function(x) res[[x]]$objective) # get the -likelihood values
  sum(val) # sum up all the -likelihood values
}
vec.get.val <- Vectorize(get.val, c("x","y"))

##Check and see if all 106 genes finished running and the .RData files are saved for the particular grid##
## 4 big grids: 1,3,5,6; each grid is splitted into 20*20 small grids indexed by x and y##
check.grid <- function(k,x,y){
  file <- paste("~/proteinevoutk20/pkg/scratch/newton/RData/SecondRun/grid.",k,".",x,".",y,sep="") #the file to chek existence
  if(!file.exists(file)){ # if it doesn't exist, print message, and return 1
    print(paste("grid.",k,".",x,".",y," does not exist! Check running time.",sep=""))
    return(1)
  }
  else # if the file exists, return 0
    return(0)
}

vcheck <- Vectorize(check.grid,c("x","y")) # vectorize x and y, (small grid indices)
###check3 <- outer(1:20,1:20,vcheck,k=3) # check on the 3rd big grid
###arrayInd(which(check3==1),c(20,20)) # the indices of small grids that didn't finish running


### This funciton get the total likelihood value for the 4 picked grids 1,3,5,6 (k)
### x and y index the small grids in the 4 big grids, range 1:20
### Similar to the function get.val
get.val.grid <- function(k,x,y){
  file <- paste("~/proteinevoutk20/pkg/scratch/newton/RData/SecondRun/grid.",k,".",x,".",y,sep="")
  if(file.exists(file)){
    load(file)
    ##  par <- sapply(1:106, function(x) res[[x]]$par)
    val <- sapply(1:106, function(x) res[[x]]$objective)
    sum(val)
  }
  else{
    print(paste("grid.",k,".",x,".",y," does not exist! Check running time.",sep=""))
    return(NA)
  }
}


##vget <- Vectorize(get.val,c("x","y")) #vectorized version of get.val, both arguments are vectorized
##likelihood.array <- outer(1:20,1:20,vget) #apply vget to recombination of all values across 1:20, return 20*20 array

vget.grid <- Vectorize(get.val.grid,c("x","y")) #vectorized version of get.val, both arguments are vectorized
### return 2 dimensional array of likelihood values for each of the 4 grids
get.likelihoods.grid <- function(grid){
  likelihood.array.grid <- outer(1:20,1:20,vget.grid,k=grid) #apply vget to recombination of all values across 1:20, return 20*20 array
  return(likelihood.array.grid)
}

get.s.grid <- function(k,x,y){
  file <- paste("~/proteinevoutk20/pkg/scratch/newton/RData/SecondRun/grid.",k,".",x,".",y,sep="")
  if(file.exists(file)){
    load(file)
    para <- sapply(1:106, function(x) res[[x]]$par)
    ## val <- sapply(1:106, function(x) res[[x]]$objective)
    ## sum(val)
    return(para)
  }
  else{
    print(paste("grid.",k,".",x,".",y," does not exist! Check running time.",sep=""))
    return(NA)
  }
}
vget.s.grid <- Vectorize(get.s.grid,c("x","y")) #vectorized version of get.s.grid, both arguments are vectorized

### A function to find all the s values for a specific gene, in one grid.
get.s.gene <- function(k,x,y,index){
  file <- paste("~/proteinevoutk20/pkg/scratch/newton/RData/SecondRun/grid.",k,".",x,".",y,sep="")
  if(file.exists(file)){
    load(file)
    para <- res[[index]]$par
    ## val <- sapply(1:106, function(x) res[[x]]$objective)
    ## sum(val)
    return(para)
  }
  else{
    print(paste("grid.",k,".",x,".",y," does not exist! Check running time.",sep=""))
    return(NA)
  }
}
### vectorize the indices of small grids in the 4 grids
vget.s.gene <- Vectorize(get.s.gene,c("x","y"))

### for a gene and a big grid, find all 400 mles for s in the small grids
s.gene.grid <- function(grid,index){
  c(outer(1:20,1:20,vget.s.gene,k=grid,index=index))
}

### Big grids
grids <- c(1,3,5,6)

### for one gene with "index" find all 1600 mles in all the small grids
s.gene <- function(index){
  c(sapply(grids, s.gene.grid,index))
}

### mean value of mle s from all the 1600 grids
mean.s.gene <- function(index){
  mean(s.gene(index))
}
