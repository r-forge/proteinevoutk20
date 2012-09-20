library(akima) #package for function "interp"
filedir <- "~/proteinevoutk20/pkg/scratch/newton/grid_log_0914/"
prefix <- "bglog"
suffix <- ".RData"
separator <- "_"

# filedir <- "~/proteinevoutk20/pkg/scratch/newton/RData/FirstRun/"
# prefix <- "bg"
# suffix <- ""
# separator <- "."

get.val <- function(x,y){
  file <- paste(filedir,prefix,separator,x,separator,y,suffix,sep="")
  #file <- paste("~/proteinevoutk20/pkg/scratch/newton/RData/FirstRun/bg.",x,".",y,sep="")
  if(file.exists(file)){
  load(file)
##  par <- sapply(1:106, function(x) res[[x]]$par)
  val <- sapply(1:106, function(x) res[[x]]$objective) # get the -likelihood values
  sum(val) # sum up all the -likelihood values
  }
  else
    return (370000)
}
vec.get.val <- Vectorize(get.val, c("x","y"))

ll <- 15
beta <- seq(-15,0,length.out=(ll+1))[-1]
gamma <- seq(-15,0,length.out=(ll+1))[-1]
z = c(outer(1:15,1:15,vec.get.val))
x = rep(beta,15)
y = rep(gamma, each=15)
bg.li = interp(x,y,z)
image(bg.li,main="contour plot on [-15,0] by [-15,0]")
contour(bg.li,add=T)

ll1 <- 26
beta1 <- seq(-5,0,length.out=ll1)
gamma1 <- seq(-7,-3,length.out=ll1)
z1 = c(outer(1:26,1:26,vec.get.val))
x1 = rep(beta1,26)
y1 = rep(gamma1, each=26)
bg1.li = interp(x1,y1,z1)
image(bg1.li,main="contour plot on [-5,0] by [-7,-3]")
contour(bg1.li,add=T)

l <- 26
beta <- seq(0,5,length.out=l)
gamma <- seq(-7,-3,length.out=l)
z2 = c(outer(1:26,1:26,vec.get.val))
x2 = rep(beta,26)
y2 = rep(gamma, each=26)
bg2.li = interp(x2,y2,z2)
image(bg2.li,main="contour plot on [0,5] by [-7,-3]")
contour(bg2.li,add=T)

X = c(x,x1)
Y = c(y,y1)
Z = c(z,z1)
Bg.li = interp(X,Y,Z,duplicate="mean",xo=seq(-14,0,length=100),yo=seq(-14,0,length=100))
image(Bg.li,main="contour plot combined")
contour(Bg.li,add=T)
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
