load("~/proteinevoutk20/pkg/scratch/lab9/bgsearch/bgsearch.RData")
load("~/proteinevoutk20/pkg/scratch/lab9/bgsearch/bgsearch_svals.RData")
source("~/proteinevoutk20/pkg/R/main.R")
source("~/proteinevoutk20/pkg/R/simulation.R")
source("~/proteinevoutk20/pkg/R/readRokas.R")
# GM weights, same for all the genes
beta = bg$par[1]
gamma = bg$par[2]
# s values for all genes under the values for beta and gamma
bg.sval = sapply(1:106,function(x) bg.s[[x]]$par)
simgenes <- vector("list",length=106)
simdata <- NULL
for(i in 1:106){
  cat("start gene",i,"\n")
  res = mllm(ROKAS_DATA[[i]],tree,s=bg.sval[i],beta=beta,gamma=gamma,Q=GTRvec) #likelihood at the ML estimates, the whole data structure
  index <- attr(res$data,"index")
  opaa <- res$ll$opaa
  protein_op <- opaa[index]
  simgenes[[i]] <- simTree(tree,protein_op,bg.sval[i],GTRvec,alpha=al,beta=beta,gamma=gamma,
                           rootseq=sample(20,length(opaa),replace=T))
  cat("finish gene",i,"\n")
  simdata <- cbind(simdata,simgenes[[i]]$data)
}
save.image(file="simRokas.RData",compress=TRUE)