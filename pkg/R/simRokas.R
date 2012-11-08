load("~/proteinevoutk20/pkg/scratch/newton/rokas_mle_newmat/rokasMLE.RData")
source("~/proteinevoutk20/pkg/R/main.R")
source("~/proteinevoutk20/pkg/R/simulation.R")
index <- attr(res_op$data,"index")
opaa <- res_op$ll$opaa
protein_op <- opaa[index]
tree <- res_op$tree
s <- res_op$s
beta <- res_op$GMweights[2]
gamma <- res_op$GMweights[3]
bfaa <- res_op$bfaa
GTRvec <- res_op$Q

sim <- vector("list",length=42)
for(i in 1:42){
  start_ind <- 1 + (i-1)*1000
  end_ind <- i*1000
  cat("the", i, "th simulation:", "\n")
  sim[[i]] <- simTree(tree,protein_op[start_ind:end_ind],s,GTRvec,alpha=al,beta=beta,gamma=gamma,bfaa=bfaa)
}

simdata <- sim[[1]]$data
for(i in 2:42){
  simdata <- cbind(simdata,sim[[i]]$data)
}
simdata.phy <- phyDat(simdata,type="AA")
sim.res <- mllm(data=simdata.phy,tree=tree,s=s,beta=beta,gamma=gamma,Q=GTRvec,bfaa=bfaa)
str(sim.res)
sim.res_op = optim.mllm(sim.res,optQ=T,optBranch=T,optsWeight=T,control=list(epsilon=1e-08,maxit=300,hmaxit=50,trace=0,htrace=1))

save(res_op,sim.res_op,file="simRokas.RData")