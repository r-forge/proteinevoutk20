source("~/proteinevoutk20/pkg/R/prune.R") 
gene = 1
prottestfile <- paste("~/proteinevoutk20/pkg/Result/Prottest/Rokas7/rokas_",gene,"_prottest.txt",sep="") 
RDatafile <- paste("~/BackupProEvo/Newton/rokas7/rootEqm/gene",gene,".RData",sep="")
imagefile <- paste("gene",gene,".RData",sep="") 
best_emp_model <- get_best_model(prottestfile) 
load(RDatafile)
#######################################################################
## simulation on one branch only
# test = FALSE 
# p2 <- prune_emp(fastafile,"Smik",ROKAS_TREE,best_emp_model$model) 
# p1 <- prune_new(fastafile,"Smik",ROKAS_TREE,ancestral="eqm") 
# save.image(RDatafile,compress=TRUE) 
#######################################################################
## simulation along the tree
datanum <- sapply(1:length(dataphy),function(x) as.numeric(factor(data[x,],tolower(AA))))
datanum <- t(datanum)
dimnames(datanum)[[1]] <- dimnames(data)[[1]]
nr <- attr(dataphy,"nr")
index <- attr(dataphy,"index")
opaa <- res_op$ll$opaa
ftny.vec <- apply(datanum,MARGIN=1,FUN=Ftny_protein,protein_op=opaa[index],s=s,DisMat=GM_cpv(GM_CPV,al,beta,gamma))
s = res_op$s
beta = res_op$GMweights[2]
gamma = res_op$GMweights[3]
# nsim <- 1
# sim <- vector(mode="list")
# sim$new <- vector(mode="list",length=nsim)
# sim$emp <- vector(mode="list",length=nsim)
# for(i in 1:nsim){
#   if(exists(".Random.seed"))
#     rm(.Random.seed)
#   ## root sequence, from the inference of ancestral states in mllm1, length equal to number of sites
#   rootseq <- sapply(1:length(index),function(x) sample(1:20,size=1,replace=TRUE,prob=res$ll$ancestral[,index[x]]))
#   sim$new[[i]] <- simTree(tree=res_op$tree,protein_op=opaa[index],matall=res_op$Qall,rootseq=rootseq)
#   sim$emp[[i]] <- sim_emp(dataphy,best_emp_model) ## best empirical model
# }
# 
# ## find the starting and ending position of all edges in the tree
# plot_trace(sim=sim$emp[[1]],plottip=FALSE)
# plot_trace(sim=sim$new[[1]])
# save.image(file=imagefile,compress=TRUE)