# rm(list=ls())
# gene = 60
# RDatafile <- paste("~/BackupProEvo/Newton/rokas/simtree/gene",gene,".RData",sep="")
# load(RDatafile)
# source("~/proteinevoutk20/pkg/R/prune.R")
# source("~/proteinevoutk20/pkg/R/plot_simtree.R")
# obs.data <- conv(filename=fastafile,type="num") #observed data sequences
# ftny.vec <- apply(obs.data,MARGIN=1,FUN=Ftny_protein,protein_op=opaa[index],s=s,DisMat=GM_cpv(GM_CPV,al,beta,gamma))
# ## find the starting and ending position of all edges in the tree
# plot_trace(sim=sim$emp[[1]],model="emp",plottip=FALSE)
# plot_trace(sim=sim$new[[1]],model="new")
#################################################
#gene <- Sys.getenv("SGE_TASK_ID")
gene <- 1
print(gene)
gene <- as.numeric(gene)
imagefile <- paste("gene",gene,".RData",sep="")
load("~/proteinevoutk20/pkg/Data/Rokas7/rokas7.RData")
source("~/proteinevoutk20/pkg/R/main.R")

tree <- read.nexus("~/proteinevoutk20/pkg/Data/Rokas7/GTR7.tre") ## tree
data <-  rokas7[,seq(AArange[gene,1],AArange[gene,2])] ##data
sort(unique(c(data)))
dataphy <- phyDat(data=data,type="AA") ## data in phyDat format

res=mllm1(dataphy,tree,s=0.1,beta=be,gamma=ga,Q=NULL) 
res$ll$loglik 
res_op = optim.mllm1(res,optQ=T,optBranch=T,optsWeight=T,optOpw=FALSE,
                     control=list(epsilon=1e-08,hmaxit=200,htrace=1,print_level=0,maxeval="50")) 
res_op$ll$loglik
save.image(file=imagefile,compress=TRUE)
