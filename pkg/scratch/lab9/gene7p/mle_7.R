source("~/proteinevoutk20/pkg/R/main.R") 
 gene7 = conv("gene7.fasta","phyDat")
tree = read.tree("gene7prune.tre")
 res=mllm(gene7,tree,0.1,be,ga,Q=NU_VEC) 
 res$ll$loglik 
 res_op = optim.mllm(res,optQ=T,optBranch=T,optsWeight=T,optOpw=FALSE,control=list(epsilon=1e-08,hmaxit=500,htrace=1,print_level=0,maxeval="50")) 
 res_op$ll$loglik 
 save.image(file="gene7.RData",compress=TRUE) 
