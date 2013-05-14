source("~/proteinevoutk20/pkg/R/prune.R") 
 test = FALSE 
 gene = 1 
 fastafile <- paste("~/proteinevoutk20/pkg/Data/Rokas/gene",gene,".fasta",sep="") 
 prottestfile <- paste("~/proteinevoutk20/pkg/Result/Prottest/Rokas/rokas_",gene,"_prottest.txt",sep="") 
 RDatafile <- paste("gene",gene,".RData",sep="") 
 best_emp_model <- get_best_model(prottestfile) 
 p2 <- prune_emp(fastafile,"Smik",ROKAS_TREE,best_emp_model$model) 
 p1 <- prune_new(fastafile,"Smik",ROKAS_TREE,ancestral="max") 
 save.image(RDatafile,compress=TRUE) 
