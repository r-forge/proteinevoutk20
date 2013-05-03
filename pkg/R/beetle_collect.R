eqmdir <- "~/BackupProEvo/Newton/beetle/tree1/rootEqm/"
l <- 8
res_max <- vector("list",length=8)

for(genect in 1:l){
  filename = paste(eqmdir, "gene_",genect,".RData",sep="")
  if(!file.exists(filename))
    cat("load RData for gene", genect,"failed, file does not exist","\n")
  else{
    load(filename)
    res_max[[genect]] <- res_op
  }
}

s_max <- sapply(1:8,function(x) res_max[[x]]$s)
loglik_max <- sapply(1:8,function(x) res_max[[x]]$ll$loglik)
GM_max <- sapply(1:8,function(x) res_max[[x]]$GMweights)
Q_max <- sapply(1:8,function(x) res_max[[x]]$Q)
br_max <- sapply(1:8,function(x) sum(res_max[[x]]$tree$edge.length))

model_comp.b <- function(i){
  protTest_file <- paste("~/proteinevoutk20/pkg/Result/Prottest/beetle/beetle",i,".txt",sep="")
  maxroot_file <- paste("~/BackupProEvo/Newton/beetle/tree1/rootMax/gene_",i,".RData",sep="")
  emproot_file <- paste("~/BackupProEvo/Newton/beetle/tree1/rootEmp/gene_",i,".RData",sep="")
  oproot_file <- paste("~/BackupProEvo/Newton/beetle/tree1/rootOp/gene_",i,".RData",sep="")
  eqmroot_file <- paste("~/BackupProEvo/Newton/beetle/tree1/rootEqm/gene_",i,".RData",sep="")
  #### protTest result on empirical models
  r1 <- read.table(protTest_file,header=TRUE,as.is=TRUE)
  names(r1)[5] <- "neg.LnL"
  r1$para <- round((r1$AIC - 2*r1$neg.LnL)/2)
  r1$AICw <- NULL
  #### root states maximized
  load(maxroot_file)
  df_bs <- length(res_op$tree$edge.length)+3+5+19 #number of free parameters, without counting op sites
  nr <- attr(res_op$data,"nr")
  r1 <- insertRow(r1,list("New+max+maxroot",0,0,round(-res_op$ll$loglik,2),df_bs+2*nr),1)
  #### root states with empirical frequencies
  load(emproot_file)
  r1 <- insertRow(r1,list("New+max+emproot",0,0,round(-res_op$ll$loglik,2),df_bs+nr),2)
  #### root states with equilibrium frequencies
  load(eqmroot_file)
  r1 <- insertRow(r1,list("New+max+eqmroot",0,0,round(-res_op$ll$loglik,2),df_bs+nr),2)
  #### root states with empirical frequencies
  load(oproot_file)
  r1 <- insertRow(r1,list("New+max+oproot",0,0,round(-res_op$ll$loglik,2),df_bs+nr),2)
  r1$AIC <- 2*(r1$neg.LnL+r1$para)
  r1$deltaAIC <- r1$AIC - min(r1$AIC)
  return(r1)
}
getres.b <- function(i){
  maxroot_file <- paste("~/BackupProEvo/Newton/beetle/tree1/rootMax/gene_",i,".RData",sep="")
  emproot_file <- paste("~/BackupProEvo/Newton/beetle/tree1/rootEmp/gene_",i,".RData",sep="")
  oproot_file <- paste("~/BackupProEvo/Newton/beetle/tree1/rootOp/gene_",i,".RData",sep="")
  eqmroot_file <- paste("~/BackupProEvo/Newton/beetle/tree1/rootEqm/gene_",i,".RData",sep="")
  res <- vector(mode="list",length=4)
  load(maxroot_file)
  res[[1]] <- res_op
  load(emproot_file)
  res[[2]] <- res_op
  load(oproot_file)
  res[[3]] <- res_op
  load(eqmroot_file)
  res[[4]] <- res_op
  return(res)
}
op.comp <- function(res){
  ind <- combinations(n=4,r=2)
  opaa <- sapply(1:4,function(x) res[[x]]$ll$opaa)
  dimnames(opaa)[[2]] <- c("max","emp","op","eqm")
  comp <- NULL
  for(i in 1:6){
    comp <- c(comp,sum(opaa[,ind[i,1]]==opaa[,ind[i,2]]))
  }
  return(list(opaa=opaa,comp=comp))
}