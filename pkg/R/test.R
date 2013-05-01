## insert newrow to existing data frame at rth place
insertRow <- function(existingDF, newrow, r) {
  existingDF[seq(r+1,nrow(existingDF)+1),] <- existingDF[seq(r,nrow(existingDF)),]
  existingDF[r,] <- newrow
  existingDF
}
# output clean model comparison for rokas data, gene by gene
model_comp <- function(i){
  protTest_file <- paste("~/proteinevoutk20/pkg/Result/Prottest/Rokas/rokas",i,".txt",sep="")
  maxroot_file <- paste("~/BackupProEvo/Newton/rokas/rootMax/gene",i,".RData",sep="")
  emproot_file <- paste("~/BackupProEvo/Newton/rokas/rootEmp/gene",i,".RData",sep="")
  oproot_file <- paste("~/BackupProEvo/Newton/rokas/rootOp/gene",i,".RData",sep="")
  eqmroot_file <- paste("~/BackupProEvo/Newton/rokas/rootEqm/gene",i,".RData",sep="")
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

# plot.grid <- function(index){
#   xyz = grid.gen(opw_mean,index,res_op=res_op,gridnum=20)
#   akima.xyz = interp(xyz$x,xyz$y,xyz$z)
#   image(akima.xyz)
#   contour(akima.xyz,add=TRUE)
#   points(xyz,pch=3,col="blue")
# }

beta=be
gamma=ga
s = 2
Q = NU_VEC
dismat = GM_cpv(GM_CPV,al,beta,gamma)
fixmatall <- fixmatAll(s,DisMat=dismat)
mumat = aa_MuMat_form(Q)
Qall = QAllaa1(fixmatall,mumat)
getres <- function(i){
  maxroot_file <- paste("~/BackupProEvo/Newton/rokas/rootMax/gene",i,".RData",sep="")
  emproot_file <- paste("~/BackupProEvo/Newton/rokas/rootEmp/gene",i,".RData",sep="")
  oproot_file <- paste("~/BackupProEvo/Newton/rokas/rootOp/gene",i,".RData",sep="")
  eqmroot_file <- paste("~/BackupProEvo/Newton/rokas/rootEqm/gene",i,".RData",sep="")
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
  comp <- NULL
  for(i in 1:6){
    comp <- c(comp,sum(opaa[,ind[i,1]]==opaa[,ind[i,2]]))
  }
  return(list(opaa=opaa,comp=comp))
}
