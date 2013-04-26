## insert newrow to existing data frame at rth place
insertRow <- function(existingDF, newrow, r) {
  existingDF[seq(r+1,nrow(existingDF)+1),] <- existingDF[seq(r,nrow(existingDF)),]
  existingDF[r,] <- newrow
  existingDF
}
# output clean model comparison for rokas data, gene by gene
model_comp <- function(i){
  protTest_file <- paste("~/proteinevoutk20/pkg/Result/Prottest/Rokas/rokas",i,".txt",sep="")
  maxroot_file <- paste("~/BackupProEvo/Newton/rokas/gene",i,".RData",sep="")
  emproot_file <- paste("~/BackupProEvo/Newton/rokas/emp_root/gene",i,"_s_weight.RData",sep="")
  r1 <- read.table(protTest_file,header=TRUE,as.is=TRUE)
  names(r1)[5] <- "neg.LnL"
  r1$para <- round((r1$AIC - 2*r1$neg.LnL)/2)
  r1$AICw <- NULL
  load(maxroot_file)
  df_bs <- length(res_op$tree$edge.length)+3+5+19 #number of free parameters, without counting op sites
  nr <- attr(res_op$data,"nr")
  r1 <- insertRow(r1,list("New+max+maxroot",0,0,round(-res_op$ll$loglik,2),df_bs+2*nr),1)
  load(emproot_file)
  r1 <- insertRow(r1,list("New+max+emproot",0,0,round(-res_op$ll$loglik,2),df_bs+nr),2)
  r1$AIC[1:2] <- 2*(r1$neg.LnL[1:2]+r1$para[1:2])
  r1$deltaAIC <- r1$AIC - min(r1$AIC)
  return(r1)
}

trimRData <- function(DataFile,obj="res_op"){
  load(DataFile)
  save(list=obj,file=DataFile,compress=TRUE)
}
# plot.grid <- function(index){
#   xyz = grid.gen(opw_mean,index,res_op=res_op,gridnum=20)
#   akima.xyz = interp(xyz$x,xyz$y,xyz$z)
#   image(akima.xyz)
#   contour(akima.xyz,add=TRUE)
#   points(xyz,pch=3,col="blue")
# }

# beta=be
# gamma=ga
# s = 2
# Q = NU_VEC
# dismat = GM_cpv(GM_CPV,al,beta,gamma)
# fixmatall <- fixmatAll(s,DisMat=dismat)
# mumat = aa_MuMat_form(Q)
# Qall = QAllaa1(fixmatall,mumat)
