##simulations using Rokas's tree information
line1 <- 'source("~/proteinevoutk20/pkg/R/hessian.R")'
# line2 <- 'tree_short <- read.nexus("~/proteinevoutk20/pkg/Data/GTR.tre")' #Rokas's tree
# line3 <- 'tree_long <- read.nexus("~/proteinevoutk20/pkg/Data/GTRlong.tre")' #Rokas's tree with longer branches
# line4 <- 'op_seq <- sample(20,1000,replace=T)'
# line5 <- 'root_seq <- sample(20,1000,replace=T)'
svalues <- c(0,0.1,0.5,0.9,1.5,2.0) #6 different s values
pbsCommands=""
queue="medium*"

for(sval in svalues){
  line2 <- paste('load("simRokas_',sval,'.RData")',sep="")
  line3 <- 'startpt <- log(c(sval,be,ga))'
  line4 <- 'mle_log_nop_short <- MLE_sw_log(startpt,tree_short,simdata_short[1:8,],20,al,mumat,trace=1)'
  line5 <- 'mle_log_nop_long <- MLE_sw_log(startpt,tree_long,simdata_long[1:8,],20,al,mumat,trace=1)'
  line6 <- 'mle_log_op_short <- MLE_sw_log(startpt,tree_short,simdata_short[1:8,],20,al,mumat,protein_op=op_seq,trace=1)'
  line7 <- 'mle_log_op_long <- MLE_sw_log(startpt,tree_long,simdata_long[1:8,],20,al,mumat,protein_op=op_seq,trace=1)'

  line8 <- paste('save.image(file="estim_Rokas_',sval,'.RData",compress=TRUE))',sep="")
  cat(line1,"\n",line2,"\n",line3,"\n",line4,"\n",line5,"\n",line6,"\n", line7,"\n",line8,"\n",
      file=paste("estimRokas_",sval,".R",sep=""))
  pbsCommands=paste('#!/bin/bash','#$ -cwd','#$ -o /dev/null','#$ -e /dev/null',sep="\n")
  pbsCommands=paste(pbsCommands,'\n#$ -q ',queue,sep="")
  pbsCommands=paste(pbsCommands,'#$ -M chaij@umail.iu.edu', '#$ -m beas', '#$ -S /bin/bash',sep="\n")
  pbsCommands=paste(pbsCommands,"\n","#$ -N EstimRokas_",sval,sep="")
  pbsCommands=paste(pbsCommands,"\n","/data/apps/R/2.14.0/bin/R CMD BATCH estimRokas_",sval,".R",sep="")
  cat(pbsCommands,file=paste('runEstimRokas_',sval,'.sh',sep=""),append=FALSE)
  system(paste('chmod u+x runEstimRokas_',sval,'.sh',sep=""))
  #system(paste('/opt/sge/bin/lx24-amd64/qsub runEstimRokas_',sval,'.sh',sep=""))
}