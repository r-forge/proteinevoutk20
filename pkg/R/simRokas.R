##simulations using Rokas's tree information
line1 <- 'source("~/proteinevoutk20/pkg/R/hessian.R")'
line2 <- 'tree_short <- read.nexus("~/proteinevoutk20/pkg/Data/GTR.tre")' #Rokas's tree
line3 <- 'tree_long <- read.nexus("~/proteinevoutk20/pkg/Data/GTRlong.tre")' #Rokas's tree with longer branches
line4 <- 'op_seq <- sample(20,1000,replace=T)'
line5 <- 'root_seq <- sample(20,1000,replace=T)'
svalues <- c(0,0.1,0.5,0.9,1.5,2.0) #6 different s values
pbsCommands=""
queue="medium*"

for(sval in svalues){
  line7 <- paste("simdata_short <- simTree(tree_short,1000,protein_op=op_seq,20,",sval,",Nu_vec,rootseq=root_seq)",sep="")
  line8 <- paste("simdata_long <- simTree(tree_long,1000,protein_op=op_seq,20,",sval,",Nu_vec,rootseq=root_seq)",sep="")
  line9 <- paste('save(simdata_short,simdata_long,file="simRokas_',sval,'.RData",compress=TRUE))',sep="")
  cat(line1,"\n",line2,"\n",line3,"\n",line4,"\n",line5,"\n",line7,"\n",line8,"\n",line9,"\n",
      file=paste("simRokas_",sval,".R",sep=""))
  pbsCommands=paste('#!/bin/bash','#$ -cwd','#$ -o /dev/null','#$ -e /dev/null',sep="\n")
  pbsCommands=paste(pbsCommands,'\n#$ -q ',queue,sep="")
  pbsCommands=paste(pbsCommands,'#$ -M chaij@umail.iu.edu', '#$ -m beas', '#$ -S /bin/bash',sep="\n")
  pbsCommands=paste(pbsCommands,"\n","#$ -N SimRokas_",sval,sep="")
  pbsCommands=paste(pbsCommands,"\n","/data/apps/R/2.14.0/bin/R CMD BATCH simRokas_",sval,".R",sep="")
  cat(pbsCommands,file=paste('runSimRokas_',sval,'.sh',sep=""),append=FALSE)
  system(paste('chmod u+x runSimRokas_',sval,'.sh',sep=""))
  #system(paste('/opt/sge/bin/lx24-amd64/qsub runSimRokas_',sval,'.sh',sep=""))
}