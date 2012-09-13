sites <- c(100,300,500,700,1000)
svalues <- c(0.1,0.3,0.5,0.7,0.9,1)
tipnums <- c(4,4,8,8,10,10,12,12,14,14,16,16,18,18)

line1 <- 'source("~/proteinevoutk20/pkg/R/hessian.R")'
##line2 <- 'trees <- read.tree("~/proteinevoutk20/pkg/R/simtrees.tre")'
pbsCommands=""
queue="short*"

for(k in 1:14){
  for(i in sites){
    for(j in svalues){
      ##load the workspace which stores the simulated data
      line3 <- paste('load("op_simu','_',k,'_',i,'_',j,'.RData")',sep="")
      ##recalculate the hessian at the mle_log_op, previous calculation used mle_log$par, which is wrong
      line4 <-  paste("hes_op <- find_hessian_log(mle_log_op$par,trees[[",k,"]],sim[1:",tipnums[k],",],al,mumat,20,protein_op=op_seq,root=start_seq)",sep="")
      line5 <- paste('save.image(file="hes_op_simu_',k,'_',i,'_',j,'.RData",compress=TRUE)',sep="")
       
      cat(line1,"\n",line3,"\n",line4,"\n",line5,"\n",file=paste("hes_op_sim_",k,"_",i,"_",j,".R",sep=""))
      
       pbsCommands=paste('#!/bin/bash','#$ -cwd','#$ -o /dev/null','#$ -e /dev/null',sep="\n")
       pbsCommands=paste(pbsCommands,'\n#$ -q ',queue,sep="")
       pbsCommands=paste(pbsCommands,'#$ -M chaij@umail.iu.edu', '#$ -m beas', '#$ -S /bin/bash',sep="\n")
       pbsCommands=paste(pbsCommands,"\n","#$ -N Op1Sim_",k,"_",i,"_",j,sep="")
       pbsCommands=paste(pbsCommands,"\n","/data/apps/R/2.14.0/bin/R CMD BATCH hes_op_sim_",k,"_",i,"_",j,".R",sep="")
       cat(pbsCommands,file=paste('runSimOp_',k,'_',i,'_',j,'.sh',sep=""),append=FALSE)
       system(paste('chmod u+x runSimOp_',k,'_',i,'_',j,'.sh',sep=""))
       system(paste('/opt/sge/bin/lx24-amd64/qsub runSimOp_',k,'_',i,'_',j,'.sh',sep=""))
     }
  }
}
