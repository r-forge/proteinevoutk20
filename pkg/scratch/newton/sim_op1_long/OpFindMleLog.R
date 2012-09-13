sites <- c(100,300,500,700,1000)
svalues <- c(0.1,0.3,0.5,0.7,0.9,1)
tipnums <- c(16,16,18,18)
##tipnums <- c(4,4,8,8,10,10,12,12,14,14,16,16,18,18)

line1 <- 'source("~/proteinevoutk20/pkg/R/hessian.R")'
line2 <- 'trees <- read.tree("~/proteinevoutk20/pkg/R/simtrees_correct.tre")'
pbsCommands=""
queue="medium*"

for(k in 11:14){
  for(i in sites){
    for(j in svalues){
      ##load the workspace which stores the simulated data
      ##line1 <- paste('load("~/proteinevoutk20/pkg/scratch/newton/RData/sim_',k+6,'_',i,'_',j,'.RData")',sep="")
      ## find the mle for s, beta and gamma, given alpha value

      ##simulate data on the list of trees, number of sites, and values of s
       line3 <- paste("start_seq=sample(20,",i,",replace=T)",sep="") #root sequence
       line4 <- paste("op_seq=rep(1,",i,")",sep="") #optimal aa sequence
       line5 <- paste("system.time(sim <- simTree(trees[[",k,"]]",",",i,",op_seq,20,",j,",Nu_vec,alpha=al,beta=be,gamma=ga,rootseq=start_seq, ancestral=TRUE))",sep="")
       ##MLE given optimal sequence and root sequence
       line6 <- paste("system.time(mle_log_op <- MLE_sw_log(c(-2,-1,-8),trees[[",k,"]],sim[1:",tipnums[k-10],",],m=20,al,mumat,protein_op=op_seq,root=start_seq,trace=1))",sep="")
       line7 <- paste("hes_op <- find_hessian_log(mle_log_op$par,trees[[",k,"]],sim[1:",tipnums[k-10],",],al,mumat,20,protein_op=op_seq,root=start_seq)",sep="")

       ## MLE without optimal and root sequences given
       line8 <- paste("system.time(mle_log_nop <- MLE_sw_log(c(-2,-1,-8),trees[[",k,"]],sim[1:",tipnums[k-10],",],m=20,al,mumat,trace=1))",sep="")
       line9 <- paste("hes_nop <- find_hessian_log(mle_log_nop$par,trees[[",k,"]],sim[1:",tipnums[k-10],",],al,mumat,20)",sep="")
       line10 <- paste('save.image(file="op_simu_',k+4,'_',i,'_',j,'.RData",compress=TRUE)',sep="")
       
       cat(line1,"\n",line2,"\n",line3,"\n",line4,"\n",line5,"\n",line6,"\n",line7,"\n",line8,"\n",line9,"\n",line10,"\n",file=paste("op_sim_",k+4,"_",i,"_",j,".R",sep=""))
      
       pbsCommands=paste('#!/bin/bash','#$ -cwd','#$ -o /dev/null','#$ -e /dev/null',sep="\n")
       pbsCommands=paste(pbsCommands,'\n#$ -q ',queue,sep="")
       pbsCommands=paste(pbsCommands,'#$ -M chaij@umail.iu.edu', '#$ -m beas', '#$ -S /bin/bash',sep="\n")
       pbsCommands=paste(pbsCommands,"\n","#$ -N OpSim_",k+4,"_",i,"_",j,sep="")
       pbsCommands=paste(pbsCommands,"\n","/data/apps/R/2.14.0/bin/R CMD BATCH op_sim_",k+4,"_",i,"_",j,".R",sep="")
       cat(pbsCommands,file=paste('runSimOp_',k+4,'_',i,'_',j,'.sh',sep=""),append=FALSE)
       system(paste('chmod u+x runSimOp_',k+4,'_',i,'_',j,'.sh',sep=""))
       system(paste('/opt/sge/bin/lx24-amd64/qsub runSimOp_',k+4,'_',i,'_',j,'.sh',sep=""))
     }
  }
}
