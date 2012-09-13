##numer of sites in the simulation
numsites <- c(100,300,500,700,1000)
##s values used in simulation
svalues <- c(0.1,0.3,0.5,0.7,0.9,1)
##tips of trees used in simulation
tip <- c(4,4,8,8,12,12,10,10,14,14,16,16,18,18)
##result to store all the data
datares <- array(list(NULL),dim=c(14,5,6))
##change this to "simMle" for information on MLE without log
fileprefix <- "op_simu_"
DataFile <- "op1_short.RData"

for(treeindex in 1:14){
  for(sitesindex in 1:5){
    for(sindex in 1:6){
      print(paste("tree ",treeindex,", ",numsites[sitesindex], " sites, s = ",svalues[sindex], sep=""))
      RDataName = paste(fileprefix,treeindex,"_",numsites[sitesindex],"_",svalues[sindex],".RData",sep="")
      if(file.exists(RDataName)){
        load(RDataName)
        datares[[treeindex,sitesindex,sindex]]$para <- log(c(svalues[sindex],be,ga))
        datares[[treeindex,sitesindex,sindex]]$mle_op <- mle_log_op$par
        datares[[treeindex,sitesindex,sindex]]$mle_nop <- mle_log_nop$par
        datares[[treeindex,sitesindex,sindex]]$hes_op <- hes_op
        datares[[treeindex,sitesindex,sindex]]$hes_nop <- hes_nop
        }
      else
        print("file doesn't exist!")
    }
  }
}
fileprefix <- "op_simu_"
for(treeindex in 15:18){
  for(sitesindex in 1:5){
    for(sindex in 1:6){
      print(paste("tree ",treeindex,", ",numsites[sitesindex], " sites, s = ",svalues[sindex], sep=""))
      RDataName = paste(fileprefix,treeindex,"_",numsites[sitesindex],"_",svalues[sindex],".RData",sep="")
      if(file.exists(RDataName)){
        print("file exists!")
        load(RDataName)
        datares[[treeindex,sitesindex,sindex]]$para <- log(c(svalues[sindex],be,ga))
        datares[[treeindex,sitesindex,sindex]]$mle_op <- mle_log_op$par
        datares[[treeindex,sitesindex,sindex]]$mle_nop <- mle_log_nop$par
        datares[[treeindex,sitesindex,sindex]]$hes_op <- hes_op
        datares[[treeindex,sitesindex,sindex]]$hes_nop <- hes_nop
      }
      else
        print("file doesn't exist!")
    }
  }
}
save(datares,file=DataFile)

