##numer of sites in the simulation
numsites <- c(100,300,500,700,1000)
##s values used in simulation
svalues <- c(0.1,0.3,0.5,0.7,0.9,1)
##tips of trees used in simulation
tip <- c(4,4,8,8,12,12,10,10,14,14,16,16,18,18)
#tip <- c(4,4,8,8,12,12)
#tip <- c(10,10,14,14,16,16,18,18)
##result to store all the data
datares14 <- array(list(NULL),dim=c(14,5,6))
##change this to "simMle" for information on MLE without log
fileprefix <- "simMleLog_"
DataFile <- "collect_data_log.RData"
##fileprefix <- "simMle"
##DataFile <- "collect_data.RData"
for(treeindex in 1:14){
  for(sitesindex in 1:5){
    for(sindex in 1:6){
      print(paste("tree ",treeindex,", ",numsites[sitesindex], " sites, s = ",svalues[sindex], sep=""))
      RDataName = paste(fileprefix,treeindex,"_",numsites[sitesindex],"_",svalues[sindex],".RData",sep="")
      if(file.exists(RDataName)){
        print("file exists!")
        load(RDataName)
        datares14[[treeindex,sitesindex,sindex]]$para <- log(c(svalues[sindex],be,ga))
        datares14[[treeindex,sitesindex,sindex]]$mle <- mle_log$par
        rm(mle)
        if(exists("hes")){
          datares14[[treeindex,sitesindex,sindex]]$hessian <- hes
          datares14[[treeindex,sitesindex,sindex]]$eigenval <- eigen(hes)$values
          datares14[[treeindex,sitesindex,sindex]]$invhes <- solve(hes)
          print("remove hessian now")
          rm(hes)
        }
        else{
          print("hessian not calculated")
##           datares14[[treeindex,sitesindex,sindex]]$hessian <- NULL
##           datares14[[treeindex,sitesindex,sindex]]$eigenval <- NULL
##           datares14[[treeindex,sitesindex,sindex]]$invhes <- NULL
        }
      }
      else
        print("file doesn't exist!")
    }
  }
}
save(datares14,file=DataFile)

