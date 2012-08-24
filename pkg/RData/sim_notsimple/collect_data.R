numsites <- c(100,300,500,700,1000)
svalues <- c(0.1,0.3,0.5,0.7,0.9,1)
tip <- c(4,4,8,8,12,12)
datares <- array(list(NULL),dim=c(6,5,6))
for(treeindex in 1:6){
  for(sitesindex in 1:5){
    for(sindex in 1:6){
      print(paste("tree ",treeindex,", ",numsites[sitesindex], " sites, s = ",svalues[sindex], sep=""))
      RDataName = paste("simMle_",treeindex,"_",numsites[sitesindex],"_",svalues[sindex],".RData",sep="")
      if(file.exists(RDataName)){
        print("file exists!")
        load(RDataName)
        datares[[treeindex,sitesindex,sindex]]$para <- c(svalues[sindex],be,ga)
        datares[[treeindex,sitesindex,sindex]]$mle <- mle$par
        rm(mle)
        if(exists("hes")){
          datares[[treeindex,sitesindex,sindex]]$hessian <- hes
          datares[[treeindex,sitesindex,sindex]]$eigenval <- eigen(hes)$values
          datares[[treeindex,sitesindex,sindex]]$invhes <- solve(hes)
          print("remove hessian now")
          rm(hes)
        }
        else{
          print("hessian not calculated")
##           datares[[treeindex,sitesindex,sindex]]$hessian <- NULL
##           datares[[treeindex,sitesindex,sindex]]$eigenval <- NULL
##           datares[[treeindex,sitesindex,sindex]]$invhes <- NULL
        }
      }
      else
        print("file doesn't exist!")
    }
  }
}
save(datares,file="collect_data.RData")
##       print(paste("tree ",treeindex,", ",sitesindex, " sites, s = ",sindex, sep=""))
##       print(" ")

##               print("mle for parameters: ")
##               print(mle$par)
##               print(" ")
##               print("Hessian:")
##               print(hes)
##               print(" ")
##               print("eigen values of Hessian matrix:")
##               print(eigen(hes)$values)
##               print(" ")
##               print("inverse of Hesssian Matrix:")
##               print(solve(hes))
