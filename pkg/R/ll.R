svalues <- c(0.1,0.3,0.5,0.7,0.9)
sites <- c(500,1000,2000,2500,3000)
for(siteval in sites){
  for(sval in svalues){
    for(op in 1:10){
      sim <- simTree(ROKAS_TREE,rep(op,siteval),sval,NU_VEC,alpha=al,beta=0.1,ga=0.003)
      AAdata <- phyDat(sim$data,type="AA")
      sw_est <- optim.s.weight(AAdata,ROKAS_TREE,sval,be,ga,trace=1,Q=NU_VEC,bfaa=rep(1/20,20))
      sw_op <- optim.s.weight(AAdata,ROKAS_TREE,sval,be,ga,trace=1,Q=NU_VEC,bfaa=rep(1/20,20),opaa=protein_op)
      sw_opw <- optim.sw.ops(AAdata,ROKAS_TREE,sval,be,ga,maxit=1000,trace=1,Q=NU_VEC,bfaa=rep(1/20,20))
      save.image(file=paste("sim_",tind,"_",siteval,"_",sval,"_",elscale,".RData",sep=""))
    }
  }
}
