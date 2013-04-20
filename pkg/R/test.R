# #for a given vector of log likelihood, and the corresponding opaa, find the smallest set of aa
# # that cover the 95% of the total likelihood
# aa.set <- function(llvec){
#   lvec = exp(llvec) #likelihood
#   ord = order(llvec,decreasing=T) #order of llvec (increasing)
#   lvec = lvec[ord] # ordered likelihood
#   tol = sum(lvec) #total likelihood
#   CumSum = cumsum(lvec) # cumulative sum 
#   ind = sum(CumSum < tol*0.95) + 1 #the index that gives > 95% totol likelihood
#   Sum = CumSum[ind] # sum of the first "ind" terms
#   #return the set of amino acids that give > 95% likelihood, percentile, and 
#   #the percentile of the likelihood given by the optimal aa (max rule)
#   return(list(aa=ord[1:ind], percentile=Sum/tol,op.percentile=lvec[1]/tol))
# }
# aa.conf = apply(X=mat,MARGIN=1,FUN=aa.set) #here mat is the llmat from result
# numaa = sapply(1:9128,function(x) length(aa.conf[[x]]$aa))
# op.lik = sapply(1:9128, function(x) aa.conf[[x]]$op.per)
# 
# plot.grid <- function(index){
#   xyz = grid.gen(opw_mean,index,res_op=res_op,gridnum=20)
#   akima.xyz = interp(xyz$x,xyz$y,xyz$z)
#   image(akima.xyz)
#   contour(akima.xyz,add=TRUE)
#   points(xyz,pch=3,col="blue")
# }

#MLE for s and Ne, using subplex method by default
#mllm1 <- function(data,tree,s=NULL,beta=be,gamma=ga,Q=NULL,dismat=NULL,fixmatall=NULL,mumat=NULL,Qall=NULL,
#                  opaa=NULL,opw=NULL,bfaa=NULL,C=2,Phi=0.5,q=4e-7,Ne=5e6)
#sample call : optim.s.weight(gene1,ROKAS_TREE,0.1,be,ga,Q=NU_VEC))
optim.s.Ne <- function(data, tree,s,Ne, method="SBPLX",maxeval="500",print_level=0,...){
  #store information from initial condition, with other parameters fixed
  res.initial = mllm1(data=data,tree=tree,s=s,beta=be,gamma=ga,Ne=Ne,...)
  #these don't change with the change of s
  mumat = res.initial$mumat
  bfaa=res.initial$bfaa
  
  ab <- c(s,Ne) ##initial value
  fn = function(ab,data,tree){
    cat("s, Ne = ",ab[1]," ",ab[2],"\n")
    result = -mllm1(data=data,tree=tree,s=ab[1],beta=be,gamma=ga, mumat=mumat,bfaa=bfaa,Ne=ab[2], ...)$ll$loglik
    cat(result,"\n")
    return(result)
    
  }
  lower <- rep(0,2)
  upper <- rep(Inf,2)
  #options for optimizer
  opts <- list("algorithm"=paste("NLOPT_LN_",method,sep=""),"maxeval"=maxeval,"xtol_rel"=1e-6,
               "ftol_rel"=.Machine$double.eps^0.5,"print_level"=print_level)
  res = nloptr(x0=ab,eval_f=fn, lb=lower,ub=upper,opts=opts,data=data,tree=tree)
  return(res)
}