## convert a string that contains many numbers separated by blanks, to a vector of numbers
str.to.num <- function(string,split=" "){
  char.vec <- strsplit(string,split=split)[[1]]
  return(as.numeric(char.vec))
}
#############################################################################
# #for a given vector of log likelihood, and the corresponding opaa, find the smallest set of aa
# # that cover the 95% of the total likelihood
aa.set <- function(llvec){
  lvec = exp(llvec) #likelihood
  ord = order(llvec,decreasing=T) #order of llvec (increasing)
  lvec = lvec[ord] # ordered likelihood
  tol = sum(lvec) #total likelihood
  CumSum = cumsum(lvec) # cumulative sum 
  ind = sum(CumSum < tol*0.95) + 1 #the index that gives > 95% totol likelihood
  Sum = CumSum[ind] # sum of the first "ind" terms
  #return the set of amino acids that give > 95% likelihood, percentile, and 
  #the percentile of the likelihood given by the optimal aa (max rule)
  return(list(aa=ord[1:ind], percentile=Sum/tol,op.percentile=lvec[1]/tol))
}
# aa.conf = apply(X=mat,MARGIN=1,FUN=aa.set) #here mat is the llmat from result
# numaa = sapply(1:length(aa.conf),function(x) length(aa.conf[[x]]$aa))
# op.lik = sapply(1:length(aa.conf), function(x) aa.conf[[x]]$op.per)
# numaa_freq <- as.numeric(table(numaa))
# numaa_freq <- numaa_freq/sum(numaa_freq)


# heatmap of probabilities that each amino acid accounts for
hmap.op <- function(llmat){
  probmat <- exp(llmat)
  lsite <- dim(probmat)[1]
  #probmat <- t(sapply(1:lsite,function(x) expmat[x,]/sum(expmat[x,])))
  #heatmap(probmat[1:lsite,],Rowv=NA,Colv=NA,col=grey.colors(256))
  heatmap(t(probmat),Rowv=NA,Colv=NA,col=grey.colors(256))
}