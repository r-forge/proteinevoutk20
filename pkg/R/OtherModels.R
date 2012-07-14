load("~/proteinevoutk20/pkg/Data/RokasAA.phyDat") #load Protein Data
## Read in protein data
##datAA <- read.nexus.data("rokasAAnoInterleave.nex")
tree <- read.nexus("gtrUnrooted.tre") # Read in unrooted tree
####################################################################
###Codon models, using codon data###
### w: dN/dS ratio; k: transition transversion ratio
### codon0: w=1 k=1
### codon1: w and k are both estimated
### codon2: k=1
### codon3: w=1
fit <- pml(tree,datCodon)
fit0 <- optim.pml(fit)
fit0 <- optim.pml(fit, optEdge=FALSE,control = pml.control(trace = 1))
fit1 <- optim.pml(fit, optEdge=FALSE,model="codon1", control=pml.control(trace=1))
fit2 <- optim.pml(fit, optEdge=FALSE,model="codon2", control=pml.control(trace=1))
fit3 <- optim.pml(fit, optEdge=FALSE,model="codon3", control=pml.control(trace=1))
####################################################################
### Amino acid models, using AA data###

mtU <- modelTest.User(datAA,tree,model=c("JTT","LG","WAG","Dayhoff","cpREV","mtmam","mtArt","Mtzoa","mtREV24"))
mt <- modelTest(datAA,tree,model=c("JTT","LG","WAG","Dayhoff","cpREV","mtmam","mtArt","Mtzoa","mtREV24"))

###fitStart=eval(get(mt$Model[which.min(mt$BIC)],env),env)

### optimization under 3 models ###
fitJTT <- pml(tree,datAA,model="JTT",k=4,inv=0.2)
fitJTT <- optim.pml(fitJTT,optInv=TRUE,optGamma=TRUE,optEdge=TRUE)
fitJTT
fitLG <- pml(tree,datAA,model="LG",k=4,inv=0.2)
fitLG <- optim.pml(fitLG,optInv=TRUE,optGamma=TRUE,optEdge=FALSE)
fitLG
fitWAG <- pml(tree,datAA,model="WAG",k=4,inv=0.2)
fitWAG <- optim.pml(fitWAG,optInv=TRUE,optGamma=TRUE,optEdge=FALSE)
fitWAG

### Other models: Dayhoff, cpREV, mtmam, mtArt, Mtzoa, mtREV24 ###

### Modified modelTest function where user can specify whether branch lengths are getting optimized ###
modelTest.User <- function (object, tree = NULL, model = c("JC", "F81", "K80","HKY", "SYM", "GTR"), optEdge=FALSE, G = TRUE, I = TRUE, k = 4, control = pml.control(epsilon = 1e-08, maxit = 3, trace = 1), multicore = FALSE) 
{
  if (class(object) == "phyDat") 
    data = object
  if (class(object) == "pml") {
    data = object$data
    if (is.null(tree)) 
      tree = object$tree
  }
  if (attr(data, "type") == "DNA") 
    type = c("JC", "F81", "K80", "HKY", "TrNe", "TrN", "TPM1", 
             "K81", "TPM1u", "TPM2", "TPM2u", "TPM3", "TPM3u", 
             "TIM1e", "TIM1", "TIM2e", "TIM2", "TIM3e", "TIM3", 
             "TVMe", "TVM", "SYM", "GTR")
  if (attr(data, "type") == "AA") 
    type = c("WAG", "JTT", "Dayhoff", "LG", "cpREV", "mtmam", 
             "mtArt", "MtZoa", "mtREV24")
  model = match.arg(model, type, TRUE)
  env = new.env()
  assign("data", data, envir = env)
  if (is.null(tree)) 
    tree = NJ(dist.hamming(data))
  trace <- control$trace
  control$trace = trace - 1
  fit = pml(tree, data)
  fit = optim.pml(fit, optEdge=optEdge, control = control)
  l = length(model)
  n = 1L + sum(I + G + (G & I))
  nseq = sum(attr(data, "weight"))
  fitPar = function(model, fit, G, I, k) {
    m = 1
    res = matrix(NA, n, 5)
    res = as.data.frame(res)
    colnames(res) = c("Model", "df", "logLik", "AIC", "BIC")
    data.frame(c("Model", "df", "logLik", "AIC", "BIC"))
    calls = vector("list", n)
    trees = vector("list", n)
    fittmp = optim.pml(fit, optEdge=optEdge, model = model, control = control)
    res[m, 1] = model
    res[m, 2] = fittmp$df
    res[m, 3] = fittmp$logLik
    res[m, 4] = AIC(fittmp)
    res[m, 5] = AIC(fittmp, k = log(nseq))
    calls[[m]] = fittmp$call
    trees[[m]] = fittmp$tree
    m = m + 1
    if (I) {
      fitI = optim.pml(fittmp, model = model, optInv = TRUE, optEdge=optEdge,
                       control = control)
      res[m, 1] = paste(model, "+I", sep = "")
      res[m, 2] = fitI$df
      res[m, 3] = fitI$logLik
      res[m, 4] = AIC(fitI)
      res[m, 5] = AIC(fitI, k = log(nseq))
      calls[[m]] = fitI$call
      trees[[m]] = fitI$tree
      m = m + 1
    }
    if (G) {
      fitG = update(fittmp, k = k)
      fitG = optim.pml(fitG, model = model, optGamma = TRUE, optEdge=optEdge,
                       control = control)
      res[m, 1] = paste(model, "+G", sep = "")
      res[m, 2] = fitG$df
      res[m, 3] = fitG$logLik
      res[m, 4] = AIC(fitG)
      res[m, 5] = AIC(fitG, k = log(nseq))
      calls[[m]] = fitG$call
      trees[[m]] = fitG$tree
      m = m + 1
    }
    if (G & I) {
      fitGI = optim.pml(fitG, model = model, optGamma = TRUE, optEdge=optEdge,
                        optInv = TRUE, control = control)
      res[m, 1] = paste(model, "+G+I", sep = "")
      res[m, 2] = fitGI$df
      res[m, 3] = fitGI$logLik
      res[m, 4] = AIC(fitGI)
      res[m, 5] = AIC(fitGI, k = log(nseq))
      calls[[m]] = fitGI$call
      trees[[m]] = fitGI$tree
      m = m + 1
    }
    list(res, trees, calls)
  }
  eval.success <- FALSE
  if (!eval.success & multicore) {
    if (!require(parallel) || .Platform$GUI != "X11") {
      warning("package 'parallel' not found or GUI is used, \n      analysis is performed in serial")
    }
    else {
      RES <- mclapply(model, fitPar, fit, G, I, k)
      eval.success <- TRUE
    }
  }
  if (!eval.success) 
    res <- RES <- lapply(model, fitPar, fit, G, I, k)
  RESULT = matrix(NA, n * l, 5)
  RESULT = as.data.frame(RESULT)
  colnames(RESULT) = c("Model", "df", "logLik", "AIC", "BIC")
  for (i in 1:l) RESULT[((i - 1) * n + 1):(n * i), ] = RES[[i]][[1]]
  for (i in 1:l) {
    for (j in 1:n) {
      mo = RES[[i]][[1]][j, 1]
      tname = paste("tree_", mo, sep = "")
      tmpmod = RES[[i]][[3]][[j]]
      tmpmod["tree"] = call(tname)
      if (!is.null(tmpmod[["k"]])) 
        tmpmod["k"] = k
      if (attr(data, "type") == "AA") 
        tmpmod["model"] = RES[[i]][[1]][1, 1]
      assign(tname, RES[[i]][[2]][[j]], envir = env)
      assign(mo, tmpmod, envir = env)
    }
  }
  attr(RESULT, "env") = env
  RESULT
}
