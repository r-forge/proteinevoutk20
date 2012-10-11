library(seqinr) #package to convert codons to amino acids
library(phangorn) #phylogenetics package
library(optimx) #optimization package
library(expm) #matrix exponentiation
#library(parallel) #parallel computing
library(multicore)
library(minqa) #optimization function bobyqa
library(mgcv) #find the unique rows in a matrix - uniquecombs
library(numDeriv) # calculate hessian of a function at a point
##################################################################
#directory of the data need to read, including gene, tree, Grantham data
datadir <- "~/proteinevoutk20/pkg/Data/"
##################################################################
##   Read in parameters and data  ##
##################################################################
##Nucleotide mutation rates estimated from Rokas's data (using PAUP):
##Model GTR, with Gamma rate distribution and invariable sites and etc
##Nuvec <-  c(2.94194,8.23705,1.39133,2.33991,14.86435,1.00000)
##Model GTR
NU_VEC <- c(1.96575,4.08655,1.39431,1.46172,6.36024,1.00000)
NU_VEC_JC <- rep(1.0,6)
##empirical base frequencies
BF <- c(0.31065,0.18365,0.20955,0.29615)
##Best tree with branch lengths output by PAUP
ROKAS_TREE <- read.nexus(paste(datadir,"GTR.tre",sep=""))
##Read the properties of amino acids, including c(composition),p(polarity) and 
##v(molecular volume), save the data.frame
GM_CPV <- read.csv(paste(datadir,"Grantham_cpv.csv",sep=""),header=TRUE, sep=",",row.names=1)
GM_CPV <- data.matrix(GM_CPV)

##################################################################
##  Grantham Matrix     ##
##################################################################
##Build the Grantham matrix with coefficients for all components as parameters
##alpha for composition,beta for polarity,gamma for molecular volume
##The mean of the distance is normalized to 1,(mean of only the upper triangle elements of the matrix)
##The parameter values in Grantham matrix are:
##alpha=1.833, beta=0.1018,gamma=0.000399
##First find the avg chemical distance, and then the weight is (1/avg_dis)^2

##this function computes the weights that are used by Grantham
## (dispite the name, they are not the average distances though, instead they are functions of average distances)
avg_dis <- function(vector){
  l <- length(vector)
  Sum <- 0
  ## find the sum of all pairwise chemical distances
  for(i in 1:(l-1)){
    Sum <- Sum + sum(abs(vector[i]-vector[(i+1):l]))
  }
  ct <- l*(l-1)/2 #total combination count
  return((ct/Sum)^2)
}
##Grantham weights
al <- avg_dis(GM_CPV[,1])
be <- avg_dis(GM_CPV[,2])
ga <- avg_dis(GM_CPV[,3])


###Given the data (20 by 3 matrix/data frame) from Grantham on all three factors, and the weights,
###Find the distance matrix for all the amino acids (20 by 20) with mean 1
###Grantham matrix has mean 100 instead of 1
GM_cpv <- function(datamatrix, alpha=al, beta=be, gamma=ga){
  A <- datamatrix
  DistanceMatrix <- matrix(0,nrow=20,ncol=20) #diagonal entries have values 0
  for(i in 1:19){
    for(j in (i+1):20)
      DistanceMatrix[i,j] <- (alpha*(A[i,1] - A[j,1])^2 + beta*(A[i,2]-A[j,2])^2+gamma*(A[i,3]-A[j,3])^2)^(1/2)
  }
  ##make a symmetric matrix
  result <- DistanceMatrix + t(DistanceMatrix)
  ##normalize the mean to 1
  m <- sum(result)/(400-20) #exclude the diagonal entries
  return(result/m)
}
##Default Grantham matrix (one used in Grantham's paper)
GM <- GM_cpv(GM_CPV,al,be,ga)


##################################################################
##   construct matrices for mutation and substitution  ##
##################################################################
##Short names for all the amino acids, in alphabetical order
AA <- s2c("ARNDCQEGHILKMFPSTWYV") #amino aicds in alphabetical order
Nu <- s2c("acgt") #nucleotide ACGT - now the order is changed !!!
stop_cd <- c("taa","tag","tga") #stop codons
nu_list <- NULL #list of all 64 codons (including stop codons)
for(i in 1:4){
  for(j in 1:4){
    for(k in 1:4){
      nu_list <- c(nu_list,c2s(c(Nu[i],Nu[j],Nu[k])))
    }
  }
}
CDS <- nu_list[!nu_list %in% stop_cd] #61 non-stop codons
rm(nu_list)
rm(stop_cd)
CDS_AA = sapply(1:61, function(x) translate(s2c(CDS[x]))) ## amino acids coded by 61 codons in order
## list of 20, each includes all the codons coding for the particular amino acid
cdlist <- list()
for(i in 1:20){
  cdlist[[i]] = CDS[CDS_AA==AA[i]]
}
names(cdlist) <- AA
##################################################################
##   Read gene data ##
##################################################################
## read in fasta file of nucleotide, convert it into amino acid data
## type = "num" -- integers ; "AA" -- amino acid names; "phyDat" --- phyDat data type used in phangorn
conv <- function(filename,type="num"){
  levels <- AA #amino acids in alphabeticla order (not the single letter names, the 3-letter names)
  data <- read.fasta(file=filename) #read fasta file including nucleotide data
  dna.to.aa <- function(x){ #given the order of the data, convert nucleotide data to amino acid data
    aas <- translate(seq=data[[x]])
    if(type=="num"){
      aas <- as.numeric(factor(aas,levels)) #convert from letters to numbers, as levels
      return(aas)
    }
    return(aas)
  }
  ## type = AA, return a list of amino acid sequences
  seqdata <- lapply(1:length(data),dna.to.aa)
  names(seqdata) = attr(data,"name")
  if(type=="AA")
    return(seqdata)
  if(type=="phyDat") ## convert amino acid sequences to phyDat type
    seqdata = phyDat(seqdata,type="AA")
  return(seqdata)
}

fasta.to.nex <- function(geneNum){
  for(i in geneNum){
    file = paste("~/proteinevoutk20/pkg/Data/gene",i,".fasta",sep="")
    gene = conv(file,"AA")
    write.nexus.data(gene,paste("~/proteinevoutk20/pkg/Data/gene",i,"AA.nex",sep=""),format="protein",interleaved=F)
  }
}
##   Read in 106 gene data  from Rokas's ##
# l <- 106
# ROKAS_DATA <- vector("list",length=l)
# ##data stores the data of all genes in rokas's data
#   for(i in 1:l){
#     ROKAS_DATA[[i]] <- conv(paste(datadir,"gene",i,".fasta",sep=""))
#   }
#############################################################################
## Remove columns with NA in the data
PruneMissing <- function(x){
  ##Find the array indices of the NA entries
  naInd <- which(is.na(x),arr.ind=T)
  dimnames(naInd) <- NULL
  ## column indices
  naCol <- naInd[,2]
  x[,-naCol]
}
#############################################################################
## find the empirical base frequencies for AMINO ACID list data (not phyDat data)
findBf <- function(datalist){
  data = unlist(datalist)
  datatb = table(data)
  datatb = datatb[AA]
  datatb = as.numeric(datatb)
  return(datatb/sum(datatb))
}
## for phyDat format
findBf2 <- function(data){
  ind = attr(data,"index")
  l = length(data)
  datalist = vector("list",length=l)
  for(i in 1:l)
    datalist[[i]] = data[[i]][ind]
  data = unlist(datalist)
  datatb = table(data)
  datatb = as.numeric(datatb)
  return(datatb/sum(datatb))
}
## Find the most frequent element of a vector
## It works for numbers and characters
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

#############################################################################
##given lower triangular part of R (i.e. Q = R %*% diag(bf)), and base frequencies, find the scaled Q
## Default: 4 by 4 matrix for Jukes-Cantor model
mat_form_lowtriQ <- function(Q=rep(1,6),bf=rep(0.25,4)){
  if(sum(bf==0)) stop("base frequencies can't all be 0!")
  l=length(bf)
  bf=bf/sum(bf)
  res = matrix(0, l, l)
  res[lower.tri(res)] = Q
  res = res+t(res) #symmetric matrix with diagonals 0
  res = res * rep(bf,each=l) #multiply cols by bf
  diag(res) = -rowSums(res) #set row sum to 0
  res2 = res * rep(bf,l) #multiply rows by bf
  diag(res2)=0 
  res = res/sum(res2) #normalize rate matrix
  return(res)
}
## given matrix Q and bf, find scaled Q
## Goal: sum(pi_i * Q_ii) = -1
## and return its eigen values, eigen vectors, and inverse of eigen vector matrix
## These are used to calcualte matrix exponentiation later.
eigQ <- function(Q,bf=NULL){
  l = dim(Q)[1]
  if(is.null(bf)) bf = rep(1/l,l)
  bf=bf/sum(bf)
  scalefactor = - sum(diag(Q)*bf) 
  Q = Q/scalefactor
  e = eigen(Q,FALSE)
  e$inv = solve(e$vectors)
  return(e)
}
#scale Q w.r.t bf, if bf is not given, use equal frequencies
# bf doesn't have sum to 1
scaleQ <- function(Q,bf=NULL){
  l = dim(Q)[1]
  if(is.null(bf)) bf = rep(1/l,l)
  bf = bf/sum(bf)
  scalef = -sum(diag(Q)*bf)
  Q = Q/scalef
  return(Q)
}
## from symmetric matrix of rates, find GTR rate matrix
sym.to.Q <- function(A,bf=NULL){
  if(!isSymmetric(A)) stop("matrix is not symmetric!")
  l = dim(A)[1] #dimension of rate matrix i.e. number of states
  if(is.null(bf)) bf = rep(1/l,l)
  if(sum(bf==0)) stop("base frequencies can't all be 0!")
  bf=bf/sum(bf)
  diag(A) = 0
  A = A * rep(bf,each=l) #multiply cols by bf
  diag(A) = -rowSums(A) #set row sum to 0
  res2 = A * rep(bf,l) #multiply rows by bf
  diag(res2)=0 
  A = A/sum(res2) #normalize rate matrix
  return(A)
}
#############################################################################
## given the base frequencies of nucleotides, find base frequencies of codons
## assuming the codon positions are all independent. Stop codons are excluded
freq_codon <- function(bf = rep(1/4,4)){
  bf=bf/sum(bf)
  freq = rep(0,61)
  for(i in 1:61){
    cd = as.numeric(factor(s2c(CDS[i]),Nu))
    freq[i] = prod(bf[cd])
  }
  freq = as.table(freq)
  freq = freq/sum(freq)
  dimnames(freq)[[1]] = CDS
  return(freq)
}
##given the nucleotide base frequencies, find the bf for amino acids
freq_aa <- function(nubf=rep(1/4,4)){
  bf = freq_codon(nubf) #codon frequencies
  freq = matrix(0,20,1)
  dimnames(freq)[[1]] = AA
  CDS_AA = sapply(1:61, function(x) translate(s2c(CDS[x])))
  for(i in 1:20)
    freq[AA[i],1] = sum(bf[CDS_AA==AA[i]])
  freq = as.table(as.vector(freq))
  dimnames(freq)[[1]] = AA
  return(freq)
}
## Given the rate matrix for 61 codons, and the base frequencies of codons, 
## find the rate matrix for 20 amino acids.Only mutation is taken into account, not selection and fixation.
AAMat <- function(CdMat){
  Q = matrix(0,20,20)
  for(i in 1:19){
    fcodons = cdlist[[i]]
    for(j in (i+1):20){
      tcodons = cdlist[[j]]
      Q[i,j] = sum(CdMat[fcodons,tcodons])
    }
  }
  #diag(Q) = -rowSums(Q)
  dimnames(Q) = list(AA,AA)
  return(Q)
}

## compare 2 codons (in character vectors), return the number of positions that differ and the different positions
cdcmp <- function(cd1,cd2){
  cmp = (cd1 != cd2)
  num = sum(cmp)
  pos = which(cmp==TRUE)
  return(list(num=num,pos=pos))
}
#Find the mutation rate matrix A (61 by 61) from the 6 rates in the symmetric generating matrix for GTR's Q
# # result is a symmetric matrix with row sum 0
cd_MuMat_form <- function(vec=rep(1,6)){
  l = 4
  Q = matrix(0, l, l)
  Q[lower.tri(Q)] = vec
  Q = Q+t(Q) #symmetric matrix with diagonals 0
  dimnames(Q) = list(Nu,Nu)
  ##mutation matrix for 61 codons
  codon_array <- array(0,dim=c(61,61),dimnames=list(CDS,CDS))
  for(m in 1:60){                       #loop through all 61 codons
    fcd <- CDS[m]
    fcd_ch <- s2c(fcd)
    for(i in (m+1):61){
      tcd = CDS[i]
      tcd_ch = s2c(tcd)
      cmp = cdcmp(fcd_ch,tcd_ch) #compare the 2 codons
      if(cmp$num==1){
        pos = cmp$pos
        codon_array[fcd,tcd] = Q[fcd_ch[pos],tcd_ch[pos]]
      }
    }
  }
  codon_array = codon_array + t(codon_array)
  diag(codon_array) = -rowSums(codon_array)
  return(codon_array)
}

#Find the mutation rate matrix A (20 by 20) from the 6 rates in GTR mutation rate matrix for nucleotides, row sum is 0
# result is a symmetric matrix with row sum 0
aa_MuMat_form <- function(vec=rep(1,6)){
  l = 4
  Q = matrix(0, l, l)
  Q[lower.tri(Q)] = vec
  Q = Q+t(Q) #symmetric matrix with diagonals 0
  dimnames(Q) = list(Nu,Nu)
  ##mutation matrix for 61 codons
  codon_array <- array(0,dim=c(61,61),dimnames=list(CDS,CDS))
  for(m in 1:60){                       #loop through all 61 codons
    fcd <- CDS[m]
    fcd_ch <- s2c(fcd)
    for(i in (m+1):61){
        tcd = CDS[i]
        tcd_ch = s2c(tcd)
        cmp = cdcmp(fcd_ch,tcd_ch) #compare the 2 codons
        if(cmp$num==1){
          pos = cmp$pos
          codon_array[fcd,tcd] = Q[fcd_ch[pos],tcd_ch[pos]]
        }
      }
    }
  codon_array = codon_array + t(codon_array)
  arr = AAMat(codon_array)
  arr = arr+ t(arr)
  diag(arr) = -rowSums(arr)
  return(arr)
}
# ##check if the transition between amino acids is time reversible
# checkSymmetric <- function(nuvec,bf){
#   bfq = freq_codon(bf)
#   bfaa = freq_aa(bf)
#   aaq = aa_MuMat_form(nuvec,bf)
#   return(isSymmetric(diag(bfaa) %*% aaq))
# }
## Amino acid mutation rate matrix for the GTR model of nucleotide
#MUMAT <- aa_MuMat_form(NU_VEC)
#MUMAT_JC <- aa_MuMat_form()

reorderPruning <- function (x, ...)
{
  parents <- as.integer(x$edge[, 1])
  child <- as.integer(x$edge[, 2])
  root <- as.integer(parents[!match(parents, child, 0)][1])  # unique out                                                                        
  if (length(root) > 2)
    stop("more than 1 root found")
  n = length(parents)
  m = max(x$edge)  # edge  parents                                                                                                               
  neworder = .C("reorder", parents, child, as.integer(n), as.integer(m), integer(n), as.integer(root-1L), DUP=FALSE, PACKAGE = "phangorn")[[5]]
  x$edge = x$edge[neworder,]
  x$edge.length = x$edge.length[neworder]
  attr(x, "order") <- "pruningwise"
  x
}

#############################################################################
#############################################################################
##find the physiochemical distance vector between two proteins, given the distance matrix
pchem_d <- function(protein1, protein2,DisMat){
  if(length(protein1)!=length(protein2)) #throw error if length of proteins are not the same
    stop("error: 2 proteins are of different lengths!")
  site_d <- function(k){
    return(DisMat[protein1[k],protein2[k]])
  }
  d <- sapply(c(1:length(protein1)),site_d,simplify=TRUE) 
  return(d)
}
#################################################################
#Functionality function, arguments: 
#d: vector of distances between the optimal and observed amino acids
#s: seleciton coeffecient, can be passed as a vector or a number if it's the same for all sites
Ftny <- function(d, s){  
  if((length(s)==1)&&(length(d)!=1)){ #if s is given as a scalar, then treat it to be the same across all sites
    s <- rep(s,length(d))
  }
  result <- length(d)/(sum(1+d*s)) #harmonic mean of ftny at all sites (F_i = 1 + d_i * s_i)
  return(result)
} #d and s are given

#Functionality of a protein given the optimal protein (aa sequence) and parameters
Ftny_protein <- function(protein,protein_op,s,DisMat){
  d <- pchem_d(protein,protein_op,DisMat)
  return(Ftny(d,s))
}

#fixation probability, since we consider only site by site, so n=1
#d1, d2: distances (vectors)
#s: selection (vector or scalar)
#Ne: effective population size
#q, Phi: constants
#C: cost function, linear function of the length of the protein
#C_n = a1 + a2*n, cost --  linear function of the length
fix <- function(d1,d2,s,C=2,Phi=0.5,q=4e-7,Ne=1.36e7){
  if((d1==d2)||(s==0)) #When the fitnesses are the same, neutral case, pure drift
    return(1/(2*Ne))
  else{
    fit_ratio <- exp(-C*Phi*q*s*(d1-d2))
    if(fit_ratio==Inf)
      return(0)
    else if(fit_ratio==1)
      return(1/(2*Ne))
    else
      return((1-fit_ratio)/(1-fit_ratio^(2*Ne)))
  }
}


#fixation probability of a mutant with initial freq 1/(2Ne) in a diploid population
#with proteins as parameteres instead of distance vectors for the proteins
fix_protein <- function(protein1, protein2, protein_op, s, DisMat,
                        C=2,Phi=0.5,q=4e-7,Ne=1.36e7){
  d1 <- pchem_d(protein1,protein_op,DisMat)
  d2 <- pchem_d(protein2,protein_op,DisMat)
  return(fix(d1,d2,s,C,Phi,q,Ne))
}

######################################################
#Function to generate the transition matrix (m*m for now) 
#given the optimal amino acid. This is only for one site (amino acid)
#parameter: aa_op: optimal amino acid at this site
#s: selection coefficient(SCALAR, since it's only one site)
#m: number of states (amino acids) considered
#DisMat: distance matrix between amino acids
#MuMat: mutation rate matrix between amino acids
mat_gen_indep <- function(aa_op,s,DisMat,MuMat,C=2, Phi=0.5,q=4e-7,Ne=1.36e7){
  m = 20
  mat <- matrix(0,nrow=m,ncol=m) #set diagonal entries to be 0 at first
  for(i in 1:(m-1)){
    for(j in (i+1):m){
      if(MuMat[i,j]!=0){
        mat[i,j] <- 2*Ne*MuMat[i,j]*fix_protein(i,j,aa_op,s,DisMat,C,Phi,q,Ne) #fixation prob -> transition rate
        mat[j,i] <- 2*Ne*MuMat[j,i]*fix_protein(j,i,aa_op,s,DisMat,C,Phi,q,Ne) #symmetric entry
      }
    }#end for j
  }#end for i
  diag(mat) <- -rowSums(mat)
  return(mat)
}
# given optimal amino acid, find the rate matrix Q and the eigen-decomposition of Q
eigOp <- function(op_aa,s,DisMat,MuMat,bf=rep(1/20,20),C=2,Phi=0.5,q=4e-07,Ne=1.36e07){
  mat = mat_gen_indep(op_aa,s,DisMat,MuMat,C,Phi,q,Ne)
  res = eigQ(t(mat),bf) #in order to use functions (C and R) from phangorn
  return(res)
}
## given s, distance matrix, and mutation rate matrix, find rate matrices for all 20 amino acids as optimal
## And then find their eigen decomposition 
eigAllaa <- function(s,DisMat,MuMat,bf=rep(1/20,20),C=2, Phi=0.5,q=4e-7,Ne=1.36e7){
  res = vector(mode="list",length=20)
  for(i in 1:20){
    #print(paste("now ", i, " is the optimal amino acid"))
    res[[i]] = eigOp(i,s,DisMat,MuMat,bf,C,Phi,q,Ne)
  }
  return(as.matrix(res))
}
# list of 20 rate matrices, one of each aa as optimal
QAllaa <- function(s,DisMat,MuMat,C=2, Phi=0.5,q=4e-7,Ne=1.36e7){
  res = lapply(1:20,function(i) {mat_gen_indep(i,s,DisMat, MuMat,C, Phi, q, Ne)})
  return(as.matrix(res))
}
# given the edge lengths el, rate matrix Q and rate g, find the probability transition matrices P
# result is of dimension length(g)*length(el)
# to parse an object from getPm, say P, do P[i,j][[1]], or P[[i,j]]
# Q is a scaled matrix
getPm <- function(el, Q, g){
  ell = length(el)
  gl = length(g)
  res = array(list(NULL),dim=c(gl,ell))
  for(i in 1:gl){
    for(j in 1:ell){
      res[[i,j]] = expm(Q*el[j]*g[i])
    }
  }
  return(res)
}

ll3 <- function (dat1, tree, bf = c(0.25, 0.25, 0.25, 0.25), g = 1,
                 Q = c(1, 1, 1, 1, 1, 1), eig = NULL, assign.dat = FALSE,
                 ...)
{
  if (is.null(attr(tree, "order")) || attr(tree, "order") ==
    "cladewise")
    tree <- reorderPruning(tree)
  q = length(tree$tip.label)
  node <- tree$edge[, 1]
  edge <- tree$edge[, 2]
  m = length(edge) + 1
  dat = vector(mode = "list", length = m)
  dat[1:q] = dat1[tree$tip.label]
  if (is.null(eig))
    eig = edQt(bf = bf, Q = Q)
  el <- tree$edge.length
  P <- getP(el, eig, g)
  nr <- as.integer(attr(dat1, "nr"))
  nc <- as.integer(attr(dat1, "nc"))
  node = as.integer(node - min(node))
  edge = as.integer(edge - 1)
  nTips = as.integer(length(tree$tip))
  mNodes = as.integer(max(node) + 1)
  contrast = attr(dat1, "contrast")
  nco = as.integer(dim(contrast)[1])
  res <- .Call("LogLik4", dat1[tree$tip.label], P, nr, nc, node, edge, nTips, mNodes, contrast, nco, PACKAGE = "phangorn")
  result = res[[2]][[1]] + log(res[[1]][[1]] %*% bf)
  if (assign.dat) {
    dat[(q + 1):m] <- res
    attr(dat, "names") = c(tree$tip.label, as.character((q +
      1):m))
    assign("asdf", dat, envir = parent.frame(n = 1))
  }
  result
}


ll3m <- function (dat1, tree, bf = rep(1/20,20), g = 1, 
                  Q, assign.dat = FALSE, 
                  ...) 
{
  if (is.null(attr(tree, "order")) || attr(tree, "order") == 
    "cladewise") 
    tree <- reorderPruning(tree)
  q = length(tree$tip.label)
  node <- tree$edge[, 1]
  edge <- tree$edge[, 2]
  m = length(edge) + 1
  dat = vector(mode = "list", length = m)
  bf = bf/sum(bf)
  Q = scaleQ(Q,bf)
  Q = t(Q) # in order to use C function in phangorn package
  el <- tree$edge.length
  P <- getPm(el,Q,g)
  nr <- as.integer(attr(dat1, "nr"))
  nc <- as.integer(attr(dat1, "nc"))
  node = as.integer(node - min(node))
  edge = as.integer(edge - 1)
  nTips = as.integer(length(tree$tip))
  mNodes = as.integer(max(node) + 1)
  contrast = attr(dat1, "contrast")
  nco = as.integer(dim(contrast)[1])
  res <- .Call("LogLik4", dat1[tree$tip.label], P, nr, nc, node, edge, nTips, mNodes, contrast, nco, PACKAGE = "phangorn")
  result = res[[2]][[1]] + log(res[[1]][[1]] %*% bf)     
  if (assign.dat) {
    dat[(q + 1):m] <- res
    attr(dat, "names") = c(tree$tip.label, as.character((q + 
      1):m))
    assign("asdf", dat, envir = parent.frame(n = 1))
  }
  result
}
# loglikelihood given optimal amino acids at each site. This funciton is actually not necessary, since it takes longer than 
# calculate all the optimal aa possibilities.
llop <- function(data,tree,op=NULL,bf=rep(1/20,20),Qall,g=1,C=2,Phi=0.2,q=4e-7,Ne=1.36e7){
  result = NULL
  weight = attr(data,"weight")
  nr = attr(data,"nr") #number of different sites
  index = attr(data,"index")
  ns = length(index) #number of sites
  opaa = NULL #optimal amino acids
  optimal_aa = "" #is optimal aa given or not?
  for(i in 1:20){ #when optimal aa is i
    llopi = ll3m(data,tree,bf=bf,g=1,Q=Qall[[i]])
    result = cbind(result,llopi)
  }
  if(!is.null(op)){
    opaa = op
    optimal_aa = "given"
    sitelik = sapply(1:ns,function(i) result[index[i],op[i]])
    loglik = sum(sitelik)
  }
  else{
    optimal_aa = "estimated using maximization"
    opaa = apply(result,1,which.max)
    sitelik=apply(result,1,max)
    loglik = sum(weight*sitelik)
  }
  return(list(loglik= loglik, optimal_aa = optimal_aa,opaa=opaa,sitelik = sitelik, llmat=result))
}


## For all the distinct sites (m), find loglikelihood, result is m * 20 matrix
## each column stores the loglikelihods when the corresponding amino acid is optimal 
llaa <- function(tree,data,eigAll,bf=rep(1/20,20),C=2,Phi=0.2,q=4e-7,Ne=1.36e7){
  result = NULL
  weight = attr(data,"weight")
  for(i in 1:20){ #when optimal aa is i
    llopi = ll3(data,tree,bf=bf,g=1,eig=eigAll[[i]])
    result = cbind(result,llopi)
  }
  opaa = apply(result,1,which.max)
  sitelik = apply(result,1,max)
  loglik = sum(weight * sitelik)
  return(list(loglik= loglik, llmat=result))
}
#this one uses expm for matrix exponentiation
llaam <- function(tree,data,QAll,bf=rep(1/20,20),C=2,Phi=0.2,q=4e-7,Ne=1.36e7){
  result = NULL
  weight = attr(data,"weight")
  for(i in 1:20){ #when optimal aa is i
    llopi = ll3m(data,tree,bf=bf,g=1,Q=QAll[[i]])
    result = cbind(result,llopi)
  }
  opaa = apply(result,1,which.max)
  sitelik = apply(result,1,max)
  loglik = sum(weight * sitelik)
  res=list(loglik= loglik, opaa = opaa, sitelik = sitelik, llmat=result)
  return(res)
}
# only return the big matrix, with number of rows equal to number of different sites, and cols equal to 20
llaam1 <- function(tree,data,QAll,bf=rep(1/20,20),C=2,Phi=0.2,q=4e-7,Ne=1.36e7){
  result = NULL
  weight = attr(data,"weight")
  for(i in 1:20){ #when optimal aa is i
    llopi = ll3m(data,tree,bf=bf,g=1,Q=QAll[[i]])
    result = cbind(result,llopi)
  }
  return(result)
}
# assume every amino acid has a weight to be the optimal one, and the weights are the same for all sites, calculate the loglikelihood
# if no weights are specified, use the aa's that maximize the likelihoods (which is the same as above)
llaaw <- function(tree,data,QAll,opw=NULL,bf=rep(1/20,20),C=2,Phi=0.2,q=4e-7,Ne=1.36e7){
  result = NULL
  weight = attr(data,"weight")
  nr = attr(data,"nr")
  opaa = NULL
  for(i in 1:20){ #when optimal aa is i
    llopi = ll3m(data,tree,bf=bf,g=1,Q=QAll[[i]])
    result = cbind(result,llopi)
  }
  if(!is.null(opw)){
    if(sum(opw)!=1) opw = opw/sum(opw)
    result = result * rep(opw,each=nr)
    sitelik = apply(result,1,sum)
    loglik = sum(weight * sitelik)
  }
  else{
    opaa = apply(result,1,which.max)
    sitelik=apply(result,1,max)
    loglik = sum(weight*sitelik)
  }
  return(list(loglik= loglik, opaa=opaa,weights=opw, sitelik = sitelik, llmat=result))
}

## find the loglikelihood given 20 by 20 rate matrix Q and base frequencies bf
llaaQm <- function(data, tree, Q, bf = rep(1/20,20),C=2, Phi=0.2,q=4e-7,Ne=1.36e7){
  result = NULL
  weight = attr(data,"weight")
  ll = ll3m(data,tree,bf=bf,g=1,Q=Q)
  loglik = sum(weight*ll)
  return(list(loglik=loglik,Q=Q,bf=bf))
}
## find the loglikelihood given Q and other paramters, here Q is the lower triangular part of the 
## nucleotide transition rate matrix of length 6
mll <- function(data,tree,s,beta,gamma,Q=NULL,bfaa=NULL,C=2,Phi=0.5,q=4e-7,Ne=1.36e7){
  call <- match.call()
  if(class(tree)!="phylo") stop("tree must be of class phylo") 
  if (is.null(attr(tree, "order")) || attr(tree, "order") == 
    "cladewise") 
    tree <- reorderPruning(tree)
  if (class(data)[1] != "phyDat") stop("data must be of class phyDat")
  if(is.null(tree$edge.length)) stop("tree must have edge weights") 
  if(any(is.na(match(tree$tip, attr(data, "names"))))) stop("tip labels are not in data")
  if(is.null(Q)) Q = rep(1,6)
  if(is.null(bfaa)) bfaa = findBf2(data)
  dismat = GM_cpv(GM_CPV,al,beta,gamma)
  mumat = aa_MuMat_form(Q)
  mumat = sym.to.Q(mumat,bfaa)
  eigall = eigAllaa(s,dismat,mumat,bf=bfaa,C,Phi,q,Ne)
  ll = llaa(tree,data,eigall,bfaa,C,Phi,q,Ne)
  result = list(ll=ll,data=data,tree=tree,s=s,GMweights=c(al,beta,gamma),Q=Q,bfaa=bfaa,call=call)
  class(result) = "mllm"
  return(result)
}
## this one uses expm for matrix exponentiation
mllm <- function(data,tree,s=0,beta=be,gamma=ga,Q=NULL,dismat=NULL,mumat=NULL,opaa=NULL,opw=NULL,bfaa=NULL,C=2,Phi=0.5,q=4e-7,Ne=1.36e7){
  call <- match.call()
  if(class(tree)!="phylo") stop("tree must be of class phylo") 
  if (is.null(attr(tree, "order")) || attr(tree, "order") == 
    "cladewise") 
    tree <- reorderPruning(tree)
  if (class(data)[1] != "phyDat") stop("data must be of class phyDat")
  if(is.null(tree$edge.length)) stop("tree must have edge weights") 
  if(any(is.na(match(tree$tip, attr(data, "names"))))) stop("tip labels are not in data")
  if(!is.null(opaa) && !is.null(opw)) warning("both optimal amino acis and weights are specified, weights are ignored")
  if(is.null(Q)) Q = rep(1,6)
  if(is.null(bfaa)) bfaa = findBf2(data) #if bfaa is not given, use the empirical bf
  if(is.null(dismat))
    dismat = GM_cpv(GM_CPV,al,beta,gamma)
  if(is.null(mumat)){
    mumat = aa_MuMat_form(Q)
    mumat = sym.to.Q(mumat,bfaa)
  }
  Qall = QAllaa(s,dismat,mumat,C,Phi,q,Ne)
  if(!is.null(opaa))
    ll = llop(data,tree,op=opaa,bf=bfaa,Qall=Qall,g=1,C=C,Phi=Phi,q=q,Ne=Ne)
  else 
    ll = llaaw(tree,data,Qall,opw,bfaa,C,Phi,q,Ne)
  result = list(ll=ll,data=data,tree=tree,s=s,GMweights=c(al,beta,gamma),Q=Q,dismat=dismat,mumat=mumat,bfaa=bfaa,call=call)
  class(result) = "mllm"
  return(result)
}
#MLE for s, beta and gamma, using Nelder-Mead method by default
#mllm <- function(data,tree,s=0,beta=be,gamma=ga,Q=NULL,dismat=NULL,mumat=NULL,
#opaa=NULL,opw=NULL,bfaa=NULL,C=2,Phi=0.5,q=4e-7,Ne=1.36e7) 
#sample call : optim.s.weight(gene1,ROKAS_TREE,0.1,be,ga,trace=1,Q=NU_VEC,bfaa=rep(1/20,20))
optim.s.weight <- function(data, tree, s,beta,gamma,method="Nelder-Mead",maxit = 500, trace=0, ...){
  ab <- log(c(s,beta,gamma))
  if(method != "nlminb"){
    fn = function(ab,data,tree, ...){
      ab = exp(ab)
      result = mllm(data=data,tree=tree,s=ab[1],beta=ab[2],gamma=ab[3], ...)$ll$loglik
      return(result)
    }
    res = optim(par=ab,fn=fn,gr=NULL,method=method,lower=-Inf,upper=Inf,
                control=list(fnscale=-1,trace=trace,maxit=maxit),data=data,tree=tree, ...)
    res$par = exp(res$par)
    return(res)
  }
  else{
    fn = function(ab,data,tree, ...){
      ab = exp(ab)
      result = mllm(data=data,tree=tree,s=ab[1],beta=ab[2],gamma=ab[3], ...)$ll$loglik
      return(-result)
    }
    res = nlminb(start=ab,objective=fn,gradient=NULL,hessian=NULL,
                 control = list(trace=trace),data=data,tree=tree, lower=-Inf, upper=Inf, ...)
    res$par = exp(res$par)
    res$objective = -res$objective
    return(res)
  }
}
## MLE for opw (weights for optimal amino acids)
## mllm(data,tree,s,beta,gamma,Q=NULL,opw=NULL,bfnu=NULL,bfaa=NULL,C=2,Phi=0.5,q=4e-7,Ne=1.36e7)
optim.opw <- function(data, tree,opw=rep(1/20,20),method="Nelder-Mead", maxit=500, trace=0, ...){
  l = length(opw)
  nenner = opw[l]
  lopw = log(opw*nenner)
  lopw = lopw[-l]
  fn = function(lopw,data,tree, ...){
    opw = exp(c(lopw,0))
    result = mllm(data=data,tree=tree,opw=opw, ...)$ll$loglik
    return(result)
  }
  res = optim(par=lopw,fn=fn,gr=NULL,method=method,lower=-Inf,upper=Inf,
              control=list(fnscale=-1,trace=trace,maxit=maxit),data=data,tree=tree, ...)
  print(res[[2]])
  opw = exp(c(res[[1]],0))
  opw = opw/sum(opw)
  res$par = opw
  return(res)
}
## mllm(data,tree,s,beta,gamma,Q=NULL,opw=NULL,bfnu=NULL,bfaa=NULL,C=2,Phi=0.5,q=4e-7,Ne=1.36e7)
optim.sw.opw <- function(data,tree,s,beta,gamma,opw=rep(1/20,20),method="Nelder-Mead",maxit=500,trace=0, ...){
  ab <- log(c(s,beta,gamma))
  l = length(opw)
  nenner = opw[l]
  lopw = log(opw*nenner)
  lopw = lopw[-l]
  fn = function(ablopw,data,tree, ...){
    ab = ablopw[1:3]
    lopw = ablopw[-(1:3)]
    ab = exp(ab)
    opw = exp(c(lopw,0))
    result = mllm(data=data,tree=tree,s=ab[1],beta=ab[2],gamma=ab[3], opw=opw, ...)$ll$loglik
    return(result)
  }
  res = optim(par=c(ab,lopw),fn=fn,gr=NULL,method=method,lower=-Inf,upper=Inf,
              control=list(fnscale=-1,trace=trace,maxit=maxit),data=data,tree=tree, ...)
  par = res$par
  ab = exp(par[1:3])
  opw = exp(c(par[-(1:3)],0))
  opw = opw/sum(opw)
  res$ab = ab
  res$opw = opw
  return(res)
}
## sample call: Q1 = optimQm(rokastree,gene1,Q=NU_VEC,subs=c(1,2,3,4,5,0),trace=1,s=0.1,beta=be,gamma=ga)
optimQm <- function(tree,data,Q=rep(1,6),subs=c(1:(length(Q)-1),0),method="Nelder-Mead",maxit=500,trace=0, ...){
  m = length(Q)
  n = max(subs)
  ab = numeric(n)
  for(i in 1:n) ab[i] = log(Q[which(subs==i)[1]])
  fn = function(ab,tree,data,m,n,subs, ...){
    Q = numeric(m)
    for(i in 1:n) Q[subs==i] = ab[i]
    result = mllm(data,tree,Q=exp(Q), ...)$ll$loglik
    return(result)
  }
  res = optim(par=ab,fn=fn,gr=NULL,method=method,lower=-Inf,upper=Inf,
              control=list(fnscale=-1,trace=trace,maxit=maxit),tree=tree,data=data,m=m,n=n,subs=subs, ...)
  Q = rep(1,m)
  for(i in 1:n) Q[subs==i] = exp(res$par[i])
  res$par = Q
  return(res)
}

optimBranch <- function(data,tree,method="Nelder-Mead",maxit=100,trace = 0, ...){
  ##s,beta,gamma,Q=rep(1,6),bfnu=rep(1/4,4),bfaa=NULL,C=2,Phi=0.5,q=4e-7,Ne=1.36e7){
  if(is.null(attr(tree,"order")) || attr(tree,"order") == "cladwise")
    tree <- reorderPruning(tree)
  el <- tree$edge.length
  tree$edge.length[el < 0] <- 1e-08
  el <- log(tree$edge.length)
  fn = function(el,data,tree, ...){
    tree$edge.length = exp(el)
    result = mllm(data,tree, ...)$ll$loglik
    return(result)
  }
  res = optim(par=el, fn=fn, gr=NULL,method=method,lower=-Inf,upper=Inf,
              control=list(fnscale=-1,trace=trace,maxit=maxit),tree=tree,data=data, ...)
  res$par = exp(res$par)
  return(res)
}
#mllm <- function(data,tree,s=0,beta=be,gamma=ga,Q=NULL,dismat=NULL,mumat=NULL,opaa=NULL,opw=NULL,bfaa=NULL)
#optim.mllm <- function(data,tree,Q=NULL,bfnu=NULL,bfaa=NULL,C=2,Phi=0.5,q=4e-7,Ne=1.36e7,){}
optim.mllm <- function(object, optQ = FALSE, optBranch = FALSE, optsWeight = TRUE, 
                       control = list(epsilon=1e-08,maxit=50,hmaxit=10,trace=0,htrace=TRUE),subs=NULL,...){
  tree = object$tree
  call = object$call
  if(class(tree)!="phylo") stop("tree must be of class phylo") 
  if(!is.rooted(tree)) stop("tree must be rooted")
  if (is.null(attr(tree, "order")) || attr(tree, "order") == "cladewise") 
    tree <- reorderPruning(tree)
  if(any(tree$edge.length < 1e-08)){
    tree$edge.length[tree$edge.length < 1e-08] <- 1e-08
    object <- update(object, tree=tree)
  }
  maxit = control$maxit
  trace = control$trace
  htrace = control$htrace #print out information about steps or not?
  data = object$data
  Q = object$Q
  if(is.null(subs)) subs = c(1:(length(Q)-1),0) #default is GTR
  bfaa = object$bfaa
  ll = object$ll$loglik
  ll1 = ll
  s = object$s
  beta = object$GMweights[2]
  gamma = object$GMweights[3]
  opti = TRUE
  rounds = 0
  while(opti){
    if(htrace)
      cat("iteration ",rounds+1,"\n")
    if(optsWeight){
      res = optim.s.weight(data,tree,s=s,beta=beta,gamma=gamma,maxit=maxit,trace=trace,
                           Q=Q,bfaa=bfaa,...)
      s = res$par[1]
      beta = res$par[2]
      gamma = res$par[3]
      if(htrace)
        cat("optimize s and Grantham weights: ", ll, "--->", res[[2]], "\n")
      ll = res[[2]]
    }
    if(optQ){
      res = optimQm(tree,data,Q=Q,subs=subs,maxit=maxit,trace=trace,s=s,beta=beta,gamma=gamma,bfaa=bfaa,...)
      Q = res$par
      if(htrace){
        cat("optimize rate matrix: ", ll, "--->", res[[2]], "\n")
      }
      ll = res[[2]]
    }
    if(optBranch){
      res = optimBranch(data,tree,maxit=maxit,trace=trace,s=s,beta=beta,gamma=gamma,Q=Q,bfaa=bfaa,...)
      if(htrace)
        cat("optimize branch lengths:", ll, "--->", res[[2]], "\n")
      tree$edge.length = res[[1]]
      ll = res[[2]]
    }
    rounds = rounds + 1
    if(rounds > control$hmaxit) opti <- FALSE
    if((ll1-ll)/ll < control$epsilon) opti <- FALSE
    ll1 = ll
  }
  object = update(object, tree=tree,data=data,s=s,beta=beta,gamma=gamma,Q=Q,bfaa=bfaa,...)
  return(object)
}
#rokasdata = read.phyDat("proteinevoutk20/pkg/Result/rokasAA",format="phylip",type="AA")
