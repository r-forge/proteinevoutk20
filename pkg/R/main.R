library(seqinr) #package to convert codons to amino acids
library(phangorn) #phylogenetics package
library(optimx) #optimization package
library(expm) #matrix exponentiation
#library(parallel) #parallel computing
library(multicore)
library(minqa) #optimization function bobyqa
library(mgcv) #find the unique rows in a matrix - uniquecombs
library(numDeriv) # calculate hessian of a function at a point
library(gtools) #package used to draw random samples from dirichlet distribution -- rdirichlet, ddirichlet
library(nloptr) #optimization
#library(scatterplot3d)
#library(akima) # gridded bivariate interpolation for irredular data
#library(Rmpfr)
#library(ppso)
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
# NU_VEC <- c(1.96575,4.08655,1.39431,1.46172,6.36024,1.00000)
# NU_VEC_JC <- rep(1.0,6)
# ##empirical base frequencies
# BF <- c(0.31065,0.18365,0.20955,0.29615)
##Best tree with branch lengths output by PAUP
#ROKAS_TREE <- read.nexus(paste(datadir,"GTR.tre",sep=""))
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
## read in fasta file of nucleotide (or nucleotide sequence), convert it into amino acid data
## type = "num" -- integers ; "AA" -- amino acid names; "phyDat" --- phyDat data type used in phangorn
conv <- function(filename,range=NULL,type="num",frame=0){
  levels <- AA #amino acids in alphabeticla order (not the single letter names, the 3-letter names)
  if(class(filename)=="character")
    data <- seqinr::read.fasta(file=filename) #read fasta file including nucleotide data
  else
    data <- filename
  if(length(range)!=0){
    for(i in 1:length(data))
      data[[i]] = data[[i]][range]
  }
  dna.to.aa <- function(x){ #given the order of the data, convert nucleotide data to amino acid data
    aas <- translate(seq=data[[x]],frame=frame)
    if(type=="num"){
      aas <- as.numeric(factor(aas,levels)) #convert from letters to numbers, as levels
      return(aas)
    }
    return(aas)
  }
  ## type = AA, return a list of amino acid sequences
  seqdata <- lapply(1:length(data),dna.to.aa)
  seqdata <- matrix(unlist(seqdata),nrow=length(data),byrow=TRUE)  
  if(type%in%c("AA","num"))
    row.names(seqdata) = attr(data,"name")
  else if(type=="phyDat"){ ## convert amino acid sequences to phyDat type
    seqdata = phyDat(seqdata,type="AA")
    names(seqdata) = attr(data,"name")
  }
  return(seqdata)
}
#convert fasta file to nexus file, for protein data
fasta.to.nex <- function(geneNum){
  for(i in geneNum){
    file = paste("~/proteinevoutk20/pkg/Data/gene",i,".fasta",sep="")
    gene = conv(file,"AA")
    write.nexus.data(gene,paste("~/proteinevoutk20/pkg/Data/gene",i,"AA.nex",sep=""),format="protein",interleaved=F)
  }
}
#convert fasta file to nexus file, for nucleotide data
fasta.to.nex.Nuc <- function(geneNum){
  for(i in geneNum){
    file = paste("~/proteinevoutk20/pkg/Data/gene",i,".fasta",sep="")
    gene = read.fasta(file=file)
    write.nexus.data(gene,paste("~/proteinevoutk20/pkg/Data/gene",i,".nex",sep=""),format="dna",interleaved=F)
  }
}
#############################################################################
## Remove columns with NA in the data
PruneMissing <- function(x,return=FALSE){
  ##Find the array indices of the NA entries
  naInd <- which(is.na(x),arr.ind=T)
  if(length(naInd)>0) {
    print("NA found!")
    dimnames(naInd) <- NULL
    ## column indices
    naCol <- naInd[,2]
    x <- x[,-naCol]
  }
  else
    print("no NA found")
  if(return)
    return(x)
}
#############################################################################
## find the empirical base frequencies for AMINO ACID list data (not phyDat data)
findBf <- function(datalist){
  data = unlist(datalist) #from list to vector
  datatb = table(data) #frequencies of each amino acid
  datatb = datatb[AA] #rearrange and order as AA
  datatb = as.numeric(datatb)
  return(datatb/sum(datatb)) #make sum 1
}
## for phyDat format
findBf2 <- function(data){
  ind = attr(data,"index")
  l = length(data)
  datalist = vector("list",length=l)
  for(i in 1:l)
    datalist[[i]] = data[[i]][ind]
  data = unlist(datalist)
  ## sometimes an amino acid does not appear in the simulated result, then the frequency equals NA
  ## this is no good
#   datatb = table(data)
#   datatb = as.numeric(datatb)[1:20]
  counts = sapply(1:20,function(i){sum(data==i)})
  return(counts/sum(counts))
}
## Find the most frequent element of a vector
## It works for numbers and characters
## "exclude" is not included in the mode (missing data)
Mode <- function(x,exclude=NULL) {
  ux <- unique(x)
  if(!is.null(exclude))
    ux <- setdiff(ux,exclude)
  ux[which.max(tabulate(match(x, ux)))]
}
## Find the most frequent amino acids for each distince pattern in a gene data
## used for optimal aa in majority rule
ModeAA <- function(phydata,exclude=23){
  if(class(phydata)!="phyDat") stop("data must be of class phyDat!")
  l = length(phydata)
  mat = NULL
  for(i in 1:l)
    mat = rbind(mat,phydata[[i]])
  modeaa = apply(mat,2,Mode,exclude=exclude)
  return(modeaa)
  }

#find the (default: 10) most frequent patterns from data of phyDat format 
# return the most frequent patterns with their counts
MostFreq <- function(data,count=10){
  weight = attr(data,"weight")
  datamat <- matrix(unlist(data),nrow=length(data),byrow=T)
  Index = head(order(weight,decreasing=T),n=count)
  return(list(mat=datamat[,Index],wt=weight[Index]))
}
#############################################################################
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
##given lower triangular part of R (i.e. Q = R %*% diag(bf)), and base frequencies, find the scaled Q
## Default: 4 by 4 matrix for Jukes-Cantor model
mat_form_lowtriQ <- function(Q=rep(1,6),bf=rep(0.25,4),byrow=FALSE){
  if(sum(bf)==0) stop("base frequencies can't all be 0!")
  l=length(bf)
  bf=bf/sum(bf)
  res = matrix(0, l, l)
  if(byrow) # if the lower triangular part of matrix is read by row instead of by column (default)
    res[upper.tri(res)] = Q
  else
    res[lower.tri(res)] = Q
  res = res+t(res) #symmetric matrix with diagonals 0
  res = res * rep(bf,each=l) #multiply cols by bf
  diag(res) = -rowSums(res) #set row sum to 0
  res = scaleQ(res,bf) #scale Q so that sum(pi_i * Q_ii) = -1
  return(res)
}


## from symmetric matrix of rates, find GTR rate matrix
sym.to.Q <- function(A,bf=NULL){
  if(!isSymmetric(A)) stop("matrix is not symmetric!")
  l = dim(A)[1] #dimension of rate matrix i.e. number of states
  if(is.null(bf)) bf = rep(1/l,l)
  if(sum(bf)==0) stop("base frequencies can't all be 0!")
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

## this is not a good way to find the frequencies of codons, nucleotide freqs
## usually are not directly related to codon freqs
freq_codon <- function(bf = rep(1/4,4)){
  bf=bf/sum(bf)
  freq = rep(0,61)
  for(i in 1:61){
    cd = as.numeric(factor(s2c(CDS[i]),Nu)) #nucleotide triplet of the codon considered
    freq[i] = prod(bf[cd]) #product of frequencies of nucleotides at all 3 positions
  }
  freq = as.table(freq)
  freq = freq/sum(freq)
  dimnames(freq)[[1]] = CDS # change the names of factors to the codon triplets
  return(freq)
}
##given the nucleotide base frequencies, find the bf for amino acids
##find the codon frequencies first,and sum up those that code for the same protein
freq_aa <- function(nubf=rep(1/4,4)){
  bf = freq_codon(nubf) #codon frequencies
  freq = matrix(0,20,1)
  dimnames(freq)[[1]] = AA
  #CDS_AA = sapply(1:61, function(x) translate(s2c(CDS[x]))) #amino acids coded by the list of codons
  for(i in 1:20)
    freq[AA[i],1] = sum(bf[CDS_AA==AA[i]])
  freq = as.table(as.vector(freq))
  dimnames(freq)[[1]] = AA
  return(freq)
}
## Given the rate matrix for 61 codons,find the rate matrix for 20 amino acids.
## Assume all the codons have the same frequencies for each amino acid
## Only mutation is taken into account, not selection and fixation.
## This is an exchange rate matrix, symmetric
AAMat <- function(CdMat){
  Q = matrix(0,20,20)
  for(i in 1:19){
    fcodons = cdlist[[i]] #codons that code for the amino acid mutate from
    nfcodons = length(fcodons)
    for(j in (i+1):20){
      tcodons = cdlist[[j]] #codons that code for the amino acid mutate to
      ntcodons = length(tcodons)
      Q[i,j] = sum(CdMat[fcodons,tcodons])/(nfcodons*ntcodons) #details in Yang's paper
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
  Q = matrix(0, l, l) # from vector to matrix
  Q[lower.tri(Q)] = vec
  Q = Q+t(Q) #symmetric matrix with diagonals 0
  dimnames(Q) = list(Nu,Nu)
  ##mutation matrix for 61 codons
  codon_array <- array(0,dim=c(61,61),dimnames=list(CDS,CDS))
  for(m in 1:60){                       #loop through all 61 codons
    fcd <- CDS[m] #codon mutated from
    fcd_ch <- s2c(fcd) #convert from string to character
    for(i in (m+1):61){ #only do the upper triangular part
      tcd = CDS[i]
      tcd_ch = s2c(tcd)
      cmp = cdcmp(fcd_ch,tcd_ch) #compare the 2 codons
      if(cmp$num==1){ # if only one position differs
        pos = cmp$pos
        #set the mutation rate to be the one between nucleotides
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
# therefore this is only an exchange rate matrix, freq of amino acids are not taken into account
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
# isSymmetric(aa_MuMat_form(runif(6)))
# checkSymmetric <- function(nuvec,bf){
#   bfq = freq_codon(bf)
#   bfaa = freq_aa(bf)
#   aaq = aa_MuMat_form(nuvec,bf)
#   return(isSymmetric(diag(bfaa) %*% aaq))
# }
## Amino acid mutation rate matrix for the GTR model of nucleotide
## MUMAT <- aa_MuMat_form(NU_VEC)

## change tree to be in pruning order, if it is not
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
    if((is.na(protein1[k])) & (!is.na(protein2[k]))){
      return(mean(DisMat[,protein2[k]]))
    }
    else if((is.na(protein2[k])) & (!is.na(protein1[k]))){
      return(mean(DisMat[protein1[k],]))
    }
    else if((is.na(protein2[k])) & (!is.na(protein1[k]))){
      return(mean(DisMat))
    }
    else
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
  result <- length(d)/(sum(1+d*s)) #harmonic mean of ftny at all sites (F_i = 1/(1 + d_i * s_i))
  return(result)
} #d and s are given

#Functionality of a protein given the optimal protein (aa sequence) and parameters
Ftny_protein <- function(protein,protein_op,s,DisMat){
  d <- pchem_d(protein,protein_op,DisMat)
  return(Ftny(d,s))
}
# harmonic.mean of a vector of numbers
harmonic.mean <- function(x){
  return(1/mean(1/x))
}
## find the equilibrium functionality, given opaa, ancestral prob matirx, s and dismat
eqm_ftny <- function(opaa,ancestral,s,DisMat,index=NULL){
  site.eqm <- function(op,prob,s,DisMat){
    ft <- sapply(1:20,FUN=Ftny_protein,protein_op=op,s=s,DisMat=DisMat)
    return(sum(ft*prob))
  }
  ancestral <- lapply(1:dim(ancestral)[2],function(x) ancestral[,x])
  fty.vec <- mapply(site.eqm,opaa,ancestral,MoreArgs=list(s=s,DisMat=DisMat))
  if(!is.null(index))
    fty.vec <- fty.vec[index]
  return(harmonic.mean(fty.vec))
}

#fixation probability, since we consider only site by site, so n=1
#d1, d2: distances (vectors)
#s: selection (vector or scalar)
#Ne: effective population size
#q, Phi: constants
#C: cost function, linear function of the length of the protein
#C_n = a1 + a2*n, cost --  linear function of the length
#mpfr indicates use of Rmpfr package, with large precision (10000 bit)
fix <- function(d1,d2,s,C=2,Phi=0.5,q=4e-7,Ne=5e6){
  if((d1==d2)||(s==0)) #When the fitnesses are the same, neutral case, pure drift
    return(1/(2*Ne))
  else{
      fit_ratio <- exp(-C*Phi*q*s*(d1-d2)) #f1/f2
     if(fit_ratio==Inf) #1 is much better than 2 (the mutant)
       return(0)
     else if(fit_ratio==1)
       return(1/(2*Ne))
     else
    return((1-fit_ratio)/(1-fit_ratio^(2*Ne)))
  }
}
#check the product of fixation probability and population size, 
# see if it's a constant
# fixNe <- function(d1,d2,s,q,Ne){
#   return(fix(d1,d2,s,C=2,Phi=0.5,q=q,Ne=Ne)*Ne)
#}
## Use the Rmpfr package for higher precision
# fix.mpfr <- function(d1,d2,s,C=2,Phi=0.5,q=4e-7,Ne=5e6){
#   if((d1==d2)||(s==0)) #When the fitnesses are the same, neutral case, pure drift
#     return(1/(2*Ne))
#   else{
#     fit_ratio <- exp(mpfr(-C*Phi*q*s*(d1-d2),prec=200)) #f1/f2
#     return((1-fit_ratio)/(1-fit_ratio^(2*Ne)))
#   }
# }
####################
#comments: in the fixation probability calculation, should the length of the protein sequence be taken
#into account? here the length is considered to be 1... which could be wrong...
#Since the production cost C is actually a linear function of length, the length n cancels
####################
#vector version of the fixation function, d1, d2 and s are vectors instead of numbers
fix2 <- function(d1,d2,s,C=2,Phi=0.5,q=4e-7,Ne=5e6){
  if(length(d1)!=length(d2)) #throw error if length of proteins are not the same
    stop("error: 2 proteins are of different lengths!")
  if(length(d1)==1){ #only one amino acid
    return(fix(d1,d2,s,C=C,Phi=Phi,q=q,Ne=Ne))
  }
  else{
    if((length(s)==1)&&(length(d1)!=1)) #if s is given as a scalar, then treat it to be the same across all sites
      s <- rep(s,length(d1))
    l = length(d1)
    cmp = cdcmp(d1,d2)
    if(cmp$num > 1) return(0) #more than 1 position differ
    else if((cmp$num ==0)) return(1/(2*Ne)) #same fitness/functionality
    else{ #exactly 1 position differs
      pos = cmp$pos 
      return(fix(d1[pos],d2[pos],s[pos],C=C,Phi=Phi,q=q,Ne=Ne))
    }
  }
}


#fixation probability of a mutant with initial freq 1/(2Ne) in a diploid population
#with proteins as parameteres instead of distance vectors for the proteins
fix_protein <- function(protein1, protein2, protein_op, s, DisMat,
                        C=2,Phi=0.5,q=4e-7,Ne=5e6){
  d1 <- pchem_d(protein1,protein_op,DisMat)
  d2 <- pchem_d(protein2,protein_op,DisMat)
  return(fix2(d1,d2,s,C,Phi,q,Ne))
}

######################################################
#Function to generate the transition matrix (m*m for now) 
#given the optimal amino acid. This is only for one site (amino acid)
#parameter: aa_op: optimal amino acid at this site
#s: selection coefficient(SCALAR, since it's only one site)
#m: number of states (amino acids) considered
#DisMat: distance matrix between amino acids
#MuMat: mutation rate matrix between amino acids
mat_gen_indep <- function(aa_op,s,DisMat,MuMat,C=2, Phi=0.5,q=4e-7,Ne=5e6){
  m = 20
  mat <- matrix(0,nrow=m,ncol=m)#set diagonal entries to be 0 at first
  for(i in 1:(m-1)){
    for(j in (i+1):m){
      mat[i,j] <- fix_protein(i,j,aa_op,s,DisMat,C,Phi,q,Ne) #fixation prob -> transition rate
      mat[j,i] <- fix_protein(j,i,aa_op,s,DisMat,C,Phi,q,Ne) #symmetric entry
    }#end for j
  }#end for i
  mat = mat*2*Ne*MuMat
  diag(mat) <- -rowSums(mat)
  return(mat)
}
#matrix of fixation probabilities, given optimal aa, sensitivity coefficient, distance matrix and other parameters
fixmat <- function(aa_op,s,DisMat,C=2, Phi=0.5,q=4e-7,Ne=5e6){
  m = 20
  mat <- matrix(0,nrow=m,ncol=m)#set diagonal entries to be 0 at first
  for(i in 1:(m-1)){
    for(j in (i+1):m){
      mat[i,j] <- fix_protein(i,j,aa_op,s,DisMat,C,Phi,q,Ne) #fixation prob -> transition rate
      mat[j,i] <- fix_protein(j,i,aa_op,s,DisMat,C,Phi,q,Ne) #symmetric entry (not the same rate!)
    }#end for j
  }#end for i
  return(mat)
}
#For all amino acids as optimal
fixmatAll <- function(s,DisMat,C=2, Phi=0.5,q=4e-7,Ne=5e6){
  res = lapply(1:20,function(i) {fixmat(i,s,DisMat,C, Phi, q, Ne)})
  return(as.matrix(res))
}
#given the fixation prob matrix and mutation matrix, 
Qmat <- function(FixMat,MuMat,Ne=5e6){
  mat = FixMat*MuMat*2*Ne
  diag(mat) = -rowSums(mat)
  return(mat)
}
# list of 20 rate matrices, one of each aa as optimal
# these matrices are not scaled, did not use the base frequencies
QAllaa <- function(s,DisMat,MuMat,C=2, Phi=0.5,q=4e-7,Ne=5e6){
  #fixall <- fixmatAll(s=s,DisMat=DisMat,C=C,Phi=Phi,q=q,Ne=Ne)
  #res <- lapply(1:20, function(x) Qmat(fixall[[x,1]],MuMat,Ne=Ne))
  res = lapply(1:20,function(i) {mat_gen_indep(i,s,DisMat, MuMat,C, Phi, q, Ne)})
  return(as.matrix(res))
}
QAllaa1 <- function(fixall,MuMat,Ne=5e6){
  res <- lapply(1:20, function(x) Qmat(fixall[[x,1]],MuMat,Ne=Ne))
  return(as.matrix(res))
}
# given the edge lengths el, rate matrix Q and rate g, find the probability transition matrices P
# result is of dimension length(g)*length(el). Therefore the transition probabilities on all branches
# To parse an object from getPm, say P, do P[i,j][[1]], or P[[i,j]]
# Q is a scaled matrix
getPm1 <- function(el, Q, g){
  eig = eigen(Q,FALSE)
  eig$inv = solve.default(eig$vec)
  res <- .Call("getPM",eig,as.integer(length(eig$values)),as.double(el),as.double(g),PACKAG="phangorn")
  attr(res,"dim") <- c(length(g),length(el))
  res
}
getPm <- function(el, Q, g){
  ell = length(el)
  gl = length(g)
  res = array(list(NULL),dim=c(gl,ell))
  for(i in 1:gl){
    for(j in 1:ell){
#       cat("branch length is ",el[j],"\n")
#       cat(g[i])
#       print(Q)
#       cat("Q*el*g:","\n")
#       print(Q*el[j]*g[i])
      res[[i,j]] = expm(Q*el[j]*g[i])
      
    }
  }
  return(res)
}
## find the equilibrium frequencies from Q(instantaneous transition rate matrix), row sums of Q are 0
## This does not use matrix exponentiation, runs faster
eqmQ <- function(Q){
  n <- dim(Q)[1]
  Q[,1] <- rep(1,n) #sum of frequencies is equal to 1
  vec <- solve(t(Q),c(1,rep(0,n-1)))
  return(vec)
}
##################################################################
##   downpass/pruning method to find probabilities of finding each
##   state at root(ancestral) node. only for one site.
##################################################################
ll_site <- function(tree,data,optimal,s,Q,alpha=al, beta=be, gamma=ga,
                    bf=NULL,C=2,Phi=0.5,q=4e-7,Ne=5e6){
  ##If the given tree is not rooted and binary, then throw error and exit
  if(!is.binary.tree(tree)|!is.rooted(tree)) stop("error: the input phylogeny is not rooted binary tree!")
  m = 20
  ##if the base frequencies are not specified, do a uniform distribution
  if(is.null(bf)) bf=rep(1/m,m)#base frequency, randomly chosen from all states
  GM = GM_cpv(GM_CPV,alpha,beta,gamma)
  mumat = aa_MuMat_form(vec=Q)
  Qmat = mat_gen_indep(optimal,s,DisMat=GM,MuMat=mumat,C,Phi,q,Ne) #transition rate matrix for the site, given the optimal aa
  Qmat = scaleQ(Qmat,bf)
  
  tree <- ape:::reorder.phylo(tree,"pruningwise") #reorder the tree in pruningwise order
  edge = tree$edge #edges
  nNodes = max(edge) #number of nodes in the tree (including tips)
  probvec = matrix(NA,nNodes,m) #probability of gettting different states at nodes that evolve to the current sequences
  parent <- as.integer(edge[, 1]) #parents of the edges
  child <- as.integer(edge[, 2]) #children of the edges
  root <- as.integer(parent[!match(parent, child, 0)][1])  
  tip <- as.integer(child[!match(child,parent,0)])
  
  init.tip <- function(x){ #initiate the vector for the tips, 1 for the tip state, 0 otherwise
    vec <- rep(0,m)
    vec[data[x]] <- 1
    vec
  }
  probvec[1:length(tip),] <- t(sapply(1:length(tip),init.tip)) #all tips
  tl = tree$edge.length #lengths of the edges
  P <- getPm(tl,Qmat,1)
  for(i in 1:tree$Nnode){ #for each interior node calculate the probability vector of observing 1 of 20 states
    from = parent[2*i] #parents
    to = child[(2*i-1):(2*i)] #direct descendents
    v.left <- P[[1,2*i-1]]
    v.right <- P[[1,2*i]]
    probvec[from,] <- as.vector((v.left%*%probvec[to[1],])*(v.right%*%probvec[to[2],])) #pruning, vector form
    check.sum <- sum(probvec[from,])
    if(check.sum==0) #probability is very very low
      warning("numerical overflow",immediate.=TRUE)
  }
  return(as.numeric(probvec[root,]%*%bfaa))
  #return(list(ll=max(probvec[root,]),root=which.max(probvec[root,]))) #with the corresponding root returned 
  #return(max(probvec[root,])) #just the value
}
####################################################################################################################################
# Q does not get scaled in the function, it reflexes the difference in evolutionary rates for different optimal amino acids
# returns loglikelihood values for all distinct patterns in the data, matrix with 1 column (or a column vector)
# This calls the internal C function in phangorn package
ll3m <- function (data, tree, scale = 1,ancestral = NULL, ancStates = NULL,Q, g = 1) 
{
  if (is.null(attr(tree, "order")) || attr(tree, "order") == "cladewise")
    tree <- reorderPruning(tree)
  node <- tree$edge[, 1]
  edge <- tree$edge[, 2]
  node = as.integer(node - min(node))
  edge = as.integer(edge - 1)
  m = length(edge) + 1
  nr <- as.integer(attr(data, "nr"))
  nc <- as.integer(attr(data, "nc"))
  nTips = as.integer(length(tree$tip))
  mNodes = as.integer(max(node) + 1)
  contrast = attr(data, "contrast")
  nco = as.integer(dim(contrast)[1])
  index = as.integer(attr(data,"index"))
  
  Q = Q*scale #scale Q by a number,default 1 (no scaling)
  stnry.prob <- eqmQ(Q) #stationary probabilities, as prior on the root states
  Q = t(Q) # in order to use C function in phangorn package
  el <- tree$edge.length
  P <- getPm(el,Q,g)
  ## get likelihood for all distinct site patterns, with every aa as root state
  res <- .Call("LogLik2", data[tree$tip.label], P, nr, nc, node, edge, nTips, 
               mNodes, contrast, nco, PACKAGE = "phangorn")
  if(!is.null(ancStates)){ 
    if(length(ancStates)==nr)#if root states are given (per distinct site pattern instead of all sites)
      result = sapply(1:nr,function(x) res[[1]][x,ancStates[x]])
    else if(length(ancStates)==length(index)) #root states are given for all sites
      result = sapply(1:length(index, function(x) res[[1]][index[x],ancStates[x]]))
    else
      stop("ancestral state sequence has wrong length!")
    result = matrix(log(result),ncol=1) #take logarithm 
    ancStates = ancStates
  }
  else{ # if root states are not given
    if(is.null(ancestral)){ # method of getting root states is not specified, max as default
      ## return the root states at all sites, TAKE INTO ACCOUNT OF PRIOR DISTRIBUTION
      ## this only affects the default option, which is "optimize root states"
      ancStates =  sapply(1:nr,function(x) which.max(res[[1]][x,]*stnry.prob)) 
      result = sapply(1:nr, function(x) res[[1]][x,ancStates[x]])
      result = matrix(log(result),ncol=1)
    }
    else{ # if ancestral is given as a matrix, with 20 rows, one column for each distince site pattern
      ancestral = as.matrix(ancestral,nrow=20)
      if(ncol(ancestral)==1) #if it's the same vector for all sites
        ancestral = sapply(1:nr,function(x) ancestral)
      ancestral <- apply(ancestral,2,function(x) x/sum(x))
      result = log(diag(res[[1]]%*%ancestral)) #updated phangorn, phangorn 1.7
    }
  }
  if(is.null(ancestral)){ # convert root states to matrix with 0-1 entries
    diag.mat = diag(20)
    ancestral = diag.mat[,ancStates]
  }
  ## result: nsite*1 matrix, ancestral: 20*nsite matrix, ancStates: inferred root states if applicable
  return(list(result=result,ancestral=ancestral,ancStates=ancStates))
}


# change the orders of the arguments in function ll3m so that mapply can be used later
ll <- function(Q,scale,ancestral,ancStates,data,tree,g){
  ll3m(data,tree,scale,ancestral,ancStates,Q,g)
}

## return a matrix with 20 columns (one for each aa as optimal aa)
## one row for each site patter in data
## with the choice of ancestral states given
## The matrix is used later, when the optimal aa's are known or the rule of finding them is given
## Now it's more of a list of 20 instead of a matrix.
## every one of the 20 list has 3 components: result(sitelik),ancestral,ancStates
llmat <- function(data,tree,Qall,scale.vec=rep(1,20),ancestral=NULL,ancStates=NULL,g=1){
  result = NULL
  nr = attr(data,"nr") #number of different sites
  if(!is.null(ancStates)){ ##root states are given
    result = mapply(ll,Qall,scale.vec,MoreArgs=list(data=data,tree=tree,ancStates=ancStates,ancestral=NULL,g=g),SIMPLIFY=FALSE)
   }
  else{ ## root states are not given, ancStates==NULL
    if(is.null(ancestral)) #default is MaxRoot
      ancestral = "max"
    if(is.vector(ancestral,mode="numeric")||is.matrix(ancestral)){ #if root freq is given as vector
      result = mapply(ll,Qall,scale.vec,MoreArgs=list(data=data,tree=tree,ancestral=ancestral,ancStates=ancStates,g=g),SIMPLIFY=FALSE)
    }
    else if(ancestral=="max"){ #MaxRoot
      result = mapply(ll,Qall,scale.vec,MoreArgs=list(data=data,tree=tree,ancStates=ancStates,ancestral=NULL,g=g),SIMPLIFY=FALSE)
    }
    else if(ancestral=="opaa"){ #OpRoot -- root is the same as optimal aa
      anc=diag(20)
      anc = lapply(1:20,function(x) anc[x,])
      result = mapply(ll,Qall,scale.vec,anc,MoreArgs=list(data=data,tree=tree,ancStates=ancStates,g=g),SIMPLIFY=FALSE)
    }
    else if(ancestral=="eqm"){ # EqmRoot
      anc=lapply(X=Qall,eqmQ) #find equilibrium frequencies for all Q's
      result = mapply(ll,Qall,scale.vec,anc,MoreArgs=list(data=data,tree=tree,ancStates=ancStates,g=g),SIMPLIFY=FALSE)
    }
    else if(ancestral=="emp"){ #EmpRoot - root frequencies from observation
      result = mapply(ll,Qall,scale.vec,MoreArgs=list(data=data,tree=tree,ancestral=findBf2(data),ancStates=ancStates,g=g),SIMPLIFY=FALSE)
    }
    else stop("ancestral options: opaa, eqm,emp,max, or a numeric vector of length 20")
  }
  return(result)
}
## optimal amino acids are given
llop <- function(op,ll_mat,data){
  weight = attr(data,"weight")
  nr = attr(data,"nr") #number of different site patterns
  index = attr(data,"index")
  if(length(op)!=nr) #opaa is given for each distinct data pattern
    stop("optimal aa sequence has wrong length!")
  
  llmat = sapply(1:20,function(x) ll_mat[[x]]$result) ## nr*20 matrix
  sitelik = sapply(1:nr,function(i) llmat[i,op[i]]) ## loglikelihood for each site
  loglik <- sum(weight*sitelik) ## total log likelihood
  get.ancestral <- function(i){ ## get the matrix for root state probabilities
    opi <- op[i]
    ll_mat[[opi]]$ancestral[,i]
  }
  ancestral <- sapply(1:nr,get.ancestral)
  return(list(loglik=loglik,opaa=op,sitelik = sitelik,ancestral=ancestral))
}
## optimal amino acids are found by max rule
llmax <- function(ll_mat,data){
  llmat = sapply(1:20,function(x) ll_mat[[x]]$result)
  opaa = apply(llmat,1,which.max) ## find the aas that maximize the likelihood
  result = llop(op=opaa,ll_mat=ll_mat,data=data) ## use the above function
  return(result)
}
# assume every amino acid has a weight to be the optimal one, and the weights are the same for all sites, calculate the loglikelihood
# if no weights are specified, use the aa's that maximize the likelihoods (which is the same as above)
### root states to return need to be corrected (TO DO)
llaaw <- function(opw,weight,ll_mat){
  llmat = sapply(1:20,function(x) ll_mat[[x]]$result)
  opw = opw/sum(opw)
  eresult = exp(llmat) #likelihood values, from exp(loglikelihood)
  sitelik = eresult %*% opw
  sitelik = log(sitelik)
  loglik = sum(weight * sitelik)
  return(list(loglik=loglik,opaa="weighted",sitelik = sitelik,ancestral=NULL))
}
## given opw (weights of aas being optimal), weights of data patterns and LIKELIHOOD (not loglikelihood) matrix: n * 20
## find the -loglikelihood for all data, return only that value, instead of a list
llaaw1 <- function(opw,weight,llmat){
  sitelik = llmat %*% opw
  sitelik = log(sitelik)
  result = -sum(sitelik*weight)
  #cat("par:",opw,"val:",result,"\n")
  return(result)
}
## gradient function of the previous function, of variables opw
llaaw_grad <- function(opw,weight,llmat){
  grad_vec <- rep(0,length(opw))
  for(i in 1:length(opw)){
    grad_vec[i] <- -weight%*%(llmat[,i]/(llmat%*%opw))
  }
  return(grad_vec)
}

## find the loglikelihood given Q and other paramters, here Q is the lower triangular part of the 
## nucleotide transition rate matrix of length 6
## this one uses expm for matrix exponentiation
#more flexibility on what is provided, good for avoiding unnecessary computations
mllm1 <- function(data,tree,s=NULL,beta=be,gamma=ga,scale.vec=rep(1,20),Q=NULL,dismat=NULL,fixmatall=NULL,mumat=NULL,Qall=NULL,
                  opaa=NULL,opw=NULL,bfaa=NULL,ancestral=NULL,ancStates=NULL,C=2,Phi=0.5,q=4e-7,Ne=5e6){
  call <- match.call()
  if(class(tree)!="phylo") stop("tree must be of class phylo") 
  if (is.null(attr(tree, "order")) || attr(tree, "order") == 
    "cladewise") 
    tree <- reorderPruning(tree)
  if (class(data)[1] != "phyDat") stop("data must be of class phyDat")
  if(is.null(tree$edge.length)) stop("tree must have branch length") 
  if(any(is.na(match(tree$tip, attr(data, "names"))))) stop("tip labels are not in data")
  if(!is.null(opaa) && !is.null(opw)) warning("both optimal amino acis and weights are specified, weights are ignored")
  
  
  ##if(is.null(bfaa)) bfaa = findBf2(data) #if bfaa is not given, use the empirical bf
  #goal: to get Qall: Q matrices for all aas as optimal
  if(is.null(Qall)){
    if(is.null(fixmatall)){
      if(is.null(dismat))
        dismat = GM_cpv(GM_CPV,al,beta,gamma)
      fixmatall <- fixmatAll(s,DisMat=dismat,C=C,Phi=Phi,q=q,Ne=Ne)
    }
    if(is.null(mumat)){
      if(is.null(Q)) {Q = rep(1,6)}
      mumat = aa_MuMat_form(Q)
    }
    Qall = QAllaa1(fixmatall,mumat,Ne=Ne)
  }
  ll_mat = llmat(data=data,tree=tree,Qall=Qall,scale.vec=scale.vec,ancestral=ancestral,ancStates=ancStates)
  if(!is.null(opaa)) #opaa specified
    ll = llop(op=opaa,ll_mat=ll_mat,data=data)
  else if(!is.null(opw)) ## opw specified or max rule for optimal aa
    ll = llaaw(opw=opw,weight=weight,llmat=ll_mat)
  else{
    ll = llmax(ll_mat,data)
  }
  result = list(ll=ll,data=data,tree=tree,s=s,beta=beta,gamma=gamma,Q=Q,scale.vec=scale.vec,
                dismat=dismat,fixmatall=fixmatall,Qall=Qall,
                mumat=mumat,opaa=opaa,opw=opw,ancestral=ancestral,call=call)
  class(result) = "mllm"
  return(result)
}

######################################################################################################
optim.all <- function(data,tree,s,beta,gamma,Q,print_level=0,method="SBPLX",maxeval="100",...){
  if(is.null(attr(tree,"order")) || attr(tree,"order") == "cladwise")
    tree <- reorderPruning(tree)
  brlen=tree$edge.length
  brCt = length(brlen) # number(count) of branches
  Q1 = Q/Q[6] #make the last entry 1
  Q = Q1[1:5] #first 5 rates
  ab = c(s,beta,gamma,brlen,Q) #all parameters as a vector, length = 3+br.num+6
  paralen <- length(ab)
  fn <- function(ab){
    s <- ab[1]
    beta <- ab[2]
    gamma <- ab[3]
    cat("s,beta,gamma:",s,beta,gamma,"\n")
    tree$edge.length <- ab[4:(4+brCt-1)]
    cat("branch lengths:", tree$edge.length,"\n")
    Q <- c(ab[(4+brCt):paralen],1)
    cat("Q:",Q,"\n")
    loglik <- mllm1(data=data,tree=tree,s=s,beta=beta,gamma=gamma,Q=Q,...)$ll$loglik
    cat("-loglik",-loglik,"\n")
    return(-loglik)
  }
  lower <- rep(10e-8,paralen)
  upper <- c(10,rep(Inf,paralen-1))
  opts <- list("algorithm"=paste("NLOPT_LN_",method,sep=""),"maxeval"=maxeval,"xtol_rel"=1e-6,
               "ftol_rel"=.Machine$double.eps^0.5,"print_level"=print_level)
  res = nloptr(x0=ab,eval_f=fn, lb=lower,ub=upper,opts=opts)
  par_optim <- res$solution
  res$s = par_optim[1]
  res$beta = par_optim[2]
  res$gamma = par_optim[3]
  res$Q = c(par_optim[(4+brCt):paralen],1)
  tree$edge.length = par_optim[4:(4+brCt-1)]
  res$tree = tree
  return(res)
}


## find the path to a particular tip in a tree
## including the branches and the nodes on the path
path_to_tip <- function(tree,i){
  br <- tree$edge #all the edges
  ntips <- length(tree$tip.label)
  parents <- as.integer(br[, 1])
  child <- as.integer(br[, 2])
  root <- as.integer(parents[!match(parents, child, 0)][1]) 
  
  done <- FALSE
  node.path <- i #nodes on the path, from root to tip
  br.path <- NULL # indices of branches on the path to the tip
  while(!done){
    if(node.path[1]==root)
      done <- TRUE
    else{
      start <- node.path[1] #start tracing, go backwards
      br.ind <- which(br[,2]==start) #find the index of branch that ends in "start"
      br.path <- c(br.ind,br.path)
      pre <- br[br.ind,1] #next node on the path when going backwards
      node.path <- c(pre,node.path)
    }
  }
  return(list(br.path=br.path,node.path=node.path))
}

## given a tree and desired depth, make the tree ultrametric without changing the interior branch lengths
ult.tree <- function(tree,depth=NULL,pathlist=NULL){
  if(is.ultrametric(tree))
    return(tree)
  else{
    nTips <- length(tree$tip.label)
    if(is.null(depth))
      depth <- max(node.depth.edgelength(tree))
    if(is.null(pathlist))
      pathlist = lapply(1:nTips,function(x) path_to_tip(tree,x)$br.path) #find list of paths to tips if not given
    for(i in 1:nTips){
      brs <- pathlist[[i]]
      l <- length(brs)
      if(l==1)
        tree$edge.length[brs] = depth
      else
        tree$edge.length[brs[l]] = depth - sum(tree$edge.length[brs[1:(l-1)]])
    }
    return(tree)
  }
}
optim.all.ultrametric <- function(data,tree,s,beta,gamma,Q,print_level=0,method="SBPLX",maxeval="100",...){
  if(is.null(attr(tree,"order")) || attr(tree,"order") == "cladwise")
    tree <- reorderPruning(tree)
  
  nTips <- length(tree$tip.label)
  pathlist = lapply(1:nTips,function(x) path_to_tip(tree,x)$br.path)
  nIntBr <- nTips - 2 #number of interior branches
  edge <- tree$edge
  int.ind <- which(!(edge[,2] %in% c(1:nTips))) #index of interior branches in tree$edge
  ext.ind <- which(edge[,2] %in% c(1:nTips))
  brCt <- length(int.ind) #interior branches
  brlen=tree$edge.length[int.ind]

  Q1 = Q/Q[6] #make the last entry 1
  Q = Q1[1:5] #first 5 rates
  ab = c(s,beta,gamma,Q,brlen,max(node.depth.edgelength(tree))) #all parameters as a vector, length = 3+br.num+6
  paralen <- length(ab)
  fn <- function(ab){
    s <- ab[1]
    beta <- ab[2]
    gamma <- ab[3]
    cat("s,beta,gamma:",s,beta,gamma,"\n")
    Q <- c(ab[4:8],1)
    cat("Q:",Q,"\n")
    
    depth <- tail(ab,1) # tree depth is the last element in the vector
    tree$edge.length[int.ind] <- ab[9:(8+brCt)] # assign internal branches new lengths
    tree <- ult.tree(tree,depth=depth,pathlist=pathlist) #calculate exterior branch lengths
    cat("tree edge lengths","\n")
    cat(tree$edge.length,"\n")
    if(any(tree$edge.length<0)) return(10^7)
    else{
      loglik <- mllm1(data=data,tree=tree,s=s,beta=beta,gamma=gamma,Q=Q,...)$ll$loglik
      cat("-loglik",-loglik,"\n")
      return(-loglik)
    }
  }
  lower <- rep(10e-8,paralen)
  upper <- c(10,rep(Inf,paralen-1))
  opts <- list("algorithm"=paste("NLOPT_LN_",method,sep=""),"maxeval"=maxeval,"xtol_rel"=1e-6,
               "ftol_rel"=.Machine$double.eps^0.5,"print_level"=print_level)
  res = nloptr(x0=ab,eval_f=fn, lb=lower,ub=upper,opts=opts)
  par_optim <- res$solution
  res$s = par_optim[1]
  res$beta = par_optim[2]
  res$gamma = par_optim[3]
  res$Q = c(par_optim[4:8],1)
  tree$edge.length[int.ind] = par_optim[9:(8+brCt)] #new interior branch lengths
  depth <- tail(par_optim,1) # tree depth is the last element in the vector
  res$tree <- ult.tree(tree,depth=depth,pathlist=pathlist) #calculate exterior branch lengths
  return(res)
}


