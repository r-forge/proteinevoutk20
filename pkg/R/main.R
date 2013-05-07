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
test = TRUE
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
# Here Q gets scaled in the function, so the given matrix doesn't have to be scaled
# returns loglikelihood values for all distinct patterns in the data, matrix with 1 column (or a column vector)
# This calls the internal C function in phangorn package
ll3m <- function (dat1, tree, bf = rep(1/20,20),ancestral = NULL, ancStates = NULL,Q, g = 1) 
{
  if (is.null(attr(tree, "order")) || attr(tree, "order") == "cladewise")
    tree <- reorderPruning(tree)
  q = length(tree$tip.label)
  node <- tree$edge[, 1]
  edge <- tree$edge[, 2]
  m = length(edge) + 1
  bf = bf/sum(bf)
  Q = scaleQ(Q,bf) #scale the substitution rate matrix
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
  res <- .Call("LogLik2", dat1[tree$tip.label], P, nr, nc, node, edge, nTips, mNodes, contrast, nco, PACKAGE = "phangorn")
#   if((!is.null(ancestral))&&(!is.null(ancStates)))
#     warning("both ancestral state frequencies and states are specified, only ancestral states are used!")
## if root state frequencies and root states are both given, root states have priority and are used
  #browser()
  if(!is.null(ancStates)){
    if(length(ancStates)!=nr) #length has to be the same as the dinstince site patterns
      stop("ancestral state sequence has wrong length!")
    result = sapply(1:nr,function(x) res[[1]][x,ancStates[x]])
    result = matrix(log(result),ncol=1) #take logarithm 
    ancStates = ancStates
  }
  else{
    if(is.null(ancestral)){
      result = sapply(1:nr,function(x) max(res[[1]][x,])) #find the root state that maximizes the likelihood
      result = matrix(log(result),ncol=1) #take logarithm 
      ancStates =  sapply(1:nr,function(x) which.max(res[[1]][x,])) # return the root states at all sites
    }
    else{
      ancestral = as.matrix(ancestral,nrow=20)
      if(ncol(ancestral)==1)
        ancestral = sapply(1:nr,function(x) ancestral)
      ancestral <- apply(ancestral,2,function(x) x/sum(x))
      result = log(diag(res[[1]]%*%ancestral)) #updated phangorn, phangorn 1.7
    }
  }
  if(is.null(ancestral)){
    diag.mat = diag(20)
    ancestral = diag.mat[,ancStates]
  }
  return(list(result=result,ancestral=ancestral,ancStates=ancStates))
  #return(result)
}
# change the orders of the arguments in function ll3m so that mapply can be used later
ll <- function(Q,ancestral,ancStates,data,tree,bf,g){
  ll3m(data,tree,bf,ancestral,ancStates,Q,g)
}

## return a matrix with 20 columns (one for each aa as optimal aa)
## one row for each site patter in data
## with the choice of ancestral states given
## The matrix is used later, when the optimal aa's are known or the rule of finding them is given
## Now it's more of a list of 20 instead of a matrix.
## every one of the 20 list has 3 components: result(sitelik),ancestral,ancStates
llmat <- function(data,tree,Qall,bf=rep(1/20,20),ancestral=NULL,ancStates=NULL,C=2,Phi=0.5,q=4e-7,Ne=5e6){
  result = NULL
  #weight = attr(data,"weight")
  nr = attr(data,"nr") #number of different sites
  #index = attr(data,"index")
  #ns = length(index) #number of sites
#   if((!is.null(ancestral))&&(!is.null(ancStates)))
#     warning("both ancestral state frequencies and states are specified, only ancestral states are used!")
  if(!is.null(ancStates)){ ##root states are given
    result = lapply(Qall,ll3m,dat1=data,tree=tree,bf=bf,ancStates=ancStates,ancestral=NULL,g=1)
   }
  else{ ## root states are not given, ancStates==NULL
    if(is.null(ancestral)) #default is MaxRoot
      ancestral = "max"
    if(is.vector(ancestral,mode="numeric")||is.matrix(ancestral)){ #if root freq is given as vector
      result = lapply(Qall,ll3m,dat1=data,tree=tree,bf=bf,ancestral=ancestral,ancStates=ancStates,g=1)
    }
    else if(ancestral=="max"){ #MaxRoot
      result = lapply(Qall,ll3m,dat1=data,tree=tree,bf=bf,ancestral=NULL,ancStates=ancStates,g=1)
    }
    else if(ancestral=="opaa"){ #OpRoot -- root is the same as optimal aa
      anc=diag(20)
      anc = lapply(1:20,function(x) anc[x,])
      result = mapply(ll,Qall,anc,MoreArgs=list(data=data,tree=tree,bf=bf,ancStates=ancStates,g=1),SIMPLIFY=FALSE)
    }
    else if(ancestral=="eqm"){ # EqmRoot
      anc=lapply(X=Qall,eqmQ) #find equilibrium frequencies for all Q's
      result = mapply(ll,Qall,anc,MoreArgs=list(data=data,tree=tree,bf=bf,ancStates=ancStates,g=1),SIMPLIFY=FALSE)
    }
    else if(ancestral=="emp"){ #EmpRoot - root frequencies from observation
      result = lapply(Qall,ll3m,dat1=data,tree=tree,bf=bf,ancestral=findBf2(data),ancStates=ancStates,g=1)
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
  loglik <- sum(weight*sitelik)
  get.ancestral <- function(i){
    opi <- op[i]
    ll_mat[[opi]]$ancestral[,i]
  }
  ancestral <- sapply(1:nr,get.ancestral)
  return(list(loglik=loglik,opaa=op,sitelik = sitelik,ancestral=ancestral))
}
## optimal amino acids are found by max rule
llmax <- function(ll_mat,data){
  llmat = sapply(1:20,function(x) ll_mat[[x]]$result)
  opaa = apply(llmat,1,which.max)
  result = llop(op=opaa,ll_mat=ll_mat,data=data)
  return(result)
}
# assume every amino acid has a weight to be the optimal one, and the weights are the same for all sites, calculate the loglikelihood
# if no weights are specified, use the aa's that maximize the likelihoods (which is the same as above)
### root states to return need to be corrected (TO DO)
llaaw <- function(opw,weight,llmat){
  optimal_aa = "weighted"
  opw = opw/sum(opw)
  eresult = exp(llmat) #likelihood values, from exp(loglikelihood)
  sitelik = eresult %*% opw
  sitelik = log(sitelik)
  loglik = sum(weight * sitelik)
  return(list(loglik=loglik,optimal_aa=optimal_aa,weights=opw, sitelik = sitelik, llmat=llmat))
}
## given opw (weights of aas being optimal), weights of data patterns and LIKELIHOOD (not loglikelihood) matrix: n * 20
## find the -loglikelihood for all data
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
mllm1 <- function(data,tree,s=NULL,beta=be,gamma=ga,Q=NULL,dismat=NULL,fixmatall=NULL,mumat=NULL,Qall=NULL,
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
  
  
  if(is.null(bfaa)) bfaa = findBf2(data) #if bfaa is not given, use the empirical bf
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
  ll_mat = llmat(data=data,tree=tree,Qall=Qall,bf=bfaa,ancestral=ancestral,ancStates=ancStates,C=C,Phi=Phi,q=q,Ne=Ne)
#   llmat = ll_mat$mat
#   root = llmat$root
  if(!is.null(opaa)) #opaa specified
    ll = llop(op=opaa,ll_mat=ll_mat,data=data)
  else if(!is.null(opw)) ## opw specified or max rule for optimal aa
    ll = llaaw(opw=opw,weight=weight,llmat=llmat)
  else{
    ll = llmax(ll_mat,data)
  }
  result = list(ll=ll,data=data,tree=tree,s=s,GMweights=c(al,beta,gamma),Q=Q,dismat=dismat,fixmatall=fixmatall,Qall=Qall,
                mumat=mumat,opaa=opaa,opw=opw,bfaa=bfaa,ancestral=ancestral,call=call)
  class(result) = "mllm"
  return(result)
}

######################################################################################################
#MLE for s, given beta and gamma and all other parameter values
#  mllm1 <- function(data,tree,s=NULL,beta=be,gamma=ga,Q=NULL,dismat=NULL,fixmatall=NULL,mumat=NULL,
#                    Qall=NULL,opaa=NULL,opw=NULL,bfaa=NULL,ancestral=NULL,C=2,Phi=0.5,q=4e-7,Ne=5e6)
#sample call : optim.s(gene1,ROKAS_TREE,print_level=1,beta=be,gamma=ga,Q=NU_VEC))
optim.s <- function(data, tree, s,method="BOBYQA",maxeval="100",print_level=0, ...){ # s is the initial
  #store information from initial condition, with other parameters fixed
  res.initial = mllm1(data=data,tree=tree,s=s,...)
  #these don't change with the change of s
  dismat = res.initial$dismat
  mumat = res.initial$mumat
  bfaa=res.initial$bfaa
  fn = function(ab,data,tree){
    result = -mllm1(data=data,tree=tree,s=ab,dismat=dismat,mumat=mumat,bfaa=bfaa, ...)$ll$loglik
    return(result)
  }
  lower <- 0 #lower bound
  upper <- Inf #upper bound
  #options for optimizer
  opts <- list("algorithm"=paste("NLOPT_LN_",method,sep=""),"maxeval"=maxeval,
               "xtol_rel"=1e-6,"ftol_rel"=.Machine$double.eps^0.5,"print_level"=print_level)
  res = nloptr(x0=s,eval_f=fn, lb=lower,ub=upper,opts=opts,data=data,tree=tree)
  return(res)
}
# find mle of s for a range of genes in rokas data, given values of beta and gamma
## only apply to rokas data, for other data sets, need to tweak the codes
optim.s.range<- function(beta,gamma,generange,tree,multicore=FALSE,method="BOBYQA",maxeval="100",print_level=0, ...){
  mle.s.one <- function(k){ ## find mle of s for one gene
    ## for now all search start with initial value "1" for s, could let user control the starting value             
    mle <- optim.s(ROKAS_DATA[[k]],tree,s=1,method=method,maxeval=maxeval,print_level=print_level,beta=beta,gamma=gamma,...)
    return(mle)
  }
  if(multicore)
    res <- mclapply(generange,mle.s.one)
  else
    res <- lapply(generange,mle.s.one)
  
  return(res)
}
######################################################################################################
# find mle for beta and gamma, that maximize the likelihood for all genes,Q is given
optim.w <- function(beta,gamma,generange,tree,multicore=FALSE,method="BOBYQA",maxeval="100",print_level=0, ...){
  ## a function of x that return the sum of -loglikelihood values for all genes with s optimized separately for different genes
  ab <- c(beta,gamma)
  fn <- function(ab,generange,tree){
    cat("beta, gamma = ", ab, "\n")
    #call the previous function to optimize s for all genes
    mle <- optim.s.range(ab[1],ab[2],generange,tree=tree,multicore=multicore,method=method,maxeval=maxeval,print_level=print_level,...) 
    mle.val <- sapply(1:length(generange),function(ind) mle[[ind]]$objective) # best -loglikelihood values
    return(sum(mle.val)) #summation of all values
  }
  
  lower <- c(0,0) #lower bound
  upper <- c(Inf,Inf) #upper bound
  #options for optimizer
  opts <- list("algorithm"="NLOPT_LN_SBPLX","maxeval"="100","xtol_rel"=1e-6,"ftol_rel"=.Machine$double.eps^0.5,"print_level"=1)
  res = nloptr(x0=ab,eval_f=fn, lb=lower,ub=upper,opts=opts,generange=generange,tree=tree)
  return(res)
}

######################################################################################################
#MLE for s, beta and gamma, using subplex method by default
#mllm <- function(data,tree,s,beta=be,gamma=ga,Q=NULL,
#             dismat=NULL,mumat=NULL,opaa=NULL,opw=NULL,bfaa=NULL,C=2,Phi=0.5,q=4e-7,Ne=5e6)
#sample call : optim.s.weight(gene1,ROKAS_TREE,0.1,be,ga,Q=NU_VEC))
optim.s.weight <- function(data, tree, s,beta,gamma, method="SBPLX",maxeval="50",print_level=0,...){
  #store information from initial condition, with other parameters fixed
  res.initial = mllm1(data=data,tree=tree,s=s,beta=beta,gamma=gamma,...)
  #these don't change with the change of s
  mumat = res.initial$mumat
  bfaa=res.initial$bfaa
  
  ab <- c(s,beta,gamma) ##initial value
#   ab[ab<=0] <- 1e-4
#   ab <- log(ab)
  fn = function(ab,data,tree){
    #ab <- exp(ab)
    #print(ab)
    result = -mllm1(data=data,tree=tree,s=ab[1],beta=ab[2],gamma=ab[3], mumat=mumat,bfaa=bfaa, ...)$ll$loglik
    return(result)
  }
  lower <- rep(0,3)
  upper <- rep(Inf,3)
  #options for optimizer
  opts <- list("algorithm"=paste("NLOPT_LN_",method,sep=""),"maxeval"=maxeval,"xtol_rel"=1e-6,
               "ftol_rel"=.Machine$double.eps^0.5,"print_level"=print_level)
  res = nloptr(x0=ab,eval_f=fn, lb=lower,ub=upper,opts=opts,data=data,tree=tree)
  return(res)
}

######################################################################################################
## MLE for opw (weights for optimal amino acids)
## supply the function, the gradient function, and the restriction function (sum=1)
#mllm <- function(data,tree,s,beta=be,gamma=ga,Q=NULL,
#             dismat=NULL,mumat=NULL,opaa=NULL,opw=NULL,bfaa=NULL,C=2,Phi=0.5,q=4e-7,Ne=5e6)
#sample call: optim.opw(data,tree,opw=rep(1,20),s=0.1,beta=beta,gamma=gamma,Q=NU_VEC...)
optim.opw <- function(data, tree,opw=NULL,print_level=0, ...){
  if(is.null(opw))
    opw = findBf2(data)
  opw = opw/sum(opw) #opw given in the function call doesn't have to sum to 1
  
  res = mllm(data=data,tree=tree,opw=opw,...) #store llmat (loglikelihood values for all opaa) 
  llmat = exp(res$ll$llmat)    #so that they don't need to be evaluated again and again
  weight = attr(data,"weight")
  #cat("opw optimization, starting loglikelihood = ", res$ll$loglik, "\n") #function value at the starting point
  
  # function to optimize on and its gradient function
  eval_f_list <- function(opw){
    return(list("objective"=llaaw1(opw,weight=weight,llmat=llmat),
                "gradient"=llaaw_grad(opw,weight=weight,llmat=llmat)))
  }
  # linear equality constraint, and its jacobian
  eval_g_list <- function(opw){
    return(list("constraints"=sum(opw)-1,"jacobian"=rep(1,20)))
  }
  lower <- rep(0,20) #lower bound
  upper <- rep(1,20) #upper bound
  local_opts <- list("algorithm"="NLOPT_LD_MMA","xtol_rel"=1e-7) #options for local optimizer
  #options for optimizer
  opts <- list("algorithm"="NLOPT_LD_AUGLAG","maxeval"="1000000","xtol_rel"=1e-7,"ftol_rel"=.Machine$double.eps,
               "local_opts"=local_opts,"print_level"=print_level)
  res = nloptr(x0=opw,eval_f=eval_f_list,eval_g_eq=eval_g_list, lb=lower,ub=upper,opts=opts)
  # res$objective: best function value found
  # res$solution: best parameter values
  return(res)
}

# find mle for beta and gamma, that maximize the likelihood for all genes,Q is given
optim.opw.range <- function(s,beta,gamma,generange,Q,tree,multicore=FALSE,print_level=0){
  if(multicore)
    res <- mclapply(ROKAS_DATA[generange],optim.opw,tree=tree,opw=rep(1,20),
                    print_level=print_level,s=s,beta=beta,gamma=gamma,Q=Q)
  else
    res <- lapply(ROKAS_DATA[generange],optim.opw,tree=tree,opw=rep(1,20),
                  print_level=print_level,s=s,beta=beta,gamma=gamma,Q=Q)
  return(res)
}
## optimize beta, gamma and s, for each gene find the best opw for them, instead of using a universal one
optim.opw.sbg <- function(s,beta,gamma,Q,tree,generange,multicore=FALSE,print_level=0,...){
  ab <- c(s,beta,gamma)
  fn <- function(ab){
    cat("s,beta,gamma = ", ab, "\n")
    mle <- optim.opw.range(ab[1],ab[2],ab[3],generange=generange,Q=Q,tree=tree,multicore=multicore,
                           print_level=print_level,...)
    mle.val <- sapply(1:length(generange), function(x) mle[[x]]$objective)
    return(sum(mle.val))
  }
  lower <- rep(0,3)
  upper <- rep(Inf,3)
  opts <- list("algorithm"="NLOPT_LN_SBPLX","maxeval"="100","xtol_rel"=1e-6,
               "ftol_rel"=.Machine$double.eps^0.5,"print_level"=1)
  res = nloptr(x0=ab,eval_f=fn, lb=lower,ub=upper,opts=opts)
  return(res)
}
optim.all <- function(s,beta,gamma,Q,tree,generange,multicore=FALSE,print_level=0,maxeval="100",...){
  if(is.null(attr(tree,"order")) || attr(tree,"order") == "cladwise")
    tree <- reorderPruning(tree)
  br=tree$edge.length
  br.num = length(br)
  Q <- Q/Q[6]
  ab = c(s,beta,gamma,br,Q[1:5]) #all parameters as a vector, length = 3+br.num+5
  ablen <- length(ab)
  fn <- function(ab){
    s <- ab[1]
    beta <- ab[2]
    gamma <- ab[3]
    cat("s,beta,gamma:",s,beta,gamma,"\n")
    tree$edge.length <- ab[4:(4+br.num-1)]
    cat("branch lengths:", tree$edge.length,"\n")
    Q <- c(ab[(4+br.num):ablen],1)
    cat("Q:",Q,"\n","\n")
    mle <- optim.opw.range(s,beta,gamma,generange=generange,Q=Q,tree=tree,multicore=multicore,print_level=print_level,...)
    mle.val <- sapply(1:length(generange), function(x) mle[[x]]$objective)
    return(sum(mle.val))
  }
  lower <- rep(0,ablen)
  upper <- rep(Inf,ablen)
  opts <- list("algorithm"="NLOPT_LN_BOBYQA","maxeval"=maxeval,"xtol_rel"=1e-6,
               "ftol_rel"=.Machine$double.eps^0.5,"print_level"=1)
  res = nloptr(x0=ab,eval_f=fn, lb=lower,ub=upper,opts=opts)
  par_optim <- res$solution
  res$s = par_optim[1]
  res$GMweights = c(al,par_optim[2:3])
  res$br = par_optim[4:(4+br.num-1)]
  res$Q = c(par_optim[(4+br.num):ablen],1)
  return(res)
}

######################################################################################################
## MLE for branch lengths, specify starting values for el, otherwise it starts with the tree supplied with branch length
#mllm <- function(data,tree,s,beta=be,gamma=ga,Q=NULL,
#             dismat=NULL,mumat=NULL,opaa=NULL,opw=NULL,bfaa=NULL,C=2,Phi=0.5,q=4e-7,Ne=5e6)
#sample call: optim.br(data,tree,el=NULL,s=0.1,beta=be,gamma=ga,Q=NU_VEC...)
optim.br <- function(data,tree,el=NULL,method="COBYLA",maxeval="100", print_level=0, ...){
  if(is.null(el))
  {el <- tree$edge.length}
  br.num <- length(el)
  
  res.initial = mllm1(data=data,tree=tree,...)
  #these don't change with the change of s
  Qall = res.initial$Qall
  bfaa = res.initial$bfaa
  fn = function(el,data,tree){
    tree$edge.length = el
    #print(el)
    result = -mllm1(data,tree, Qall=Qall,bfaa=bfaa,...)$ll$loglik
    return(result)
  }
  lower=rep(0,br.num)
  upper=rep(Inf,br.num)
  #options for optimizer
  opts <- list("algorithm"=paste("NLOPT_LN_",method,sep=""),"maxeval"= maxeval,"xtol_rel"=1e-6,"ftol_rel"=.Machine$double.eps,
               "stopval"=-Inf,"print_level"=print_level)
  res = nloptr(x0=el,eval_f=fn, lb=lower,ub=upper,opts=opts,data=data,tree=tree)
  #print(res)
  return(res)
}

######################################################################################################

## optimize the mutation rates for nucleotides, this is only for GTR model, with 5 parameters
## can be easily changed to account for other models, e.g. Jukes-Cantor, etc.
#mllm <- function(data,tree,s,beta=be,gamma=ga,Q=NULL,
#             dismat=NULL,mumat=NULL,opaa=NULL,opw=NULL,bfaa=NULL,C=2,Phi=0.5,q=4e-7,Ne=5e6)
optimQ <- function(tree,data,Q=rep(1,6),method="SBPLX",maxeval="100",print_level=0, ...){
  Q = Q/Q[6] #make last rate 1
  #store information from initial condition, with other parameters fixed
  res.initial = mllm1(data=data,tree=tree,Q=Q,...)
  #these don't change with the change of s
  fixmatall = res.initial$fixmatall
  bfaa = res.initial$bfaa
  
  ab <- Q[1:5] # optimize on the first 5 rates
  fn = function(ab,tree,data){
    #print(ab)
    result = -mllm1(data,tree,Q=c(ab,1),fixmatall=fixmatall,bfaa=bfaa, ...)$ll$loglik
    return(result)
  }
  lower=rep(0,5)
  upper=rep(Inf,5)
  #options for optimizer
  opts <- list("algorithm"=paste("NLOPT_LN_",method,sep=""),"maxeval"= maxeval,"xtol_rel"=1e-6,
               "ftol_rel"=.Machine$double.eps^0.5,"print_level"=print_level)
  res = nloptr(x0=ab,eval_f=fn, lb=lower,ub=upper,opts=opts,data=data,tree=tree)
  res$solution = c(res$solution,1) # append the last rate (1) to the rate vector
  return(res)
}
######################################################################################################
## optimize parameters
## first optimize all parameters except optimal weights
## then optimize the weights of amino acids being optimal
optim.mllm1 <- function(object, optQ = FALSE, optBranch = FALSE, optsWeight = TRUE, optOpw = FALSE,
                       control = list(epsilon=1e-08,hmaxit=10,htrace=TRUE,print_level=0,maxeval="300"),...){
  tree = object$tree
  ## these are automatically satisfied if object comes from mllm1
#   if(class(tree)!="phylo") stop("tree must be of class phylo") 
#   if(!is.rooted(tree)) stop("tree must be rooted")
#   if (is.null(attr(tree, "order")) || attr(tree, "order") == "cladewise") 
#     tree <- reorderPruning(tree)
  if(any(tree$edge.length < 1e-08)){
    tree$edge.length[tree$edge.length < 1e-08] <- 1e-08
    #object <- update(object, tree=tree)
  }
  call = object$call
  #maxit = control$maxit #maximum number of iterations for each sub optimizer
  htrace = control$htrace #print out information about steps or not?
  print_level=control$print_level
  maxeval = control$maxeval
  data = object$data
  Q = object$Q
  #if(is.null(subs)) subs = c(1:(length(Q)-1),0) #default is GTR
  bfaa = object$bfaa #this is going to be the same, no matter from data or given -- empirical frequencies
  opw = NULL
  ll = object$ll$loglik
  ll1 = ll
  s = object$s
  beta = object$GMweights[2]
  gamma = object$GMweights[3]
  opti = TRUE # continue optimizing or not
  rounds = 0 #index of iterations
  while(opti){
    if(htrace){
      cat("\n","Round ",rounds+1,"\n")
      cat("opw = ", opw, "\n")
    }
    if(optsWeight){
      res = optim.s.weight(data,tree,maxeval=maxeval,print_level=print_level,s=s,beta=beta,gamma=gamma,Q=Q,opw=opw,...)
      s = res$solution[1]
      beta = res$solution[2]
      gamma = res$solution[3]
      if(htrace){
        cat("optimize s and Grantham weights: ", ll, "--->", -res$objective, "\n")
        cat("s, beta, gamma are now: ", res$solution, "\n")
      }
      ll = -res$objective
    }
    if(optQ){
      res = optimQ(tree,data,Q=Q,maxeval=maxeval,print_level=print_level,s=s,beta=beta,gamma=gamma,opw=opw, ...)
      Q = res$solution
      if(htrace){
        cat("optimize rate matrix: ", ll, "--->", -res$objective, "\n")
        cat("Q is now: ", res$solution, "\n")
      }
      ll = -res$objective
    }
    if(optBranch){
      res = optim.br(data,tree,maxeval=maxeval,print_level=print_level,s=s,beta=beta,gamma=gamma,Q=Q,opw=opw, ...)
      if(htrace){
        cat("optimize branch lengths:", ll, "--->", -res$objective, "\n")
        cat("branch lengths are now: ", res$solution, "\n")
      }
      tree$edge.length = res$solution
      ll =-res$objective
    }
    if(optOpw){
      res = optim.opw(data,tree,opw=opw,print_level=print_level,s=s,beta=beta,gamma=gamma,Q=Q,bfaa=bfaa,...) #new optimizer using nloptr
      if(htrace){
        ##notice that the loglikelihood will decrease, from maximizing rule to weighted rule
        cat("optimize weights of optimal aa:", ll, "--->", -res$objective, "\n")
        cat("weights are now: ", res$solution, "\n")
      }
      opw = res$solution
      ll = -res$objective
    }
    rounds = rounds + 1
    if(rounds >= control$hmaxit) opti <- FALSE
    if(((ll1-ll)<=0) && ((ll1-ll)/ll < control$epsilon)) opti <- FALSE
    ll1 = ll
  }
  
  object = update(object, tree=tree,data=data,s=s,beta=beta,gamma=gamma,Q=Q,bfaa=bfaa,opw=opw,...)
  return(object)
}
#MLE for s and Ne, using subplex method by default
#mllm1 <- function(data,tree,s=NULL,beta=be,gamma=ga,Q=NULL,dismat=NULL,fixmatall=NULL,mumat=NULL,Qall=NULL,
#                  opaa=NULL,opw=NULL,bfaa=NULL,C=2,Phi=0.5,q=4e-7,Ne=5e6)
#sample call : optim.s.weight(gene1,ROKAS_TREE,0.1,5000,beta=be,gamma=ga,Q=NU_VEC,...)
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

optim.Ne <- function(data, tree,s,Ne, method="SBPLX",maxeval="500",print_level=0,...){
  #store information from initial condition, with other parameters fixed
  res.initial = mllm1(data=data,tree=tree,s=s,beta=be,gamma=ga,Ne=Ne,...)
  #these don't change with the change of s
  mumat = res.initial$mumat
  bfaa=res.initial$bfaa
  
  ab <- Ne ##initial value
  fn = function(ab,data,tree){
    cat("Ne = ",ab,"\n")
    result = -mllm1(data=data,tree=tree,s=s,beta=be,gamma=ga, mumat=mumat,bfaa=bfaa,Ne=ab, ...)$ll$loglik
    cat(result,"\n")
    return(result)
    
  }
  lower <- 0
  upper <- Inf
  #options for optimizer
  opts <- list("algorithm"=paste("NLOPT_LN_",method,sep=""),"maxeval"=maxeval,"xtol_rel"=1e-6,
               "ftol_rel"=.Machine$double.eps^0.5,"print_level"=print_level)
  res = nloptr(x0=ab,eval_f=fn, lb=lower,ub=upper,opts=opts,data=data,tree=tree)
  return(res)
}
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
# numaa = sapply(1:9128,function(x) length(aa.conf[[x]]$aa))
# op.lik = sapply(1:9128, function(x) aa.conf[[x]]$op.per)

# heatmap of probabilities that each amino acid accounts for
hmap.op <- function(llmat){
  probmat <- exp(llmat)
  lsite <- dim(probmat)[1]
  #probmat <- t(sapply(1:lsite,function(x) expmat[x,]/sum(expmat[x,])))
  #heatmap(probmat[1:lsite,],Rowv=NA,Colv=NA,col=grey.colors(256))
  heatmap(t(probmat),Rowv=NA,Colv=NA,col=grey.colors(256))
}

