#tester on nucleotide evolution
ll <- function(tree,data,Q){
    #browser()
    bf = rep(1/4,4) #base frequency, randomly chosen from all states
    tree <- ape:::reorder.phylo(tree,"p") #reorder the tree in pruningwise order
    edge = tree$edge #edges
    nNodes = max(edge) #number of nodes in the tree
    probvec = matrix(NA,nNodes,4) #probability of gettting different states at nodes that evolve to the current sequences
    parent <- as.integer(edge[, 1]) #parents of the edges
    child <- as.integer(edge[, 2]) #children of the edges
    root <- as.integer(parent[!match(parent, child, 0)][1])  
    tip <- as.integer(child[!match(child,parent,0)])
    for(i in 1:length(tip)){
      probvec[i,] = 0
      probvec[i,data[i]] = 1
    }
    #or use this to find tips:
    #tip <- c(1:length(tree$tip.label))
    tl = tree$edge.length #lengths of the edges
    for(i in 1:tree$Nnode){
        from = parent[2*i] 
        to = child[(2*i-1):(2*i)]
        t_left = tl[2*i-1]
        t_right = tl[2*i]
        for(j in 1:m){
          #find the likelihood of a node when likelihood of both its descendents is found --- pruning
          probvec[from,j] = (expm(Q*t_left)[j,] %*% probvec[to[1],])*(expm(Q*t_right)[j,] %*% probvec[to[2],]) 
        }
    }
    probvec[root,] %*% bf
}


parse_array <- function(Array,sites){
  ParsedArray <- array(dim=c(dim(Array)[1],dim(Array)[2]-sites+1))
  l <- dim(ParsedArray)[2]
  for(i in 1:dim(Array)[1]){
    ParsedArray[i,1] <- paste(Array[i,seq(1,sites,by=1)],collapse="_")
  }
  ParsedArray[,seq(2,l,by=1)] <- Array[,-seq(1,sites,by=1)]
  ParsedArray
}
simOut0.1_parsed <- parse_array(sim0.1,2)
data <- simOut0.1_parsed[,c(1,5)]
data <- data[order(data[,1]),]
data_unique <- unique(data)
data_level <- table(data[,1])
counts <- as.numeric(data_level)
freq <- counts/sum(counts)
Ne <- 1.37*10^3
nu <- 2*Ne-1
fit_vec <- as.numeric(data_unique[,2])
freq_fit <- fit_vec^nu/sum(fit_vec^nu)

ftny_vec <- NULL
for(i in 1:3){
  distance <- pchem_d(i,2)
  ftny_vec <- c(ftny_vec,Ftny(distance,0.1))
}
fit_vec <- exp(-q*Phi*C_n/ftny_vec)

ftny_vec <- NULL
for(i in 1:dim(arr)[1]){
  distance <- pchem_d(arr[i,],c(2,2))
  ftny_vec <- c(ftny_vec,Ftny(distance,c(0.1,0.1)))
}
fit_vec <- exp(-q*Phi*C_n/ftny_vec)

sim <- simOutInd1[-seq(1,1000,by=1),]
sim_level <- table(sim[,9])
fit_vec <- as.numeric(dimnames(sim_level)[[1]])
counts <- as.vector(sim_level)
freq <- counts/sum(counts)
Ne <- 1.37*10^7
log_fit_vec <- log(fit_vec)
Fit1 <- (fit_vec/fit_vec[310])^(2*Ne-1)
k <- 310
i <- 100
exp((log_fit_vec[k]-log_fit_vec[i])*Ne)

parse_simOut0.1 <- parse_array(simOut0.1,5)
nu <- 2*Ne-1
exp(-q*Phi*C*nu*(1/fit_vec[k]-1/fit_vec[i]))

y <- as.vector(table(simOutInd1[-seq(1,10000,1),2]))
y <- y/sum(y)
mat <- mat_gen_indep(2,1)
x <- expm(mat*10^11)[1,]   
x-y

sim_simple <- function(start_state,matrix,t){
  browser()
  t_now <- 0
  P <- matrix
  path <- array(c(start_state,0),dim=c(1,2))
  colnames(path) <- colnames(path,do.NULL=FALSE)
  colnames(path)[1] <- "State"
  colnames(path)[2] <- "Time Now"
  while(t_now < t){
    lambda <- -P[start_state,start_state]
    t_wait <- rexp(1,lambda)
    t_now <- t_now + t_wait
    if(t_now < t){
      rates <- P[start_state,-start_state]/(-P[start_state,start_state])
      states <- c(1:9)[-start_state]
      start_state <- states[mkv(runif(1),rates)]
      path <- rbind(path,c(start_state,t_now))
    }
  }
  path
}

data_simple <- sim_simple(2,P,100000)
data_ext <- sim_extract(data_simple,10)
counts <- as.numeric(table(data_ext[,1]))
freq <- counts/sum(counts)

###function to find the average time of m transitions in a simulation###
avg_time <- function(sim,m){
  #time <- sim[,"Time Now"] #extract the times when the transitions happen
  time <- sim[,2]
  l <- length(time)
  time <- head(time,l-1) #in the last step, there is no transition
  time_start <- head(time, l-1-m)
  time_end <- tail(time, l-1-m)
  return(mean(time_end-time_start))
}

##Extract the states every m transitions###
sim_extract <- function(sim,m){
  time_trans <- avg_time(sim,m) #average time of m transitions
  #time <- sim[,"Time Now"] #times of transition happening
  time <- sim[,2]
  time_mod <- time %/% time_trans #
  pick_time_vec <- time_mod[-length(time_mod)]-time_mod[-1]
  return(sim[pick_time_vec < 0,])
}

sim0.1 <- sim_extract(simOut0.1,10)
sim0.1 <- sim0.1[,1]
vec <- as.numeric(table(sim0.1))


f_tsn <- function(t){
  R <- 10
  return(0.25 - 0.5*exp((-((2*R-1)/(R+1))*t)) + 0.25*exp((-(2)/(R+1))*t))
}

f_tsvn <- function(t){
  R <- 10
  return(0.5-0.5*exp((-2/(R+1))*t))
}
a <- seq(0,10,0.001)
tsn_rate <- as.numeric(lapply(a,f_tsn))
plot(a,tsn_rate,type="l",ylim=c(0,0.6))
tsvn_rate <- as.numeric(lapply(a,f_tsvn))
points(a,tsvn_rate,type="l",col='red')


TREE <- "((1:500,2:500):500,(3:500,4:500):500);" #symmetric topology
TREE <- "(((1:500,2:500):500,3:500):500,4:500);" #comb topology
TREE <- "(((1,2),(3,4)),((5,6),(7,8)));"
length <- rep(500,14)
length <- rep(50,6)
length <- rep(10,6)
length <- rep(1000,6)
length <- rep(10000,6)
length <- rep(100,10000,10000,200,10000,10000)
length <- c(50,200,200,50,50,50)
length <- c(200,50,200,200,100,100)
length <- c(200,50,100,50,100,50)
tree <- read.tree(text=TREE)
tree$edge.length <- length

string <- paste("fourtip/comb",1:9,".RData",sep="")
pdf(file="combtree.pdf",paper="us")
par(mfrow=c(3,3))
for(index in 1:9){
  load(string[index])
  length <- tree$edge.length
  tle <- paste(length[1],length[2],length[3],length[4],length[5],length[6], sep=",")
  plot(tree,show.tip.label=F)
  nodelabels()
  tiplabels()
  title(tle)
}
dev.off()


simTree(tree,l=10,protein_op=rep(1,10),10,0.1,indep=T,rootseq=c(1:10),ancestral=T)
sim <- simTree(tree,l=100,protein_op=rep(1,100),10,0.01,indep=T,rootseq=sample(10,20,replace=TRUE),ancestral=T)
simulation(c(1:10),rep(1,10),1000,10,0.1,indep=T,model=2,record.fitness=F)

sim1.1000 <- simTree(tree,1,rep=1000,protein_op=1,m=10,s=0.1,rootseq=NULL,ancestral=T)
sim1.1000.tip <- sim1.1000[1:4,]
tree$edge.length
sim0.1_1000[,1:27]
sim0.01_1000[,1:27]
sim0_1000[,1:27]
sim0.1_1000[,1:27]
sim1_1000[,1:27]
sim10_1000[,1:27]
mle0.1
mle1
mle10

rm(list=ls())
load("eighttip/comb8tipSites6.Rdata")
mle_0 <- c(mle0_200$position,mle0_400$position,mle0_800$position,mle0$position)
mle_0.01 <- c(mle0.01_200$position,mle0.01_400$position,mle0.01_800$position,mle0.01$position)
mle_0.1 <- c(mle0.1_200$position,mle0.1_400$position,mle0.1_800$position,mle0.1$position)
mle_1 <- c(mle1_200$position,mle1_400$position,mle1_800$position,mle1$position)
mle_0 <- abs(mle_0-0)
mle_0.01 <- abs((mle_0.01-0.01)/0.01)
mle_0.1 <- abs((mle_0.1-0.1)/0.1)
mle_1 <- abs((mle_1 - 1)/1)
x <- c(200, 400, 800, 1000)
pdf(file="comb8tipSites.pdf")
par(mfrow=c(2,2))
plot(x,mle_0,type='o',main="s=0")
plot(x,mle_0.01,type='o',main="s=0.01")
plot(x,mle_0.1,type='o',main="s=0.1")
plot(x,mle_1,type='o',main="s=1")
dev.off()

load("comb_eqlm.RData")
mle_0.01 <- c(mle0.01_200$position,mle0.01_400$position,mle0.01_800$position)
mle_0.1 <- c(mle0.1_200$position,mle0.1_400$position,mle0.1_800$position)
mle_0.01 <- abs((mle_0.01-0.01)/0.01)
mle_0.1 <- abs((mle_0.1-0.1)/0.1)
x <- c(200, 400, 800)
pdf(file="comb_eqlm.pdf")
par(mfrow=c(1,2))
plot(x1,mle_0.01,type='o',main="s=0.01")
plot(x1,mle_0.1,type='o',main="s=0.1")
dev.off()

GranthamMatrix <- read.csv("GranthamMatrix.csv",
                           header=TRUE, sep=",",row.names=1)
#View(GranthamMatrix) #check the matrix
GM1 <- GranthamMatrix #store it in a different matrix and manipulate it
#complete the upper triangle matrix to whole, with diagonals equal to 0
GM1[1,1] <- 0
for(i in 2:20){
  GM1[i,i] <- 0
  for(j in 1:(i-1)){
    GM1[i,j] <- GM1[j,i] #symmetric matrix
  }
}

avg_dis <- function(vector){
  l <- length(vector)
  Sum <- 0
  for(i in 1:(l-1)){
    Sum <- Sum + sum(abs(vector[i]-vector[(i+1):l]))
  }
  ct <- l*(l-1)/2
  return((ct/Sum)^2)
}

C <- matrix(result,nrow=20)
Clim <- range(B)
Clen <- Clim[2]-Clim[1]+1
colorlut <- terrain.colors(clen,alpha=0)
colorlut <- terrain.colors(Clen,alpha=0)
col <- colorlut[C-Clim[1]+1]
open3d()
rgl.surface(A,B,C,color=col,alpha=0.75,back="lines")
colorlut <- heat.colors(Clen,alpha=1)
col <- colorlut[C-Clim[1]+1]
rgl.surface(A,B,matrix(1580,nrow(C),ncol(C)),color=col,back="fill")

for(i in 2:19){
  rdata <- paste("contour",i,".RData",sep="")
  load(rdata)
  C <- c(C,result)
}

#order the data frame d by row1, then row 2, then row 3 and row 4 at last
d[with(d,order(d[1,],d[2,],d[3,],d[4,]))]
#Question: how to do the same with data frame with arbitrary number or rows?

data <- read.fasta("RokasGene1.fasta")
l <- length(data[[1]]) #number of nucleotides in this gene
l_aa <- l/3 #number of amino acids in the gene
aa_order <- sample(l_aa,l_aa,replace=T) #orders of amino acids in the new sample
data_resample <- data
for(i in 1:length(data)){
  d <- data[[i]] #original data for i-th species
  resample <- vector(mode="character",length=1701)
  for(j in 1:l_aa){
    resample[((j-1)*3+1):(j*3)] <- d[((aa_order[j]-1)*3+1):(aa_order[j]*3)]
  }
  data_resample[[i]] <- resample
}

library("seqinr")
conv <- function(filename){
  levels <- levels(factor(s2c("SRLPTAVGIFYCHQNKDEMW")))
  data <- read.fasta(file=filename)
  seqdata <- vector()
  for(i in 1:length(data)){
    aas <- translate(seq=data[[i]])
    aan <- as.numeric(factor(aas,levels))
    seqdata <- rbind(seqdata,aan)
  }
  dimnames(seqdata)[[1]]<- attr(data,"name")
  return(seqdata)
}

Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

optimal <- apply(data,2,Mode)

#all the names of the nucleotides and codons that'll be needed 
aa <- s2c("SRLPTAVGIFYCHQNKDEMW")
aa <- levels(factor(aa)) #all the amino acids
Nu <- s2c("tcag") #nucleotide TCAG
stop_cd <- c("taa","tag","tga") #stop codons
nu_list <- NULL #list of all 64 codons (including stop codons)
for(i in 1:4){
  for(j in 1:4){
    for(k in 1:4){
      nu_list <- c(nu_list,c2s(c(Nu[i],Nu[j],Nu[k])))
    }
  }
}
cds <- nu_list[!nu_list %in% stop_cd] #61 non-stop codons

#Given a vector of length 6 (number of parameters for GTR matrix)
#get the GTR rate matrix Q
mat_form <- function(vec){
  rate_mat <- matrix(0,nrow=4,ncol=4)
  rate_mat[1,2] <- vec[1]
  rate_mat[1,3] <- vec[2]
  rate_mat[1,4] <- vec[3]
  rate_mat[2,3] <- vec[4]
  rate_mat[2,4] <- vec[5]
  rate_mat[3,4] <- vec[6]
  rate_mat[lower.tri(rate_mat)] <- rate_mat[upper.tri(rate_mat)]
  for(i in 1:4){ #row sums equal to 0
    rate_mat[i,i] <- -sum(rate_mat[i,-i])
  }
  return(rate_mat)
}
vec <- runif(6)

#Find the mutation rate matrix (20 by 20) from the 6 rates in GTR mutation rate
#matrix for nucleotides
aa_MuMat_form <- function(vec){
  Q <- mat_form(vec)
  dimnames(Q) <- list(Nu, Nu) #change the dimnames to be the nucleotides
  #sfreq <- expm(Q*1000)[1,] #stationary distribution
  
  #mutation matrix for 61 codons
  codon_array <- array(0,dim=c(61,61),dimnames=list(cds,cds))
  for(m in 1:61){ #loop through all 61 codons
    cd_str <- cds[m]
    cd <- s2c(cd_str) #one codon, convert from string to char vector
    #print(paste("starting from",cd_str,":" )) 
    for(i in 1:3){ #every nucleotide in the codon can mutate
      ngb <- cd #neighboring codon, set it equal to the starting codon for now
      #prob <- prod(sfreq[cd]) #equilibrium frequency of this codon
      for(j in 1:4){
        if(Nu[j]!=cd[i]){ #nucleotide has to change to a different one
          ngb[i] <- Nu[j] 
          a_new <- translate(ngb) #new amino acid
          ngb_str <- c2s(ngb)
          #print(paste("neighbor of",cd_str,":",ngb_str))
          if(a_new != "*")
            codon_array[cd_str,ngb_str] <- codon_array[cd_str,ngb_str] + Q[cd[i],Nu[j]]
          #else
          #print("stop codon!")
        }
      }
    }
  }
  #diagonal equals negative of sum of off-diagonal entries
  for(i in 1:61){ #row sums equal to 0
    codon_array[i,i] <- -sum(codon_array[i,-i])
  }
  #equilibrium frequenciew for 61 codons
  cd_eqlm <- expm(codon_array*10000)[1,]
  # aas <- vector()
  # for(i in 1:61){
  #   aas <- c(aas,translate(s2c(cds[i])))
  # }
  
  #20*20 array: mutations rates between amino acids
  arr <- array(0,dim=c(20,20),dimnames=list(aa,aa)) #mutation rate matrix for amino acids
  for(i in 1:61){
    for(j in 1:61){
      from_codon <- cds[i]
      to_codon <- cds[j]
      from_aa <- translate(s2c(from_codon))
      to_aa <- translate(s2c(to_codon))
      arr[from_aa,to_aa] <- arr[from_aa,to_aa] + cd_eqlm[from_codon]*codon_array[from_codon,to_codon]
    }
  }
  return(arr)
}