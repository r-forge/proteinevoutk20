#Functionality function, arguments: 
#d: distance between the optimal and observed amino acids, a vector of length n
#s: seleciton strength coeffecient vector, comes from a probability distribution, e.g. exponential or gamma
#n: length of the protein
#dis: distribution function from which s is drawn. 1: exponential, 2: gamma, 3: uniform(0,1)
#para: parameter for the distribution. Rate for exponential, shape for gamma
#model: model number. 1: exponential 2: parabolic
Ftny <- function(d, s, model=1, dis=0, para=1){  #default: n=15, exponential distribution with rate 1
  #browser()
  #different distributions where s comes from
  n <- length(d) #length of the protein as an amino acid sequence
  if(dis==0){
    s <- rep(s,n)
  }
  else if(dis==1)
    {s <- rexp(n,para)} #random generation of n numbers from exponential distribution
  else if(dis==2) #gamma distribution
    {s <- rgamma(n,para)}
  else if(dis==3) # uniform distribution
    {s <- runif(n)}

  #different models for functionality function
  if(model==1){ #model 1
    result <- prod(exp(-d*s)) #the function
  }
  else if(model==2){ #model 2
    result <- prod(1/(1+d*s)) 
  }
  return(result)
}

#fitness function
#q, Phi: constants
#C: cost function, linearly function of the length of the protein
#ftn: functionality value
#C_n = a1 + a2*n, cost --  linear function of the length
#n: length of the protein

#fixation probability function from protein 1 to protein 2,
#which is also the instantaneous rate from protein1 to protein2
#parameters: fit1, fit2: fitnesses of protein 1 and 2
# n: length of protein
# a1, a2: parameters in cost function C_n
#Ne: effective population size
#formula number: 1: Sela & Hirsh's formula (default); 2: classic formula

#####
#try to put it in a further calculated form, with same terms for ftny's pulle out
#as common factors. To do this, have to recognize which ONE is the one that's different
#in the two proteins
#####
fix <- function(ftny1, ftny2, n, a1, a2, Ne, formula=1){
  C_n <- a1 + a2*n
  fit_ratio <- exp(-q*Phi*C_n*((1/ftny2)-(1/ftny1))) #w2/w1
  if(formula==1){
    result <- (1-fit_ratio)/(1-fit_ratio^(2*Ne))
  }
  else if(formula==2){
    se <- fit_ratio - 1 # s=(w2-w1)/w1
    result <- (1-exp(-2*se))/(1-exp(-4*Ne*se))
  }
  else{
    print("Error! Please choose 1 or 2 as formula number. 1: Hirsh's; 2: classic pop gen's")
  }
  return(result)
}

#transition rate from protein 1 to protein 2
#Q <- Ne*mu*fix

#Assign values to different parameters
##################################
#physiochemical distance matrix#
##################################
#notice that the proteins in matrix rows and columns are not in the same order
GranthamMatrix <- read.csv("GranthamMatrix.csv",
                           header=TRUE, sep=",",row.names=1)
#View(GranthamMatrix) #check the matrix
GM <- GranthamMatrix #store it in a different matrix and manipulate it
#complete the upper triangle matrix to whole, with diagonals equal to 0
GM[1,1] <- 0
for(i in 2:20){
  GM[i,i] <- 0
  for(j in 1:(i-1)){
    GM[i,j] <- GM[j,i]
  }
}
GM <- GM/100 #make the mean 1 instead of 100


#Function to generate the transition matrix (20*20 for now) 
#given the optimal amino acid
#parameter: op(aa_op): optimal amino acid at this site
mat_gen <- function(op){
  mat <- matrix(nrow=20,ncol=20) #set diagonal entries to be 0 at first
  for(i in 1:19){
    mat[i,i] <- 0
    for(j in (i+1):20){
      d1 <- GM[i,op] #distance between aa_i and aa_op
      d2 <- GM[j,op]
      fn1 <- Ftny(d1,s,n=1,model=1,dis=0) #first pass, use model 1 and constant s
      fn2 <- Ftny(d2,s,n=1,model=1,dis=0)
      mat[i,j] <- fix(fn1,fn2,n=1,a1,a2,Ne,formula=1) #fixation prob -> transition rate
      fit_ratio <- exp(-q*Phi*C_n*((1/fn2)-(1/fn1))) #w2/w1
      mat[j,i] <- mat[i,j]*(fit_ratio)^(2*Ne-1) #symmetric entry
    }#end for j
    mat[i,i] <- -(sum(mat[i,])) #diagonal entry; sum of each row is 0
  }#end for i
  return(mat)
}

matgen_indep <- function(arr,op){
  m <- dim(arr)[1]
  n <- dim(arr)[2]
  mat <- matrix(NA,dim=c(m,m))
  for(i in 1:(m-1)){
    mat[i,i] <- 0
    for(j in (i+1):m){
      dis <- sum(abs(arr[i,] - arr[j])) #number of different sites
      if(dis>1){
        mat[i,j] <- 0
      }
      else{
        d1 <- GM[i,op]
        d2 <- GM[j,op]
        fn1 <- Ftny(d1,s,m,model=1,dis=0)
        fn2 <- Ftny(d2,s,m,model=1,dis=0)
        mat[i,j] <- fix(fn1,fn2,n=1,a1,a2,Ne,formula=1)
        fit_ratio <- exp(-q*Phi*C_n*((1/fn2)-(1/fn1))) #w2/w1
        mat[j,i] <- mat[i,j]*(fit_ratio)^(2*Ne-1) #symmetric entry
      }
    }
    mat[i,i] <- -sum(mat[i,])
  }
  return(mat)
}
#Given 2 amino acid sequences, find the transition rate matrix using parsimony path
#1.Find all the sequences on the path
#2.Calculate the entries between 2 sequences that differ at only one site
#3.If more than 1 site differ, the rate is 0
#4.What about the optimal sequence? (Right now the optimal amino acid sequence is given
#If the number of sites that are differ is m, then first step, there are m choices to make, 
# second step, there are m-1 choices, then m-2 choices, hence the total number of pasimony paths
# is m!. How many sequences are possible at each step? 1st: m; 2nd: C(m,2); 3rd: C(m,3)..
# Hence in total the number of sequences need to consider is 2^m. Use binary numbers to do the
# count. From 0000..0 to 111..1


#Given two protein sequences of the same lengths, calculated the probability (likelihood) of 
#going from protein1 to protein2 using site independent method.
#Calculate the probability of going from one amino acid to another, then
#multiply them together to get the total probability.
#Parameter: branch length t. There is only one branch and two nodes in this simplest case
prob <- function(protein1, protein2, protein_op, t){
  l <- length(protein1) #length of protein
  
}
  #parameter values
  Ne <- 10000
  q <- 10^(-6)
  Phi <- 100
  mu <- 0.0001
  a1 <- 4 #reasonable values or not??
  a2 <- 2
  n <- 8
########################################################################
########################################################################
#   pr1 <- sample(1:20,n,replace=T) #amino acids in protein1
#   pr2 <- pr1 #amino acids in protein2
#   pos <- sample(1:n,1) # position where pr2 has different aa from pr1
#   pr2[pos] <- sample(1:20,1) #assign a different amino acid at this position
#   pr0 <- sample(1:20,n,replace=T) #amino acids in optimal protein
#   d1 <- vector(length=n) #distance vector from pr1 to pr0
#   d2 <- d1 #distance vector from pr2 to pr0
#   for(i in 1:n){
#     d1[i] <- GM[pr1[i],pr0[i]]
#     d2[i] <- GM[pr2[i],pr0[i]]
#   }

  #Plot the fixation probability of going from one protein to another (fixed)
  #as a function of selection strength s. s ranges from 0.01 to 1.
fixplot <- function(n,model){
  pr1 <- sample(1:20,n,replace=T) #amino acids in protein1
  pr0 <- sample(1:20,n,replace=T) #amino acids in optimal protein
  d1 <- vector(length=n) #distance vector from pr1 to pr0
  for(i in 1:n){
    d1[i] <- GM[pr1[i],pr0[i]]
  }
  d2 <- d1
  pos <- sample(1:n,1) # position where pr2 has different aa from pr1
  s0 <- seq(0.01,1,by=0.01)
  l <- length(s0)
  d0 <- seq(-0.99,0.99,length.out=l)
  fix_1 <- matrix(nrow=l,ncol=l)
  fix_2 <- matrix(nrow=l,ncol=l)
  for(i in 1:l){
    for(j in 1:l){
      s <- rep(s0[i],n) #selection strength
      d2[pos] <- d2[pos] + d0[j]
      ftn1 <- Ftn(d1,s,n,model,dis=0) #functionality
      ftn2 <- Ftn(d2,s,n,model,dis=0)
      fit1 <- Fit(q,Phi,ftn1,a1,a2,n) #fitness
      fit2 <- Fit(q,Phi,ftn2,a1,a2,n)
      fix_1[i,j] <- fix1(fit1,fit2,Ne) #fixation probability
      fix_2[i,j] <- fix2(fit1,fit2,Ne)
    }
  }
  par(mfrow=c(1,2))
  persp(s0,d0,fix_1,col="red",scale=F)
  persp(s0,d0,fix_2,col="blue",scale=F)
  #plot(s0,fix_1,type="l",col="red",ylab="fixation prob")
  #points(s0,fix_2,type="l",col="blue")
}

fixplot_s <- function(n,model,d0){
  pr1 <- sample(1:20,n,replace=T) #amino acids in protein1
  pr0 <- sample(1:20,n,replace=T) #amino acids in optimal protein
  d1 <- vector(length=n) #distance vector from pr1 to pr0
  for(i in 1:n){
    d1[i] <- GM[pr1[i],pr0[i]]
  }
  d2 <- d1
  pos <- sample(1:n,1) # position where pr2 has different aa from pr1
  s0 <- seq(0.01,1,by=0.01)
  l <- length(s0)
  fix_1 <- vector(length=l)
  fix_2 <- vector(length=l)
  for(i in 1:l){
      s <- rep(s0[i],n) #selection strength
      d2[pos] <- d2[pos] + d0
      ftn1 <- Ftn(d1,s,n,model,dis=0) #functionality
      ftn2 <- Ftn(d2,s,n,model,dis=0)
      fit1 <- Fit(q,Phi,ftn1,a1,a2,n) #fitness
      fit2 <- Fit(q,Phi,ftn2,a1,a2,n)
      fix_1[i] <- fix1(fit1,fit2,Ne) #fixation probability
      fix_2[i] <- fix2(fit1,fit2,Ne)
    }
  par(mfrow=c(1,2))
  plot(s0,fix_1,type="l",col="red",ylab="fixation prob")
  points(s0,fix_2,type="l",col="blue")
  plot(s0,fix_1-fix_2,type="l",col="green")
}

fixplot_d <- function(n,model,s0){
  pr1 <- sample(1:20,n,replace=T) #amino acids in protein1
  pr0 <- sample(1:20,n,replace=T) #amino acids in optimal protein
  d1 <- vector(length=n) #distance vector from pr1 to pr0
  for(i in 1:n){
    d1[i] <- GM[pr1[i],pr0[i]]
  }
  d2 <- d1
  pos <- sample(1:n,1) # position where pr2 has different aa from pr1
  d0 <- seq(-0.99,0.99,by=0.01)
  l <- length(d0)
  fix_1 <- vector(length=l)
  fix_2 <- vector(length=l)
  for(i in 1:l){
      s <- rep(s0,n) #selection strength
      d2[pos] <- d2[pos] + d0[i]
      ftn1 <- Ftn(d1,s,n,model,dis=0) #functionality
      ftn2 <- Ftn(d2,s,n,model,dis=0)
      fit1 <- Fit(q,Phi,ftn1,a1,a2,n) #fitness
      fit2 <- Fit(q,Phi,ftn2,a1,a2,n)
      fix_1[i] <- fix1(fit1,fit2,Ne) #fixation probability
      fix_2[i] <- fix2(fit1,fit2,Ne)
    }
  par(mfrow=c(1,3))
  plot(d0,fix_1,type="l",col="red",ylab="fixation prob")
  plot(d0,fix_2,type="l",col="blue",ylab="fixation prob")
  plot(d0,fix_1-fix_2,type="l",col="blue")
}

# N <- 10000
# ff <- seq(1,10,by=0.1)
# ll <- length(ff)
# ffix <- vector(length=ll)
# for(i in 1:ll){
#   ffix[i] <- (1-ff[i])/(1-(ff[i])^(2*N))
# }
# plot(ff,ffix,type='l')