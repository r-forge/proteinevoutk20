####################################################################################
Ftny_sim <- function(d, s, model=1, dis=0, para=1){
  #browser()
  #different distributions where s comes from
  n <- length(d) #length of the protein as an amino acid sequence
  if(dis==0){
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
} #This is used to do simulation, not the calculation
#############################################################################
#Functions to find the means for categories of discrete gamma distribution
#x*dgamma(x)
xgamma <- function(x,shape,scale){
  return(x*dgamma(x,shape,rate=scale))
}

#Given the shape and scale of the gamma distribution,
#the number of categories of discretization.
#Find the mean of portion of the gamma distribution falling in 
#each category
#THIS FUNCTION WILL RETURN ERROR WHEN THE SHAPE PARAMETER IS TOO SMALL
#BECAUSE THE INTEGRATION VALUE WILL BE NON-FINITE!!!!
dis_gamma <- function(shape,scale,category=4){
  #browser()
  z <- vector(length=(category+1))
  x <- vector(length=category)
  z[1] <- 0
  z[category+1] <- Inf
  for(i in 2:category){
    z[i] <- qgamma((i-1)/category,shape=shape,rate=scale)
  }
  for(j in 1:category){
    x[j] <- integrate(xgamma,shape=shape,scale=scale,lower=z[j],upper=z[j+1])$value/(integrate(dgamma,shape=shape,rate=scale,lower=z[j],upper=z[j+1])$value)
  }
  return(x)
}

#means <- dis_gamma(sh,sc)
####################################################################################################
#Given two protein sequences of the same lengths, calculated the probability (-loglikelihood) of 
#going from protein1 to protein2 using site independent method.
#Calculate the probability of going from one amino acid to another, then
#multiply them together to get the total probability.
#Parameter: branch length t. There is only one branch and two nodes in this simplest case
prob_indep_gamma <- function(protein1, protein2, protein_op, t,shape,category=4,formula=1,model=1,a1=2, a2=1, Phi=0.5,q=4*10^(-7),Ne=1.37*10^7){
  #browser()
  l <- length(protein1) #length of protein
  means <- dis_gamma(shape,shape,category) #means for discrete gamma distribution
  pr <- 0
  for(i in 1:l){ # for each amino acid in the protein, calculate the probability of transition
    p_t <- vector() #store the probabilities for all categories
    for(j in 1:category){
        ma <- mat_gen_indep(protein_op[i],means[j],formula, model,a1,a2,Phi,q,Ne) #transition matrix Q
        Pt <- expm(ma*t) #P(t) = exp(Qt)
        p_t <- c(p_t,Pt[protein1[i],protein2[i]])#the entry corresponding to the amino acids in the sequences
    }
    p <- mean(p_t) #f(x)=(1/k)*sum(f(x|means))
    pr <- pr-log(p)
  }  
  return(pr)
}

prob_dep_gamma <- function(protein1,protein2,protein_op,t,shape,category=4,formula=1,model=1,a1=2, a2=1, Phi=0.5,q=4*10^(-7),Ne=1.37*10^7){
  Arr <- path(protein1,protein2) # all the sequences lying on the path from protein 1 to protein 2
  m <- dim(Arr)[1]
  l <- length(protein1) #number of amino acids
  means <- dis_gamma(shape,shape,category) #means in discrete gamma distribution for each category
  means_index <- integer.base.b(seq(0,category^l-1),category)+1 #all possible combinations of (discrete) gamma means
  index <- dim(means_index)[1] #number of assignments of means to each site (s)
  pr <- vector(length=index)
  for(i in 1:index){
    s <- means[means_index[i,]]
    ma <- mat_gen_dep(Arr,protein_op,s,formula,model,a1,a2,Phi,q,Ne)
    Pt <- expm(ma*t)
    pr <- c(pr,Pt[1,m]) #first row last column entry of P(t)
  }
  p <- mean(pr)
  return(-log(p))
}


prob_gamma <- function(protein1, protein2, protein_op, t, shape, indep=TRUE,category=4,formula=1,model=1,a1=2, a2=1, Phi=0.5,q=4*10^(-7),Ne=1.37*10^7){
  if(indep){ #site independent case
    prob_indep_gamma(protein1,protein2,protein_op,t,shape,category,formula,model,a1,a2,Phi,q,Ne)
  }
  else{ #site dependent case
    prob_dep_gamma(protein1,protein2,protein_op,t,shape,category,formula,model,a1,a2,Phi,q,Ne)
  }
  return(pr)
}
###################################################################################################
#Fixed rate for s -> likelihood
#Independent-site version
prob_indep_fix <- function(protein1,protein2,protein_op,t,s=1,formula=1,model=1,a1=2, a2=1, Phi=0.5,q=4*10^(-3),Ne=1.37*10^3){
  l <- length(protein1) #length of protein
  pr <- 0
  for(i in 1:l){ # for each amino acid in the protein, calculate the probability of transition
    ma <- mat_gen_indep(protein_op[i],s,formula, model,a1,a2,Phi,q,Ne) #transition matrix Q
    Pt <- expm(ma*t) #P(t) = exp(Qt)
    p <- Pt[protein1[i],protein2[i]]#the entry corresponding to the amino acids in the sequences
    pr <- pr-log(p)
  }  
  return(pr)
}


#Dependent-site version
prob_dep_fix <- function(protein1,protein2,protein_op,t,s=1,formula=1,model=1,a1=2, a2=1, Phi=0.5,q=4*10^(-3),Ne=1.37*10^3){
  Arr <- path(protein1,protein2) # all the sequences lying on the path from protein 1 to protein 2
  m <- dim(Arr)[1] # number of amino acid sequences on the parsimony path
  l <- length(protein1) #number of amino acids in the protein
  sv <- rep(s,l) #s is the same for all sites
  ma <- mat_gen_dep(Arr,protein_op,sv,formula,model,a1,a2,Phi,q,Ne) #transition rate matrix for the path Q
  Pt <- expm(ma*t) # P(t) = Exp(Q*t)
  p <- Pt[1,m] #first row last column entry of P(t)
  return(-log(p))
}

prob_fix <- function(para, protein1, protein2, protein_op, indep=TRUE,formula=1,model=1,a1=2, a2=1, Phi=0.5,q=4*10^(-3),Ne=1.37*10^3){
  t <- para[1]
  s <- para[2]
  if(indep){ #site independent case
    prob_indep_fix(protein1,protein2,protein_op,t,s,formula,model,a1,a2,Phi,q,Ne)
  }
  else{ #site dependent case
    prob_dep_fix(protein1,protein2,protein_op,t,s,formula,model,a1,a2,Phi,q,Ne)
  }
}
#####################################################
#Find the rates of moving to other amino acids for only one site, with other sites fixed either at the optimal or non-optimal amino acid
#change_site is the site that is evolving, default is the first one
rate_move_1site <- function(protein,protein_op,m,change_site=1,s,indep=T,formula=1,model=1,a1=2,a2=1,Phi=0.5,q=4*10^(-3),Ne=1.37*10^3,mu=1/(2*Ne)){
  #browser()
  rates <- vector()
  d1 <- pchem_d(protein,protein_op)
  d2 <- d1
  aa <- seq(1:m)[-protein[change_site]]
  for(i in 1:(m-1)){
    d2[change_site] <- GM[protein_op[change_site],aa[i]]
    rates <- c(rates,2*Ne*mu*fix(d1,d2,s,indep,formula,model,a1,a2,Phi,q,Ne))
  }
  return(rates)
}

#simulation when only one site is evolving, with other sites fixed
#finish writing this function
simulation_1site <- function(protein,protein_op,t,m,change_site=1, s,indep=FALSE,formula=1,model=1,a1=2, a2=1, Phi=0.5,q=4*10^(-3),Ne=1.37*10^3,mu=1/(2*Ne),record.fitness=TRUE){
  #browser()
  l <- length(protein) #number of sites
  C <- a1 + a2*l #protein production cost
  t_now <- 0 #time until the current step of simulation
  path <- array(c(protein,0,0),dim=c(1,l+2)) #the array that stores the protein sequences
  colnames(path) <- colnames(path,do.NULL = FALSE, prefix = "Site.") #colomn names
  colnames(path)[l+1] <- "Time Now"
  colnames(path)[l+2] <- "Waiting Time"
  while(t_now < t){ #when current time is less than t
    rates <- rate_move_1site(protein,protein_op,m,change_site,s,indep,formula,model,a1,a2,Phi,q,Ne,mu) #moving rates of a protein to its neighbors
    lambda <- sum(rates) # rate of staying at the same state
    t_wait <- rexp(1,lambda) #waiting time, coming from exponential distribution with rate lambda
    t_now <- t_now + t_wait
    if(t_now < t){
      rates <- rates/lambda #probability vector of moving
      index <- mkv(runif(1),rates) #index of the protein it moves to
      protein[change_site] <- seq(1:m)[-protein[change_site]][index] #change the state of the site that is evolving
      path <- rbind(path,c(protein,t_now,t_wait)) #record this moving step
    }
    else{
      path <- rbind(path,c(protein,t,NA)) #record this moving step
    }
  }
  if(record.fitness){
    ftny_vec <- NULL
    fitness_vec <- NULL
    for(i in 1:dim(path)[1]){
      ftny_vec <- c(ftny_vec,Ftny(pchem_d(path[i,seq(1:l)],protein_op),s,model)) #functionality
    }
    fitness_vec <- exp(-q*Phi*C/ftny_vec) #fitness
    path <- cbind(path,ftny_vec)
    path <- cbind(path,fitness_vec)
    colnames(path)[l+3] <- "Functionality"
    colnames(path)[l+4] <- "Fitness"
  }
  path
}
#####################################################
#Given a protein, the optimal protein, and the selection vector (or constant selection)
#find the vector of rates of moving from the protein to all its neighbors
#Here we fix one particular site and let other sites evolve
#If the length of the protein is l, then it has (m-1)*(l-1) neighbors
rate_move_fixsite <- function(protein, protein_op, m, fix_site,s,indep=TRUE,formula=1,model=1,a1=2, a2=1, Phi=0.5,q=4*10^(-3),Ne=1.37*10^3,mu=1/(2*Ne)){
  l <- length(protein) #number of sites
  rates <- vector() # rates of moving to neighboring proteins, vector to return
  d1 <- pchem_d(protein,protein_op) #distance vector between given protein and optimal protein
  for(i in 1:l){ #loop over every position
    if(i!=fix_site){ #the fixed site does not change, we don't calculate fixation probabilities for these
    d2 <- d1 #set the distance vectors to be the same, later on, change the specific one with change
    aa <- seq(1:m)[-protein[i]] #All the m-1 aa's the site can change to
    for(j in 1:(m-1)){ #loop over the m-1 amino acids
        d2[i] <- GM[protein_op[i],aa[j]]
        rates <- c(rates,2*Ne*mu*fix(d1,d2,s,indep,formula,model,a1,a2,Phi,q,Ne))
    }
  }
#     else{
#       rates <- c(rates,rep(0,m-1))
#     }
  }
  rates
}


#Simulation. Given the starting protein and the optimal protein
#t: running time of the chain
simulation_fixsite <- function(protein,protein_op,t,m,fix_site,s,indep=FALSE,formula=1,model=1,a1=2, a2=1, Phi=0.5,q=4*10^(-3),Ne=1.37*10^3,mu=1/(2*Ne),record.fitness=TRUE){
  #browser()
  l <- length(protein) #number of sites
  C <- a1 + a2*l #protein production cost
  t_now <- 0 #time until the current step of simulation
  path <- array(c(protein,0,0),dim=c(1,l+2)) #the array that stores the protein sequences
  colnames(path) <- colnames(path,do.NULL = FALSE, prefix = "Site.") #colomn names
  colnames(path)[l+1] <- "Time Now"
  colnames(path)[l+2] <- "Waiting Time"
  while(t_now < t){ #when current time is less than t
    rates <- rate_move_fixsite(protein,protein_op,m,fix_site,s,indep,formula,model,a1,a2,Phi,q,Ne,mu) #moving rates of a protein to its neighbors
    lambda <- sum(rates) # rate of staying at the same state
    t_wait <- rexp(1,lambda) #waiting time, coming from exponential distribution with rate lambda
    t_now <- t_now + t_wait
    if(t_now < t){
      rates <- rates/lambda #probability vector of moving
      index <- mkv(runif(1),rates) #index of the protein it moves to
      pos_nofixsite <- seq(1:l)[-fix_site]
      pos <- (index-1) %/% (m-1) + 1 #index of the site that changes
      rmd <- (index-1) %% (m-1) + 1 # the amino acid that site changes to
      protein[pos_nofixsite[pos]] <- seq(1:m)[-protein[pos_nofixsite[pos]]][rmd] #new protein
      path <- rbind(path,c(protein,t_now,t_wait)) #record this moving step
    }
    else{
      path <- rbind(path,c(protein,t,NA)) #record this moving step
    }
  }
  if(record.fitness){
    ftny_vec <- NULL
    fitness_vec <- NULL
    for(i in 1:dim(path)[1]){
      ftny_vec <- c(ftny_vec,Ftny(pchem_d(path[i,seq(1:l)],protein_op),s,model)) #functionality
    }
    fitness_vec <- exp(-q*Phi*C/ftny_vec) #fitness
    path <- cbind(path,ftny_vec)
    path <- cbind(path,fitness_vec)
    colnames(path)[l+3] <- "Functionality"
    colnames(path)[l+4] <- "Fitness"
  }
  path
}

#####################################################
#simulate a big number(num) of chains up to time t, and record the last protein state in each chain.
#the frequencies of each protein in this set should serve as a good approximation to the stationary probability
simulation_large <- function(protein,protein_op,t,m,s,indep=FALSE,num,indep=TRUE,formula=1,model=1,a1=2, a2=1, Phi=0.5,q=4*10^(-3),Ne=1.37*10^3,mu=1/(2*Ne)){
  result <- array()
  for(i in 1:num){
    sim <- simulation(protein,protein_op,t,m,s,indep,formula=1,model=1,a1=2, a2=1, Phi=0.5,q=4*10^(-3),Ne=1.37*10^3,mu=1/(2*Ne),record.fitness=FALSE)
    result <- rbind(result,tail(sim,1)[c(1,2)])
  }
  return(result[-1,])
}
########################################################################
#Given an array whose rows are proteins, find the substitution rate matrix between them
#arr: the protein sequences; each row represents a protein 
#protein_op: optimal protein
#s: selection coefficients 
#########################################################################3
mat_gen_arr<- function(arr,protein_op,s,indep=FALSE,formula=1,model=1,a1=2, a2=1, Phi=0.5,q=4*10^(-3),Ne=1.37*10^3, mu=1/(2*Ne)){
  n1 <- dim(arr)[1] #number of proteins to consider
  #n <- dim(arr)[2] #number of sites in each protein
  mat <- matrix(0,nrow=n1,ncol=n1) #matrix to return, every entry is initiated to be 0
  for(i in 1:(n1-1)){
    for(j in (i+1):n1){
      mat[i,j] <- 2*Ne*mu*fix_protein(arr[i,],arr[j,],protein_op,s,indep,formula,model,a1,a2,Phi,q,Ne) #fixation prob -> transition rate
      mat[j,i] <- 2*Ne*mu*fix_protein(arr[j,],arr[i,],protein_op,s,indep,formula,model,a1,a2,Phi,q,Ne) #symmetric entry
    }
    mat[i,i] <- -sum(mat[i,])
  }
  mat[n1,n1] <- -sum(mat[n1,])
  return(mat)
}
#######################################################################################################


#######################################################################################################
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
#Convert a decimal integer to a binary number
binary <- function(x,k){
  #ndigits <- (floor(logb(x, base=2))+1)
  Base <- rep(0,k)
  for(i in 1:k){#i <- 1
    Base[k-i+1] <- (x %% 2)
  	x <- (x %/% 2)
	}
	return(Base)
}

integer.base.b <-function(x, b=2){
  xi <- as.integer(x)
  if(any(is.na(xi) | ((x-xi)!=0)))
  	print(list(ERROR="x not integer", x=x))
	N <- length(x)
	xMax <- max(x)	
	ndigits <- (floor(logb(xMax, base=b))+1)
	Base.b <- array(NA, dim=c(N, ndigits))
	for(i in 1:ndigits){#i <- 1
		Base.b[, ndigits-i+1] <- (x %% b)
		x <- (x %/% b)
	}
	if(N ==1) Base.b[1, ] else Base.b
}
#integer.base.b(seq(0,4^3),b=4) + 1
########################################################
#Given two sequences, find all the sequences on
#the possible paths
path <- function(protein1,protein2){
  l <- length(protein1)
  U <- seq(l) # (1,2,3,4,...,l)
  diff <- protein1 - protein2
  S <- U[diff==0] #set of same sites
  D <- U[diff!=0] #set of different sites
  d <- length(D) #number of different sites
  m <- 2^(length(D))#number of sequences
  paths <- array(NA,dim=c(m,l))
  for(i in 0:(m-1)){
    pos <- binary(i,d)
    paths[i+1,] <- protein1
    paths[i+1,D[pos==1]] <- protein2[D[pos==1]]
  }
  return(paths)
}
# p1 <- c(4,3,6,8,3)
# p2 <- c(2,3,5,8,4)
# path(p1,p2)


################################################################################

#Function to choose a state given the probabilities of changing to 
#all states (probability vector)
#CONSIDER USING "SAMPLE" INSTEAD OF THIS FUNCTION
mkv <- function(odds,vec){
  vec <- vec/sum(vec)
  vec_order <- order(vec,decreasing=T)
  vec <- vec[vec_order]
  if(odds<vec[1]){
    pick <- 1
  }
  else{
    k <- 1
    while(odds > sum(vec[1:k])){
      k <- k + 1
      if(odds < sum(vec[1:k])){
        pick <- k
      }
    }
  }
  return(vec_order[pick])
}
########################################################################

###function to find the average time of m transitions in a simulation###
avg_time <- function(sim,m){
  time <- sim[,"Time Now"] #extract the times when the transitions happen
  #time <- as.numeric(sim[,2])
  l <- length(time)
  time <- head(time,l-1) #in the last step, there is no transition
  time_start <- head(time, l-1-m)
  time_end <- tail(time, l-1-m)
  return(mean(time_end-time_start))
}

##Extract the states every m transitions###
sim_extract <- function(sim,m){
  time_trans <- avg_time(sim,m) #average time of m transitions
  time <- sim[,"Time Now"] #times of transition happening
  #time <- as.numeric(sim[,2])
  time_mod <- time %/% time_trans #
  pick_time_vec <- time_mod[-length(time_mod)]-time_mod[-1]
  return(sim[pick_time_vec < 0,])
}

#Put the sites together as one column so that we can order them and get the levels of them
parse_array <- function(Array,sites){
  ParsedArray <- array(dim=c(dim(Array)[1],dim(Array)[2]-sites+1)) #the array to return
  l <- dim(ParsedArray)[2] #number of records
  for(i in 1:dim(Array)[1]){
    ParsedArray[i,1] <- paste(Array[i,seq(1,sites,by=1)],collapse="_") #collapse the sites to one column
  }
  for(j in 2:l){
    ParsedArray[,j] <- as.numeric(Array[,sites+j-1])
  }
  ParsedArray
}



sim1 <- sim_extract(sim_ind,10)
sim1_state <- sim1[,1]
vec <- as.numeric(table(sim1_state))
table(sim1_state)
freq_sim <- vec/sum(vec)

q <- 4e-3
Ne <- 1.37e3
Phi <- 0.5
C <- 3
ftny_vec <- NULL
for(i in 1:m){
  distance <- pchem_d(i,2)
  ftny_vec <- c(ftny_vec,Ftny(distance,0.1))
}
fit_vec <- exp(-q*Phi*C/ftny_vec)
freq <- (fit_vec/max(fit_vec))^(2*Ne-1) #this is a better way to evaluate S-H formula
freq_SH <- freq/sum(freq)

find_freq <- function(sim, sites,transitions,burnin,Ne=1.37e3){
  #browser()
  sim <- sim[-c(1:burnin),]
  sim_ex <- sim_extract(sim,transitions)
  sim_parsed <- parse_array(sim_ex,sites)
  data <- sim_parsed[,1]
  data <- data[order(data)]
  data_unique <- unique(data)
  data_level <- table(data)
  counts <- as.numeric(data_level)
  freq_sim <- counts/sum(counts)
  ##fit_vec <- as.numeric(data_unique[,2])
  ##freq <- (fit_vec/max(fit_vec))^(2*Ne-1) #this is a better way to evaluate S-H formula
  ##freq_SH <- freq/sum(freq)
  ##freq <- list(freq_sim=freq_sim, freq_SH=freq_SH)
  return(freq_sim)
}

#Expected value of the fitness ratio between 2 states at a certain sites
#This is a simple case where there are only 3 sites and 3 states for each site, which means, there are 27 proteins in total
#protein_op is the optimal AA sequence
#freq_sta is the stationary probability vector for all 27 AA sequences
exp_fratio <- function(state1, state2,site,protein_op,s,freq_sta){
  exp_fit_ratio <- 0
  for(i in 1:3){
    for(j in 1:3){
      if(site==1){
        index <- (i-1)*3 + j
        freq <- sum(freq_sta[index],freq_sta[index+9],freq_sta[index+18])
        exp_fit_ratio <- exp_fit_ratio + (Ftny_protein(c(state1,i,j),protein_op,s)/Ftny_protein(c(state2,i,j),protein_op,s))*freq
      }
      else if(site==2){
        index <- index <- (i-1)*9 + j
        freq <- sum(freq_sta[index],freq_sta[index+3],freq_sta[index+6])
        exp_fit_ratio <- exp_fit_ratio + (Ftny_protein(c(i,state1,j),protein_op,s)/Ftny_protein(c(i,state2,j),protein_op,s))*freq
      }
      else if(site==3){
        index <- index <- (i-1)*9 + (j-1)*3 + 1
        freq <- sum(freq_sta[index],freq_sta[index+1],freq_sta[index+2])
        exp_fit_ratio <- exp_fit_ratio + (Ftny_protein(c(i,j,state1),protein_op,s)/Ftny_protein(c(i,j,state2),protein_op,s))*freq
      }
    }
  }
  exp_fit_ratio
}
sim0.1_parsed <- parse_array(sim0.1,2)
data <- simOut0.1_parsed[,c(1,5)]
data <- data[order(data[,1]),]
data_unique <- unique(data)
data_level <- table(data[,1])
counts <- as.numeric(data_level)
freq <- counts/sum(counts)
fit_vec <- as.numeric(data_unique[,2])

ssim <- sim_simple(4,W,10^7)
ssim_ex <- sim_extract(ssim,5)
ssim_state <- ssim_ex[,1]
vec <- as.numeric(table(ssim_state))
freq <- vec/sum(vec)

ftny <- vector("numeric",length=20)
for(i in 1:20){
  ftny[i] <- Ftny_protein(i,1,0.1,GM)
}

sites <- vector("numeric",20)
for(i in 1:20){
  sites[i] <- length(which(A2[,i]!=0))-1
}

i=3
file <- paste("rokas/valley",i,".RData",sep="")
load(file)
op <- as.vector(apply(fungi_data[1:8,],2,Mode))
op[op!=i]

l <- 100
Root <- sample(20,l,replace=T)
Protein_op <- sample(20,l,replace=T)
system.time(fungi_data <- simTree(tree,l,protein_op=Protein_op,m=20,s=0.1,
                                  Nu_vec,scale,al,be,ga,q=4e-7,Ne=1.36e7,root=Root,ancestral=TRUE))

gene1tree <- read.nexus(file="rokas/gene1.tre")
########################################################################
## convert between file formats for data
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
########################################################################
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