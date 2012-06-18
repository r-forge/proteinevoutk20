worker.init <- function(packages){
  for(p in packages){
    library(p,character.only=TRUE) #interpret the argument as a character variable so that library() won't attempt to load a package named p repeatedly.
  }
  NULL #don't accidentally send unnecessary data transfer back to the master
}
clusterCall(cl,worker.init(c("MASS","boot")))
clusterCall(cl,function(x) Sys.info()[c("nodename","machine")])

## clusterApplyLB() lets the workers pull tasks as needed, not in a "one at a time" fashion. This can be more efficient if some tasks take longer than others, or if some cluster workers are slower.
set.seed(7777442)
sleeptime <- abs(rnorm(10,10,10))
tm <- snow.time(clusterApplyLB(cl,sleeptime,Sys.sleep))
plot(tm)

## parLapply() splits up x into a list of subvectors, and process those subvectors on the cluster workers using lapply(). It is prescheduling the work by dividing the tasks into as many chunks as there are workers in the cluster. clusterApply() gives you more control over what gets sent to who, while parLapply() provides a convenient way to efficiently divide the work among the cluster workers. parLapply() is much more efficient than clusterApply() if you have many more tasks than workers, and one or more large, additional arguments to pass to parLapply().
clusterSplit(cl, 1:30)
parVapply <- function(cl, x, fun, ...){
  do.call("c",clusterApply(cl, clusterSplit(cl,x),fun, ...))
}
parVapply(cl,1:10,"^",1/3) #cubic root of numbers from 1 to 10

##Load balancing Redux
##First send the function and the fixed arguments to the cluster workers using clusterCall(), which saves them in the global environment, and then send the varying argument values using clusterApplyLB(), specifying a function that will execute the user-supplied funciton along with the full collection of arguments.
parLapplyLB <- function(cl,x,fun, ...){
  clusterCall(cl, LB.init, fun, ...)
  r <- clusterApplyLB(cl, x, LB.worker)
  clusterEvalQ(cl, rm('LB.fun','.LB.args',pos=globalenv()))
  r
}
LB.init <- function(fun, ...){
  assign('.LB.fun',fun, pos=globalenv())
  assign('.LB.args',list(...),pos=globalenv)
  NULL
}
LB.worker <- function(x){
  do.call('.LB.fun',c(list(x), .LB.args))
}

##compare clusterApplyLB() and parLapplyLB()
bigsleep <- function(sleeptime,mat) Sys.sleep(sleeptime)
bigmatrix <- matrix(0,2000,2000)
sleeptime <- rep(1,100)
tm <- snow.time(clusterApplyLB(cl, sleeptime, bigsleep, bigmatrix))
plot(tm)

tm <- snow.time(parLapplyLB(cl, sleeptime, bigsleep, bigmatrix))
plot(tm)

serialize()
unserialize()
##Name space environments are serialized by name, not by value.

##Random Number Generation. snow provides support for the 'rlecuyer' and 'rsprng' packages

##use rlecuyer:
clusterSetupRNG(cl, type="RNGstream")
clusterSetupRNG(cl, type="RNGstream",seed=c(1,22,333,444,55,6)) #have to specify the random seed explicitly using the seed argument
##use rsprng:
clusterSetupRNG(cl, type="SPRNG")

##We can get reproducible results using clusterApply(), but not with clusterApplyLB() because clusterApply() always uses the same task scheduling, while clusterApplyLB() does not. (use the same random seeds to reproduce results)

setDefaultClusterOptions() #change a default configuration option during an R session

##Rmpi package
##mpi.R
library(snow)
library(Rmpi)
cl <- makeMPICluster(mpi.universe.size()-1)
r <- clusterEvalQ(cl, R.version.string)
print(unlist(r))
stopCluster(cl)
mpi.quit()
##command to run this script
% orterun -H localhost, n1,n2,n3,n4 -n 1 R --slave -f mpi.R
% orterun -wdir /tmp -H localhost, n1, n2,n3, n4 -n 1 R --slave -f ~/mpi.R

##Batch execution
## batchmpi.sh
#!/bin/sh
#PBS -N SNOWMPI
#PBS -j oe
cd $PBS_O_WORKDIR
orterun -n 1 /usr/bin/R --slave -f mpi.R > mpi-$PBS_JOBID.out 2>&1

% qsub -q devel -l nodes=2:ppn=4 batchmpi.sh

## manual mode
cl <- makeCluster(2,type="SOCK",manual=TRUE,outfile="")

## use ping to see if computers can be reached and commmunicate

##############################################################
##   multicore package ##
##############################################################
mcapply()
pvec()
parallel()
collect()

##2 workers instead of 2 cores despite what it looks like
unique(unlist(mclapply(1:100,function(i) Sys.getpid(),mc.cores=2)))
#or
options(cores=3)
unique(unlist(mclapply(1:100,function(i) Sys.getpid())))
##seed each of the workers to a different value after they have been created.
mclapply(1:3,function(i) rnorm(3),mc.cores=3,mc.set.seed=TRUE)

##Load Balancing with mclapply - set mc.preschedule to FALSE
set.seed(93564990)
sleeptime <- abs(rnorm(10,10,10))
system.time(mclapply(sleeptime,Sys.sleep,mc.cores=4))
system.time(mclapply(sleeptime,Sys.sleep,mc.cores=4,mc.preschedule=FALSE))

##pvec() : high level function used to execute vector functions in parallel
x <- 1:10
pvec(x,"^",1/3)

library(multicore)
fun1 <- function() {Sys.sleep(10); 1}
fun2 <- function() {Sys.sleep(5); 2}
fun3 <- function() {Sys.sleep(1); 3}

f1 <- parallel(fun1())
f2 <- parallel(fun2())
f3 <- parallel(fun3())
collect(list(f1,f2,f3))

collect(list(f1,f2,f3), wait=FALSE)
##parallel() createds a new process using fork() to evaluate an expression in parallel with the calling process. It returns a parallelJob object which is passed to the collect() function to retrieve the result of the computation. collect() can be called with either a single parallelJob object, or a list of parallelJob objects. One can think of parallel() as a submit operation, and collect() as a wait operation, similar to batch queueing commands. 

## Use initSprngNode() in snow to initialize each of the workers to use parallel random numbers at the start of the task.
library(snow)
nw <- 3
seed <- 7777442
kind <- 0
para <- 0
f1 <- parallel({
  initSprngNode(0,nw.seed,kind,para)
  rnorm(1)
})
f2 <- parallel({
  initSprngNode(1,nw.seed,kind,para)
  rnorm(1)
})
f3 <- parallel({
  initSprngNode(2,nw.seed,kind,para)
  rnorm(1)
})
unlist(collect(list(f1,f2,f3)),use.names=FALSE)

cl <- makeCluster(3,type="SOCK")
seed <- 7777442
clusterSetupSPNG(cl, seed=seed)
unlist(clusterEvalQ(cl, rnorm(1)), use.names=FALSE)
stopCluster(cl)


##############################################################
##   parallel  package ##
##############################################################
detectCores() #detect number or cores in the machine
options(mc.cores=detectCores())

## FOCK vs. PSOCK clusters
type <- if(exists("mcfork",mode="function")) "FORK" else "PSOCK"

RNGkind("L'Ecuyer-CMRG")
clusterSetRNGStream(cl,7777442)
.Random.seed <<- nextRNGSubStream(.Random.seed)

##############################################################
##   MapReduce and Hadoop  ##
##############################################################
forqlift http://www.forqlift.net/
