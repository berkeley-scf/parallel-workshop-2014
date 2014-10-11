## @knitr RlinAlg

require(RhpcBLASctl)
x <- matrix(rnorm(5000^2), 5000)

omp_set_num_threads(4)
system.time({
x <- crossprod(x)
U <- chol(x)
})

omp_set_num_threads(1)
system.time({
x <- crossprod(x)
U <- chol(x)
})

## @knitr foreach

require(parallel) # one of the core R packages
require(doParallel)
# require(multicore); require(doMC) # alternative to parallel/doParallel
# require(Rmpi); require(doMPI) # to use Rmpi as the back-end
library(foreach)
library(iterators)

taskFun <- function(){
	mn <- mean(rnorm(10000000))
	return(mn)
}
nCores <- 8  # manually for non-cluster machines
#nCores <- as.numeric(Sys.getenv('NSLOTS')) # for use on cluster
registerDoParallel(nCores) 
# registerDoMC(nCores) # alternative to registerDoParallel
# cl <- startMPIcluster(nCores); registerDoMPI(cl) # when using Rmpi as the back-end

out <- foreach(i = 1:100) %dopar% {
	cat('Starting ', i, 'th job.\n', sep = '')
	outSub <- taskFun()
	cat('Finishing ', i, 'th job.\n', sep = '')
	outSub # this will become part of the out object
}

## @knitr parallelApply

require(parallel)
nCores <- 8  # manually for non-cluster machines
# nCores <- as.numeric(Sys.getenv('NSLOTS')) # for use on cluster

#############################
# using sockets
#############################

# ?clusterApply
cl <- makeCluster(nCores) # by default this uses the PSOCK 
#  mechanism as in the SNOW package - starting new jobs via Rscript 
#  and communicating via sockets
nSims <- 60
input <- seq_len(nSims) # same as 1:nSims but more robust
testFun <- function(i){
	mn <- mean(rnorm(1000000))
	return(mn)
}
# clusterExport(cl, c('x', 'y')) # if the processes need objects 
#   (x and y, here) from the master's workspace
system.time(
	res <- parSapply(cl, input, testFun)
)
system.time(
	res2 <- sapply(input, testFun)
)
res <- parLapply(cl, input, testFun)

################################
# using forking
################################

system.time(
	res <- mclapply(input, testFun, mc.cores = nCores) 
)


## @knitr pvec

require(parallel)
nCores <- 8
# nCores <- as.numeric(Sys.getenv('NSLOTS')) # for use on cluster
cl <- makeCluster(nCores) 
library(fields)
ds <- runif(6000000, .1, 10)
ds_exp <- pvec(ds, exp, mc.cores = nCores)
# here's a more computationally intensive function
system.time(
	corVals <- pvec(ds, Matern, .1, 2, mc.cores = nCores)
)
system.time(
	corVals <- Matern(ds, .1, 2)
)

## @knitr mcparallel

library(parallel)
n <- 10000000
system.time({
	p <- mcparallel(mean(rnorm(n)))
	q <- mcparallel(mean(rgamma(n, shape = 1)))
	res <- mccollect(list(p,q))
})
system.time({
	p <- mean(rnorm(n))
	q <- mean(rgamma(n, shape = 1))
})

## @knitr RNG-apply

require(parallel)
require(rlecuyer)
nSims <- 250
testFun <- function(i){
	val <- runif(1)
	return(val)
}

nSlots <- 4
RNGkind()
cl <- makeCluster(nSlots)
iseed <- 0
# ?clusterSetRNGStream
clusterSetRNGStream(cl = cl, iseed = iseed)
RNGkind() # clusterSetRNGStream sets RNGkind as L'Ecuyer-CMRG
# but it doesn't show up here on the master
res <- parSapply(cl, 1:nSims, testFun)
clusterSetRNGStream(cl = cl, iseed = iseed)
res2 <- parSapply(cl, 1:nSims, testFun)
identical(res,res2)
stopCluster(cl)

## @knitr RNGstream

RNGkind("L'Ecuyer-CMRG") 
seed <- 0
set.seed(seed) ## now start M workers 
s <- .Random.seed 
for (i in 1:M) { 
	s <- nextRNGStream(s) 
	# send s to worker i as .Random.seed 
} 

## @knitr RNG-mclapply

require(parallel)
require(rlecuyer)
RNGkind("L'Ecuyer-CMRG")
res <- mclapply(seq_len(nSims), testFun, mc.cores = nSlots, 
    mc.set.seed = TRUE) 
# this also seems to reset the seed when it is run
res2 <- mclapply(seq_len(nSims), testFun, mc.cores = nSlots, 
    mc.set.seed = TRUE) 
identical(res,res2)

## @knitr RNG-doMPI

nslaves <- 4
library(doMPI, quietly = TRUE)
cl <- startMPIcluster(nslaves)
registerDoMPI(cl) 
result <- foreach(i = 1:20, .options.mpi = list(seed = 0)) %dopar% { 
	out <- mean(rnorm(1000)) 
}
result2 <- foreach(i = 1:20, .options.mpi = list(seed = 0)) %dopar% { 
	out <- mean(rnorm(1000)) 
}
identical(result, result2)

## @knitr RNG-doRNG

rm(result, result2)
nCores <- 4
library(doRNG, quietly = TRUE)
library(doParallel)
registerDoParallel(nCores) 
result <- foreach(i = 1:20, .options.RNG = 0) %dorng% { 
	out <- mean(rnorm(1000)) 
}
result2 <- foreach(i = 1:20, .options.RNG = 0) %dorng% { 
	out <- mean(rnorm(1000)) 
}
identical(result, result2)

## @knitr RNG-doRNG2

rm(result, result2)
library(doRNG, quietly = TRUE)
library(doParallel)
registerDoParallel(nCores)
registerDoRNG(seed = 0) 
result <- foreach(i = 1:20) %dopar% { 
	out <- mean(rnorm(1000)) 
}
registerDoRNG(seed = 0) 
result2 <- foreach(i = 1:20) %dopar% { 
	out <- mean(rnorm(1000)) 
}
identical(result,result2)



## @knitr fork

library(fork)
{ # this set of braces is REQUIRED, unless you pass a function 
  # to the slave argument of fork()
	pid <- fork(slave = NULL) 
	if(pid==0) {
		cat("Starting child process execution.\n") 
		tmpChild <- mean(rnorm(10000000))
		cat("Result is ", tmpChild, "\n", sep = "")
		save(tmpChild, file = 'child.RData')
		cat("Finishing child process execution.\n")
		exit() 
	} else {
		cat("Starting parent process execution.\n")
		tmpParent <- mean(rnorm(10000000))
		cat("Finishing parent process execution.\n")
		wait(pid)  # wait til child is finished so can read in 
                   # updated child.RData below
	} 
} 
load('child.RData')
print(c(tmpParent, tmpChild))


## @knitr Rmpi-foreach

# start R as usual: "R" or via a batch job
library(Rmpi)
library(doMPI)
nslaves = 4
cl = startMPIcluster(nslaves)
registerDoMPI(cl)
clusterSize(cl) # just to check

result <- foreach(i = 1:20) %dopar% {
  out = mean(rnorm(1e7))
}

closeCluster(cl)

# you can also start as
#"mpirun -np 1 R --no-save"
#"mpirun -np 1 R CMD BATCH --no-save example.R example.out"

# if you do this, you should quit via:
#mpi.quit()

## @knitr Rmpi-usingMPIsyntax

# example syntax of standard MPI functions

library(Rmpi)
mpi.spawn.Rslaves(nslaves = 4)

n = 5
mpi.bcast.Robj2slave(n)
mpi.bcast.cmd(id <- mpi.comm.rank())
mpi.bcast.cmd(x <- rnorm(id))

mpi.remote.exec(ls(.GlobalEnv), ret = TRUE)

mpi.bcast.cmd(y <- 2 * x)
mpi.remote.exec(print(y))

objs <- c('y', 'z')
# next command sends value of objs on _master_ as argument to rm
mpi.remote.exec(rm, objs)  
mpi.remote.exec(print(z))

# collect results back via send/recv
mpi.remote.exec(mpi.send.Robj(x, dest = 0, tag = 1))
results = list()
for(i in 1:(mpi.comm.size()-1)){
  results[[i]] = mpi.recv.Robj(source = i, tag = 1)
}
  
print(results)

## @knitr sockets

# example with PSOCK cluster

library(parallel)

machine = "arwen"
cl = makeCluster(rep(machine, 4))
# cl = makeCluster(4) # this would do it by default on your local machine
n = 1e7
clusterExport(cl, c('n'))
fun = function(i)
  out = mean(rnorm(n))
  
result <- parSapply(cl, 1:20, fun)

stopCluster(cl) # not strictly necessary

## @knitr Rmpi-foreach-multipleNodes

# start R as either:
# mpirun -hostfile .hosts -np 1 R --no-save
# mpirun -hostfile .hosts -np 1 R CMD BATCH --no-save

library(Rmpi)
library(doMPI)
nslaves = 2
# 2 slaves since my host file specifies 3 slots,
# and one will be used for master
cl = startMPIcluster(nslaves)

registerDoMPI(cl)
clusterSize(cl) # just to check

results <- foreach(i = 1:200) %dopar% {
  out = mean(rnorm(1e7))
}

if(FALSE){
  foreach(i = 1:20) %dopar% {
    out = chol(crossprod(matrix(rnorm(3000^2),3000)))[1,2]
  }
}

closeCluster(cl)

mpi.quit()

## @knitr sockets-multipleNodes

# multinode example with PSOCK cluster

library(parallel)

machineVec = c(rep("arwen", 4), rep("treebeard", 2), rep("beren", 2))
cl = makeCluster(machineVec)

n = 1e7
clusterExport(cl, c('n'))
fun = function(i)
  out = mean(rnorm(n))
  
result <- parSapply(cl, 1:20, fun)

stopCluster(cl) # not strictly necessary

