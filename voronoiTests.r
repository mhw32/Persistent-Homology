# A test script to play around with Voronoi simulations and measuring hypothesis tests between them.
source('NHST.r')
source('twosample.r')
source("voronoi3dfct.r")
source('euler.r')
source('distance.r')

# First, create two sets of voronoi diagrams.
# Simulation settings. Small diagrams for computational purposes.
Boxlim <- c(0,10)
Xlim <- Boxlim
Ylim <- Boxlim
Zlim <- Boxlim
resolution <- 0.5  ## grid space for approximating the voronoi cells.
perturb <- 1 ## variance around the filaments
N <- 1000   ## number of particles
groupN <- 5 ## each "set" has 15 diagrams.

# Generate group 1.
set1 <- sapply(seq(1:groupN), function(i) {
  percFil1 <- 0.9
  vf1 <- voronoi3d(Boxlim, resolution, perturb, Ncells=64, N, percClutter=0, percWall=1-0.02-percFil, percFil=percFil1, percClust=0.02)
  diag1 <- gridDiag(vf1, dtm, lim=cbind(Xlim,Ylim,Zlim), by=resolution, sublevel=T, printProgress=T, m0=0.001)
  return(diag1)
})

# Generate group 2.
set2 <- sapply(seq(1:groupN), function(i) {
  percFil2 <- 0.1
  vf2 <- voronoi3d(Boxlim, resolution, perturb, Ncells=64, N, percClutter=0, percWall=1-0.02-percFil, percFil=percFil2, percClust=0.02)
  diag2 <- gridDiag(vf2, dtm, lim=cbind(Xlim,Ylim,Zlim), by=resolution, sublevel=T, printProgress=T, m0=0.001)
  return(diag2)
})

# The datasets are fully generated. They may be too sparse and the box limits may also be too restrictive but I think they should be operational.

X <- cbind(set1, set2)
L <- cbind(rep(0, groupN), rep(1, groupN))
permN <- 100

# ---------------------------------------------------------------------

# 1. NHST test.
# Should check that nhst is working as expected.
distfunc <- bottleneckDist
# distfunc <- wassersteinDist
proba <- nhst(X, L, permN, bottleneckDist)
# This works but it is so slow. 
# system.time(bottleneck(X[[1,1]], X[[2,1]], dimension=1))
# 1.335 * 225 * 2 * 1000 (group size 15 --> 15^2, 2 groups, 1000 iterations.)
# Even 100 iterations = 16 hours...
 
# system.time(bottleneck(X[[1,1]], X[[2,1]], dimension=2))
# user  system elapsed
# 0.212  0.024  0.237

# Is the 2-wasserstein faster (dimension 1)? This is so slow.
# Crashes. 160.656 seconds.

# 2-wasserstein faster (dimension 2)
# user  system elapsed
# 5.760   0.014   5.796

# No this is much slower. How do I make these raw functions faster? Otherwise, tests for as big as the Universe will take decades. I don't think this is a viable path.

# ---------------------------------------------------------------------
# In fact, tests 2-4 all require some form of a distance measurement. This is going to be very problematic since everything takes such a long time. Perhaps we must move to landscapes directly. How did you run these? 

# 2. Gaussian kernel (max likelihood h).
distfunc <- bottleneckDist
kernel <- gaussianKernel
# If we want to find the best value for parameter h:
best <- findBestParam(X, L)
proba <- permutationTest(permN, X, L, kernelStat, distfunc, best)

# ---------------------------------------------------------------------

# 3. Energy kernel.
distfunc <- bottleneckDist
proba <- permutationTest(permN, X, L, energyStat, distfunc, 1)

# ---------------------------------------------------------------------
 
# 4. Rossenbaum kernel.
# No need for permutation tests for still takes the annoying distance func.
distfunc <- bottleneckDist
proba <- rosenbaumStat(X, L, distfunc)

# ---------------------------------------------------------------------

# 5. Euler's constant. 
proba <- permutationTest(permN, X, L, eulerStat)

