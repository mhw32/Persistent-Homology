# A test script to play around with Voronoi simulations and measuring hypothesis tests between them.
source('NHST.r')
source('twosample.r')
source("voronoi3dfct.r")
source('euler.r')
source('distance')

# First, create two sets of voronoi diagrams.
# Simulation settings. Small diagrams for computational purposes.
Boxlim <- c(0,10)
Xlim <- Boxlim
Ylim <- Boxlim
Zlim <- Boxlim
resolution <- 0.5  ## grid space for approximating the voronoi cells.
perturb <- 1 ## variance around the filaments
N <- 1000   ## number of particles
groupN <- 15 ## each "set" has 15 diagrams.

# Generate group 1.
set1 <- sapply(seq(1:groupN), function(i) {
  percFil1 <- 0.9
  vf1 <- voronoi3d(Boxlim, resolution, perturb, Ncells=64, N, percClutter=0, percWall=1-0.02-percFil, percFil=percFil1, percClust=0.02)
  diag1 <- gridDiag(vf1, dtm, lim=cbind(Xlim,Ylim,Zlim), by=resolution, sublevel=T, printProgress=T, m0=0.001)
  return(diag1$diagram)
})

# Generate group 2.
set2 <- sapply(seq(1:groupN), function(i) {
  percFil2 <- 0.1
  vf2 <- voronoi3d(Boxlim, resolution, perturb, Ncells=64, N, percClutter=0, percWall=1-0.02-percFil, percFil=percFil2, percClust=0.02)
  diag2 <- gridDiag(vf2, dtm, lim=cbind(Xlim,Ylim,Zlim), by=resolution, sublevel=T, printProgress=T, m0=0.001)
  return(diag2$diagram)
})

# The datasets are fully generated. They may be too sparse and the box limits may also be too restrictive but I think they should be operational.

# 1. NHST test.
X <- cbind(set1, set2)
L <- cbind(rep(0, groupN), rep(1, groupN))
# Should check that nhst is working as expected.
distfunc <- bottleneckDist
# distfunc <- wassersteinDist
proba <- nhst(X, L, 1000, bottleneckDist)

# 2. Gaussian kernel (max likelihood h).
# 3. Energy kernel.
# 4. Rossenbaum kernel.
# 5, Euler's constant. 
# 6. Wu's Custom kernels.