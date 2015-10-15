# ***************************************************************************
# Title  : Voronoi Hypothesis Testing 
# Author : Mike Wu
# Description : A test script to play around with Voronoi simulations and measuring hypothesis tests between them. In designing, there is a preference for speed. Each of the datasets contains a sample from some high dimensional estimation of the Universe. Because this is a random variable, maximum matching doesn't really make sense since the exact locations are variable. So algorithms using distance metrics are most likely meaningless. Instead, the focus here is on distributional and summary algorithms.
# ***************************************************************************

source('VoronoiFoam.r')
source('euler.r')
source('summarize.r')
source('distribution.r')
source('tools.r')
library(abind)
library(Hotelling)
library(ks)

# voronoi_compilation() 
foam <- readRDS('./voronoifoam.rds')

# Everything will be compared to baseline of 0.1
setnum <- length(foam)
colnum <- length(foam[[1]])

# =====================================================================
# Euler Characteristic [Baseline]
# =====================================================================
eulerMat <- gridOperation(foam, eulerIntegration)
# Do a T-test on the two sets of AUC.
# Set 0.1 as the baseline and calculate on top of it.
eulerProba <- rep(0, setnum)
for (i in 1:setnum) {
  currproba <- t.test(eulerMat[,1], eulerMat[,i])
  eulerProba[i] <- currproba$p.value
}
# Plot the probabilities.
plot(seq(1:setnum), proba)
lines(seq(1:setnum), proba, type="b")

# [1] 1.00000000 0.47366485 0.18146536 0.05285394 0.20642641 0.45542569 
# [7] 0.51701674 0.44226854 0.90792663

# =====================================================================
# Try landscapes: Individual dimension at a time.
# =====================================================================
landDimProba <- matrix(NA, nrow=setnum, ncol=3)
for (i in 0:2) {
  landDimFxn <- dimWrapper(i, landscapeAUC)
  landDimMat <- gridOperation(foam, landDimFxn)
  # Calculate probabilities with t-test
  for (j in 1:setnum) {
    currproba <- t.test(landDimMat[,1], landDimMat[,j])
    landDimProba[j, i+1] <- currproba$p.value
  }
}

 #              [,1]      [,2]      [,3]
 # [1,] 1.0000000000 1.0000000 1.0000000
 # [2,] 0.9496798198 0.2026717 0.4859109
 # [3,] 0.5152820473 0.9918285 0.2281313
 # [4,] 0.4357253549 0.1242701 0.2797075
 # [5,] 0.6706838754 0.5950084 0.2499642
 # [6,] 0.1485794949 0.3256215 0.4520231
 # [7,] 0.0794191314 0.2294548 0.4101071
 # [8,] 0.0239369177 0.6921783 0.2995093
 # [9,] 0.0002065439 0.1358652 0.5459708

# =====================================================================
# Try landscapes: Combining all 3 dimensions (multi-D test).
# =====================================================================
# landMat <- gridOperation(foam, landscapeAUC)
landfxn <- allWrapper(landscapeAUC)
landMat <- gridOperation(foam, landfxn)
# Calculate probabilities through multi-D t-test
landProba <- rep(0, setnum)
for (i in 1:setnum) {
  currproba <- hotelling.test(t(landMat[,,1]), t(landMat[,,i]))
  landProba[i] <- currproba$pval
}

# [1] 1.00000000 0.58836091 0.48818559 0.88391109 0.85846633 0.23193217 
# [7] 0.22234285 0.52068225 0.07237922

# =====================================================================
# Try Silhouettes: Individual dimension at a time.
# =====================================================================
silhDimProba <- matrix(NA, nrow=setnum, ncol=3)
for (i in 0:2) {
  silhDimFxn <- dimWrapper(i, silhouetteAUC)
  silhDimMat <- gridOperation(foam, silhDimFxn)
  # Calculate probabilities with t-test
  for (j in 1:setnum) {
    currproba <- t.test(silhDimMat[,1], silhDimMat[,j])
    silhDimProba[j, i+1] <- currproba$p.value
  }
}

 #              [,1]       [,2]      [,3]
 # [1,] 1.0000000000 1.00000000 1.0000000
 # [2,] 0.9633687841 0.55760612 0.5388922
 # [3,] 0.4171907642 0.63043277 0.3856845
 # [4,] 0.2246211250 0.31691032 0.7358537
 # [5,] 0.6896973898 0.13374005 0.2879359
 # [6,] 0.0765590148 0.27137786 0.5141326
 # [7,] 0.0581606799 0.43319966 0.2095858
 # [8,] 0.0136696725 0.04339062 0.4384976
 # [9,] 0.0003270311 0.39483547 0.9780893

# =====================================================================
# Try silhouttes: Combining all 3 dimensions (multi-D test)
# Silhouttes are so freaking fast!
# =====================================================================
# silhMat <- gridOperation(foam, silhouetteAUC)
silhfxn <- allWrapper(silhouetteAUC)
silhMat <- gridOperation(foam, silhfxn)
# Calculate probabilities through multi-D t-test
silhProba <- rep(0, setnum)
for (i in 1:setnum) {
  currproba <- hotelling.test(t(silhMat[,,1]), t(silhMat[,,i]))
  silhProba[i] <- currproba$pval
}

# [1] 1.000000000 0.890041611 0.730164341 0.387986322 0.302669180 0.137406704
# [7] 0.117997710 0.011137029 0.003031294

# =====================================================================
# Kernel density based local two-sample comparison test 
# http://www.mvstat.net/tduong/research/publications/duong-2013-jns.pdf
# =====================================================================
# Ask Jessi Cisewski. Not sure how to handle grid-based p-values.

# =====================================================================
# Distribution Testing (Separate): Given diagrams, split it into 0, 1, 2 dimensions and compare their distributions via nonparametric test.
# This is by far the fastest method compared.
# =====================================================================

# Loop through dimensions.
distrDimProba <- matrix(NA, nrow=setnum, ncol=3)
for (d in 0:2) {
  # Loop through the set and do the distribDimStat for each.
  distrDimList <- vector("list", setnum)
  for (i in 1:setnum) { distrDimList[[i]] <- distribDimStat(foam[[i]], d) }
  # Do a wilcox test for each with the baseline being the 0.1.
  baseline <- distrDimList[[1]] 
  for (i in 1:setnum) {
    currproba <- wilcox.test(baseline, distrDimList[[i]])
    distrDimProba[i, d+1] <- currproba$p.value
  }
}

 #              [,1]         [,2]         [,3]
 # [1,] 1.000000e+00 1.000000e+00 1.000000e+00
 # [2,] 2.672220e-01 3.999023e-01 9.388415e-01
 # [3,] 8.332263e-01 3.902581e-09 8.312350e-03
 # [4,] 1.016991e-02 9.380618e-13 1.880157e-06
 # [5,] 8.363584e-03 1.034085e-15 8.307920e-09
 # [6,] 6.978095e-05 4.624225e-23 2.154619e-13
 # [7,] 4.408826e-08 7.069825e-34 1.605872e-16
 # [8,] 4.942691e-11 8.549844e-39 1.015916e-17
 # [9,] 3.880575e-17 9.091459e-52 4.047969e-27

# =====================================================================
# Distribution Testing (Combined): Given diagrams, Create 3 dimensional vectors (using parametrization) and just use a Hotelling T2. (is there a non-normal one?)
# =====================================================================



