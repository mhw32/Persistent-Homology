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

# We made a giant foam generation. Combine the foam[[10]] as the base line. So the last one is 0.1 percFil (but an independent sample).
foam <- readRDS('./voronoifoamfull.rds')
baseline <- readRDS('./voronoibaseline.rds')
foam[[10]] = baseline
foam <- cleanFoam(foam) # Remove the first 0th homology point.

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
  eulerProba[i] <- log(currproba$p.value)
}

# Trial 1 (size 5 grid)
# [1] 1.00000000 0.47366485 0.18146536 0.05285394 0.20642641 0.45542569 
# [7] 0.51701674 0.44226854 0.90792663

# Trial 2 (size 20 grid)
# [1] 1.000000e+00 7.139153e-01 3.722447e-01 1.938558e-01 8.116045e-03
# [6] 8.876126e-06 1.377133e-02 1.326899e-06 2.326814e-05

# Trial 3 (size 50 grid, 10k particles, removed outlier.)
# [1] 1.000000e+00 6.669889e-02 3.506863e-05 2.097664e-07 9.100366e-09
# [6] 4.485769e-15 8.174813e-18 3.680558e-15 5.486470e-22 2.476125e-01

# =====================================================================
# Try landscapes: Individual dimension at a time.
# Landscapes are the slowest.
# =====================================================================
landDimProba <- matrix(NA, nrow=setnum, ncol=3)
for (i in 0:2) {
  landDimFxn <- dimWrapper(i, landscapeAUC)
  landDimMat <- gridOperation(foam, landDimFxn)
  # Calculate probabilities with t-test
  for (j in 1:setnum) {
    currproba <- t.test(landDimMat[,1], landDimMat[,j])
    landDimProba[j, i+1] <- log(currproba$p.value)
  }
}

# Trial 1 (size 5 grid)
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

# Trial 2 (size 20 grid)
#              [,1]         [,2]       [,3]
# [1,] 1.000000e+00 1.000000e+00 1.00000000
# [2,] 1.804851e-01 1.711003e-01 0.79257006
# [3,] 4.241058e-01 1.061118e-01 0.31537189
# [4,] 9.417385e-01 4.288224e-03 0.58457744
# [5,] 4.192151e-01 5.232015e-06 0.25161676
# [6,] 9.204761e-02 8.105954e-06 0.04813845
# [7,] 2.409869e-02 1.529855e-05 0.20765584
# [8,] 1.536369e-03 6.857390e-07 0.15432086
# [9,] 4.880993e-09 7.164622e-07 0.10884305


# =====================================================================
# Try landscapes: Combining all 3 dimensions (multi-D test).
# Landscapes are the slowest.
# =====================================================================
# landMat <- gridOperation(foam, landscapeAUC)
landfxn <- allWrapper(landscapeAUC)
landMat <- gridOperation(foam, landfxn)
# Calculate probabilities through multi-D t-test
landProba <- rep(0, setnum)
for (i in 1:setnum) {
  currproba <- hotelling.test(t(landMat[,,1]), t(landMat[,,i]))
  landProba[i] <- log(currproba$pval)
}

# Trial 1 (size 5 grid)
# [1] 1.00000000 0.58836091 0.48818559 0.88391109 0.85846633 0.23193217 
# [7] 0.22234285 0.52068225 0.07237922

# Trial 2 (size 20 grid)
# [1] 1.000000e+00 2.948253e-01 2.651703e-01 3.475018e-02 1.182310e-05
# [6] 1.594319e-05 5.215882e-06 2.463658e-10 3.824774e-09]

# Trial 3 (size 50 grid)
#  [1] 1.000000e+00 1.231139e-01 6.526770e-08 5.028286e-07 5.781315e-10
 # [6] 1.571498e-11 1.112888e-12 3.730349e-14 2.886580e-15 1.215550e-01

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
    silhDimProba[j, i+1] <- log(currproba$p.value)
  }
}

# Trial 1 (size 5 grid)
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

# Trial 2 (size 20 grid)
#              [,1]         [,2]      [,3]
# [1,] 1.000000e+00 1.000000e+00 1.0000000
# [2,] 1.476883e-01 2.660581e-01 0.5933302
# [3,] 2.898063e-01 7.610234e-01 0.9669927
# [4,] 5.774533e-01 1.903724e-01 0.4685740
# [5,] 1.603077e-01 1.456238e-01 0.7165078
# [6,] 1.070692e-02 7.336313e-02 0.1997858
# [7,] 4.386880e-03 2.694533e-05 0.4295506
# [8,] 4.353056e-04 2.393692e-04 0.5926335
# [9,] 5.002512e-10 3.674042e-06 0.9206319

# Trial 3 (size 50 grid, 10k particles, removed outlier.)
#             [,1]         [,2]         [,3]
# [1,] 1.000000e+00 1.000000e+00 1.0000000000
# [2,] 7.240694e-02 4.362882e-05 0.0739629423
# [3,] 1.776881e-06 5.975643e-11 0.6462000624
# [4,] 9.516119e-06 2.987558e-12 0.7049165940
# [5,] 3.245846e-06 1.563172e-12 0.8349860776
# [6,] 4.165985e-08 2.841770e-11 0.0850386598
# [7,] 1.519206e-05 6.158313e-16 0.0008663968
# [8,] 3.698749e-04 4.714918e-14 0.0003421091
# [9,] 1.645924e-04 4.468696e-13 0.0003544215
# [10,] 3.527798e-01 2.723179e-01 0.3958034100


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
  silhProba[i] <- log(currproba$pval)
}

# Trial 1 (size 5 grid)
# [1] 1.000000000 0.890041611 0.730164341 0.387986322 0.302669180 0.137406704
# [7] 0.117997710 0.011137029 0.003031294

# Trial 2 (size 20 grid)
# [1] 1.000000e+00 1.295513e-01 7.601015e-01 5.830569e-01 5.307828e-02
# [6] 1.557875e-02 1.587168e-06 6.446052e-07 1.758393e-09

# Trial 3 (size 50 grid, 10k particles, removed outlier.)
# [1] 1.000000e+00 2.047602e-04 2.468581e-11 1.041389e-13 0.000000e+00
# [6] 1.776357e-15 0.000000e+00 0.000000e+00 0.000000e+00 4.322806e-01

# =====================================================================
# Kernel density based local two-sample comparison test 
# http://www.mvstat.net/tduong/research/publications/duong-2013-jns.pdf
# =====================================================================
test <- kde.local.test(foam[[10]][[1]], foam[[10]][[2]])
testmat <- matrix(test$pvalue, ncol=ncol(test$pvalue), nrow=nrow(test$pvalue))
test2 <- kde.local.test(foam[[10]][[1]], foam[[9]][[2]])
testmat2 <- matrix(test2$pvalue, ncol=ncol(test2$pvalue), nrow=nrow(test2$pvalue))

visualize <- function(df) {
  image(matrix(df, ncol=ncol(df), nrow=nrow(df)))
}

# =====================================================================
# Distribution Testing (Separate): Given diagrams, split it into 0, 1, 2 dimensions and compare their distributions via nonparametric test.
# This is by far the fastest method compared.
# Runs in < 1 second.
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
    distrDimProba[i, d+1] <- log(currproba$p.value)
  }
}

# Trial 1 (size 5 grid)
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

# Trial 2 (size 20 grid)
# [1,]  1.000000e+00  1.000000e+00 1.000000e+00
# [2,]  1.085204e-03  3.666638e-01 3.122177e-01
# [3,]  6.876741e-12  1.148448e-08 2.550250e-02
# [4,]  7.728523e-34  9.421209e-32 2.420532e-11
# [5,]  1.019956e-61  2.417327e-58 3.389703e-07
# [6,] 1.114050e-141 1.725675e-113 8.755756e-22
# [7,] 2.522842e-200 8.872042e-142 5.617762e-23
# [8,]  0.000000e+00 1.739559e-227 5.568554e-51
# [9,]  0.000000e+00  0.000000e+00 6.389747e-60

# Trial 3 (size 50 grid, 10k particles, removed outlier.)
#                [,1]          [,2]          [,3]
#  [1,]  1.000000e+00  1.000000e+00  1.000000e+00
#  [2,]  2.850133e-59  6.630248e-01  3.716201e-04
#  [3,] 1.627787e-198  2.127742e-01  2.475116e-07
#  [4,]  0.000000e+00  9.812215e-04  1.036986e-19
#  [5,]  0.000000e+00  7.630445e-09  9.161538e-70
#  [6,]  0.000000e+00  2.987253e-01  4.401864e-87
#  [7,]  0.000000e+00  1.530538e-19 1.212260e-105
#  [8,]  0.000000e+00  5.854761e-60 3.747478e-165
#  [9,]  0.000000e+00 2.828989e-114 8.928295e-155
# [10,]  2.072744e-02  1.619691e-01  3.522129e-01

# Combining these into one is a little harder. Concatenating across dimensions here makes much less sense. And I would have to assume normality which is just false (to do a hotelling test). But another way to combine things but lose even more information is just compare distribution of all points regardless of dimension (treat them as equivalent features).
distrProba <- rep(0, setnum)
distrList <- vector("list", setnum)
for (i in 1:setnum) { 
  distrList[[i]] <- c(distribDimStat(foam[[i]], 0), distribDimStat(foam[[i]], 1), distribDimStat(foam[[i]], 2))
}
baseline <- distrList[[1]] 
for (i in 1:setnum) {
  currproba <- wilcox.test(baseline, distrList[[i]])
  distrProba[i] <- log(currproba$p.value)
}

# [1]  1.000000e+00  2.257676e-26  8.315331e-57 8.321309e-122 3.797138e-231
# [6]  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  1.145516e-03
# Useless... not good.


##########
# KS test
##########
distrDimProba <- matrix(NA, nrow=setnum, ncol=3)
for (d in 0:2) {
  # Loop through the set and do the distribDimStat for each.
  distrDimList <- vector("list", setnum)
  for (i in 1:setnum) { distrDimList[[i]] <- distribDimStat(foam[[i]], d) }
  # Do a wilcox test for each with the baseline being the 0.1.
  baseline <- distrDimList[[1]] 
  for (i in 1:setnum) {
    currproba <- ks.test(baseline, distrDimList[[i]])
    distrDimProba[i, d+1] <- log(currproba$p.value)
  }
}

#               [,1]         [,2]         [,3]
#  [1,] 1.0000000000 1.000000e+00 1.000000e+00
#  [2,] 0.0000000000 8.306134e-12 1.743079e-08
#  [3,] 0.0000000000 0.000000e+00 0.000000e+00
#  [4,] 0.0000000000 0.000000e+00 0.000000e+00
#  [5,] 0.0000000000 0.000000e+00 0.000000e+00
#  [6,] 0.0000000000 0.000000e+00 0.000000e+00
#  [7,] 0.0000000000 0.000000e+00 0.000000e+00
#  [8,] 0.0000000000 0.000000e+00 0.000000e+00
#  [9,] 0.0000000000 0.000000e+00 0.000000e+00
# [10,] 0.0005928798 1.016147e-01 1.117919e-01

