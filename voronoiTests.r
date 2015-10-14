# A test script to play around with Voronoi simulations and measuring hypothesis tests between them.
source('VoronoiFoam.r')
source('euler.r')
source('summarize.r')
source('tools.r')
library(abind)
library(corpcor)
library(Hotelling)

# voronoi_compilation() 
foam <- readRDS('./voronoifoam.rds')

# Each of the datasets contains a sample from some high dimensional estimation of the universe. Because this is a random variable, maximum matching doesn't really mean much. So algorithms using distance metrics are most likely meaningless.

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

# =====================================================================
# Distribution Testing [Hypothesis]
# Given diagrams, split it into 0, 1, 2 dimensions and compare their distributions via a Wilcoxon / Manning-Whitney U. 
# The question is how to compare both the birth and death times.
# =====================================================================


# =====================================================================
# Try landscapes: How do the tests on these compare?
# =====================================================================
# landMat <- gridOperation(foam, landscapeAUC)
landfxn <- dimWrapper(landscapeAUC)
landMat <- gridOperation(foam, landfxn)
# Calculate probabilities through multi-D t-test
landProba <- rep(0, setnum)
for (i in 1:setnum) {
  currproba <- hotelling.test(t(landMat[,,1]), t(landMat[,,i]))
  landProba[i] <- currproba$pval
}

# =====================================================================
# Try silhouttes: How do the tests on these compare?
# Silhouttes are so freaking fast!
# =====================================================================
# silhMat <- gridOperation(foam, silhouetteAUC)
silhfxn <- dimWrapper(silhouetteAUC)
silhMat <- gridOperation(foam, silhfxn)
# Calculate probabilities through multi-D t-test
silhProba <- rep(0, setnum)
for (i in 1:setnum) {
  currproba <- hotelling.test(t(silhMat[,,1]), t(silhMat[,,i]))
  silhProba[i] <- currproba$pval
}

# =====================================================================
# Blobbing: Representations of high dimensional space.
# =====================================================================


