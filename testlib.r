source('VoronoiFoam.r')
source('euler.r')
source('summarize.r')
source('distribution.r')
source('tools.r')
library(abind)
library(Hotelling)

voronoi_tests <- function(foam, baseline, ticker=1:6) {
  # Pre-setup on baseline.
  foam[[length(foam)]] = baseline
  foam <- cleanFoam(foam)
  setnum <- length(foam)
  colnum <- length(foam[[1]])
  # Euler Characteristic.
  euler_test <- function() {
    eulerMat <- gridOperation(foam, eulerIntegration)
    eulerProba <- rep(0, setnum)
    for (i in 1:setnum) {
      currproba <- t.test(eulerMat[,1], eulerMat[,i])
      eulerProba[i] <- log(currproba$p.value)
    }
    return(eulerProba)
  }
  # Individual Landscape Test.
  land_indiv_test <- function() {
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
    return(landDimProba)
  }
  # Combined Landscape Test.
  land_all_test <- function() {
    landfxn <- allWrapper(landscapeAUC)
    landMat <- gridOperation(foam, landfxn)
    # Calculate probabilities through multi-D t-test
    landProba <- rep(0, setnum)
    for (i in 1:setnum) {
      currproba <- hotelling.test(t(landMat[,,1]), t(landMat[,,i]))
      landProba[i] <- log(currproba$pval)
    }
    return(landProba)
  }
  # Individual Silhouette Test.
  silh_indiv_test <- function() {
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
    return(silhDimProba)
  }
  # Combined Silhouette Test.
  silh_all_test <- function() {
    silhfxn <- allWrapper(silhouetteAUC)
    silhMat <- gridOperation(foam, silhfxn)
    # Calculate probabilities through multi-D t-test
    silhProba <- rep(0, setnum)
    for (i in 1:setnum) {
      currproba <- hotelling.test(t(silhMat[,,1]), t(silhMat[,,i]))
      silhProba[i] <- log(currproba$pval)
    }
    return(silhProba)
  }
  # Distribution Test.
  distr_test <- function(test='wilcox') {
    # Loop through dimensions.
    distrDimProba <- matrix(NA, nrow=setnum, ncol=3)
    for (d in 0:2) {
      # Loop through the set and do the distribDimStat for each.
      distrDimList <- vector("list", setnum)
      for (i in 1:setnum) { distrDimList[[i]] <- distribDimStat(foam[[i]], d) }
      # Do a wilcox test for each with the baseline being the 0.1.
      baseline <- distrDimList[[1]] 
      for (i in 1:setnum) {
        if (test == 'wilcox') {
          currproba <- wilcox.test(baseline, distrDimList[[i]])
        } else {
          currproba <- ks.test(baseline, distrDimList[[i]])
        } 
        distrDimProba[i, d+1] <- log(currproba$p.value)
      }
    }
    return(distrDimProba)
  }
  # Run the tests.
  keys <- c('euler', 'indiv-land', 'all-land', 'indiv_silh', 'all-silh', 'distr-test')
  tests <- c(euler_test, land_indiv_test, land_all_test, silh_indiv_test, silh_all_test, distr_test)
  # Create an empty list and initialize it.
  tests_results <- vector(mode="list", length=length(keys))
  names(tests_results) <- keys
  # Loop through them.
  for (i in ticker) { 
    tests_results[[keys[i]]] <- tests[i]() 
  }
  return(tests_results)
}



