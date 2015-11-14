source('VoronoiFoam.r')
source('euler.r')
source('summarize.r')
source('distribution.r')
source('tools.r')
library(FNN)
library(abind)
library(Hotelling)

voronoi_tests <- function(foam, baseline) {
  # Pre-setup on baseline.
  foam[[length(foam)+1]] = baseline
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
  distr_test <- function() {
    # Loop through dimensions.
    distrDimProba <- matrix(NA, nrow=setnum, ncol=3)
    for (d in 0:2) {
      # Loop through the set and do the distribDimStat for each.
      distrDimList <- vector("list", setnum)
      for (i in 1:setnum) { distrDimList[[i]] <- distribDimStat(foam[[i]], d) }
      # Do a wilcox test for each with the baseline being the 0.1.
      baseline <- distrDimList[[1]] # first index.
      for (i in 1:setnum) {
        counter <- 0
        for (j in 1:colnum) {
          counter <- counter + ks.test(baseline[[j]], distrDimList[[i]][[j]])$p.value
        }
        distrDimProba[i, d+1] <- log(counter / colnum)
      }
    }
    return(distrDimProba)
  }
  # Contour Test.
  contour_test <- function() {
    contourDimProba <- matrix(NA, nrow=setnum, ncol=3)
    for (d in 0:2) {
      # Loop through the set and grab the stat for each.
      contourDimList <- vector("list", setnum)
      for (i in 1:setnum) {
        contourDimList[[i]] <- contourDimStat(foam[[i]], d)
      }
      baseline <- contourDimList[[1]]
      # Loop through each set and each diag, calculate the differences.
      for (i in 1:setnum) {
        counter <- 0
        for (j in 1:colnum) {
          counter <- counter + max(KL.divergence(baseline[[j]], contourDimList[[i]][[j]]))
        }
        contourDimProba[i, d+1] <- counter / colnum
      }
    }
    return(contourDimProba)
  }
  # Return the tests.
  keys <- c('euler', 'indiv-land', 'all-land', 'indiv_silh', 'all-silh', 'distr', 'contour')
  tests <- c(euler_test, land_indiv_test, land_all_test, silh_indiv_test, silh_all_test, distr_test, contour_test)
  testsfxns <- vector(mode="list", length=length(keys))
  names(testsfxns) <- keys
  for (i in 1:length(keys)) {
    testsfxns[keys[i]] <- tests[i]
  }
  return(testsfxns)
}

test_wrapper <- function(foam, base1, base9) {
  # 'indiv-land', 'all-land' not included.
  keys <- c('euler', 'indiv_silh', 'all-silh', 'distr', 'contour')
  # Direct output to a file.
  sink("./results.txt", append=FALSE, split=FALSE)
  print("Starting Tests with Baseline 0.1")
  print("--------------------------------")
  t1 <- voronoi_tests(foam, base1)
  for (i in keys) {
    print(paste("Test for", i, ":"))
    currfxn <- t1[[i]]
    response <- currfxn()
    print(format_response(response))
  }
  print("\n")
  print("Starting Tests with Baseline 0.9")
  print("--------------------------------")
  t2 <- voronoi_tests(foam, base9)
  for (i in keys) {
    print(paste("Test for", i, ":"))
    currfxn <- t2[[i]]
    response <- currfxn()
    print(format_response(response))
  }
  sink()
}

format_response <- function(mat) {
  response <- ""
  if (dim(mat) == NULL) {
    for (i in mat) { response <- paste(response, i) }
  } else {
    for (i in mat) {
      for (j in i) { response <- paste(response, j) }
      response <- paste(response, i, sep="\n")
    }
  }
}
