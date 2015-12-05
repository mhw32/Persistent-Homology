source('VoronoiFoam.r')
source('euler.r')
source('summarize.r')
source('twosample.r')
source('distribution.r')
source('tools.r')
library(FNN)
library(abind)
library(Hotelling)
library(ks)

voronoi_tests <- function(foam, baseline, norm=FALSE) {
  # Pre-setup on baseline.
  foam[[length(foam)+1]] = baseline
  # Must do prior to cleaning.
  if (norm=TRUE) { foam <- normFoam(foam) }
  # Now clean
  foam <- cleanFoam(foam)
  setnum <- length(foam)
  colnum <- length(foam[[1]])
  basenum <- setnum # This represents the added baseline.
  # Euler Characteristic.
  euler_test <- function() {
    eulerMat <- gridOperation(foam, eulerIntegration)
    eulerProba <- rep(0, setnum)
    for (i in 1:setnum) {
      currproba <- t.test(eulerMat[,basenum], eulerMat[,i])
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
        currproba <- t.test(landDimMat[,basenum], landDimMat[,j])
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
      currproba <- hotelling.test(t(landMat[,,basenum]), t(landMat[,,i]))
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
        currproba <- t.test(silhDimMat[,basenum], silhDimMat[,j])
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
      currproba <- hotelling.test(t(silhMat[,,basenum]), t(silhMat[,,i]))
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
      baseline <- distrDimList[[basenum]] # baseline index.
      for (i in 1:setnum) {
        counter <- 0
        for (j in 1:colnum)
          counter <- counter + ks.test(baseline[[j]], distrDimList[[i]][[j]])$p.value
        distrDimProba[i, d+1] <- log(counter / colnum)
      }
    }
    return(distrDimProba)
  }
  # Contour Test.
  contour_test <- function() {
    contourDimProba <- matrix(NA, nrow=setnum, ncol=3)
    baseline <- foam[[basenum]]
    for (d in 0:2) {
      # Perform a permutation test for each thing.
      for (i in 1:setnum) {
        # Artificial labels: 0 (base), 1 (percFil)
        L <- c(rep(0, colnum), rep(1, colnum))
        # Pass a contatenated vector of contours maps.
        X <- c(contourDimStat(baseline, d), contourDimStat(foam[[i]], d))
        contourDimProba[i, d+1] <- log(1 - permutationTest(500, X, L, 1))
      }
    }
    # Here, we have 2 sample t-test results comparing each 2D shape with the baseline.
    return(contourDimProba)
  }
  # global kde test. Similar to contour but a chi^2 instead.
  global_kde_test <- function() {
    globalDimProba <- matrix(NA, nrow=setnum, ncol=3)
    baseline <- foam[[basenum]]
    for (d in 0:2) {
      for (i in 1:setnum) {
        counter <- 0
        base_edit <- globalDimStat(baseline, d)
        foam_edit <- globalDimStat(foam[[i]], d)
        for (j in 1:colnum)
          counter <- counter + ks::kde.test(base_edit[[j]], foam_edit[[j]])$pvalue
        globalDimProba[i, d+1] <- log(counter / colnum)
      }
    }
    return(globalDimProba)
  }

  # Return the tests.
  keys <- c('euler', 'indiv-land', 'all-land', 'indiv_silh', 'all-silh', 'distr', 'contour', 'global-kde')
  tests <- c(euler_test, land_indiv_test, land_all_test, silh_indiv_test, silh_all_test, distr_test, contour_test, global_kde_test)
  testsfxns <- vector(mode="list", length=length(keys))
  names(testsfxns) <- keys
  for (i in 1:length(keys)) {
    testsfxns[keys[i]] <- tests[i]
  }
  return(testsfxns)
}

test_wrapper <- function(foam, base, ext) {
  # 'indiv-land', 'all-land' not included.
  keys <- c('euler', 'indiv_silh', 'all-silh', 'distr', 'contour', 'global-kde')
  # Direct output to a file.
  sink(paste("./saved_states/results-", ext, ".txt", sep=""), append=FALSE, split=FALSE)
  print("--------------------------------")
  t <- voronoi_tests(foam, base, norm=TRUE)
  for (i in keys) {
    print(paste("Test for", i, ":"))
    currfxn <- t[[i]]
    response <- currfxn()
    print(response)
    print("")
  }
  print("--------------------------------")
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
