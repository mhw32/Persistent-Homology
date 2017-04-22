library(TDA)
library(FNN)
library(Hotelling)
library(abind)
source("euler.r")
source("tools.r")

silhouetteAUC <- function(diagram, p=1, dim=1) {
  tseq <- seq(min(diagram[,2:3]), max(diagram[,2:3]), length=1000)
  silh <- silhouette(diagram, p=p, dimension=dim, tseq)
  return(integrate(tseq, silh))
}

sileuler <- function(diagram, p=1, length=1000) {
  tseq <- seq(min(diagram[,2:3]), max(diagram[,2:3]), length=length)

  # Calculate silhouette for each dimension.
  s0 <- silhouette(diagram, p=p, dimension=0, tseq)
  s1 <- silhouette(diagram, p=p, dimension=1, tseq)
  s2 <- silhouette(diagram, p=p, dimension=2, tseq)

  # Calculate the alternating 'betti'.
  seuler <- rep(0, length)
  for (i in seq(length)) {
    seuler[i] = s0[i] - s1[i] + s2[i]
  }

  # Integrate the absolute value.
  score <- integrate(tseq, abs(seuler))
  return(score)
}

gridOperation <- function(foam, fxn) {
  # Get dimension of foams.
  setnum <- length(foam)
  colnum <- length(foam[[1]])
  # Do something to the first column.
  matrix <- sapply(seq(1:colnum), function(j) {
    fxn(foam[[1]][[j]])
  })
  maxdim <- if (is.null(dim(matrix))) 1 else length(dim(matrix))
  # Do something to the rest of the columns.
  for (i in 2:setnum) {
    vector <- sapply(seq(1:colnum), function(j) {
      fxn(foam[[i]][[j]])
    })
    matrix <- abind(matrix, vector, along=maxdim+1)
  }
  colnames(matrix) <- NULL
  return(matrix)
}

dimWrapper <- function(p, dim, fxn) {
  inner <- function(diagram) { return(fxn(diagram, p, dim)) }
  return(inner)
}

allWrapper <- function(p, fxn) {
  inner <- function(diagram, mindim=0, maxdim=2) {
    areas <- sapply(seq(from=mindim, to=maxdim, by=1), function(i) {
      fxn(diagram, p, i)
    })
    return(areas)
  }
  return(inner)
}

customWrapper <- function(p, fxn) {
  inner <- function(diagram) { return(fxn(diagram, p=p)) }
  return(inner)
}

silh_tests <- function(foam, baseline, p) {
  # Pre-setup on baseline.
  foam[[length(foam)+1]] = baseline
  # Must do prior to cleaning.
  # if (norm == TRUE) { foam <- normFoam(foam) }
  # Now clean
  foam <- cleanFoam(foam)
  setnum <- length(foam)
  colnum <- length(foam[[1]])
  basenum <- setnum # This represents the added baseline.

  silh_indiv_test <- function(p) {
    silhDimProba <- matrix(NA, nrow=setnum, ncol=3)
    for (i in 0:2) {
      silhDimFxn <- dimWrapper(p, i, silhouetteAUC)
      silhDimMat <- gridOperation(foam, silhDimFxn)
      # Calculate probabilities with t-test
      for (j in 1:setnum) {
        currproba <- t.test(silhDimMat[,basenum], silhDimMat[,j], paired=TRUE)
        silhDimProba[j, i+1] <- log(currproba$p.value)
      }
    }
    return(silhDimProba)
  }

  # Combined Silhouette Test.
  silh_all_test <- function(p) {
    silhfxn <- allWrapper(p, silhouetteAUC)
    silhMat <- gridOperation(foam, silhfxn)
    # Calculate probabilities through multi-D t-test
    silhProba <- rep(0, setnum)
    for (i in 1:setnum) {
      currproba <- T2.test(t(silhMat[,,basenum]) - t(silhMat[,,i]))
      silhProba[i] <- log(hotelling_pval(currproba))
    }
    return(silhProba)
  }

  silh_euler_test <- function(p) {
    silhEulerFxn <- customWrapper(p, sileuler)
    silhEulerMat <- gridOperation(foam, silhEulerFxn)

    silhEulerProba <- rep(0, setnum)
    for (i in 1:setnum) {
      currproba <- t.test(
        silhEulerMat[,basenum], 
        silhEulerMat[,i],
        paired=TRUE
      )
      silhEulerProba[i] <- log(currproba$p.value)
    }
    return(silhEulerProba)
  }

  keys <- c('indiv_silh', 'all-silh', 'silh-euler')
  tests <- c(silh_indiv_test, silh_all_test, silh_euler_test)
  testsfxns <- vector(mode="list", length=length(keys))
  names(testsfxns) <- keys

  for (i in 1:length(keys)) {
    testsfxns[keys[i]] <- tests[i]
  }

  return(testsfxns)
}


silh_wrapper <- function() {
  for (i in 1:5) {
    for (p in seq(0,1,0.1)) {
      print(paste("Processing for tuning parameter ", p, " and iteration ", i))
      foam <- readRDS(paste('./saved_states/test_set/foam', i, '.rds', sep=''))
      base <- readRDS(paste('./saved_states/test_set/baseline', i, '-0.1.rds', sep=''))

      keys <- c('silh-euler')
      sink(paste("./saved_states/silh_euler_results/results-iter-", i, "-tune-", p, sep=""), append=FALSE, split=FALSE)

      print("--------------------------------")
      t <- silh_tests(foam, base, p)

      for (k in keys) {
        print(paste("Test for", k, ":"))
        currfxn <- t[[k]]
        response <- currfxn(p)
        print(response)
        print("")
      }
      print("--------------------------------")
      sink()
    }
  }
}

silh_wrapper()

