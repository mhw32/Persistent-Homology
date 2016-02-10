library(TDA)
library(FNN)
library(Hotelling)
library(abind)

silhouetteAUC <- function(diagram, p=1, dim=1) {
  tseq <- seq(min(diagram[,2:3]), max(diagram[,2:3]), length=1000)
  silh <- silhouette(diagram, p=p, dimension=dim, tseq)
  return(integrate(tseq, silh))
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
  inner <- function(diagram) { return(fxn(diagram, p=p, dim=dim)) }
  return(inner)
}

allWrapper <- function(p, fxn) {
  inner <- function(diagram, mindim=0, maxdim=2) {
    areas <- sapply(seq(from=mindim, to=maxdim, by=1), function(i) {
      fxn(diagram, p=p, dim=i)
    })
    return(areas)
  }
  return(inner)
}

silh_tests <- function(foam, baseline, p) {
  # Pre-setup on baseline.
  foam[[length(foam)+1]] = baseline
  # Must do prior to cleaning.
  if (norm == TRUE) { foam <- normFoam(foam) }
  # Now clean
  foam <- cleanFoam(foam)
  setnum <- length(foam)
  colnum <- length(foam[[1]])
  basenum <- setnum # This represents the added baseline.

  silh_indiv_test <- function(p) {
    silhDimProba <- matrix(NA, nrow=setnum, ncol=3)
    for (i in 0:2) {
      silhDimFxn <- dimWrapper(p, i, silhouetteAUC())
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
  silh_all_test <- function(p) {
    silhfxn <- allWrapper(p, silhouetteAUC())
    silhMat <- gridOperation(foam, silhfxn)
    # Calculate probabilities through multi-D t-test
    silhProba <- rep(0, setnum)
    for (i in 1:setnum) {
      currproba <- hotelling.test(t(silhMat[,,basenum]), t(silhMat[,,i]))
      silhProba[i] <- log(currproba$pval)
    }
    return(silhProba)
  }

  keys <- c('indiv_silh', 'all-silh')
  tests <- c(silh_indiv_test, silh_all_test)
  testsfxns <- vector(mode="list", length=length(keys))
  names(testsfxns) <- keys

  for (i in 1:length(keys)) {
    testsfxns[keys[i]] <- tests[i]
  }

  return(testsfxns)
}

for (i in 1:5) {
  for (p in seq(1,5,0.1)) {
    print(paste("Processing for tuning parameter ", p, " and iteration ", i))

    foam <- readRDS(paste('./saved_states/test_set/foam', i, '.rds', sep=''))
    base <- readRDS(paste('./saved_states/test_set/baseline', i, '-0.1.rds', sep=''))

    keys <- c('indiv_silh', 'all-silh')
    sink(paste("./saved_states/silh_results/results-iter-", i, "-tune-", p, sep=""), append=FALSE, split=FALSE)

    print("--------------------------------")
    t <- silh_tests(foam, base)

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
}

