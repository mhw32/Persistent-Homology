library(rrcov)

corrTestSuite <- function(baseStats, foamStats) {
	# define some parameters
	dims     <- dim(foamStats)
	nreps    <- dims[1]
	npercfil <- dims[2]
	nsamples <- dims[3]
	ndims    <- dims[4]

	hotelling_pval <- function(p) {
		pval <- pf(1/p$statistic[2],p$parameter[2], p$parameter[1])
		return(pval)
	}

	corrDimTest <- function(baseStats, foamStats) {
		corrDimProba <- array(NA, dim=c(nreps, npercfil, ndims))
		for (i in 1:nreps) {
			for (j in 1:npercfil) {
				for (k in 1:ndims) {
					pval <- t.test(baseStats[i,,k], foamStats[i,j,,k], paired=TRUE)
					corrDimProba[i, j, k] <- log(pval$p.value)
				}
			}
		}
		return(corrDimTest)
	}

	corrParallelTest <- function(baseStats, foamStats) {
		corrProba <- array(NA, dim=c(nreps, npercfil))
		for (i in 1:nreps) {
			for (j in 1:npercfil) {
				pval <- T2.test(baseStats[i,,] - foamStats[i,j,,])
				corrProba[i, j] <- log(hotelling_pval(pval))
			}
		}
		return(corrParallelTest)
	}

	# prep the functions and return them 
	# in a dictionary like structure.

	keys  <- c('indiv-corr', 'all-corr')
  	tests <- c(corrDimTest, corrParallelTest)
  
  	fxns <- vector(mode="list", length=length(keys))
  	names(fxns) <- keys
	for (i in 1:length(keys)) {
		fxns[keys[i]] <- fxns[i]
	}
  
  	return(fxns)
}
