library(rrcov)

corrVoronoiTestSuite <- function(baseStats, foamStats) {
	# define some parameters
	dims     <- dim(foamStats)
	nreps    <- dims[1]
	npercfil <- dims[2]
	nsamples <- dims[3]

	hotelling_pval <- function(p) {
		pval <- pf(1/p$statistic[2],p$parameter[2], p$parameter[1])
		return(pval)
	}

	corrDimTest <- function(baseStats, foamStats) {
		corrDimProba <- array(NA, dim=c(nreps, npercfil))
		for (i in 1:nreps) {
			for (j in 1:npercfil) {
				pval <- t.test(baseStats[i,], foamStats[i,j,], paired=TRUE)
				corrDimProba[i, j] <- log(pval$p.value)
			}
		}
		return(corrDimProba)
	}

	# prep the functions and return them 
	# in a dictionary like structure.

	keys  <- c('indiv-corr')
  	tests <- c(corrDimTest)
  
  	fxns <- vector(mode="list", length=length(keys))
  	names(fxns) <- keys
	for (i in 1:length(keys)) {
		fxns[keys[i]] <- tests[i]
	}
  
  	return(fxns)
}

corrSimuTestSuite <- function(cdmStats, wdmStats) {
	# define some parameters
	dims     <- dim(wdmStats)
	nreps    <- dims[1]
	nsamples <- dims[2]

	hotelling_pval <- function(p) {
		pval <- pf(1/p$statistic[2],p$parameter[2], p$parameter[1])
		return(pval)
	}

	corrDimTest <- function(cdmStats, wdmStats) {
		corrDimProba <- array(NA, dim=c(nreps))
		for (i in 2:nreps) {
			num_to_keep <- i^3
			pval <- t.test(cdmStats[i,1:num_to_keep], wdmStats[i,1:num_to_keep], paired=TRUE)
			corrDimProba[i] <- log(pval$p.value)
		}
		return(corrDimProba)
	}

	# prep the functions and return them 
	# in a dictionary like structure.

	keys  <- c('indiv-corr')
  	tests <- c(corrDimTest)
  
  	fxns <- vector(mode="list", length=length(keys))
  	names(fxns) <- keys
	for (i in 1:length(keys)) {
		fxns[keys[i]] <- tests[i]
	}
  
  	return(fxns)
}

