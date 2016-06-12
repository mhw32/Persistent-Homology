library(rrcov)

corrVoronoiTestSuite <- function(baseStats, foamStats) {
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
		return(corrDimProba)
	}

	corrParallelTest <- function(baseStats, foamStats) {
		corrProba <- array(NA, dim=c(nreps, npercfil))
		for (i in 1:nreps) {
			for (j in 1:npercfil) {
				pval <- T2.test(baseStats[i,,] - foamStats[i,j,,])
				corrProba[i, j] <- log(hotelling_pval(pval))
			}
		}
		return(corrProba)
	}

	# prep the functions and return them 
	# in a dictionary like structure.

	keys  <- c('indiv-corr', 'all-corr')
  	tests <- c(corrDimTest, corrParallelTest)
  
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
	ndims    <- dims[3]

	hotelling_pval <- function(p) {
		pval <- pf(1/p$statistic[2],p$parameter[2], p$parameter[1])
		return(pval)
	}

	corrDimTest <- function(cdmStats, wdmStats) {
		corrDimProba <- array(NA, dim=c(nreps, ndims))
		for (i in 1:nreps) {
			num_to_keep <- i^3
			for (j in 1:ndims) {
				pval <- t.test(cdmStats[i,1:num_to_keep,j], wdmStats[i,1:num_to_keep,j], paired=TRUE)
				corrDimProba[i, j] <- log(pval$p.value)
			}
		}
		return(corrDimProba)
	}

	corrParallelTest <- function(cdmStats, wdmStats) {
		corrProba <- array(NA, dim=c(nreps))
		for (i in 1:nreps) {
			num_to_keep <- i^3
			pval <- T2.test(baseStats[i,1:num_to_keep,] - foamStats[i,1:num_to_keep,])
			corrProba[i] <- log(hotelling_pval(pval))
		}
		return(corrProba)
	}

	# prep the functions and return them 
	# in a dictionary like structure.

	keys  <- c('indiv-corr', 'all-corr')
  	tests <- c(corrDimTest, corrParallelTest)
  
  	fxns <- vector(mode="list", length=length(keys))
  	names(fxns) <- keys
	for (i in 1:length(keys)) {
		fxns[keys[i]] <- tests[i]
	}
  
  	return(fxns)
}


applyCorrVoronoiTestSuite <- function() {
	load('output/base_corr.gzip')
	load('output/foam_corr.gzip')
	
	tlibs <- corrTestSuite(base, foam)
	corrDimTest <- tlibs$`indiv-corr`
	corrParallelTest <- tlibs$`all-corr`

	corrDimProba <- corrDimTest(base, foam)
	corrParallelProba <- corrParallelTest(base, foam)

	saveRDS(corrDimProba, file='output/indiv_proba.rds')
	saveRDS(corrParallelProba, file='output/parallel_proba.rds')
}

applyCorrSimuTestSuite <- function() {
	load('output/cdm_corr.gzip')
	load('output/wdm_corr.gzip')
	
	tlibs <- corrTestSuite(cdm, wdm)
	corrDimTest <- tlibs$`indiv-corr`
	corrParallelTest <- tlibs$`all-corr`

	corrDimProba <- corrDimTest(cdm, wdm)
	corrParallelProba <- corrParallelTest(cdm, wdm)

	saveRDS(corrDimProba, file='output/simu_indiv_proba.rds')
	saveRDS(corrParallelProba, file='output/simu_parallel_proba.rds')
}


