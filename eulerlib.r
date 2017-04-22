source('euler.r')
source('tools.r')

euler_tests <- function(foam, baseline, norm=FALSE, run_base=TRUE) {
  foam[[length(foam)+1]] = baseline
  # normalize foam if needed
  if (norm == TRUE) { foam <- normFoam(foam) }
  # no need to remove vestigial for Euler tests
  setnum <- length(foam)
  colnum <- length(foam[[1]])
  basenum <- setnum
  euler_length <- 1000
  # if we don't need to compare base with base
  if (run_base == FALSE) { setnum <- setnum - 1 }

  baseline_storage <- array(0, dim=c(colnum, euler_length))
  foam_storage <- array(0,dim=c(setnum, colnum, euler_length))

  for (i in 1:colnum) {  # loop over each col
    baseline_euler <- eulerVanillaFunction(baseline[[i]])
    baseline_storage[i,] <- baseline_euler
    for (t in 1:setnum) {  # loop over each param
      foam_euler <- eulerVanillaFunction(foam[[t]][[i]])
      foam_storage[t,i,] <- foam_euler
    }
  }

  return(list(
    "baseline"=baseline_storage, 
    "foam"=foam_storage,
  ))
}

# apply euler_tests to an entire directory
euler_directory_test <- function(directory, num=100, norm=FALSE) {
  for (i in 1:100) {
    foam <- readRDS(paste(directory, '/foam', i, '.rds', sep=''))
    base <- readRDS(paste(directory, '/baseline', i, '-0.1.rds', sep=''))
    response <- euler_tests(foam, base, norm=norm)
    if (i == 1) {
      foam_storage <- response["foam"]
      baseline_storage <- response["baseline"]
    } else {
      foam_storage <- rbind(foam_storage, response["foam"])
      baseline_storage <- rbind(baseline_storage, response["baseline"])
    }
  }
  eulerComparePlot(baseline_storage, foam_storage)
}

# plot the average EC function (point wise average) and a 1 sd band around the 
# average against the same thing for the null model.
eulerComparePlot <- function(baselineMat, foamMat) {
  # @param baselineMat : 1500 
}